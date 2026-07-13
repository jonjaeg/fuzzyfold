//! [`PseudoEnergyModel`] implementation for [`ViennaRNA`].

use ff_structure::NAIDX;

use crate::{Base, EnergyError, EnergyModel, NearestNeighborLoop};
use crate::pseudoknots::{
    Loop, LoopType, ClosingDescriptor, LocationStatus,
    PseudoEnergyModel, PseudoloopContext,
    single_pair, double_pair, collect_single_branches,
};
use super::ViennaRNA;

impl PseudoEnergyModel for ViennaRNA {
    fn energy_of_pseudo_loop(
        &self,
        sequence: &[Base],
        lp: &Loop,
    ) -> Result<i32, EnergyError> {
        match lp.loop_type {
            // ── Step 6: SpanBand stacks / interior loops use scaled Turner energies ──
            LoopType::Stack | LoopType::Interior | LoopType::Bulge => {
                let (i, j) = single_pair(lp.closing)?;
                let (k, l) = lp.inner.ok_or(EnergyError::InvalidClosingPair)?;
                let raw = self.energy_of_loop(sequence, &NearestNeighborLoop::Interior {
                    closing: (i as NAIDX, j as NAIDX),
                    inner:   (k as NAIDX, l as NAIDX),
                })?;
                if lp.location == LocationStatus::SpanBand {
                    if let Some(pk) = &self.pk_params {
                        let scale = if lp.loop_type == LoopType::Stack {
                            pk.e_stp
                        } else {
                            pk.e_intp
                        };
                        return Ok((raw as f64 * scale).round() as i32);
                    }
                }
                Ok(raw)
            }

            LoopType::Hairpin => {
                let (i, j) = single_pair(lp.closing)?;
                self.energy_of_loop(sequence, &NearestNeighborLoop::Hairpin {
                    closing: (i as NAIDX, j as NAIDX),
                })
            }

            // ── Step 5: SpanBand multiloops use ap/bp/cp instead of Turner ──
            LoopType::Multiloop => {
                if lp.location == LocationStatus::SpanBand {
                    if let Some(pk) = &self.pk_params {
                        let n_branches = lp.children.len();
                        let n_unpaired = lp.unpaired_5p + lp.unpaired_3p;
                        return Ok(
                            pk.ap
                            + pk.bp * n_branches as i32
                            + pk.cp * n_unpaired  as i32,
                        );
                    }
                }
                let (i, j) = single_pair(lp.closing)?;
                let branches: Vec<(NAIDX, NAIDX)> = collect_single_branches(&lp.children)?
                    .into_iter()
                    .map(|(a, b)| (a as NAIDX, b as NAIDX))
                    .collect();
                self.energy_of_loop(sequence, &NearestNeighborLoop::Multibranch {
                    closing: (i as NAIDX, j as NAIDX),
                    branches,
                })
            }

            LoopType::External => {
                // Each branch (Single or Double descriptor) contributes independently:
                // AU/GU terminal penalty + dangling ends based on adjacent sequence context.
                // Crossing pairs from a Double descriptor cannot share a single slice pass,
                // so every pair is evaluated by slicing [..=i] and [j..] individually.
                let mut en = 0;
                for cd in &lp.children {
                    match cd {
                        ClosingDescriptor::Single((a, b)) => {
                            en += self.exterior(&[&sequence[..=*a], &sequence[*b..]])?;
                        }
                        ClosingDescriptor::Double((a1, b1), (a2, b2)) => {
                            en += self.exterior(&[&sequence[..=*a1], &sequence[*b1..]])?;
                            en += self.exterior(&[&sequence[..=*a2], &sequence[*b2..]])?;
                        }
                    }
                }
                Ok(en)
            }

            // ── Steps 1-4: Full 6-term D&P Pseudoloop formula ──
            LoopType::Pseudoloop => {
                let _ = double_pair(lp.closing)?; // validate descriptor
                if let Some(pk) = &self.pk_params {
                    // E = init(ctx) + pb×n_bands + pup×(n_loop1+n_loop2) + pps×n_nested
                    let init = match lp.pk_context.unwrap_or(PseudoloopContext::External) {
                        PseudoloopContext::External   => pk.init_external,
                        PseudoloopContext::Multiloop  => pk.init_multiloop,
                        PseudoloopContext::Pseudoloop => pk.init_pseudoloop,
                    };
                    let n_unpaired = lp.n_loop1 + lp.n_loop2;
                    Ok(init
                        + pk.pb  * lp.n_bands  as i32
                        + pk.pup * n_unpaired  as i32
                        + pk.pps * lp.n_nested as i32)
                } else {
                    // Fallback: simplified multiloop approximation.
                    let n_stems = lp.children.len() + 2;
                    Ok(self.ml_closing + self.ml_intern * n_stems as i32)
                }
            }
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::{NucleotideVec, parse_structure};
    use crate::parameters::RNA_EXTENDED;

    fn model() -> ViennaRNA {
        ViennaRNA::from_thermo_params(&RNA_EXTENDED, 37.0)
    }

    /// A non-pseudoknotted structure should give the same total energy via
    /// `energy_of_pseudoknotted_structure` and the standard `energy_of_structure`.
    #[test]
    fn test_non_pk_matches_standard() {
        use ff_structure::PairTable;

        let seq = NucleotideVec::try_from_rna("GGGAAACCC").unwrap();
        let dot = "(((...)))";
        let pt    = PairTable::try_from(dot).unwrap();
        let loops = parse_structure(dot).unwrap();

        let m = model();
        let e_standard = m.energy_of_structure(&seq, &pt).unwrap();
        let e_pseudo   = m.energy_of_pseudoknotted_structure(&seq, &loops).unwrap();
        assert_eq!(e_standard, e_pseudo,
            "non-PK energy mismatch: standard={e_standard}, pseudo={e_pseudo}");
    }

    /// An H-type pseudoknot should evaluate without error.
    #[test]
    fn test_h_type_pseudoknot_runs() {
        // GGGCCCCCCAAAGGG pairs as ((([[[)))...]]]
        let seq  = NucleotideVec::try_from_rna("GGGCCCCCCAAAGGG").unwrap();
        let loops = parse_structure("((([[[)))...]]]").unwrap();

        let m = model();
        let result = m.energy_of_pseudoknotted_structure(&seq, &loops);
        println!("Pseudoknot energy: {:?}", result);
        assert!(result.is_ok(), "pseudoknot energy evaluation failed: {:?}", result);
    }

    /// Smoke test: dp03 junction formula is pure arithmetic, no Turner tables needed.
    ///
    /// Structure: GCGAUUUCUGACCGCUUUUUUGUCAG / [[[....(((((]]]......]))))
    ///   n_bands=2, n_loop1=4, n_loop2=6, n_nested=0
    ///   init_ext=960 + pb=20×2 + pup=10×(4+6) + pps=10×0 = 960+40+100 = 1100
    #[test]
    fn test_dp03_junction_formula() {
        use crate::parameters::{RNA_TURNER_2004, RNA_DP03};

        let seq   = NucleotideVec::try_from_rna("GCGAUUUCUGACCGCUUUUUUGUCAG").unwrap();
        let loops = parse_structure("[[[....(((((]]]......)))))").unwrap();
        let m = ViennaRNA::from_thermo_params(&RNA_TURNER_2004, 37.0)
            .with_pseudoknot_params(RNA_DP03);
        let pk_loop = loops.iter().find(|l| l.loop_type == LoopType::Pseudoloop).unwrap();

        assert_eq!(m.energy_of_pseudo_loop(&seq, pk_loop).unwrap(), 1100);
    }

    /// Smoke test: dp09 junction formula is pure arithmetic, no Turner tables needed.
    ///
    /// Same structure as above.
    ///   init_ext=−138 + pb=246×2 + pup=6×(4+6) + pps=96×0 = −138+492+60 = 414
    #[test]
    fn test_dp09_junction_formula() {
        use crate::parameters::{RNA_TURNER_2004, RNA_DP09};

        let seq   = NucleotideVec::try_from_rna("GCGAUUUCUGACCGCUUUUUUGUCAG").unwrap();
        let loops = parse_structure("[[[....(((((]]]......)))))").unwrap();
        let m = ViennaRNA::from_thermo_params(&RNA_TURNER_2004, 37.0)
            .with_pseudoknot_params(RNA_DP09);
        let pk_loop = loops.iter().find(|l| l.loop_type == LoopType::Pseudoloop).unwrap();

        assert_eq!(m.energy_of_pseudo_loop(&seq, pk_loop).unwrap(), 414);
    }

    /// Full dp03 total energy for an H-type pseudoknot.
    ///
    /// HotKnots v2 reference: −2.42 kcal/mol = −242 dcal/mol.
    /// Tolerance ±1 dcal/mol for f64→i32 rounding in stack scaling.
    
    #[test]
    fn test_dp03_total_energy() {
        use crate::parameters::{RNA_TURNER_2004 , RNA_DP03};

        let seq   = NucleotideVec::try_from_rna("GCGAUUUCUGACCGCUUUUUUGUCAG").unwrap();
        let loops = parse_structure("(((....[[[[[)))......]]]]]").unwrap();
        
        let m = ViennaRNA::from_thermo_params(&RNA_TURNER_2004, 37.0)
            .with_pseudoknot_params(RNA_DP03);
        let e = m.energy_of_pseudoknotted_structure(&seq, &loops).unwrap();
        assert!((e - (-242)).abs() <= 1,
            "dp03 total energy: expected ~−242 dcal/mol, got {e}");
    }

    /// Full dp09+mt09 total energy for an H-type pseudoknot.
    ///
    /// HotKnots v2 reference (with mt09 NN + dp09 PK): −5.11 kcal/mol = −511 dcal/mol.
    /// dp09 PK parameters were trained jointly with the mt09 NN parameters; using them
    /// with Turner 2004 tables instead gives systematically wrong results (~−3.43 kcal/mol).
    /// Tolerance ±1 dcal/mol for f64→i32 rounding in stack scaling.
    #[test]
    fn test_dp09_total_energy() {
        use crate::parameters::{RNA_MT09, RNA_DP09};

        let seq   = NucleotideVec::try_from_rna("GCGAUUUCUGACCGCUUUUUUGUCAG").unwrap();
        let loops = parse_structure("[[[....(((((]]]......)))))").unwrap();
        let m = ViennaRNA::from_andrunescu_params(&RNA_MT09)
            .with_pseudoknot_params(RNA_DP09);
        let e = m.energy_of_pseudoknotted_structure(&seq, &loops).unwrap();
        assert!((e - (-511)).abs() <= 1,
            "dp09+mt09 total energy: expected ~−511 dcal/mol, got {e}");
    }
}
