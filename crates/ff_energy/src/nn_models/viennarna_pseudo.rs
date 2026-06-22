//! [`PseudoEnergyModel`] implementation for [`ViennaRNA`].

use ff_structure::NAIDX;

use crate::{Base, EnergyError, EnergyModel, NearestNeighborLoop};
use crate::pseudoknots::{
    Loop, LoopType, ClosingDescriptor,
    PseudoEnergyModel, single_pair, double_pair, collect_single_branches,
};
use super::ViennaRNA;

impl PseudoEnergyModel for ViennaRNA {
    fn energy_of_pseudo_loop(
        &self,
        sequence: &[Base],
        lp: &Loop,
    ) -> Result<i32, EnergyError> {
        match lp.loop_type {
            LoopType::Stack | LoopType::Interior | LoopType::Bulge => {
                let (i, j) = single_pair(lp.closing)?;
                let (k, l) = lp.inner.ok_or(EnergyError::InvalidClosingPair)?;
                self.energy_of_loop(sequence, &NearestNeighborLoop::Interior {
                    closing: (i as NAIDX, j as NAIDX),
                    inner:   (k as NAIDX, l as NAIDX),
                })
            }

            LoopType::Hairpin => {
                let (i, j) = single_pair(lp.closing)?;
                self.energy_of_loop(sequence, &NearestNeighborLoop::Hairpin {
                    closing: (i as NAIDX, j as NAIDX),
                })
            }

            LoopType::Multiloop => {
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
                let n = sequence.len();
                // Double descriptors (top-level crossing pairs) cannot be sliced as
                // monotone branches — skip them here.
                // TODO: add terminal AU/GU penalties for Double-descriptor children.
                let branches: Vec<(NAIDX, NAIDX)> = lp.children.iter()
                    .filter_map(|cd| match cd {
                        ClosingDescriptor::Single((a, b)) => Some((*a as NAIDX, *b as NAIDX)),
                        ClosingDescriptor::Double(..) => None,
                    })
                    .collect();
                self.energy_of_loop(sequence, &NearestNeighborLoop::Exterior {
                    ends: (0, (n - 1) as NAIDX),
                    branches,
                })
            }

            LoopType::Pseudoloop => {
                // No standard NN parameter exists for crossing-pair loops.
                // Approximation: treat as a multiloop with two closing pairs.
                //   E ≈ ml_closing + ml_intern × (n_children + 2)
                // Unpaired-base (ml_base) and mismatch terms are omitted.
                // TODO: add dedicated pseudoloop parameters (Strategy C).
                let _ = double_pair(lp.closing)?; // validate descriptor
                let n_stems = lp.children.len() + 2;
                Ok(self.ml_closing + self.ml_intern * n_stems as i32)
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
        assert!(result.is_ok(), "pseudoknot energy evaluation failed: {:?}", result);
    }
}
