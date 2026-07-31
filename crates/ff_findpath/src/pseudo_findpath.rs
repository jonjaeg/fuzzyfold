//! PK-aware findpath: beam search from the empty structure to a pseudoknotted target.
//!
//! Unlike the standard [`findpath`](crate::findpath::findpath), this variant:
//! - starts from the fully-unpaired structure (not an arbitrary S1)
//! - allows crossing base-pair insertions (pseudoknots)
//! - evaluates energy via [`energy_of_pseudoknotted_structure`] — a full
//!   closed-region-tree recalculation on every step
//!
//! The algorithm is a directed beam search: only moves that add pairs present
//! in the target structure are considered. Since the start is fully unpaired,
//! all moves are insertions; the search explores orderings of inserting the
//! P target pairs and returns the ordering with the lowest saddle energy.

use std::cmp::Ordering;
use std::collections::HashSet;

use ff_energy::{
    pair_table_to_dot_bracket, parse_structure, EnergyError, NucleotideVec,
    PseudoEnergyModel, ViennaRNA,
    extended_dot_bracket_to_pair_table,
};
use ff_structure::{ExtendedDotBracketVec, PairTable, StructureError, NAIDX};

use crate::utils::{Move, PathStats, PathStep, compare_structures};

// ---------------------------------------------------------------------------
// Internal beam state
// ---------------------------------------------------------------------------

struct PseudoIntermediate {
    pt:              PairTable,
    saddle_energy:   f64,
    current_energy:  f64,
    remaining_moves: Vec<Move>,
    path:            Vec<Move>,
}

// ---------------------------------------------------------------------------
// Energy helper
// ---------------------------------------------------------------------------

/// Evaluate the full pseudoknot energy of `pt` in kcal/mol.
///
/// Converts the pair table to an extended dot-bracket string, builds the
/// closed-region tree via [`parse_structure`], and sums loop energies.
fn eval_energy(model: &ViennaRNA, seq: &[ff_energy::Base], pt: &PairTable) -> Result<f64, String> {
    let dot = pair_table_to_dot_bracket(pt)
        .map_err(|e: StructureError| format!("dot-bracket conversion: {e}"))?;
    let loops = parse_structure(&dot)
        .map_err(|e: StructureError| format!("parse_structure: {e}"))?;
    model.energy_of_pseudoknotted_structure(seq, &loops)
        .map(|e| e as f64 / 100.0)
        .map_err(|e: EnergyError| format!("energy: {e:?}"))
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Find the minimum-barrier folding path from the empty structure to `target`.
///
/// # Arguments
/// - `model`       — ViennaRNA energy model (must have PK parameters attached via
///                   [`ViennaRNA::with_pseudoknot_params`] for PK-accurate energies)
/// - `sequence`    — RNA sequence string (`ACGU` alphabet)
/// - `target`      — target structure in extended dot-bracket notation (e.g. `"((([[[)))...]]]"`)
/// - `beam_width`  — number of beam members kept per step (`1` = greedy)
/// - `max_energy`  — optional energy ceiling in kcal/mol; intermediates above
///                   this threshold are pruned
///
/// # Returns
/// `Ok((path, stats))` where `path` is the full step-by-step trajectory and
/// `stats` summarises the saddle energy and barrier.
///
/// # Note on bracket families
/// Structure strings in the returned path are produced by [`pair_table_to_dot_bracket`],
/// which applies a greedy graph-colouring. The family labels (`()`, `[]`, …) may differ
/// from those in `target` even when the pair tables are identical.
pub fn findpath_pseudo(
    model: &ViennaRNA,
    sequence: &str,
    target: &str,
    beam_width: usize,
    mut max_energy: Option<f64>,
) -> Result<(Vec<PathStep>, PathStats), String> {
    // ── parse inputs ─────────────────────────────────────────────────────────
    let seq_vec = NucleotideVec::try_from_rna(sequence)
        .map_err(|_| "Failed to parse RNA sequence".to_string())?;
    let seq: &[ff_energy::Base] = &seq_vec;

    let target_edbv = ExtendedDotBracketVec::try_from(target)
        .map_err(|e: StructureError| format!("Invalid target structure: {e}"))?;
    let target_pt = extended_dot_bracket_to_pair_table(&target_edbv)
        .map_err(|e: StructureError| format!("Target pair table: {e}"))?;

    let n = target_pt.len();
    if n != seq.len() {
        return Err(format!(
            "Sequence length ({}) != structure length ({})", seq.len(), n
        ));
    }

    // ── start state ───────────────────────────────────────────────────────────
    let start_pt = PairTable::new(n);
    let start_energy = eval_energy(model, seq, &start_pt)?;

    // All moves are insertions of the target pairs (start is fully unpaired).
    let diff = compare_structures(&start_pt, &target_pt);
    let total_steps = diff.move_list.len();

    if total_steps == 0 {
        // Target is the empty structure — trivial path.
        let stats = PathStats {
            saddle_energy: start_energy,
            barrier_energy: 0.0,
            start_energy,
            end_energy: start_energy,
        };
        return Ok((
            vec![PathStep {
                structure: ".".repeat(n),
                move_applied: None,
                energy: start_energy,
                step_index: 0,
            }],
            stats,
        ));
    }

    // ── beam initialisation ───────────────────────────────────────────────────
    let mut beam = vec![PseudoIntermediate {
        pt:              start_pt,
        saddle_energy:   start_energy,
        current_energy:  start_energy,
        remaining_moves: diff.move_list,
        path:            Vec::new(),
    }];

    // ── beam search ───────────────────────────────────────────────────────────
    for _step in 0..total_steps {
        let mut candidates: Vec<PseudoIntermediate> = Vec::new();

        for parent in &beam {
            for (idx, mv) in parent.remaining_moves.iter().enumerate() {
                let i = mv.i as usize;
                let j = mv.j as usize;

                // Validity check.
                // Insertions: both positions must currently be unpaired.
                // Crossing pairs (pseudoknots) are allowed — no LoopTable check.
                if mv.is_insertion {
                    if parent.pt[i].is_some() || parent.pt[j].is_some() {
                        continue;
                    }
                } else {
                    // Deletion: the pair must exist (shouldn't happen for start=empty,
                    // but kept for correctness if used with a non-empty start later).
                    if parent.pt[i] != Some(mv.j) {
                        continue;
                    }
                }

                // Apply move.
                let mut new_pt = parent.pt.clone();
                if mv.is_insertion {
                    new_pt[i] = Some(mv.j);
                    new_pt[j] = Some(mv.i);
                } else {
                    new_pt[i] = None;
                    new_pt[j] = None;
                }

                // Full PK energy evaluation.
                let energy = match eval_energy(model, seq, &new_pt) {
                    Ok(e) => e,
                    Err(_) => continue,
                };

                // Energy ceiling filter.
                if max_energy.is_some_and(|cap| energy > cap) {
                    continue;
                }

                let new_saddle = f64::max(parent.saddle_energy, energy);

                let mut new_remaining = parent.remaining_moves.clone();
                new_remaining.remove(idx);

                let mut new_path = parent.path.clone();
                new_path.push(mv.clone());

                candidates.push(PseudoIntermediate {
                    pt:              new_pt,
                    saddle_energy:   new_saddle,
                    current_energy:  energy,
                    remaining_moves: new_remaining,
                    path:            new_path,
                });
            }
        }

        if candidates.is_empty() {
            return Err(
                "Search stuck: no valid moves remain (energy ceiling too tight?)".to_string()
            );
        }

        // Sort: lowest saddle first; break ties by current energy.
        candidates.sort_by(|a, b| {
            a.saddle_energy
                .partial_cmp(&b.saddle_energy)
                .unwrap_or(Ordering::Equal)
                .then_with(|| {
                    a.current_energy
                        .partial_cmp(&b.current_energy)
                        .unwrap_or(Ordering::Equal)
                })
        });

        // Deduplicate by pair table (first occurrence = best path to that structure).
        let mut seen: HashSet<PairTable> = HashSet::new();
        candidates.retain(|c| seen.insert(c.pt.clone()));

        // Tighten the energy ceiling to the best saddle found so far.
        if let Some(best_saddle) = candidates.first().map(|c| c.saddle_energy) {
            max_energy = Some(match max_energy {
                Some(cap) => cap.min(best_saddle),
                None      => best_saddle,
            });
        }

        // Prune to beam width.
        candidates.truncate(beam_width);
        beam = candidates;
    }

    // ── reconstruct full trajectory ───────────────────────────────────────────
    let winner = beam.into_iter().next().ok_or("No beam survivors")?;
    reconstruct_path(model, seq, n, &winner.path, start_energy)
}

// ---------------------------------------------------------------------------
// Path reconstruction
// ---------------------------------------------------------------------------

fn reconstruct_path(
    model:        &ViennaRNA,
    seq:          &[ff_energy::Base],
    n:            usize,
    moves:        &[Move],
    start_energy: f64,
) -> Result<(Vec<PathStep>, PathStats), String> {
    let mut trajectory = Vec::with_capacity(moves.len() + 1);
    let mut pt = PairTable::new(n);

    trajectory.push(PathStep {
        structure:     ".".repeat(n),
        move_applied:  None,
        energy:        start_energy,
        step_index:    0,
    });

    let mut saddle = start_energy;

    for (step_idx, mv) in moves.iter().enumerate() {
        let i = mv.i as usize;
        let j = mv.j as usize;
        if mv.is_insertion {
            pt[i] = Some(mv.j);
            pt[j] = Some(mv.i);
        } else {
            pt[i] = None;
            pt[j] = None;
        }

        let en = eval_energy(model, seq, &pt)?;
        if en > saddle {
            saddle = en;
        }

        let dot = pair_table_to_dot_bracket(&pt)
            .unwrap_or_else(|_| "?".repeat(n));

        trajectory.push(PathStep {
            structure:    dot,
            move_applied: Some(mv.clone()),
            energy:       en,
            step_index:   step_idx + 1,
        });
    }

    let stats = PathStats {
        saddle_energy:  saddle,
        barrier_energy: saddle - start_energy,
        start_energy,
        end_energy:     trajectory.last().unwrap().energy,
    };

    Ok((trajectory, stats))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use ff_energy::{ViennaRNA, parameters::{RNA_MT09, RNA_DP09}};

    fn model() -> ViennaRNA {
        ViennaRNA::from_andrunescu_params(&RNA_MT09).with_pseudoknot_params(RNA_DP09)
    }

    /// Empty target — trivial path, no steps.
    #[test]
    fn test_empty_target() {
        let m = model();
        let seq = "GCGAAACGC";
        let target = ".........";
        let (path, stats) = findpath_pseudo(&m, seq, target, 1, None).unwrap();
        assert_eq!(path.len(), 1);
        assert_eq!(stats.barrier_energy, 0.0);
    }

    /// Non-pseudoknotted target: result must equal a standard hairpin energy.
    #[test]
    fn test_non_pk_target() {
        use ff_energy::{EnergyModel, NucleotideVec};
        use ff_structure::PairTable;

        let m = model();
        let seq = "GGGAAACCC";
        let target = "(((...)))";
        let (path, _stats) = findpath_pseudo(&m, seq, target, 4, None).unwrap();

        // Final structure must be the target.
        assert_eq!(path.last().unwrap().structure, target);

        // Final energy must match energy_of_pseudoknotted_structure on the target.
        let seq_vec = NucleotideVec::try_from_rna(seq).unwrap();
        let loops = parse_structure(target).unwrap();
        let expected = m.energy_of_pseudoknotted_structure(&seq_vec, &loops).unwrap() as f64 / 100.0;
        assert!(
            (path.last().unwrap().energy - expected).abs() < 0.01,
            "energy mismatch: got {}, expected {expected}", path.last().unwrap().energy
        );
    }

    /// H-type pseudoknot target: path must reach target and all steps must be valid.
    #[test]
    fn test_h_type_pk_target() {
        let m = model();
        // GCGAUUUCUGACCGCUUUUUUGUCAG / [[[....(((((]]]......)))))
        let seq    = "GCGAUUUCUGACCGCUUUUUUGUCAG";
        let target = "[[[....(((((]]]......)))))";
        let (path, stats) = findpath_pseudo(&m, seq, target, 40, None).unwrap();

        // pair_table_to_dot_bracket may assign different bracket families than the input,
        // so compare pair tables rather than strings.
        let final_pt = extended_dot_bracket_to_pair_table(
            &ExtendedDotBracketVec::try_from(path.last().unwrap().structure.as_str()).unwrap()
        ).unwrap();
        let target_pt = extended_dot_bracket_to_pair_table(
            &ExtendedDotBracketVec::try_from(target).unwrap()
        ).unwrap();
        assert_eq!(final_pt, target_pt, "final pair table must match target");
        println!("H-type PK path ({} steps), saddle = {:.2} kcal/mol:",
            path.len() - 1, stats.saddle_energy);
        for step in &path {
            println!("  [{:2}] {} {:.2} kcal/mol",
                step.step_index, step.structure, step.energy);
        }
    }
}
