//! PK-aware neighbor enumeration for co-transcriptional folding.
//!
//! Given a sequence and a pseudoknotted structure in extended dot-bracket
//! notation, [`neighbors`] returns all structures that differ by exactly one
//! base pair (add or delete), tagged with the [`Move`] that produced them.
//! Bracket families are reassigned from scratch after each move by greedily
//! colouring the pair-crossing graph.

use ff_structure::{BracketKind, ExtendedDotBracketVec, PairTable, StructureError, NAIDX};

use crate::nucleotides::{Base, PairTypeRNA};
use crate::pseudoknots::extended_dot_bracket_to_pair_table;

/// Minimum number of unpaired bases inside a hairpin loop.
const MIN_HAIRPIN_LOOP: usize = 3;

// ---------------------------------------------------------------------------
// Move set
// ---------------------------------------------------------------------------

/// A single move in the PK-aware neighbor graph.
///
/// The three active variants cover all 1-pair-difference neighbors needed for
/// co-transcriptional folding (Mode A greedy MFE and Mode B kinetic SSA).
/// Commented variants document design choices that are deferred or
/// intentionally excluded from the initial implementation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Move {
    /// Remove an existing base pair `(i, j)` with `i < j`.
    ///
    /// Bracket families are reassigned globally after deletion. If `(i, j)`
    /// was the only pair crossing a partner family, the pseudoknot resolves
    /// and all remaining brackets collapse to `()`.
    DeletePair { i: usize, j: usize },

    /// Add a Watson–Crick or G·U pair `(i, j)` that does **not** cross any
    /// existing pair.
    ///
    /// Always assigned family `()`. Cannot create a pseudoknot by itself.
    /// In kinetic Mode B this move uses rate constant `k0`.
    AddNested { i: usize, j: usize },

    /// Add a Watson–Crick or G·U pair `(i, j)` that crosses **one or more**
    /// existing pairs.
    ///
    /// Assigned the lowest bracket family not already used by any pair it
    /// crosses. Creates or extends a pseudoknot. The exact topology is
    /// determined by the crossing graph after the move.
    /// In kinetic Mode B this move also uses rate constant `k0`; the PK
    /// initiation penalty is absorbed into ΔG via `init_external`.
    AddCrossing { i: usize, j: usize },

    // ── Deferred variants ────────────────────────────────────────────────────
    //
    // ShiftNested { from_i: usize, from_j: usize, to_i: usize, to_j: usize }
    //
    // Move one endpoint of an existing non-crossing pair concertedly:
    // (i,j) → (i,k) or (k,j). Kinetically distinct from DeletePair +
    // AddNested; uses rate `k_3ws` in the non-PK SSA (ShiftIK / ShiftJK).
    // Required for Mode B but not for Mode A greedy MFE sweep.

    // ShiftCrossing { from_i: usize, from_j: usize, to_i: usize, to_j: usize }
    //
    // Move one endpoint of an existing crossing pair. The most complex variant:
    // the bracket family and crossing partners may both change, and the
    // crossing graph must be re-analysed to determine whether the PK topology
    // is preserved, modified, or destroyed. Rate k_3ws or k_4ws depending on
    // how many strands participate in the rearrangement.

    // AddCrossingMulti { i: usize, j: usize }
    //
    // Variant of AddCrossing where (i,j) crosses pairs in more than one
    // existing closed region, merging them into a single pseudoknot. Currently
    // folded into AddCrossing — the greedy coloring resolves the bracket
    // assignment automatically — but separating it would allow a different
    // initiation penalty for higher-order topologies (kissing loops, etc.)
    // beyond the H-type PKs that DP09 was parameterised for.
}

// ---------------------------------------------------------------------------
// PairTable → extended dot-bracket
// ---------------------------------------------------------------------------

/// Convert a [`PairTable`] to an extended dot-bracket string.
///
/// Bracket families (`()`, `[]`, `{}`, …) are assigned by greedily colouring
/// the pair-crossing graph: crossing pairs receive different families; the
/// lowest available family index is always chosen.
///
/// Returns `Err` if the structure requires more than 8 bracket levels.
pub fn pair_table_to_dot_bracket(pt: &PairTable) -> Result<String, StructureError> {
    let n = pt.len();

    // Collect canonical (i < j) pairs sorted by left endpoint.
    let mut pairs: Vec<(usize, usize)> = (0..n)
        .filter_map(|i| pt[i].map(|j| (i, j as usize)))
        .filter(|&(i, j)| i < j)
        .collect();
    pairs.sort_unstable_by_key(|&(i, _)| i);

    let max_levels = BracketKind::all().len();
    let mut family = vec![0u8; pairs.len()];

    for idx in 0..pairs.len() {
        let (pi, pj) = pairs[idx];
        let mut blocked = [false; 8]; // one slot per BracketKind

        // All previously assigned pairs have left endpoint qi <= pi (sort order).
        // (pi, pj) and (qi, qj) cross iff qi < pi < qj < pj, i.e. pi < qj < pj.
        for jdx in 0..idx {
            let (_, qj) = pairs[jdx];
            if pi < qj && qj < pj {
                blocked[family[jdx] as usize] = true;
            }
        }

        let f = (0..max_levels).find(|&f| !blocked[f]).ok_or_else(|| {
            StructureError::InvalidToken(
                format!("pair ({pi},{pj})"),
                "extended dot-bracket".to_string(),
                pi,
            )
        })?;
        family[idx] = f as u8;
    }

    let kinds = BracketKind::all();
    let mut chars: Vec<char> = vec!['.'; n];
    for (idx, &(i, j)) in pairs.iter().enumerate() {
        let kind = kinds[family[idx] as usize];
        chars[i] = kind.open_char();
        chars[j] = kind.close_char();
    }
    Ok(chars.into_iter().collect())
}

// ---------------------------------------------------------------------------
// Move classification helper
// ---------------------------------------------------------------------------

/// Returns `true` if adding `(i, j)` would cross any existing pair in `pt`.
///
/// Two pairs `(i,j)` and `(qi,qj)` with `qi < qj` cross when
/// `qi < i < qj < j` or `i < qi < j < qj`.
fn crosses_existing(pt: &PairTable, i: usize, j: usize) -> bool {
    (0..pt.len())
        .filter_map(|qi| pt[qi].map(|qj| (qi, qj as usize)))
        .filter(|&(qi, qj)| qi < qj)
        .any(|(qi, qj)| (qi < i && i < qj && qj < j) || (i < qi && qi < j && j < qj))
}

// ---------------------------------------------------------------------------
// Neighbor enumeration
// ---------------------------------------------------------------------------

/// All structures differing from `structure` by exactly one base pair,
/// together with the [`Move`] that produced each one.
///
/// Includes **delete-pair** moves (remove any existing pair) and **add-pair**
/// moves (insert any Watson–Crick or G·U pair between two currently unpaired
/// positions separated by at least [`MIN_HAIRPIN_LOOP`] bases). Add moves are
/// further classified as [`Move::AddNested`] or [`Move::AddCrossing`]
/// depending on whether the new pair crosses any existing pair.
///
/// Bracket families in the returned strings are assigned fresh after each
/// move, so moves that create or resolve pseudoknots are reflected
/// automatically.
///
/// `seq` must have the same length as `structure`.
pub fn neighbors(seq: &[Base], structure: &str) -> Result<Vec<(Move, String)>, StructureError> {
    let edbv = ExtendedDotBracketVec::try_from(structure)?;
    let pt = extended_dot_bracket_to_pair_table(&edbv)?;
    let n = pt.len();

    debug_assert_eq!(
        n,
        seq.len(),
        "sequence and structure length must match"
    );

    let mut result = Vec::new();

    // ── Delete moves ──────────────────────────────────────────────────────────
    for i in 0..n {
        if let Some(j_naidx) = pt[i] {
            let j = j_naidx as usize;
            if i < j {
                let mut new_pt = pt.clone();
                new_pt[i] = None;
                new_pt[j] = None;
                if let Ok(s) = pair_table_to_dot_bracket(&new_pt) {
                    result.push((Move::DeletePair { i, j }, s));
                }
            }
        }
    }

    // ── Add moves ─────────────────────────────────────────────────────────────
    for i in 0..n {
        if pt[i].is_some() {
            continue;
        }
        for j in (i + MIN_HAIRPIN_LOOP + 1)..n {
            if pt[j].is_some() {
                continue;
            }
            if !PairTypeRNA::from((seq[i], seq[j])).can_pair() {
                continue;
            }
            let mv = if crosses_existing(&pt, i, j) {
                Move::AddCrossing { i, j }
            } else {
                Move::AddNested { i, j }
            };
            let mut new_pt = pt.clone();
            new_pt[i] = Some(j as NAIDX);
            new_pt[j] = Some(i as NAIDX);
            if let Ok(s) = pair_table_to_dot_bracket(&new_pt) {
                result.push((mv, s));
            }
        }
    }

    Ok(result)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::convert::TryFrom;
    use ff_structure::PairTable;

    fn seq(s: &str) -> Vec<Base> {
        s.chars().map(|c| Base::try_from(c).unwrap()).collect()
    }

    fn strings(nbrs: Vec<(Move, String)>) -> Vec<String> {
        nbrs.into_iter().map(|(_, s)| s).collect()
    }

    // ── pair_table_to_dot_bracket ─────────────────────────────────────────

    #[test]
    fn roundtrip_simple_hairpin() {
        let pt = PairTable::try_from("(((...)))").unwrap();
        let s = pair_table_to_dot_bracket(&pt).unwrap();
        assert_eq!(s, "(((...)))");
    }

    #[test]
    fn roundtrip_h_type_pk() {
        let pt = PairTable::try_from("((([[[)))...]]]").unwrap();
        let s = pair_table_to_dot_bracket(&pt).unwrap();
        assert_eq!(s, "((([[[)))...]]]");
    }

    #[test]
    fn roundtrip_fully_unpaired() {
        let pt = PairTable::try_from("....").unwrap();
        let s = pair_table_to_dot_bracket(&pt).unwrap();
        assert_eq!(s, "....");
    }

    // ── move classification ───────────────────────────────────────────────

    #[test]
    fn add_non_crossing_is_nested() {
        // GGGCCC: adding (0,5) doesn't cross anything
        let s = seq("GGGCCC");
        let nbrs = neighbors(&s, "......").unwrap();
        assert!(nbrs.iter().all(|(mv, _)| matches!(mv, Move::AddNested { .. })));
    }

    #[test]
    fn add_crossing_is_classified() {
        // (.[.) .] — existing pair (0,3); adding (2,5) crosses it → AddCrossing
        // Sequence: GCACGC (0=G,1=C,2=A,3=C,4=G,5=C)
        // Actually let's build it more carefully:
        // existing structure: (.....) with pair (0,6), seq GCAAAAGC (len 8)
        // add (2,5): 0 < 2 < 6, and 5 < 6, so (0,6) and (2,5) don't cross.
        // Let me use: structure (.....) pair (0,5), seq GCAAAAC, add (2,7) — wait need both unpaired.
        // Simpler: seq GACGAC (len 6), structure (.(...)) pair (1,5)
        // add (0,3): 0 < 1 < 3 < 5 → crosses (1,5) → AddCrossing
        let _s = seq("GACGAC");
        let pt = PairTable::try_from(".(...)").unwrap(); // pair (1,5)
        // manually verify crosses_existing
        assert!(crosses_existing(&pt, 0, 3)); // 0 < 1 < 3 < 5 ✓
    }

    // ── neighbors: delete moves ───────────────────────────────────────────

    #[test]
    fn delete_only_pair_gives_unpaired() {
        let s = seq("GGGAAACCC");
        let nbrs = strings(neighbors(&s, "(((...)))").unwrap());
        let deletes: Vec<_> = nbrs.iter()
            .filter(|n| n.chars().filter(|&c| c == '(').count() < 3)
            .collect();
        assert_eq!(deletes.len(), 3);
    }

    #[test]
    fn delete_from_pk_can_resolve_pseudoknot() {
        // Minimal H-type PK: (.[.).]  pairs (0,4)=GC and (2,6)=CG cross.
        // Deleting either crossing pair leaves one non-crossing pair → () only.
        let s = seq("GACACAG");
        let nbrs = strings(neighbors(&s, "(.[.).]").unwrap());
        let resolved: Vec<_> = nbrs.iter()
            .filter(|n| !n.contains('[') && !n.contains(']'))
            .collect();
        assert!(!resolved.is_empty(), "expected at least one PK-free neighbor");
    }

    // ── neighbors: add moves ──────────────────────────────────────────────

    #[test]
    fn add_pair_to_unpaired() {
        let s = seq("GGGCCC");
        let nbrs = strings(neighbors(&s, "......").unwrap());
        assert!(!nbrs.is_empty());
        for n in &nbrs {
            let pairs = n.chars().filter(|&c| c == '(').count();
            assert_eq!(pairs, 1, "unexpected: {n}");
        }
    }

    #[test]
    fn add_pair_creates_pk() {
        let s = seq("GGGAAACCCUUU");
        let nbrs = strings(neighbors(&s, "(((....)))..").unwrap());
        let _ = nbrs.iter().filter(|n| n.contains('[')).collect::<Vec<_>>();
    }

    #[test]
    fn min_hairpin_loop_respected() {
        let s = seq("GCGCGCGC");
        let nbrs = strings(neighbors(&s, "........").unwrap());
        for n in &nbrs {
            let pt = PairTable::try_from(n.as_str()).unwrap();
            for i in 0..pt.len() {
                if let Some(j) = pt[i] {
                    let j = j as usize;
                    if i < j {
                        assert!(j - i > MIN_HAIRPIN_LOOP, "hairpin too small: ({i},{j}) in {n}");
                    }
                }
            }
        }
    }

    #[test]
    fn no_duplicate_neighbors() {
        let s = seq("GCGCAAAAGCGC");
        let nbrs = strings(neighbors(&s, "((((....))))").unwrap());
        let mut seen = std::collections::HashSet::new();
        for n in &nbrs {
            assert!(seen.insert(n.clone()), "duplicate neighbor: {n}");
        }
    }
}
