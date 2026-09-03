//! Band finding algorithm by Condon et al.
//! main enumeration file

use std::collections::HashSet;
use ff_structure::{PairTable, StructureError, ExtendedDotBracketVec};


use crate::pseudoknots::{
    ClosedRegion,
    RegionTree,
    LocationStatus,
    ClosingDescriptor,
    LoopType,
    Loop,
    PseudoloopContext,
    extended_dot_bracket_to_pair_table,
    is_pseudo,
    collect_bands,
    closing_descriptor,
    closing_pairs,
    location_status,
    interior_loop_type,
    build_closed_regions_tree,
    nested_pairs,
};


/// Enumerates the SPAN_BAND loops for `region`: the interior loops (or
/// multiloops, if other children sit in the gaps) between consecutive
/// rungs of each band found by `collect_bands`.
pub fn enumerate_band_spanning_loops(tree: &RegionTree, region: &ClosedRegion, pt: &PairTable, loops: &mut Vec<Loop>) {
    let bands = collect_bands(tree, region, pt);

    for chain in &bands {
        for idx in 0..chain.len().saturating_sub(1) {
            let outer_left = chain[idx];
            let inner_left = chain[idx + 1];
            let outer_right = pt[outer_left].expect("band rung must be paired") as usize;
            let inner_right = pt[inner_left].expect("band rung must be paired") as usize;

            let closing_bp = (outer_left, outer_right);
            let inner_bp = (inner_left, inner_right);

            let mut nest = nested_pairs(tree, &region.children, pt, outer_left, inner_left);
            nest.extend(nested_pairs(tree, &region.children, pt, inner_right, outer_right));

            if nest.is_empty() {
                let (ltype, n5, n3) = interior_loop_type(closing_bp, inner_bp);
                loops.push(
                    Loop::new(ltype, LocationStatus::SpanBand)
                        .with_closing(ClosingDescriptor::Single(closing_bp))
                        .with_inner(inner_bp)
                        .with_unpaired(n5, n3),
                );
            } else {
                let n_unpaired_5p = count_unpaired(pt, outer_left + 1, inner_left);
                let n_unpaired_3p = count_unpaired(pt, inner_right + 1, outer_right);
                let mut children = vec![ClosingDescriptor::Single(inner_bp)];
                children.extend(nest.into_iter().map(ClosingDescriptor::Single));

                loops.push(
                    Loop::new(LoopType::Multiloop, LocationStatus::SpanBand)
                        .with_closing(ClosingDescriptor::Single(closing_bp))
                        .with_children(children)
                        .with_unpaired(n_unpaired_5p, n_unpaired_3p),
                );
            }
        }
    }
}



// Main loop enumeration function

/// Builds the list of `Loop`s for a structure, given its closed-regions tree.
/// Mirrors the Python `enumerate_loops(root, pt)`, with the `region.is_root`
/// branch replaced by a final pass over `tree.root_children`.
pub fn enumerate_loops(tree: &RegionTree, pt: &PairTable) -> Vec<Loop> {
    let mut loops = Vec::new();

    for &idx in &tree.root_children {
        visit(tree, idx, pt, &mut loops, PseudoloopContext::External);
    }

    // External loop: in Python this was the `region.is_root` branch of
    // `visit`, using `root.children`. Here `root_children` plays that role.
    let child_bps: Vec<ClosingDescriptor> = tree.root_children.iter()
        .map(|&idx| closing_descriptor(&tree.nodes[idx], pt))
        .collect();

    loops.push(
        Loop::new(LoopType::External, LocationStatus::Standard)
            .with_children(child_bps),
    );

    loops
}

/// Post-order visit of `idx` and its descendants, appending one `Loop` per
/// region (plus any SPAN_BAND loops from `enumerate_band_spanning_loops`).
///
/// `pk_context` is the context this node sees — i.e. what type of loop
/// encloses this node. Passed top-down so `Pseudoloop` loops can record the
/// correct initiation-penalty context.
fn visit(tree: &RegionTree, idx: usize, pt: &PairTable, loops: &mut Vec<Loop>, pk_context: PseudoloopContext) {
    let region = &tree.nodes[idx];
    let n_children = region.children.len();
    let pseudo = is_pseudo(region, pt);

    // Determine the context that children of THIS node will inherit.
    let is_multiloop = !pseudo && (
        n_children > 1 ||
        (n_children == 1 && is_pseudo(&tree.nodes[region.children[0]], pt))
    );
    let child_context = if pseudo {
        PseudoloopContext::Pseudoloop
    } else if is_multiloop {
        PseudoloopContext::Multiloop
    } else {
        pk_context
    };

    for &child_idx in &region.children {
        visit(tree, child_idx, pt, loops, child_context);
    }

    let closing = closing_pairs(region, pt);
    let loc = location_status(region, region.parent.map(|p| &tree.nodes[p]), pt);

    if pseudo {
        let closing_set: HashSet<(usize, usize)> = closing.iter().copied().collect();
        let mut child_pairs: Vec<ClosingDescriptor> = Vec::new();

        let bands: Vec<Vec<usize>> = collect_bands(tree, region, pt);

        for chain in &bands {
            let outer_left_arm = chain[0];
            let inner_left_arm = *chain.last().unwrap();
            let outer_bp = (outer_left_arm, pt[outer_left_arm].expect("band rung must be paired") as usize);
            let inner_bp = (inner_left_arm, pt[inner_left_arm].expect("band rung must be paired") as usize);

            if !closing_set.contains(&outer_bp) {
                child_pairs.push(ClosingDescriptor::Single(outer_bp));
            }
            child_pairs.push(ClosingDescriptor::Single(inner_bp));
        }

        for &child_idx in &region.children {
            child_pairs.push(closing_descriptor(&tree.nodes[child_idx], pt));
        }

        child_pairs.sort_by_key(|cd| match cd {
            ClosingDescriptor::Single((i, _)) => *i,
            ClosingDescriptor::Double((i, _), _) => *i,
        });

        // Compute free linker bases in the pseudoloop gaps (the `pup` term).
        //
        // Mirrors HotKnots' pseudoEnergyDP directly: sort all arm segments by
        // start position, then sum the 2n−1 gaps between consecutive segments,
        // subtracting the span of any nested child region in each gap.
        // No special case for n=2 — the formula is the same for all densities.
        // DP09 was trained on 2-band structures; applying it to n>2 is an
        // extrapolation, but the counting method is correct in all cases.
        if bands.len() > 2 {
            eprintln!(
                "Warning: pseudoloop with {} bands. \
                 DP09 parameters are extrapolated beyond their 2-band training domain.",
                bands.len()
            );
        }
        let mut arm_segs: Vec<(usize, usize)> = bands.iter().flat_map(|chain| {
            let outer    = chain[0];
            let inner    = *chain.last().unwrap();
            let inner_3p = pt[inner].expect("band rung must be paired") as usize;
            let outer_3p = pt[outer].expect("band rung must be paired") as usize;
            [(outer, inner), (inner_3p, outer_3p)]
        }).collect();
        arm_segs.sort_by_key(|&(s, _)| s);
        let n_pup: usize = arm_segs.windows(2)
            .map(|w| count_pup_linkers(tree, &region.children, pt, w[0].1 + 1, w[1].0))
            .sum();

        loops.push(
            Loop::new(LoopType::Pseudoloop, loc)
                .with_closing(ClosingDescriptor::Double(closing[0], closing[1]))
                .with_children(child_pairs)
                .with_pup(n_pup)
                .with_bands(bands.len())
                .with_nested(n_children)
                .with_pk_context(pk_context),
        );
    } else if n_children == 0 {
        let (ci, cj) = closing[0];
        loops.push(
            Loop::new(LoopType::Hairpin, loc)
                .with_closing(ClosingDescriptor::Single(closing[0]))
                .with_unpaired(cj - ci - 1, 0),
        );
    } else if n_children == 1 && !is_pseudo(&tree.nodes[region.children[0]], pt) {
        let inner_pair = closing_pairs(&tree.nodes[region.children[0]], pt)[0];
        let (ltype, n5, n3) = interior_loop_type(closing[0], inner_pair);
        loops.push(
            Loop::new(ltype, loc)
                .with_closing(ClosingDescriptor::Single(closing[0]))
                .with_inner(inner_pair)
                .with_unpaired(n5, n3),
        );
    } else {
        let child_bps: Vec<ClosingDescriptor> = region.children.iter()
            .map(|&c| closing_descriptor(&tree.nodes[c], pt))
            .collect();
        loops.push(
            Loop::new(LoopType::Multiloop, loc)
                .with_closing(ClosingDescriptor::Single(closing[0]))
                .with_children(child_bps),
        );
    }

    enumerate_band_spanning_loops(tree, region, pt, loops);
}

/// Count unpaired positions in the half-open range `[start, end)`.
fn count_unpaired(pt: &PairTable, start: usize, end: usize) -> usize {
    if start >= end {
        return 0;
    }
    (start..end).filter(|&p| pt[p].is_none()).count()
}

/// Count free linker positions in the pseudoloop gap `[start, end)`.
///
/// Mirrors HotKnots' `pseudoEnergyDP`: takes the raw positional span of the
/// gap and subtracts the full span of every nested child region that overlaps
/// with it.  This matches the physical intent of the DP09 `pup` parameter —
/// only genuinely single-stranded linker nucleotides between structured
/// elements are counted, not bases that belong to a child hairpin or stem.
fn count_pup_linkers(
    tree: &RegionTree,
    children: &[usize],
    pt: &PairTable,
    start: usize,
    end: usize,
) -> usize {
    if start >= end {
        return 0;
    }
    let raw = end - start;
    let child_overlap: usize = children
        .iter()
        .filter_map(|&ci| {
            let pairs = closing_pairs(&tree.nodes[ci], pt);
            if pairs.is_empty() {
                return None;
            }
            let cb = pairs.iter().map(|&(i, _)| i).min().unwrap();
            let ce = pairs.iter().map(|&(_, j)| j).max().unwrap();
            // Intersection of child span [cb, ce] with gap [start, end-1].
            let lo = cb.max(start);
            let hi = ce.min(end - 1);
            if lo <= hi { Some(hi - lo + 1) } else { None }
        })
        .sum();
    raw - child_overlap
}

/// Builds the list of `Loop`s for a structure directly from its `PairTable`,
/// skipping the string round-trip.
///
/// Prefer this over [`parse_structure`] whenever the pair table is already
/// available (e.g. during beam search): it avoids the O(P²)
/// `pair_table_to_dot_bracket` step and the redundant re-parse back to a
/// `PairTable` that `parse_structure` performs internally.
pub fn parse_loops_from_pt(pt: &PairTable) -> Vec<Loop> {
    let tree = build_closed_regions_tree(pt);
    enumerate_loops(&tree, pt)
}

/// Main entry point for the whole enumeration process.
/// Takes a dot-bracket string, builds the pair table and closed-regions tree,
/// and returns the list of `Loop`s.
pub fn parse_structure(s: &str) -> Result<Vec<Loop>, StructureError> {
    let edbv = ExtendedDotBracketVec::try_from(s)?;
    let pt = extended_dot_bracket_to_pair_table(&edbv)?;
    Ok(parse_loops_from_pt(&pt))
}








// TESTS

#[cfg(test)]
mod band_spanning_loop_tests {
    use super::*;
    use std::convert::TryFrom;

    #[test]
    fn test_h_type_pseudoknot_band_spanning_loops() {
        // ((([[[)))...]]] -> bands [[0,1,2],[3,4,5]], no children
        // -> 4 stacked-pair (Stack) loops, all SPAN_BAND
        let pt = PairTable::try_from("((([[[)))...]]]").unwrap();
        let tree = build_closed_regions_tree(&pt);
        let region = &tree.nodes[tree.root_children[0]];

        let mut loops = Vec::new();
        enumerate_band_spanning_loops(&tree, region, &pt, &mut loops);

        assert_eq!(loops.len(), 4);
        for lp in &loops {
            assert_eq!(lp.loop_type, LoopType::Stack);
            assert_eq!(lp.location, LocationStatus::SpanBand);
            assert_eq!((lp.unpaired_5p, lp.unpaired_3p), (0, 0));
        }

        assert_eq!(loops[0].closing, Some(ClosingDescriptor::Single((0, 8))));
        assert_eq!(loops[0].inner, Some((1, 7)));

        assert_eq!(loops[1].closing, Some(ClosingDescriptor::Single((1, 7))));
        assert_eq!(loops[1].inner, Some((2, 6)));

        assert_eq!(loops[2].closing, Some(ClosingDescriptor::Single((3, 14))));
        assert_eq!(loops[2].inner, Some((4, 13)));

        assert_eq!(loops[3].closing, Some(ClosingDescriptor::Single((4, 13))));
        assert_eq!(loops[3].inner, Some((5, 12)));
    }

    #[test]
    fn test_loop_display() {
        let lp = Loop::new(LoopType::Stack, LocationStatus::SpanBand)
            .with_closing(ClosingDescriptor::Single((0, 8)))
            .with_inner((1, 7));

        assert_eq!(lp.to_string(), "Loop(Stack, SpanBand, closing=(0, 8), inner=(1, 7))");
    }
}

#[cfg(test)]
mod enumerate_loops_tests {
    use super::*;

    #[test]
    fn test_simple_hairpin() {
        // (((...)))
        let loops = parse_structure("(((...)))").unwrap();
        assert_eq!(loops.len(), 4);

        assert_eq!(
            loops[0],
            Loop::new(LoopType::Hairpin, LocationStatus::Standard)
                .with_closing(ClosingDescriptor::Single((2, 6)))
                .with_unpaired(3, 0)
        );

        assert_eq!(
            loops[1],
            Loop::new(LoopType::Stack, LocationStatus::Standard)
                .with_closing(ClosingDescriptor::Single((1, 7)))
                .with_inner((2, 6))
        );

        assert_eq!(
            loops[2],
            Loop::new(LoopType::Stack, LocationStatus::Standard)
                .with_closing(ClosingDescriptor::Single((0, 8)))
                .with_inner((1, 7))
        );

        assert_eq!(
            loops[3],
            Loop::new(LoopType::External, LocationStatus::Standard)
                .with_children(vec![ClosingDescriptor::Single((0, 8))])
        );
    }

    #[test]
    fn test_h_type_pseudoknot() {
        // ((([[[)))...]]]
        let loops = parse_structure("((([[[)))...]]]").unwrap();
        assert_eq!(loops.len(), 6);

        assert_eq!(
            loops[0],
            Loop::new(LoopType::Pseudoloop, LocationStatus::Standard)
                .with_closing(ClosingDescriptor::Double((0, 8), (3, 14)))
                .with_children(vec![
                    ClosingDescriptor::Single((2, 6)),
                    ClosingDescriptor::Single((5, 12)),
                ])
                .with_pup(3)
                .with_bands(2)
                .with_nested(0)
                .with_pk_context(PseudoloopContext::External)
        );

        for lp in &loops[1..5] {
            assert_eq!(lp.loop_type, LoopType::Stack);
            assert_eq!(lp.location, LocationStatus::SpanBand);
        }

        assert_eq!(
            loops[5],
            Loop::new(LoopType::External, LocationStatus::Standard)
                .with_children(vec![ClosingDescriptor::Double((0, 8), (3, 14))])
        );
    }
}
