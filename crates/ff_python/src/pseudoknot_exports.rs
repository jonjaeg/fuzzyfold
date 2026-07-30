use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

use ff_structure::PairTable;
use ff_energy::BuildEvent;
use ff_energy::build_closed_regions_tree_steps;
use ff_energy::build_closed_regions_tree;
use ff_energy::NucleotideVec;
use ff_energy::PseudoEnergyModel;
use ff_energy::parse_structure;
use ff_energy::parameters::{RNA_MT09, RNA_TURNER_2004, DPParams, RNA_DP03, RNA_DP09};
use ff_energy::ViennaRNA;

/// Python-facing snapshot of one step in Algorithm 1 (Rastegari & Condon).
///
/// Fields:
/// - `lam`: the position processed this step (0-based).
/// - `event`: `"unpaired"`, `"opening"`, or `"closing"`.
/// - `merged`: (i,j) of regions absorbed during pseudoknot merging this step.
/// - `stack`: (i,j) of all open (not yet completed) regions on the construction stack.
/// - `top_level`: (i,j) of all completed top-level regions so far.
/// - `nodes`: (i,j) of every region in the arena (index-stable across steps).
/// - `parents`: parent arena index for each node (None = top-level).
/// - `children`: child arena indices for each node.
/// - `completed`: (i,j) of the region completed this step, or None.
#[pyclass(get_all)]
pub struct TreeStep {
    pub lam: usize,
    pub event: String,
    pub merged: Vec<(usize, usize)>,
    pub stack: Vec<(usize, usize)>,
    pub top_level: Vec<(usize, usize)>,
    pub nodes: Vec<(usize, usize)>,
    pub parents: Vec<Option<usize>>,
    pub children: Vec<Vec<usize>>,
    pub completed: Option<(usize, usize)>,
}

/// Parse `structure` (dot-bracket string, pseudoknots allowed) and return
/// `(pair_table, steps)` where:
/// - `pair_table[i]` is `j` if position `i` pairs with `j`, else `-1`.
/// - `steps` is one `TreeStep` per position in the structure, describing the
///   state of the closed-region-tree builder after processing that position.
#[pyfunction]
pub fn region_tree_steps(structure: &str) -> PyResult<(Vec<i64>, Vec<TreeStep>)> {
    let pt = PairTable::try_from(structure)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    let pair_table: Vec<i64> = (0..pt.len())
        .map(|i| pt[i].map(|j| j as i64).unwrap_or(-1))
        .collect();

    let steps = build_closed_regions_tree_steps(&pt)
        .into_iter()
        .map(|s| {
            let (event_str, merged) = match s.event {
                BuildEvent::Unpaired => ("unpaired".to_string(), vec![]),
                BuildEvent::Opening  => ("opening".to_string(),  vec![]),
                BuildEvent::Closing { merged } => ("closing".to_string(), merged),
            };
            TreeStep {
                lam: s.lam,
                event: event_str,
                merged,
                stack:     s.stack.iter().map(|&i| (s.nodes[i].i, s.nodes[i].j)).collect(),
                top_level: s.top_level.iter().map(|&i| (s.nodes[i].i, s.nodes[i].j)).collect(),
                nodes:     s.nodes.iter().map(|r| (r.i, r.j)).collect(),
                parents:   s.nodes.iter().map(|r| r.parent).collect(),
                children:  s.nodes.iter().map(|r| r.children.clone()).collect(),
                completed: s.completed.map(|i| (s.nodes[i].i, s.nodes[i].j)),
            }
        })
        .collect();

    Ok((pair_table, steps))
}

/// One node in the final closed-region tree, listed in postfix (bottom-up) order.
/// The root is always the last entry (`is_root == true`).
///
/// Fields:
/// - `i`, `j`: 0-based endpoints of this closed region.
/// - `children`: postfix-list indices of this node's children (left-to-right).
/// - `is_root`: true only for the virtual root (last entry).
#[pyclass(get_all)]
pub struct RegionNode {
    pub i: usize,
    pub j: usize,
    pub children: Vec<usize>,
    pub is_root: bool,
}

fn postorder_collect(nodes: &[ff_energy::ClosedRegion], arena_idx: usize, result: &mut Vec<usize>) {
    for &child in &nodes[arena_idx].children {
        postorder_collect(nodes, child, result);
    }
    result.push(arena_idx);
}

/// Parse `structure` and return `(pair_table, region_nodes)` where:
/// - `pair_table[i]` is `j` if position `i` pairs with `j`, else `-1`.
/// - `region_nodes` lists every closed region in **postfix order** (leaves first,
///   root last). The root entry (`is_root == True`) represents the implicit root
///   with `i=0, j=n`; all other entries are real closed regions.
///
/// Use this to walk the BL-construction algorithm: iterate `region_nodes[:-1]`
/// in order and for each node extract from L all positions in `[node.i, node.j]`.
#[pyfunction]
pub fn region_tree(structure: &str) -> PyResult<(Vec<i64>, Vec<RegionNode>)> {
    let pt = PairTable::try_from(structure)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    let pair_table: Vec<i64> = (0..pt.len())
        .map(|i| pt[i].map(|j| j as i64).unwrap_or(-1))
        .collect();

    let tree = build_closed_regions_tree(&pt);

    // Build postfix ordering over the arena.
    let mut postfix: Vec<usize> = Vec::with_capacity(tree.nodes.len());
    for &tl in &tree.top_level {
        postorder_collect(&tree.nodes, tl, &mut postfix);
    }

    // Map arena index → position in the postfix list.
    let mut arena_to_post = vec![0usize; tree.nodes.len()];
    for (pos, &arena_idx) in postfix.iter().enumerate() {
        arena_to_post[arena_idx] = pos;
    }

    // Build RegionNode list in postfix order.
    let mut region_nodes: Vec<RegionNode> = postfix
        .iter()
        .map(|&arena_idx| {
            let cr = &tree.nodes[arena_idx];
            RegionNode {
                i: cr.i,
                j: cr.j,
                children: cr.children.iter().map(|&c| arena_to_post[c]).collect(),
                is_root: false,
            }
        })
        .collect();

    // Append virtual root (children = top-level nodes in postfix positions).
    let root_children: Vec<usize> = tree.top_level.iter().map(|&tl| arena_to_post[tl]).collect();
    region_nodes.push(RegionNode {
        i: 0,
        j: tree.n,
        children: root_children,
        is_root: true,
    });

    Ok((pair_table, region_nodes))
}

/// Evaluate the free energy (kcal/mol) of a pseudoknotted RNA structure using
/// the Dirks-Pierce model with caller-supplied DP parameters.
///
/// Parameters
/// ----------
/// sequence : str
///     RNA sequence (uppercase A/U/G/C).
/// structure : str
///     Dot-bracket structure with pseudoknot notation (`()[]{}<>ABCDabcd`).
/// nn_params : str
///     Nearest-neighbor parameter set: `"mt09"` (recommended), `"turner2004"`,
///     `"dp03"`, or `"dp09"`.
/// init_external, init_multiloop, init_pseudoloop, pb, pup, pps, ap, bp, cp : float
///     Integer-valued DP penalties in **kcal/mol** (converted internally to dcal/mol).
/// e_stp, e_intp : float
///     Multiplicative scale factors (dimensionless).
///
/// Returns
/// -------
/// float
///     Total free energy in kcal/mol.
#[pyfunction]
#[pyo3(signature = (
    sequence,
    structure,
    nn_params = "mt09",
    init_external   = 1.38,
    init_multiloop  = 10.07,
    init_pseudoloop = 15.00,
    pb   = 2.46,
    pup  = 0.06,
    pps  = 0.96,
    ap   = 3.41,
    bp   = 0.56,
    cp   = 0.12,
    e_stp  = 0.89,
    e_intp = 0.74,
))]
pub fn pseudo_energy(
    sequence:        &str,
    structure:       &str,
    nn_params:       &str,
    init_external:   f64,
    init_multiloop:  f64,
    init_pseudoloop: f64,
    pb:    f64,
    pup:   f64,
    pps:   f64,
    ap:    f64,
    bp:    f64,
    cp:    f64,
    e_stp:  f64,
    e_intp: f64,
) -> PyResult<f64> {
    let dp = DPParams {
        init_external:   (init_external   * 100.0).round() as i32,
        init_multiloop:  (init_multiloop  * 100.0).round() as i32,
        init_pseudoloop: (init_pseudoloop * 100.0).round() as i32,
        pb:  (pb  * 100.0).round() as i32,
        pup: (pup * 100.0).round() as i32,
        pps: (pps * 100.0).round() as i32,
        ap:  (ap  * 100.0).round() as i32,
        bp:  (bp  * 100.0).round() as i32,
        cp:  (cp  * 100.0).round() as i32,
        e_stp,
        e_intp,
    };

    let model = match nn_params {
        "mt09" => ViennaRNA::from_andrunescu_params(&RNA_MT09).with_pseudoknot_params(dp),
        "turner2004" => ViennaRNA::from_thermo_params(&RNA_TURNER_2004, 37.0).with_pseudoknot_params(dp),
        "dp03" => ViennaRNA::from_thermo_params(&RNA_TURNER_2004, 37.0).with_pseudoknot_params(RNA_DP03),
        "dp09" => ViennaRNA::from_andrunescu_params(&RNA_MT09).with_pseudoknot_params(RNA_DP09),
        other => return Err(PyValueError::new_err(format!(
            "Unknown nn_params '{}'. Valid: 'mt09', 'turner2004', 'dp03', 'dp09'.", other
        ))),
    };

    let seq   = NucleotideVec::try_from_rna(sequence)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let loops = parse_structure(structure)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    let energy = model.energy_of_pseudoknotted_structure(&seq, &loops)
        .map_err(|e| PyValueError::new_err(format!("{:?}", e)))?;

    Ok(energy as f64 / 100.0)
}
