//! Energy parameters for pseudoknotted structures.
//!
//! All energy values are in **dcal/mol** (1/100 kcal/mol), consistent with
//! the rest of `ff_energy`.
//!
//! The full Dirks-Pierce model has 11 features (Andronescu et al. 2010, Table 6).
//! See `docs/pseudoknots/DPModelImplementation.md` for the complete description.
//!
//! The Pseudoloop energy formula:
//!
//! ```text
//! E = init(ctx) + pb × n_bands + pup × n_pup + pps × n_nested
//! ```
//!
//! SpanBand multiloop formula:
//!
//! ```text
//! E = ap + bp × n_branches + cp × n_unpaired
//! ```
//!
//! Stacks / interior loops *inside* a band use scaled Turner energies:
//! - stacked pair: `round(e_stp × turner_stack)`
//! - interior loop / bulge: `round(e_intp × turner_interior)`

/// Full Dirks-Pierce (DP) energy parameter set for pseudoknots.
///
/// Covers all 11 features described in Andronescu, Pop & Condon (2010), Table 6.
/// Two static constants are provided: [`RNA_DP03`] (initial 2003 values) and
/// [`RNA_DP09`] (Andronescu 2010 trained values, identical to HotKnots v2 / Spark).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DPParams {
    /// Initiation when the pseudoloop is in an exterior context (dcal/mol).
    pub init_external: i32,
    /// Initiation when the pseudoloop is nested in a standard multiloop (dcal/mol).
    pub init_multiloop: i32,
    /// Initiation when the pseudoloop is nested in another pseudoloop (dcal/mol).
    pub init_pseudoloop: i32,
    /// Per-band (per crossing helix) penalty (dcal/mol).
    pub pb: i32,
    /// Per-unpaired-base in the junction gaps (dcal/mol).
    pub pup: i32,
    /// Per-nested-closed-region inside the pseudoloop (dcal/mol).
    pub pps: i32,
    /// SpanBand multiloop: initiation penalty (dcal/mol).
    pub ap: i32,
    /// SpanBand multiloop: per-helix penalty (dcal/mol).
    pub bp: i32,
    /// SpanBand multiloop: per-unpaired-base penalty (dcal/mol).
    pub cp: i32,
    /// Scale factor applied to Turner stacking energies inside a band.
    pub e_stp: f64,
    /// Scale factor applied to Turner interior-loop / bulge energies inside a band.
    pub e_intp: f64,
}

/// Initial Dirks-Pierce 2003 pseudoknot parameters for RNA at 37 °C.
///
/// Source: Andronescu, Pop & Condon, *RNA* **16** (2010), Table 6, column dp03.
/// These are the parameter values used in the original HotKnots / NUPACK software.
pub static RNA_DP03: DPParams = DPParams {
    init_external: 960,    //  9.60 kcal/mol
    init_multiloop: 1500,  // 15.00 kcal/mol
    init_pseudoloop: 1500, // 15.00 kcal/mol
    pb: 20,                //  0.20 kcal/mol per band
    pup: 10,               //  0.10 kcal/mol per unpaired base
    pps: 10,               //  0.10 kcal/mol per nested closed region
    ap: 340,               //  3.40 kcal/mol SpanBand multiloop initiation
    bp: 40,                //  0.40 kcal/mol SpanBand per helix
    cp: 0,                 //  0.00 kcal/mol SpanBand per unpaired base
    e_stp: 0.83,
    e_intp: 0.83,
};

/// Trained Dirks-Pierce 2009 pseudoknot parameters for RNA at 37 °C.
///
/// Source: Andronescu, Pop & Condon, *RNA* **16** (2010), Table 6, column dp09.
/// These are the best-performing DP parameters (F-measure 79% vs. 68% for dp03).
/// They are identical to the HotKnots v2 and Spark `PK_globals.hh` values.
///
/// Key differences from dp03:
/// - External initiation flips from +960 → −138 (exterior pseudoknots become rewarded)
/// - Per-band penalty rises from 20 → 246 (discourages over-complex pseudoknots)
pub static RNA_DP09: DPParams = DPParams {
    init_external: -138,  // -1.38 kcal/mol (slight reward — exterior PK is favorable)
    init_multiloop: 1007, // 10.07 kcal/mol
    init_pseudoloop: 1500, // 15.00 kcal/mol
    pb: 246,              //  2.46 kcal/mol per band
    pup: 6,               //  0.06 kcal/mol per unpaired base
    pps: 96,              //  0.96 kcal/mol per nested closed region
    ap: 341,              //  3.41 kcal/mol SpanBand multiloop initiation
    bp: 56,               //  0.56 kcal/mol SpanBand per helix
    cp: 12,               //  0.12 kcal/mol SpanBand per unpaired base
    e_stp: 0.89,
    e_intp: 0.74,
};

// Cao & Chen 2006 pseudoknot parameters for RNA at 37 °C.
//
// The actual Cao & Chen 2006 model (`pkmodelCC2006.dat` in HotKnots v2) uses a
// 2D lookup table indexed by *both* loop length and stem size, not a linear formula.
// `DPParams` cannot represent it. A proper implementation requires a dedicated 2D
// table struct. Deferred.
//pub struct CaoChen2006Note;

// ---------------------------------------------------------------------------
// Why CC09 is not implemented (and DP09 is the right choice)
// ---------------------------------------------------------------------------
//
// Andronescu, Pop & Condon (2010) also trained a CC09 parameter set for the
// Cao-Chen (CC) energy model. CC09 = MT09 + dp09 + cc09, with 920 total
// features (vs. 374 for DP09). The 546 additional cc09 parameters are:
//   - 258 loop entropy features for the L1/L2 junction loops of H-type PKs
//   - 288 coaxial stacking features for H-type PKs (Tyagi & Mathews 2007)
// The CC features apply only to H-type PKs where Cao & Chen provide values;
// non-H-type PKs and unparameterised H-type cases fall back to DP.
//
// Two fundamental reasons CC09 is not implemented here:
//
// 1. Non-linear energy function.
//    The CC model takes min(coaxial_stack, dangle) at each applicable site,
//    producing quadratic terms in the parameters. This is incompatible with
//    the linear loop-sum architecture used throughout ff_energy.
//
// 2. Loop entropies are ΔS°, not ΔG°.
//    Cao & Chen treat enthalpy as zero for junction loops, so their 258
//    loop entropy parameters are entropy changes that scale as −TΔS°. They
//    require separate temperature handling, unlike all other parameters.
//
// Crucially, CC09 does not outperform DP09. On short pseudoknotted structures
// (ShPK, the relevant regime for e.g. Zika xrRNA1 at ~70 nt), F-measure is
// 0.817 for DP09 vs. 0.742 for CC09 (Table 3/4, Andronescu 2010). The paper
// attributes this to insufficient training data for the large CC feature set.
// DP09 remains the recommended parameter set for this codebase.
