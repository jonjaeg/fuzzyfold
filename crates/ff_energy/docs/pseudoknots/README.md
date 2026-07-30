# Documentation overview of pseudoknots directory

The directory is managed to have linear dependency.
All modules depend on the ff_structure crate.

## Documents

- [01_HISTORY.md](01_HISTORY.md) — Chronological log of implementation phases and pending work
- [02_Pseudoknots.md](02_Pseudoknots.md) — Full technical reference for the loop-enumeration pipeline (`parser.rs`, `closed_region_tree.rs`, `loops.rs`, `enumerate.rs`)
- [03_LocationStatus.md](03_LocationStatus.md) — Detailed explanation of the `LocationStatus` enum and arc-diagram visualization
- [04_PseudoloopParams.md](04_PseudoloopParams.md) — `DPParams` struct, `RNA_DP03` / `RNA_DP09` constants, energy formulas, and notes on Cao & Chen 2006
- [05_PseudoEnergyModel.md](05_PseudoEnergyModel.md) — `PseudoEnergyModel` trait and `ViennaRNA` implementation for pseudoknot free energy evaluation
- [06_DPModelImplementation.md](06_DPModelImplementation.md) — Dirks-Pierce model roadmap: all 11 features, dp03/dp09 parameter values, validation against HotKnots v2, implementation status
- [07_SparkEnergyModel.md](07_SparkEnergyModel.md) — HotKnots v2 / Spark DP matrix structure and parameter mapping to fuzzyfold
- [08_how_hotknots_parses_energy_params.md](08_how_hotknots_parses_energy_params.md) — How HotKnots locates and reads its parameter files

## Recommended reading order

**User / API consumer** — wants to evaluate pseudoknot energies:

1. [01_HISTORY.md](01_HISTORY.md) — what is implemented and what is still pending
2. [02_Pseudoknots.md](02_Pseudoknots.md) — what `parse_structure` produces
3. [04_PseudoloopParams.md](04_PseudoloopParams.md) — which parameter set to use and what the formulas are
4. [05_PseudoEnergyModel.md](05_PseudoEnergyModel.md) — how energy is computed per loop type

**Developer / extending the implementation:**

1. [01_HISTORY.md](01_HISTORY.md) — know the phases before touching anything
2. [02_Pseudoknots.md](02_Pseudoknots.md) — foundational pipeline (`RegionTree`, `Loop`, `LoopType`)
3. [03_LocationStatus.md](03_LocationStatus.md) — `Band` / `SpanBand` / `InsideBand` distinction (prerequisite for SpanBand scaling)
4. [04_PseudoloopParams.md](04_PseudoloopParams.md) — `DPParams` fields, formulas, dp09/mt09 warning, Cao-Chen limitation
5. [05_PseudoEnergyModel.md](05_PseudoEnergyModel.md) — full loop-type dispatch table
6. [06_DPModelImplementation.md](06_DPModelImplementation.md) — step-by-step implementation record with validation numbers
7. [07_SparkEnergyModel.md](07_SparkEnergyModel.md) — Spark DP matrix reference; read this if adding `RNA_MT09` or implementing Cao-Chen
8. [08_how_hotknots_parses_energy_params.md](08_how_hotknots_parses_energy_params.md) — how HotKnots locates and reads its parameter files

---

The `extended_dot_bracket_to_pair_table` in the `parser.rs` module is a strict
parser for single RNA molecules — strand break characters (`+`/`&`) are rejected,
in contrast to the permissive `PairTable::try_from(&str)` in `ff_structure`.
```
&str → DotBracketVec    → PairTable   (simple)
&str → ExtendedDotBracketVec → PairTable   (pseudoknots, strict)
```

````mermaid
flowchart TD
    Z(ff_structure) --> A[parser.rs]
    Z(ff_structure) --> B[closed_region_tree.rs]
    Z(ff_structure) --> C[loops.rs]
    Z(ff_structure) --> D[enumerate.rs]
    A --> D
    B --> C
    B --> D
    C --> D

````

