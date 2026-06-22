# Documentation overview of pseudoknots directory

The directory is managed to have linear dependency.
All modules depend on the ff_structure crate.

## Documents

- [Pseudoknots.md](Pseudoknots.md) — Full technical reference for the loop-enumeration pipeline (`parser.rs`, `closed_region_tree.rs`, `loops.rs`, `enumerate.rs`)
- [LocationStatus.md](LocationStatus.md) — Detailed explanation of the `LocationStatus` enum and arc-diagram visualization
- [PseudoEnergyModel.md](PseudoEnergyModel.md) — `PseudoEnergyModel` trait and `ViennaRNA` implementation for pseudoknot free energy evaluation

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

