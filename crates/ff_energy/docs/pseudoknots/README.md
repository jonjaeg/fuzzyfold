# Documentation overview of pseudoknots directory

The directory is managed to have linear dependency.
All modules depend on the ff_structure crate.

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

