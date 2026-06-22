# PseudoEnergyModel — Technical Reference

> **Crate:** `ff_energy`  
> **Files:** `pseudoknots/pseudo_energy_model.rs`, `nn_models/viennarna_pseudo.rs`  
> **Reference:** Turner 2004 nearest-neighbor model; Rastegari & Condon loop classification

---

## Table of Contents

1. [Overview](#overview)
2. [Architecture](#architecture)
3. [Module: `pseudo_energy_model.rs`](#module-pseudo_energy_modelrs)
4. [Module: `viennarna_pseudo.rs`](#module-viennarna_pseudors)
5. [Loop Type Dispatch](#loop-type-dispatch)
6. [The Pseudoloop Approximation](#the-pseudoloop-approximation)
7. [Known Limitations](#known-limitations)
8. [End-to-end Example](#end-to-end-example)

---

## Overview

`PseudoEnergyModel` is a Rust trait that extends the existing nearest-neighbor
energy framework to **pseudoknotted RNA structures**. Given the `Vec<Loop>`
produced by `parse_structure`, it computes a scalar free energy in
**dcal/mol** (1/100 kcal/mol).

The design deliberately keeps pseudoknot energy evaluation separate from the
standard `EnergyModel` trait. `PseudoEnergyModel` is its own independent trait;
`ViennaRNA` happens to implement both, but neither requires the other.

The full pipeline from string to energy is:

```
&str
  │
  ▼  parse_structure  (pseudoknots/enumerate.rs)
Vec<Loop>
  │
  ▼  energy_of_pseudoknotted_structure  (nn_models/viennarna_pseudo.rs)
i32  (dcal/mol)
```

---

## Architecture

### Relationship to `EnergyModel`

`EnergyModel` (in `energy_model.rs`) evaluates non-pseudoknotted structures
via `LoopDecomposition::for_each_loop`, which emits `NearestNeighborLoop`
variants. This pathway cannot handle pseudoknots because `NearestNeighborLoop`
has no crossing-pair type, and `PairTable::try_from` rejects pseudoknotted
dot-bracket strings.

`PseudoEnergyModel` works from the orthogonal `Vec<Loop>` representation
produced by the pseudoknot pipeline. It does **not** extend `EnergyModel` —
the two traits are parallel, not hierarchical:

```
EnergyModel  ──── impl ──── ViennaRNA  ──── impl ──── PseudoEnergyModel
     │                                                       │
  NearestNeighborLoop                                   Vec<Loop>
  (standard nested structures)                   (pseudoknotted structures)
```

### File layout

| File | Contents |
|---|---|
| `pseudoknots/pseudo_energy_model.rs` | `PseudoEnergyModel` trait; helper functions `single_pair`, `double_pair`, `collect_single_branches` |
| `nn_models/viennarna_pseudo.rs` | `impl PseudoEnergyModel for ViennaRNA` |

The impl is separated from `viennarna.rs` to keep ViennaRNA's core
nearest-neighbor logic self-contained. It accesses three `pub(crate)` fields
from `ViennaRNA`: `ml_closing`, `ml_intern`, `ml_base`.

---

## Module: `pseudo_energy_model.rs`

### `PseudoEnergyModel` — trait

```rust
pub trait PseudoEnergyModel {
    fn energy_of_pseudo_loop(
        &self,
        sequence: &[Base],
        lp: &Loop,
    ) -> Result<i32, EnergyError>;

    fn energy_of_pseudoknotted_structure(
        &self,
        sequence: &[Base],
        loops: &[Loop],
    ) -> Result<i32, EnergyError> {
        // default: sum energy_of_pseudo_loop over all loops
    }
}
```

`energy_of_pseudo_loop` is the required method. The default
`energy_of_pseudoknotted_structure` folds over the slice, short-circuiting on
the first error.

### Helper functions

These are public utilities for implementing `energy_of_pseudo_loop`. They
extract pairs from `ClosingDescriptor` values and map errors to `EnergyError`.

#### `single_pair`

```rust
pub fn single_pair(cd: Option<ClosingDescriptor>) -> Result<(usize, usize), EnergyError>
```

Extracts the pair from a `ClosingDescriptor::Single`. Returns
`Err(InvalidClosingPair)` for `None` or `Double`.

#### `double_pair`

```rust
pub fn double_pair(
    cd: Option<ClosingDescriptor>,
) -> Result<((usize, usize), (usize, usize)), EnergyError>
```

Extracts both pairs from a `ClosingDescriptor::Double`. Returns
`Err(InvalidClosingPair)` for `None` or `Single`.

#### `collect_single_branches`

```rust
pub fn collect_single_branches(
    children: &[ClosingDescriptor],
) -> Result<Vec<(usize, usize)>, EnergyError>
```

Maps a slice of `ClosingDescriptor`s to their pairs, failing if any element is
a `Double`. Used for `Multiloop` children, which are guaranteed `Single` by the
enumeration algorithm.

---

## Module: `viennarna_pseudo.rs`

Implements `PseudoEnergyModel` for `ViennaRNA`. The method
`energy_of_pseudo_loop` dispatches on `lp.loop_type`:

- For all standard loop types (`Stack`, `Hairpin`, `Interior`, `Bulge`,
  `Multiloop`, `External`): constructs the corresponding `NearestNeighborLoop`
  and delegates to the existing `ViennaRNA::energy_of_loop`.
- For `Pseudoloop`: computes a simplified multiloop approximation directly
  using `ViennaRNA`'s multiloop parameters.

The three `pub(crate)` fields accessed from `ViennaRNA` are:

| Field | Turner parameter | Role |
|---|---|---|
| `ml_closing` | `ML_closing` | Per-loop initiation penalty |
| `ml_intern` | `ML_intern` | Per-stem penalty |
| `ml_base` | `ML_base` | Per-unpaired-base penalty (reserved for future use) |

---

## Loop Type Dispatch

The table below shows how each `LoopType` is handled in
`energy_of_pseudo_loop`:

| `LoopType` | `NearestNeighborLoop` constructed | Notes |
|---|---|---|
| `Stack` | `Interior { closing, inner }` | 0 unpaired bases on both sides |
| `Interior` | `Interior { closing, inner }` | Standard interior loop |
| `Bulge` | `Interior { closing, inner }` | 1 unpaired base on one side |
| `Hairpin` | `Hairpin { closing }` | Standard hairpin |
| `Multiloop` | `Multibranch { closing, branches }` | Children are always `Single` |
| `External` | `Exterior { ends: (0, n-1), branches }` | `Double` children skipped (see below) |
| `Pseudoloop` | *(direct formula — no `NearestNeighborLoop`)* | See next section |

**`Stack`, `Interior`, `Bulge`** all map to `NearestNeighborLoop::Interior`
because `ViennaRNA::energy_of_loop` dispatches on pair geometry, not loop name.
The `LoopType` distinction matters for the pseudoknot pipeline's loop
enumeration but not for energy evaluation.

**`External` with `Double` children**: when a pseudoknot sits at the top level
of the structure, the `External` loop's children list contains one or more
`ClosingDescriptor::Double` entries. Crossing pairs cannot be sliced as
monotone sequence segments (which `Exterior` evaluation requires), so `Double`
entries are skipped and contribute 0 to the exterior loop energy. The stacking
energy of the pseudoknot's helices is already captured by the `SpanBand` Stack
loops. The terminal AU/GU penalties for the outermost crossing pairs are
currently omitted (see [Known Limitations](#known-limitations)).

---

## The Pseudoloop Approximation

A `Pseudoloop` loop has two crossing closing pairs and no standard
nearest-neighbor parameter. The current implementation uses a **simplified
multiloop model**:

```
E_pseudoloop ≈ ml_closing + ml_intern × (n_children + 2)
```

where:
- `ml_closing` is the per-loop initiation penalty
- `ml_intern` is the per-stem penalty
- `n_children` is the number of children (band-tip pairs in the loop's `children` list)
- the `+2` accounts for the two crossing closing pairs, each counted as one stem

This is the same linear term used in the Turner multiloop model, applied to the
pseudoloop treated as a generalized junction. The unpaired-base penalty
(`ml_base × n_unpaired`) is omitted because the pseudoloop body is topologically
non-planar and counting unpaired bases in it requires access to the pair table.

### Basis

Most published pseudoknot energy models (pknotsRG — Reeder & Giegerich 2004;
NUPACK — Dirks & Pierce 2003) treat the pseudoloop as a multiloop for the
purpose of free energy calculation. The linear approximation recovers the
dominant contribution (initiation + branch penalties) while avoiding
non-planar topology issues.

### Future extension (Strategy C)

A more accurate model would fit dedicated pseudoloop parameters
(`pseudoloop_closing`, `pseudoloop_intern`, `pseudoloop_base`) against
experimental pseudoknot thermodynamic data (e.g. Cao & Chen 2009). Adding this
requires only:

1. New fields in `RNAThermoParams` and each parameter file.
2. Replacing the formula in `energy_of_pseudo_loop` with the dedicated
   parameters.

The architecture does not need to change.

---

## Known Limitations

| Limitation | Location | Workaround / future fix |
|---|---|---|
| Unpaired-base term (`ml_base × n`) missing for pseudoloops | `viennarna_pseudo.rs`, `Pseudoloop` arm | Pass `&PairTable` to `energy_of_pseudo_loop` to count `n_unpaired` from the structure |
| Terminal AU/GU penalties missing for top-level crossing pairs | `viennarna_pseudo.rs`, `External` arm | Iterate `Double` children and add `terminal_ru`/`terminal_ap` per pair |
| No dedicated pseudoloop parameters | `pseudo_energy_model.rs` | Add `pseudoloop_*` fields to `RNAThermoParams` (Strategy C) |
| `energy_of_pseudoknotted_structure` panics on empty sequence | implicit | Guard at call site |

---

## End-to-end Example

### Input

```
Sequence:  G G G C C C C C C A  A  A  G  G  G
Position:  0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
Structure: ( ( ( [ [ [ ) ) ) .  .  .  ]  ]  ]
```

`parse_structure("((([[[)))...]]]")` produces 6 loops (see `Pseudoknots.md`).

### Energy evaluation

`ViennaRNA::energy_of_pseudoknotted_structure` sums `energy_of_pseudo_loop`
for each loop:

| Loop | Type | Evaluation |
|---|---|---|
| `Pseudoloop, closing=((0,8),(3,14)), children=[(2,6),(5,12)]` | `Pseudoloop` | `ml_closing + ml_intern × (2+2)` |
| `Stack, SpanBand, closing=(0,8), inner=(1,7)` | `Stack` → `Interior` | `energy_of_loop(Interior{(0,8),(1,7)})` |
| `Stack, SpanBand, closing=(1,7), inner=(2,6)` | `Stack` → `Interior` | `energy_of_loop(Interior{(1,7),(2,6)})` |
| `Stack, SpanBand, closing=(3,14), inner=(4,13)` | `Stack` → `Interior` | `energy_of_loop(Interior{(3,14),(4,13)})` |
| `Stack, SpanBand, closing=(4,13), inner=(5,12)` | `Stack` → `Interior` | `energy_of_loop(Interior{(4,13),(5,12)})` |
| `External, children=[Double((0,8),(3,14))]` | `External` | 0 (Double child skipped) |

The four `Interior` calls each look up `stack[GC][GC]` (a G-C/G-C stack energy)
from the Turner 2004 parameters. The `Pseudoloop` term uses the multiloop
linear approximation.

Total energy = sum of all six contributions, in dcal/mol.
