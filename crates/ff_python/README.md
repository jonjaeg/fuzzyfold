# ff_python

Python bindings for **fuzzyfold** (energy evaluation + kinetic simulation) and a
collection of developer tooling for parameter generation, verification, and validation.

---

## Python bindings

The crate is a `cdylib` built with [PyO3](https://pyo3.rs) and
[maturin](https://www.maturin.rs). It exports:

| Class | Description |
|---|---|
| `ViennaRNA` | Nearest-neighbor free energy evaluation |
| `Simulator` | Gillespie SSA kinetic simulation |

### Build and install (development)

```bash
maturin develop --release
```

See `example.py` for basic usage.

---

## Directory layout

```
generate/       scripts that produce artifacts
  fm363_to_rust.py      convert HotKnots fm363 → Rust parameter source (AndronescuParams)
  fm363_to_rust.md      documentation for the converter
  gen_notebook.py       generate notebooks/region_tree_explorer.ipynb
  gen_bl_notebook.py    generate notebooks/bl_construction_explorer.ipynb

audit/          scripts that verify Rust constants match their source files
  compare_nn_params.py  diff an fm363 file against a Rust parameter directory
  compare_pk_params.py  diff HotKnots PK params against pseudoknot_params.rs

validate/       scripts that evaluate energies and compare to references
  validate_pseudo.py      T-Train PK (22 pseudoknotted structures)
  validate_non_pseudo.py  T-Train PKfree (277 single-strand non-PK structures)

data/           reference data (committed)
  RNA-thermo-db_v1.4.xml        sequences and structures
  energies_DP03_T-Train.txt     HotKnots v2 DP03 energies for T-Train
  energies_DP09_T-Train.txt     HotKnots v2 DP09 energies for T-Train

results/        validation output (git-ignored)
notebooks/      Jupyter notebooks
```

---

## Parameter tooling (`generate/`, `audit/`)

```bash
make generate-mt09   # regenerate rna_mt09/ from HotKnots DP09 fm363
make compare         # verify all Rust constants match source files
make regen           # generate-mt09 + compare in one step
```

Both `validate_pseudo.py` and `validate_non_pseudo.py` accept
`--pk-params {dp03,dp09,mt09}`, `--binary`, `--xml`, `--results`, `--tol`, `--celsius`.
Output columns: `exp`, `HotKnots`, `ff-calc`, `Δ(HK−exp)`, `Δ(ff−exp)`, `Δ(ff−HK)`.

---

## Validation (`validate/`)

```bash
make validate        # runs all four combinations; saves to results/
make summary         # print stats from existing results without re-running
```

See `summary.md` for a full write-up of results and error analysis.

---

## Makefile targets

```
make help            # list all targets
make all             # build + test + validate
make build           # cargo build --release (ff-calc-pseudo)
make test            # cargo test -p ff_energy
make generate-mt09   # regenerate rna_mt09/ from HotKnots fm363
make compare         # verify all Rust constants match source files
make regen           # generate-mt09 + compare
make validate        # run all validations
make summary         # print stats from existing results/
make clean-results   # remove results/
```
