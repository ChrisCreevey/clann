# pyclann — Python interface to Clann

`pyclann` is a Python package that lets you drive the
[Clann](https://github.com/ChrisCreevey/clann) phylogenomics tool from Python
scripts, Jupyter notebooks, Snakemake workflows, and any other Python
environment.

---

## Requirements

* Python ≥ 3.9
* The `clann` binary must be on your `PATH`, *or* you must call
  `pyclann.set_clann_path(...)` before use.

---

## Installation

From inside the Clann repository root:

```bash
pip install ./pyclann
```

Or, for an editable development install:

```bash
pip install -e ./pyclann
```

---

## Quick start

```python
import pyclann

# Point pyclann at the clann binary if it is not on PATH
# pyclann.set_clann_path("/path/to/clann")

# Heuristic supertree search
result = pyclann.hs("examples/tutorial_single.ph", criterion="dfit", nreps=10)
print(result.best_tree)   # Newick string
print(result.score)       # float

# Neighbour-joining tree
result = pyclann.nj("examples/tutorial_single.ph")
print(result.best_tree)

# Majority-rule consensus
result = pyclann.consensus("examples/tutorial_single.ph")
print(result.best_tree)

# Score user-supplied topologies
result = pyclann.usertrees(
    "examples/tutorial_single.ph",
    "examples/tutorial_candidates.ph",
    criterion="ml",
)
for tree, score in zip(result.trees, result.scores):
    print(f"score={score:.4f}  {tree}")

# Bootstrap analysis
result = pyclann.bootstrap("examples/tutorial_single.ph", nreps=100, nthreads=4)
print(result.best_tree)   # majority-rule consensus Newick
```

---

## API reference

### `pyclann.hs(treefile, *, criterion=None, nreps=None, nthreads=None, nbest=None, start=None, seed=None, savetrees=None, **kwargs) → ClannResult`

Heuristic supertree search.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `treefile` | path | — | Source-tree input file |
| `criterion` | str | `"dfit"` | `"dfit"`, `"ml"`, `"rf"`, `"sfit"`, `"qfit"`, `"avcon"` |
| `nreps` | int | 10 | Number of search replicates |
| `nthreads` | int | all CPUs | OpenMP thread count |
| `nbest` | int | 1 | Number of best trees to retain |
| `start` | str | `"nj"` | Starting tree: `"nj"` or `"memory"` |
| `seed` | int | random | RNG seed for reproducibility |
| `savetrees` | str | `"Heuristic_result.txt"` | Output filename override |

Extra keyword arguments are passed directly as `key=value` tokens to Clann.

---

### `pyclann.nj(treefile, *, missing=None, savetrees=None, **kwargs) → ClannResult`

Neighbour-joining supertree from average-consensus distances.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `missing` | `"4point"` | Missing-data method: `"4point"` or `"ultrametric"` |
| `savetrees` | `"NJ-tree.ph"` | Output filename |

---

### `pyclann.consensus(treefile, *, percentage=None, savetrees=None, **kwargs) → ClannResult`

Compute a majority-rule (or strict) consensus tree.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `percentage` | 0.5 | Minimum clade support (0.5–1.0) |
| `savetrees` | `"consensus.ph"` | Output filename |

---

### `pyclann.alltrees(treefile, *, criterion=None, savetrees=None, seed=None, **kwargs) → ClannResult`

Exhaustively score all possible topologies (only feasible for ≤ 8 taxa).

---

### `pyclann.usertrees(treefile, candidates, *, criterion=None, tests=None, nboot=None, savetrees=None, **kwargs) → ClannResult`

Score a set of user-supplied candidate topologies.

| Parameter | Description |
|-----------|-------------|
| `candidates` | File of candidate Newick trees to score |
| `tests` | `True` to run ML topology tests (KH, SH, WS; requires `criterion="ml"`) |
| `nboot` | Bootstrap replicates for topology tests |

---

### `pyclann.bootstrap(treefile, *, criterion=None, nreps=None, nthreads=None, seed=None, savetrees=None, **kwargs) → ClannResult`

Bootstrap analysis: heuristic search on resampled source trees followed by a
majority-rule consensus.

---

### `pyclann.run(command, treefile=None, *, **kwargs) → ClannResult`

Low-level: run any Clann command string.  No output-file parsing is performed;
inspect `.stdout` for results.

```python
result = pyclann.run("rfdists", "examples/tutorial_single.ph")
print(result.stdout)
```

---

### `ClannResult`

All functions return a `ClannResult` dataclass:

| Attribute | Type | Description |
|-----------|------|-------------|
| `best_tree` | `str \| None` | Best Newick tree string |
| `trees` | `list[str]` | All returned trees |
| `score` | `float \| None` | Best score |
| `scores` | `list[float]` | Scores for all trees |
| `stdout` | `str` | Full captured stdout |
| `stderr` | `str` | Full captured stderr |
| `returncode` | `int` | Process exit code |
| `command` | `str` | Executed command string |

---

## Error handling

```python
import pyclann
from pyclann import ClannError

try:
    result = pyclann.hs("trees.ph", nreps=50, timeout=300)
except ClannError as exc:
    print(f"Clann error: {exc}")
except subprocess.TimeoutExpired:
    print("Clann timed out")

if result.returncode != 0:
    print("Clann returned non-zero exit code")
    print(result.stderr)
```

---

## Integration with workflow managers

### Snakemake

```python
import pyclann, os

rule supertree:
    input: "genes/{gene}.ph"
    output: "results/{gene}_supertree.txt"
    run:
        result = pyclann.hs(input[0], criterion="ml", nreps=20, nthreads=4)
        with open(output[0], "w") as fh:
            fh.write(f"{result.best_tree}\t{result.score}\n")
```

### Nextflow (PyScript process)

```python
import pyclann
result = pyclann.hs("$treefile", criterion="dfit", nreps=10)
open("supertree.ph", "w").write(result.best_tree + "\n")
```

---

## Architecture note

This package uses the **subprocess approach**: each call spawns a fresh
`clann` process in an isolated temporary directory, captures its output, and
parses the result files it writes to disk.  This is the safest integration
layer and requires no changes to the Clann C source code.

A future `libclann.so`-based in-process layer (without disk overhead) is
described in `clann_api.h` and `clann_api.c` in the Clann source tree.  When
built, it will be used automatically as a drop-in backend for these same
functions.
