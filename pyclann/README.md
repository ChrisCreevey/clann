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

---

## Building Clann for Python use

There are two distinct build modes, depending on how you want to use Clann
from Python.

---

### Mode 1 — Standard binary (required for pyclann subprocess approach)

This is the normal build and is all you need for `pyclann` as described above.
The `pyclann` package calls the `clann` command-line binary via `subprocess`.

```bash
# Prerequisites (Debian/Ubuntu)
sudo apt install build-essential libgomp1

# Optional but recommended: readline support for the interactive REPL
sudo apt install libreadline-dev

# Build
./configure
make

# The binary is now at ./clann
# Run directly or install system-wide
sudo make install          # installs to /usr/local/bin/clann (or PREFIX/bin)
```

On **macOS** with Homebrew:

```bash
brew install gcc libomp readline
./configure CC=gcc-14      # use Homebrew GCC for OpenMP
make
```

After building, tell pyclann where the binary is:

```python
import pyclann
pyclann.set_clann_path("/path/to/clann")   # or let it find clann on PATH
```

The standard binary supports OpenMP parallelism, GNU readline, and all
interactive commands.  It is the recommended backend for `pyclann`.

---

### Mode 2 — Shared library (`libclann.so`) for in-process use

An alternative build produces a shared library that can be loaded directly
into Python (or any other host process) via `ctypes` / `cffi`, **without
spawning a subprocess for each call**.  This eliminates the process-creation
and disk-I/O overhead that the subprocess approach incurs, which matters when
you need to score thousands of topologies in a tight loop.

**Key differences from the standard build:**

| Aspect | Standard binary | `libclann.so` |
|--------|----------------|---------------|
| Compilation flag | *(none)* | `-DCLANN_LIBRARY_MODE` |
| Entry point | `main()` (excluded by flag) | `clann_init()` / `clann_run_command()` |
| Fatal errors | `_exit()` terminates process | `longjmp()` returns error code to caller |
| Output | written to stdout / log file | forwarded to a registered C callback |
| Signal handlers | registered for SIGINT etc. | skipped (host process owns signals) |
| OpenMP | works normally | works, but GIL must be released around parallel sections |

**Build the shared library:**

```bash
# After running ./configure (needed to generate config.h)
./configure
make libclann.so
```

This compiles all sources with `-DCLANN_LIBRARY_MODE -fPIC` into a separate
set of object files under `libclann_objs/` and links them into `libclann.so`.
The normal `clann` binary is left untouched.

**Minimal ctypes example:**

```python
import ctypes, os

lib = ctypes.CDLL("./libclann.so")

# Declare function signatures
lib.clann_init.restype = ctypes.c_int
lib.clann_load_trees.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
lib.clann_load_trees.restype  = ctypes.c_int
lib.clann_run_command.argtypes = [ctypes.c_char_p]
lib.clann_run_command.restype  = ctypes.c_int
lib.clann_reset.restype = None

# Register an output callback so output goes to Python instead of stdout
OUTPUT_CB = ctypes.CFUNCTYPE(None, ctypes.c_char_p, ctypes.c_void_p)

output_lines = []

@OUTPUT_CB
def capture(msg, _userdata):
    output_lines.append(msg.decode())

lib.clann_set_output_fn(capture, None)

# Initialise, load trees, run a command
assert lib.clann_init() == 0
assert lib.clann_load_trees(b"examples/tutorial_single.ph", b"") == 0
rc = lib.clann_run_command(b"hs nreps=10 criterion=dfit nthreads=1")
print("return code:", rc)
print("captured lines:", len(output_lines))

# Reset state before the next analysis
lib.clann_reset()
```

The public C API (`clann_init`, `clann_set_output_fn`, `clann_load_trees`,
`clann_run_command`, `clann_reset`) is fully documented in `clann_api.h`.

**When to use each mode:**

* **Use the standard binary + pyclann subprocess** for all normal scripting
  needs.  It is the simplest, most robust approach and works out of the box
  after a normal `make`.

* **Use `libclann.so`** only when you need in-process performance (e.g.
  scoring millions of topologies from Python without subprocess overhead, or
  embedding Clann in a larger C/C++ application).  It requires more careful
  lifecycle management (init / reset / error handling via return codes).

