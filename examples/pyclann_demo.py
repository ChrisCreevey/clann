#!/usr/bin/env python3
"""pyclann_demo.py — Demonstration of the pyclann Python API for Clann.

Run from the repository root after building and installing pyclann:

    # 1. Build Clann (if not already done)
    ./configure && make

    # 2. Install the Python package
    pip install ./pyclann

    # 3. Run this script
    python3 examples/pyclann_demo.py

The script uses the example source-tree files in this directory so no
external data is required.
"""

from __future__ import annotations

import sys
import os
from pathlib import Path

# ---------------------------------------------------------------------------
# Locate the repository root and set the clann binary path so the demo works
# whether or not clann is on PATH.
# ---------------------------------------------------------------------------

HERE = Path(__file__).resolve().parent          # examples/
REPO_ROOT = HERE.parent                          # repository root
CLANN_BIN = REPO_ROOT / "clann"

# Make pyclann importable when running directly from the repository without
# a system-wide install (the package lives at <repo>/pyclann/).
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

try:
    import pyclann
except ModuleNotFoundError:
    sys.exit(
        "pyclann is not installed.  Run:\n"
        "    pip install ./pyclann\n"
        "from the repository root and then re-run this script."
    )

if CLANN_BIN.exists():
    pyclann.set_clann_path(str(CLANN_BIN))
else:
    try:
        pyclann.get_clann_path()   # will raise ClannError if not on PATH
    except pyclann.ClannError:
        sys.exit(
            "clann binary not found.  Build it first:\n"
            "    ./configure && make\n"
            "or install clann and make sure it is on your PATH."
        )

# Input files bundled with the repository
SINGLE_PH     = str(HERE / "tutorial_single.ph")      # 28 primate+rodent gene trees
TREES_PH      = str(HERE / "tutorial_trees.ph")       # 6-taxon demo trees
CANDIDATES_PH = str(HERE / "tutorial_candidates.ph")  # candidate topologies to score


def section(title: str) -> None:
    """Print a section header."""
    print()
    print("=" * 60)
    print(f"  {title}")
    print("=" * 60)


# ---------------------------------------------------------------------------
# 1. Heuristic supertree search (hs)
# ---------------------------------------------------------------------------

section("1. Heuristic supertree search  (pyclann.hs)")

print(
    "\npyclann.hs() runs Clann's SPR/TBR heuristic search to find the\n"
    "supertree that best fits all source trees under a chosen criterion.\n"
)

result = pyclann.hs(
    SINGLE_PH,
    criterion="dfit",   # distance-fit criterion (default)
    nreps=5,            # 5 independent search replicates
    nthreads=1,         # single thread for reproducible demo output
)

print(f"Command executed : {result.command}")
print(f"Return code      : {result.returncode}")
print(f"Best score       : {result.score:.4f}")
print(f"Best supertree   : {result.best_tree}")
print(f"Trees retained   : {len(result.trees)}")

# hs() with ML criterion
print("\n--- Same search with criterion=ml ---")
result_ml = pyclann.hs(
    SINGLE_PH,
    criterion="ml",
    nreps=3,
    nthreads=1,
)
print(f"ML best score  : {result_ml.score:.4f}")
print(f"ML best tree   : {result_ml.best_tree}")


# ---------------------------------------------------------------------------
# 2. Neighbour-joining supertree (nj)
# ---------------------------------------------------------------------------

section("2. Neighbour-joining supertree  (pyclann.nj)")

print(
    "\npyclann.nj() builds a neighbour-joining supertree from pairwise\n"
    "average-consensus distances between source trees.  The result\n"
    "includes branch lengths.\n"
)

result = pyclann.nj(SINGLE_PH)

print(f"Command  : {result.command}")
print(f"NJ tree  : {result.best_tree}")


# ---------------------------------------------------------------------------
# 3. Majority-rule consensus (consensus)
# ---------------------------------------------------------------------------

section("3. Majority-rule consensus  (pyclann.consensus)")

print(
    "\npyclann.consensus() computes a majority-rule (or strict) consensus\n"
    "directly from the input source trees.  Clade support values appear\n"
    "as branch labels in the Newick output.\n"
)

result = pyclann.consensus(SINGLE_PH, percentage=0.5)

print(f"Command        : {result.command}")
print(f"Consensus tree : {result.best_tree}")


# ---------------------------------------------------------------------------
# 4. Exhaustive topology search (alltrees)
# ---------------------------------------------------------------------------

section("4. Exhaustive topology search  (pyclann.alltrees)")

print(
    "\npyclann.alltrees() scores every possible unrooted topology.\n"
    "Only feasible for small taxon sets (the tutorial_trees.ph example\n"
    "has 6 taxa → 105 unrooted topologies).\n"
)

result = pyclann.alltrees(TREES_PH, criterion="dfit")

print(f"Command    : {result.command}")
print(f"Best score : {result.score}")
print(f"Best tree  : {result.best_tree}")
print(f"Total topologies evaluated : {len(result.trees)}")


# ---------------------------------------------------------------------------
# 5. Score user-supplied topologies (usertrees)
# ---------------------------------------------------------------------------

section("5. Score user-supplied topologies  (pyclann.usertrees)")

print(
    "\npyclann.usertrees() scores a set of candidate topologies that you\n"
    "supply, rather than searching for the optimum itself.  Useful for\n"
    "comparing a handful of hypotheses or for topology tests.\n"
)

result = pyclann.usertrees(
    SINGLE_PH,
    CANDIDATES_PH,
    criterion="dfit",
)

print(f"Command           : {result.command}")
print(f"Candidates scored : {len(result.trees)}")
print(f"Best score        : {result.score:.4f}")
print(f"Best tree         : {result.best_tree}")
print()
print("All candidate scores:")
for i, (tree, score) in enumerate(zip(result.trees, result.scores), 1):
    print(f"  [{i:2d}] score={score:.4f}  {tree}")


# ---------------------------------------------------------------------------
# 6. Bootstrap analysis (bootstrap)
# ---------------------------------------------------------------------------

section("6. Bootstrap analysis  (pyclann.bootstrap)")

print(
    "\npyclann.bootstrap() runs a full bootstrap analysis: each replicate\n"
    "resamples the source trees, runs a heuristic search, and the set of\n"
    "bootstrap trees is summarised as a majority-rule consensus with clade\n"
    "support proportions as branch labels.\n"
)

result = pyclann.bootstrap(
    SINGLE_PH,
    nreps=10,       # 10 bootstrap replicates (use ≥100 for real analyses)
    nthreads=1,
)

print(f"Command             : {result.command}")
print(f"Bootstrap consensus : {result.best_tree}")
print(f"Bootstrap trees     : {len(result.trees)}")


# ---------------------------------------------------------------------------
# 7. Low-level run() — any Clann command
# ---------------------------------------------------------------------------

section("7. Low-level access  (pyclann.run)")

print(
    "\npyclann.run() lets you execute any Clann command and inspect the\n"
    "raw stdout.  No output-file parsing is performed.\n"
)

result = pyclann.run("nj", SINGLE_PH)
print(f"Command     : {result.command}")
print(f"Return code : {result.returncode}")
# Print the first few lines of stdout so the demo stays brief
first_lines = result.stdout.splitlines()[:10]
print("stdout (first 10 lines):")
for line in first_lines:
    print(f"  {line}")


# ---------------------------------------------------------------------------
# 8. Accessing captured stdout for custom parsing
# ---------------------------------------------------------------------------

section("8. Accessing captured stdout")

print(
    "\nEvery ClannResult exposes the full stdout captured from the clann\n"
    "process.  This lets you extract information that pyclann doesn't yet\n"
    "parse automatically.\n"
)

result = pyclann.hs(SINGLE_PH, nreps=2, nthreads=1)
print("Lines mentioning 'score' from hs stdout:")
for line in result.stdout.splitlines():
    if "score" in line.lower():
        print(f"  {line}")


# ---------------------------------------------------------------------------
# 9. Error handling
# ---------------------------------------------------------------------------

section("9. Error handling")

print(
    "\nClannError is raised for configuration problems (missing binary,\n"
    "missing input file).  subprocess.TimeoutExpired is raised if the\n"
    "optional timeout parameter is exceeded.\n"
)

from pyclann import ClannError  # noqa: E402

# Missing input file → ClannError
try:
    pyclann.hs("/no/such/file.ph")
except ClannError as exc:
    print(f"ClannError caught (expected): {exc}")

# Missing binary path → ClannError
try:
    pyclann.hs(SINGLE_PH, clann_path="/nonexistent/clann")
except ClannError as exc:
    print(f"ClannError caught (expected): {exc}")


# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------

section("Done")
print("\nAll pyclann demo steps completed successfully.")
