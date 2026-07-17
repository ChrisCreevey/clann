#!/usr/bin/env python3
"""Step 0.1b (PLAN_web_client.md): confirm clann_reset() returns the engine to a
clean baseline — a fresh session after reset reproduces results deterministically.

In ONE process:
  init -> run scenario (A) -> clann_reset() -> init -> run scenario (B)
Asserts the stable result lines of run B are identical to run A (which is itself
a fresh-init baseline), so ending a session and starting a new one in the same
process leaves no state behind. This underpins Model A worker reuse (§2.1/§3.1).

The seed is pinned with `set seed=42` (hs has no seed= option — only the `set
seed=` command calls srand(), main.c:4335), so both runs drive the RNG from an
identical state before hs. Only timing-independent RESULT lines are compared
(final scores, dup+loss total, per-tree reconstruct scores) — the search-progress
lines carry elapsed time and are deliberately excluded.

Build first:  make libclann.so
Run:          python3 tools/session_reset_check.py
  (Apple Silicon + Homebrew-gcc x86_64 build: run under
   `arch -x86_64 /usr/bin/python3 tools/session_reset_check.py`)
"""
import ctypes
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIB = os.path.join(REPO, "libclann.so")

if not os.path.exists(LIB):
    sys.exit(f"{LIB} not found — build it first with `make libclann.so`")

os.chdir(REPO)
try:
    lib = ctypes.CDLL(LIB)
except OSError as e:
    sys.exit(f"dlopen failed ({e}).\nIf this is an architecture mismatch, run under "
             f"a matching interpreter, e.g.:\n  arch -x86_64 /usr/bin/python3 {__file__}")

lib.clann_init.restype = ctypes.c_int
lib.clann_load_trees.restype = ctypes.c_int
lib.clann_load_trees.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
lib.clann_run_command.restype = ctypes.c_int
lib.clann_run_command.argtypes = [ctypes.c_char_p]
lib.clann_reset.restype = None

OUTFN = ctypes.CFUNCTYPE(None, ctypes.c_char_p, ctypes.c_void_p)
buf = []
cb = OUTFN(lambda msg, _ud: buf.append(msg.decode(errors="replace")))
lib.clann_set_output_fn.argtypes = [OUTFN, ctypes.c_void_p]
lib.clann_set_output_fn(cb, None)

# Stable, timing-independent result markers.
KEEP = ("source trees met", "total dup+loss score", "Tree number:")


def signature(text):
    """Extract deterministic result lines (final scores, dup/loss, per-tree)."""
    out = []
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        if "Supertree" in s and "score =" in s:      # final best-tree score
            out.append(s)
        elif any(k in s for k in KEEP):
            out.append(s)
    return out


def scenario():
    """Run the full user workflow; return its captured text."""
    start = len(buf)
    lib.clann_load_trees(b"examples/tutorial_multicopy.ph", b"")
    lib.clann_run_command(b"set seed=42")
    lib.clann_run_command(b"set criterion=recon")
    lib.clann_run_command(b"nj open=no")
    lib.clann_run_command(b"hs nreps=5 nthreads=1")
    lib.clann_run_command(b"showtrees display=no")
    lib.clann_run_command(b"reconstruct speciestree=memory open=no")
    return "".join(buf[start:])


print("clann_init ->", lib.clann_init())
textA = scenario()
sigA = signature(textA)

print("clann_reset()")
lib.clann_reset()
# clann_reset() calls clann_init() internally; the callback survives it.
textB = scenario()
sigB = signature(textB)

print(f"\nrun A: {len(sigA)} result lines")
for l in sigA:
    print("   A:", l)
print(f"\nrun B (post-reset): {len(sigB)} result lines")
for l in sigB:
    print("   B:", l)

assert sigA, "run A produced no result lines — scenario broken"
assert "32.0000" in textA and "32.0000" in textB, "recon 32.0000 regression missing"
assert sigA == sigB, (
    "MISMATCH: post-reset run differs from fresh baseline — state leaked across "
    "clann_reset(). Diff:\n"
    + "\n".join(f"  A: {a!r}\n  B: {b!r}"
                for a, b in zip(sigA, sigB) if a != b)
)
print("\nPASS: clann_reset() gives a clean baseline; post-reset run reproduces the "
      "fresh run exactly (result lines identical, recon 32.0000 both times).")
