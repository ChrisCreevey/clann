#!/usr/bin/env python3
"""Step 0.1 (PLAN_web_client.md): prove one live libclann engine runs a
persistent multi-command session — state survives across commands, no reload.

Simulates the interactive user scenario in a SINGLE process:
  init -> load trees -> set criterion -> nj -> hs -> showtrees -> reconstruct
Each command sees the state the previous ones left behind (unlike the one-shot
pyclann subprocess model).

Build the library first:   make libclann.so
Run:                       python3 tools/session_persistence_check.py
  (On Apple Silicon the Homebrew-gcc build is x86_64; if you hit an
   "incompatible architecture" dlopen error, re-run under a matching
   interpreter, e.g.:  arch -x86_64 /usr/bin/python3 tools/session_persistence_check.py )

Step 0.1b will extend this with a clann_reset() round-trip for determinism.
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

# Capture all Clann output through the library callback.
OUTFN = ctypes.CFUNCTYPE(None, ctypes.c_char_p, ctypes.c_void_p)
buf = []
cb = OUTFN(lambda msg, _ud: buf.append(msg.decode(errors="replace")))
lib.clann_set_output_fn.argtypes = [OUTFN, ctypes.c_void_p]
lib.clann_set_output_fn(cb, None)


def run(label, fn):
    start = len(buf)
    fn()
    text = "".join(buf[start:])
    sig = [l.strip() for l in text.splitlines()
           if any(k in l for k in ("source trees", "score =", "17.0000",
                                   "Neighbor-joining", "dup+loss"))]
    print(f"\n### {label}")
    for l in sig[:4]:
        print("   ", l)
    return text


print("clann_init ->", lib.clann_init())

run("exe tutorial_multicopy.ph",
    lambda: lib.clann_load_trees(b"examples/tutorial_multicopy.ph", b""))
run("set criterion=recon", lambda: lib.clann_run_command(b"set criterion=recon"))
run("nj", lambda: lib.clann_run_command(b"nj open=no"))
run("hs nreps=5 seed=42",
    lambda: lib.clann_run_command(b"hs nreps=5 nthreads=1 seed=42"))
# showtrees proves the ORIGINAL 8 source trees are still in memory after nj+hs.
st = run("showtrees display=no", lambda: lib.clann_run_command(b"showtrees display=no"))
run("reconstruct speciestree=memory",
    lambda: lib.clann_run_command(b"reconstruct speciestree=memory open=no"))

alltext = "".join(buf)
assert "8 source trees" in st, "source trees NOT still in memory after nj+hs!"
assert "17.0000" in alltext or "score = 17" in alltext, "hs recon score missing"
print("\nPASS: one process handled exe->nj->hs->showtrees->reconstruct;")
print("      the 8 source trees survived across all commands (no reload).")
