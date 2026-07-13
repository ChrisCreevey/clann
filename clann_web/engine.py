"""In-process Clann engine wrapper for the web server (PLAN_web_client.md Step 1.1).

Loads ``libclann-server.so`` (the CLANN_SERVER_MODE build — no shell/system()
surface, Step 0.2) via ctypes and exposes a small, serialised API:

    eng = ClannEngine()
    eng.load("examples/tutorial_multicopy.ph")
    eng.run("set criterion=recon")
    eng.run("hs nreps=5")
    eng.state()          # {"input_file", "num_taxa", "num_source_trees", ...}

Because Clann keeps all analysis state in process globals — including OpenMP
``threadprivate`` scratch — it is NOT re-entrant AND must always run on the SAME
OS thread. Calling it from a pool of HTTP handler threads segfaults. So the
engine owns one **dedicated worker thread**: the library is loaded, initialised,
and every subsequent call is executed there; public methods marshal work onto it
and block for the result. This serialises all analysis (Model B in the plan —
one engine per server for now; Step 3.3 moves to one worker process per session)
and gives a natural home for the eventual per-session workers. ``reset()``
returns the engine to a clean baseline (proven deterministic in Step 0.1b).
"""

from __future__ import annotations

import ctypes
import os
import queue
import threading

# Criterion integer -> name (matches Clann's `set criterion` values).
_CRITERION_NAMES = {
    0: "dfit", 1: "mrp", 2: "sfit", 3: "qfit",
    5: "recon", 6: "rf", 7: "ml",
}


def _default_lib_path() -> str:
    env = os.environ.get("CLANN_LIB")
    if env:
        return env
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(repo, "libclann-server.so")


_OUTPUT_FN = ctypes.CFUNCTYPE(None, ctypes.c_char_p, ctypes.c_void_p)


class ClannError(RuntimeError):
    pass


class ClannEngine:
    """One live libclann engine, pinned to a dedicated worker thread."""

    def __init__(self, lib_path: str | None = None, workdir: str | None = None):
        self.lib_path = lib_path or _default_lib_path()
        if not os.path.exists(self.lib_path):
            raise ClannError(
                f"{self.lib_path} not found — build it with `make libclann-server.so`"
            )
        # Clann resolves relative paths (examples, output files) against cwd.
        self.workdir = workdir or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

        self._buf: list[str] = []
        self._input_file: str | None = None
        self._lib = None
        self._cb = None  # strong ref to the ctypes callback

        # Dedicated worker thread that owns the library.
        self._tasks: "queue.Queue" = queue.Queue()
        self._init_error: Exception | None = None
        self._ready = threading.Event()
        self._worker = threading.Thread(target=self._worker_loop, daemon=True,
                                        name="clann-engine")
        self._worker.start()
        self._ready.wait()
        if self._init_error is not None:
            raise self._init_error

    # -- worker thread -----------------------------------------------------
    def _worker_loop(self) -> None:
        try:
            self._load_and_init()      # library bound to THIS thread
        except Exception as e:         # noqa: BLE001 - surfaced to constructor
            self._init_error = e
            self._ready.set()
            return
        self._ready.set()
        while True:
            item = self._tasks.get()
            if item is None:
                break
            func, result = item
            try:
                result["value"] = func()
            except Exception as e:     # noqa: BLE001 - surfaced to caller
                result["error"] = e
            result["event"].set()

    def _submit(self, func):
        result = {"event": threading.Event()}
        self._tasks.put((func, result))
        result["event"].wait()
        if "error" in result:
            raise result["error"]
        return result.get("value")

    # -- internals (all run on the worker thread) --------------------------
    def _on_output(self, msg: bytes, _userdata) -> None:
        self._buf.append(msg.decode(errors="replace"))

    def _load_and_init(self) -> None:
        try:
            self._lib = ctypes.CDLL(self.lib_path)
        except OSError as e:
            raise ClannError(
                f"could not load {self.lib_path} ({e}). On Apple Silicon the "
                f"Homebrew-gcc build is x86_64 — run the server under a matching "
                f"interpreter (e.g. `arch -x86_64 /usr/bin/python3`)."
            ) from e

        self._lib.clann_init.restype = ctypes.c_int
        self._lib.clann_load_trees.restype = ctypes.c_int
        self._lib.clann_load_trees.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
        self._lib.clann_run_command.restype = ctypes.c_int
        self._lib.clann_run_command.argtypes = [ctypes.c_char_p]
        self._lib.clann_reset.restype = None
        self._lib.clann_set_output_fn.argtypes = [_OUTPUT_FN, ctypes.c_void_p]

        self._cb = _OUTPUT_FN(self._on_output)
        os.chdir(self.workdir)
        if self._lib.clann_init() != 0:
            raise ClannError("clann_init() failed (allocation error)")
        self._lib.clann_set_output_fn(self._cb, None)

    def _drain(self) -> str:
        text = "".join(self._buf)
        self._buf.clear()
        return text

    def _int_global(self, name: str, default: int = -1) -> int:
        try:
            return ctypes.c_int.in_dll(self._lib, name).value
        except ValueError:
            return default

    def _do_load(self, filename: str, parse_opts: str) -> str:
        os.chdir(self.workdir)
        self._buf.clear()
        self._lib.clann_load_trees(filename.encode(), parse_opts.encode())
        self._input_file = filename
        return self._drain()

    def _do_run(self, command: str) -> str:
        os.chdir(self.workdir)
        self._buf.clear()
        self._lib.clann_run_command(command.encode())
        return self._drain()

    def _do_reset(self) -> None:
        os.chdir(self.workdir)
        self._lib.clann_reset()                     # frees + re-inits globals
        self._lib.clann_set_output_fn(self._cb, None)  # re-arm callback
        self._buf.clear()
        self._input_file = None

    def _do_state(self) -> dict:
        crit = self._int_global("criterion", 0)
        return {
            "input_file": self._input_file,
            "num_taxa": self._int_global("number_of_taxa", 0),
            "num_source_trees": self._int_global("Total_fund_trees", 0),
            "criterion": _CRITERION_NAMES.get(crit, str(crit)),
            "trees_in_memory": self._int_global("trees_in_memory", 0),
        }

    # -- public API (marshalled onto the worker thread) --------------------
    def load(self, filename: str, parse_opts: str = "") -> str:
        """Load source trees (equivalent to `exe <filename>`). Returns log text."""
        return self._submit(lambda: self._do_load(filename, parse_opts))

    def run(self, command: str) -> str:
        """Run one Clann command against the live session. Returns log text."""
        return self._submit(lambda: self._do_run(command))

    def reset(self) -> None:
        """Return the engine to a clean baseline (new session)."""
        return self._submit(self._do_reset)

    def state(self) -> dict:
        """Current session state read from Clann's globals."""
        return self._submit(self._do_state)

    def set_workdir(self, path: str) -> None:
        """Point Clann's working directory (where relative paths resolve and
        output files land) at `path` — used to confine a session to its sandbox.
        Applied on the next load/run/reset, all of which chdir to self.workdir."""
        self.workdir = path


# One process can hold only ONE Clann engine: the library keeps all analysis
# state in globals, so a second ClannEngine would share (and corrupt) the first's
# state. This process-wide singleton enforces that — every server in a process
# shares it. (True per-session isolation comes with per-process workers, Step 3.3.)
_SHARED: ClannEngine | None = None


def get_shared_engine(lib_path: str | None = None) -> ClannEngine:
    global _SHARED
    if _SHARED is None:
        _SHARED = ClannEngine(lib_path=lib_path)
    return _SHARED
