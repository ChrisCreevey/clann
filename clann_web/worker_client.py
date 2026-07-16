"""Server-side proxy for a per-session Clann worker process (Step 3.3, Model A).

:class:`WorkerEngine` spawns ``python -m clann_web.worker`` and speaks the
newline-JSON protocol in :mod:`clann_web.worker`. It exposes the SAME public
surface as the in-process :class:`~clann_web.engine.ClannEngine`
(``load``/``run``/``reset``/``state``/``set_workdir``) so the server code is
unchanged, PLUS two out-of-process powers the in-process engine can't offer:

    * :meth:`terminate` — hard-kill the worker; the guaranteed recovery from a
      runaway/wedged search or the known repeated-``hs`` hang (the server then
      respawns a fresh worker). Session analysis state is lost, which is inherent
      to a kill. Clann's own ``controlc`` SIGINT handlers are interactive (they
      prompt Y/N on stdin, i.e. this protocol pipe), so a graceful in-process
      stop isn't usable in server mode — killing is the clean path.

Requests for one worker are serialised by a lock; the worker itself runs one
analysis at a time. Cross-session parallelism comes from each session owning a
*separate* WorkerEngine (hence a separate process and a separate copy of Clann's
globals).
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import threading

from .engine import ClannError


def _worker_command() -> list[str]:
    """Command to launch a worker interpreter that can load the x86_64 lib.

    Overridable via ``CLANN_WORKER_CMD`` (space-split). On macOS the
    Homebrew-gcc ``libclann-server.so`` is x86_64, so re-invoking the same
    interpreter must be forced under Rosetta (``arch -x86_64``); a bare
    ``/usr/bin/python3`` would launch arm64 and fail the dlopen.
    """
    override = os.environ.get("CLANN_WORKER_CMD")
    if override:
        return override.split() + ["-m", "clann_web.worker"]
    if sys.platform == "darwin":
        return ["arch", "-x86_64", sys.executable, "-m", "clann_web.worker"]
    return [sys.executable, "-m", "clann_web.worker"]


class WorkerEngine:
    """Drives one Clann worker subprocess. API-compatible with ClannEngine."""

    def __init__(self, workdir: str | None = None):
        self._lock = threading.Lock()
        self._workdir = workdir
        self._proc: subprocess.Popen | None = None
        self._start()
        if workdir:
            self.set_workdir(workdir)

    # -- process lifecycle -------------------------------------------------
    def _start(self) -> None:
        repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        env = dict(os.environ)
        # Make `clann_web` importable in the child regardless of how we launched.
        env["PYTHONPATH"] = repo_root + os.pathsep + env.get("PYTHONPATH", "")
        self._proc = subprocess.Popen(
            _worker_command(),
            stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            cwd=repo_root, env=env, text=True, bufsize=1,
        )
        # Wait for the worker to finish loading + initialising the library.
        msg = self._read_msg()
        if msg is None or msg.get("type") != "ready":
            err = (msg or {}).get("error", "worker did not start")
            self.terminate()
            raise ClannError(f"worker startup failed: {err}")

    def _read_msg(self) -> dict | None:
        assert self._proc is not None
        line = self._proc.stdout.readline()
        if not line:            # EOF — worker died/was killed
            return None
        try:
            return json.loads(line)
        except json.JSONDecodeError:
            return None

    def _request(self, req: dict, on_line=None):
        """Send one request, pump output, return the op's value (or raise)."""
        with self._lock:
            proc = self._proc
            if proc is None or proc.poll() is not None:
                raise ClannError("worker is not running")
            proc.stdin.write(json.dumps(req) + "\n")
            proc.stdin.flush()
            while True:
                msg = self._read_msg()
                if msg is None:
                    raise ClannError("worker exited during request "
                                     f"(op={req.get('op')})")
                kind = msg.get("type")
                if kind == "line":
                    if on_line is not None:
                        try:
                            on_line(msg.get("text", ""))
                        except Exception:   # a broken listener must not break us
                            pass
                elif kind == "result":
                    return msg.get("value")
                elif kind == "error":
                    raise ClannError(msg.get("error", "worker error"))

    # -- public API (mirrors ClannEngine) ----------------------------------
    def load(self, filename: str, parse_opts: str = "") -> str:
        return self._request({"op": "load", "filename": filename,
                              "parse_opts": parse_opts})

    def run(self, command: str, on_line=None) -> str:
        return self._request({"op": "run", "command": command}, on_line=on_line)

    def reset(self) -> None:
        return self._request({"op": "reset"})

    def state(self) -> dict:
        return self._request({"op": "state"})

    def set_workdir(self, path: str) -> None:
        self._workdir = path
        return self._request({"op": "set_workdir", "path": path})

    # -- cancellation (Step 3.2) -------------------------------------------
    def terminate(self) -> None:
        """Hard-kill the worker process (guaranteed recovery from a hang)."""
        proc = self._proc
        self._proc = None
        if proc is None:
            return
        try:
            proc.kill()
        except ProcessLookupError:
            pass
        try:
            proc.wait(timeout=5)
        except Exception:       # noqa: BLE001
            pass

    def is_alive(self) -> bool:
        return self._proc is not None and self._proc.poll() is None
