"""Per-session Clann worker subprocess (PLAN_web_client.md Step 3.3, Model A).

Run as ``python -m clann_web.worker``. Owns ONE in-process :class:`ClannEngine`
(one independent copy of Clann's global state) and speaks a tiny newline-delimited
JSON protocol on stdin/stdout so the web server can drive it out-of-process:

    request  (server -> worker, one JSON object per line):
        {"op": "load", "filename": "...", "parse_opts": ""}
        {"op": "run",  "command": "hs nreps=5"}
        {"op": "reset"}
        {"op": "state"}
        {"op": "set_workdir", "path": "/tmp/clannweb-xxxx"}
        {"op": "quit"}

    response (worker -> server, one JSON object per line):
        {"type": "ready"}                        # emitted once at startup
        {"type": "line", "text": "..."}          # live output chunk during a run
        {"type": "result", "value": <any>}       # the op finished; value is its return
        {"type": "error", "error": "..."}        # the op raised

Because each session gets its own worker *process*, sessions are isolated: a
crash, a runaway ``hs``, or the known repeated-high-``nreps`` hang kills only that
worker, and the server recovers by killing + respawning it (Step 3.2 cancel).
"""

from __future__ import annotations

import json
import sys

from .engine import ClannEngine, ClannError


def _emit(obj: dict) -> None:
    sys.stdout.write(json.dumps(obj) + "\n")
    sys.stdout.flush()


def main() -> int:
    try:
        eng = ClannEngine()
    except Exception as e:  # noqa: BLE001 - report to the parent and exit
        _emit({"type": "error", "error": f"engine init failed: {e}"})
        return 1
    _emit({"type": "ready"})

    for raw in sys.stdin:
        raw = raw.strip()
        if not raw:
            continue
        try:
            req = json.loads(raw)
        except json.JSONDecodeError as e:
            _emit({"type": "error", "error": f"bad request: {e}"})
            continue

        op = req.get("op")
        try:
            if op == "quit":
                break
            elif op == "load":
                value = eng.load(req["filename"], req.get("parse_opts", ""))
            elif op == "run":
                value = eng.run(req["command"],
                                on_line=lambda t: _emit({"type": "line", "text": t}))
            elif op == "reset":
                eng.reset()
                value = None
            elif op == "state":
                value = eng.state()
            elif op == "set_workdir":
                eng.set_workdir(req["path"])
                value = None
            else:
                _emit({"type": "error", "error": f"unknown op: {op!r}"})
                continue
            _emit({"type": "result", "value": value})
        except (ClannError, KeyError, Exception) as e:  # noqa: BLE001
            _emit({"type": "error", "error": str(e)})
    return 0


if __name__ == "__main__":
    sys.exit(main())
