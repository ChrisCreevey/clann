"""Minimal HTTP server exposing the Clann engine (PLAN_web_client.md Step 1.1).

Zero external dependencies (stdlib ``http.server``) so it runs under the same
x86_64 interpreter that can load ``libclann-server.so``. The HTTP layer is
deliberately thin — a later step swaps it for FastAPI without touching
``engine.py``.

Endpoints (JSON in/out), single shared engine (Model B):
    POST /api/session            -> reset engine + fresh sandbox {session_id, state}
    POST /api/files?name=FILE    -> upload a file (raw body) into the sandbox
    POST /api/load  {file}       -> load an uploaded file, return {ok, log, state}
    POST /api/run   {command}    -> run one command, return {ok, log, state}
    GET  /api/session            -> current {session_id, state, files}

Every session is confined to its own temp sandbox: uploads land there, the engine
runs there, and all file paths in commands are confined to it (Step 1.2).
Binds 127.0.0.1 only.
"""

from __future__ import annotations

import json
import os
import threading
import time
import uuid
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs, unquote

from .engine import ClannError
from .worker_client import WorkerEngine
from .sandbox import (Sandbox, UnsafePath, safe_basename,
                      sanitize_command, sanitize_options)
from .results import (RESULT_JSON, is_tree_command, build_results, parse_tree_counts,
                      parse_tree_file, looks_like_tree_file, build_tree_view_document)
from .commands import list_commands, command_schema
from .viewer import build_viewer_html, placeholder_html

STATIC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "static")

# Request-size caps (Step 4.1): reject oversized bodies up front rather than
# reading an attacker-chosen amount into memory. Tree files are small.
MAX_UPLOAD_BYTES = 64 * 1024 * 1024   # 64 MiB per uploaded file
MAX_JSON_BYTES = 256 * 1024           # 256 KiB per JSON request body


class RequestTooLarge(ValueError):
    """Raised when a request body exceeds its configured cap."""


class Job:
    """One asynchronous command run: accumulates log lines live and, on
    completion, a structured result. Read from the HTTP thread without touching
    the engine, so status/log polling never blocks behind a running search."""

    def __init__(self, job_id: str, command: str):
        self.id = job_id
        self.command = command
        self.status = "running"        # running | done | error | cancelled
        self.cancel_requested = False
        self._lines: list[str] = []
        self._lock = threading.Lock()
        self.result: dict | None = None
        self.error: str | None = None
        self.done = threading.Event()

    def append(self, text: str) -> None:
        with self._lock:
            self._lines.append(text)

    def log(self) -> str:
        with self._lock:
            return "".join(self._lines)

    def to_dict(self, include_log: bool = True) -> dict:
        d = {"id": self.id, "command": self.command, "status": self.status}
        if include_log:
            d["log"] = self.log()
        if self.result is not None:
            d.update(self.result)
        if self.error:
            d["error"] = self.error
        return d


class _App:
    """One session: its own worker process (Model A) + sandbox.

    Each session owns a separate :class:`WorkerEngine` — a distinct process with
    its own copy of Clann's globals — so sessions are isolated and a runaway or
    wedged search in one can be killed without touching the others.
    """

    def __init__(self, engine: WorkerEngine | None = None):
        self.sandbox = Sandbox()
        self.engine = engine or WorkerEngine(workdir=self.sandbox.dir)
        self.engine.set_workdir(self.sandbox.dir)
        self.last_result_json = None   # raw JSON of the most recent tree result
        self.tree_counts: dict = {}     # single-/multi-copy counts parsed from exe
        self.jobs: dict[str, Job] = {}
        self.current_job: Job | None = None
        self._job_lock = threading.Lock()
        self.session_id = uuid.uuid4().hex
        self.command_log: list[str] = []   # commands sent, for a reproducible script
        self.initial_seed: int | None = None

    def record(self, cmd: str) -> None:
        """Append a command to the reproducibility log. The first call snapshots
        the session's seed (before any `set seed=` has run) so the generated
        script can re-seed the RNG identically on replay."""
        if self.initial_seed is None:
            try:
                self.initial_seed = self.engine.state().get("seed")
            except Exception:      # noqa: BLE001
                self.initial_seed = None
        self.command_log.append(cmd.strip())

    def commands_script(self) -> str:
        """Assemble a replayable Clann command script for the session."""
        import datetime
        out = [
            "# Clann session commands — generated "
            + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "# Replay:  clann -c <this-file>",
            "#   (put the uploaded input tree file(s) in the same directory).",
            "# Reproducibility: the seed below reproduces single-threaded runs;",
            "#   multi-threaded random-start searches may still vary "
            "(known issue — run with nthreads=1 for full determinism).",
            "",
        ]
        if self.initial_seed is not None:
            out.append(f"set seed={self.initial_seed}")
        out.extend(self.command_log)
        if not self.command_log or self.command_log[-1].split()[:1] not in (["quit"], ["exit"]):
            out.append("quit")
        return "\n".join(out) + "\n"

    def state(self) -> dict:
        """Engine state augmented with the single-/multi-copy counts Clann
        reports at load time (not held in engine globals, so tracked here)."""
        s = self.engine.state()
        s.update(self.tree_counts)
        return s

    def output_files(self) -> list:
        """Files in the session sandbox the user may want to download (every
        file except the reserved result-JSON scratch file)."""
        return [f for f in self.sandbox.list() if f != RESULT_JSON]

    def busy(self) -> bool:
        """True if a command is still running on this session's engine. The
        engine runs one thing at a time, so a load issued now would block behind
        it — callers should reject instead of hanging."""
        with self._job_lock:
            return (self.current_job is not None
                    and not self.current_job.done.is_set())

    def _respawn_engine(self) -> None:
        """Start a fresh worker on the same sandbox and kill the old one.

        Used to cancel a search (Step 3.2). Session analysis state (loaded trees,
        `set` options) is lost — inherent to a kill. The new engine is installed
        BEFORE the old is killed so `self.engine` is never a dead/None handle.
        """
        old = self.engine
        self.engine = WorkerEngine(workdir=self.sandbox.dir)
        self.last_result_json = None
        try:
            old.terminate()   # in-flight run unblocks via EOF; job -> cancelled
        except Exception:      # noqa: BLE001
            pass

    def new_session(self) -> dict:
        old_engine = self.engine
        old_sandbox = self.sandbox
        self.sandbox = Sandbox()
        self.engine = WorkerEngine(workdir=self.sandbox.dir)
        try:
            old_engine.terminate()
        except Exception:      # noqa: BLE001
            pass
        old_sandbox.destroy()
        self.last_result_json = None
        self.tree_counts = {}
        self.command_log = []
        self.initial_seed = None
        with self._job_lock:
            self.jobs.clear()
            self.current_job = None
        self.session_id = uuid.uuid4().hex
        return {"session_id": self.session_id, "state": self.state(),
                "files": self.sandbox.list()}

    def cancel_current(self, job_id: str | None = None) -> Job | None:
        """Cancel the running search by killing its worker and respawning a fresh
        one. Clann's own SIGINT handlers are interactive (they prompt Y/N on
        stdin, which in server mode is the worker's protocol pipe), so a graceful
        in-process stop isn't possible — a killable worker process is the only
        clean recovery, which is exactly why each session owns one. This also
        recovers the known repeated-high-`nreps` `hs` hang."""
        with self._job_lock:
            job = self.current_job
            if job is None or job.done.is_set():
                return None
            if job_id is not None and job.id != job_id:
                return None
            job.cancel_requested = True
        self._respawn_engine()
        return job

    def start_job(self, safe_cmd: str) -> Job:
        """Begin an async run of `safe_cmd` (already sandbox-sanitised).

        Raises RuntimeError('busy') if a job is still running — the engine can
        only run one analysis at a time. Returns the new Job immediately; the
        actual work happens on a background thread.
        """
        with self._job_lock:
            if self.current_job is not None and not self.current_job.done.is_set():
                raise RuntimeError("busy")
            job = Job(uuid.uuid4().hex, safe_cmd)
            self.jobs[job.id] = job
            self.current_job = job
        threading.Thread(target=self._run_job, args=(job,), daemon=True).start()
        return job

    def _run_job(self, job: Job) -> None:
        # Inject resultjson= for tree-producing commands (controlled bare name).
        jpath = os.path.join(self.sandbox.dir, RESULT_JSON)
        if os.path.exists(jpath):
            os.remove(jpath)
        run_cmd = job.command
        toks = job.command.split()
        # showtrees: default display=no in the web app — the ASCII trees are
        # slow to write for large sets and the interactive viewer shows them
        # anyway. Only inject when the user didn't specify display=.
        if toks and toks[0] == "showtrees" and "display=" not in job.command:
            run_cmd = f"{run_cmd} display=no"
        if is_tree_command(job.command) and "resultjson=" not in job.command:
            run_cmd = f"{run_cmd} resultjson={RESULT_JSON}"
        engine = self.engine          # bind now; a cancel may swap self.engine
        try:
            engine.run(run_cmd, on_line=job.append)        # blocks this thread
            result = {"state": self.state()}
            if os.path.exists(jpath):
                try:
                    with open(jpath, encoding="utf-8") as jf:
                        self.last_result_json = jf.read()
                    res = build_results(jpath, job.log())
                    result["trees"] = res["trees"]
                    result["scores"] = res["scores"]
                    result["result_type"] = res["type"]
                    result["has_viewer"] = True
                finally:
                    os.remove(jpath)
            # surface any output files the command wrote into the sandbox
            result["output_files"] = self.output_files()
            job.result = result
            job.status = "cancelled" if job.cancel_requested else "done"
        except Exception as e:                              # noqa: BLE001
            if job.cancel_requested:
                job.status = "cancelled"
                job.error = "search cancelled; session was reset"
            else:
                job.error = str(e)
                job.status = "error"
        finally:
            job.done.set()


def make_handler(app: _App):
    class Handler(BaseHTTPRequestHandler):
        server_version = "clann-web/0.1"

        # -- helpers -------------------------------------------------------
        def _send(self, status: int, payload: dict) -> None:
            body = json.dumps(payload).encode()
            self.send_response(status)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def _read_json(self) -> dict:
            length = int(self.headers.get("Content-Length", 0) or 0)
            if length > MAX_JSON_BYTES:
                raise RequestTooLarge(
                    f"request body {length} bytes exceeds {MAX_JSON_BYTES}")
            if length == 0:
                return {}
            raw = self.rfile.read(length)
            try:
                return json.loads(raw or b"{}")
            except json.JSONDecodeError:
                return {}

        def log_message(self, *args):  # keep test output quiet
            pass

        def _send_file(self, relpath: str, content_type: str) -> None:
            full = os.path.normpath(os.path.join(STATIC_DIR, relpath))
            if not full.startswith(STATIC_DIR) or not os.path.isfile(full):
                self._send(404, {"ok": False, "error": "not found"})
                return
            with open(full, "rb") as fh:
                body = fh.read()
            self.send_response(200)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def _send_manual(self) -> None:
            """Serve USER_MANUAL.md (raw markdown) for the in-app help modal."""
            repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            path = os.path.join(repo, "USER_MANUAL.md")
            try:
                with open(path, "rb") as fh:
                    body = fh.read()
            except OSError:
                body = b"The user manual could not be found."
            self.send_response(200)
            self.send_header("Content-Type", "text/markdown; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def _download(self, raw_name: str) -> None:
            """Serve a file from the session sandbox as a download."""
            try:
                base = safe_basename(unquote(raw_name))   # rejects traversal
            except UnsafePath:
                self._send(400, {"ok": False, "error": "bad filename"})
                return
            full = os.path.join(app.sandbox.dir, base)
            if not os.path.isfile(full):
                self._send(404, {"ok": False, "error": "no such file"})
                return
            with open(full, "rb") as fh:
                body = fh.read()
            self.send_response(200)
            self.send_header("Content-Type", "application/octet-stream")
            self.send_header("Content-Disposition",
                             f'attachment; filename="{base}"')
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def _sandbox_text(self, raw_name: str):
            """Read a sandbox output file as text, or (None, error-status)."""
            try:
                base = safe_basename(unquote(raw_name))
            except UnsafePath:
                return None, 400
            full = os.path.join(app.sandbox.dir, base)
            if not os.path.isfile(full):
                return None, 404
            with open(full, "r", encoding="utf-8", errors="replace") as fh:
                return fh.read(), base

        def _view_file(self, raw_name: str) -> None:
            """Classify an output file: tree files go to the tree viewer, everything
            else to the text viewer. Returns {kind, [text]}."""
            text, base = self._sandbox_text(raw_name)
            if text is None:
                self._send(base, {"ok": False, "error": "no such file"})
                return
            if looks_like_tree_file(text) and parse_tree_file(text):
                self._send(200, {"kind": "tree", "name": base})
            else:
                self._send(200, {"kind": "text", "name": base, "text": text})

        def _view_tree(self, raw_name: str) -> None:
            """Serve the interactive viewer for a tree file from the sandbox."""
            text, base = self._sandbox_text(raw_name)
            if text is None:
                self._send(base, {"ok": False, "error": "no such file"})
                return
            trees = parse_tree_file(text)
            if not trees:
                self._send(415, {"ok": False, "error": "not a tree file"})
                return
            html = build_viewer_html(build_tree_view_document(trees, base))
            body = html.encode()
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        # -- routes --------------------------------------------------------
        def do_GET(self):
            p = urlparse(self.path).path
            if p in ("/", "/index.html"):
                self._send_file("index.html", "text/html; charset=utf-8")
            elif p == "/api/session":
                self._send(200, {"session_id": app.session_id,
                                 "state": app.state(),
                                 "files": app.sandbox.list(),
                                 "output_files": app.output_files(),
                                 "workdir": app.sandbox.dir})
            elif p.startswith("/api/download/"):
                self._download(p[len("/api/download/"):])
            elif p == "/api/session-commands":
                body = app.commands_script().encode()
                disp = ("attachment" if parse_qs(urlparse(self.path).query).get("download")
                        else "inline")
                self.send_response(200)
                self.send_header("Content-Type", "text/plain; charset=utf-8")
                self.send_header("Content-Disposition",
                                 f'{disp}; filename="clann_session_commands.txt"')
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
            elif p == "/api/viewfile":
                self._view_file(parse_qs(urlparse(self.path).query).get("name", [""])[0])
            elif p == "/api/viewtree":
                self._view_tree(parse_qs(urlparse(self.path).query).get("name", [""])[0])
            elif p == "/api/manual":
                self._send_manual()
            elif p == "/api/commands":
                self._send(200, {"commands": list_commands()})
            elif p.startswith("/api/commands/") and p.endswith("/schema"):
                name = p[len("/api/commands/"):-len("/schema")]
                self._send(200, {"name": name, **command_schema(name)})
            elif p == "/api/viewer":
                html = (build_viewer_html(app.last_result_json)
                        if app.last_result_json else placeholder_html())
                body = html.encode()
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
            elif p.startswith("/api/jobs/") and p.endswith("/stream"):
                job = app.jobs.get(p[len("/api/jobs/"):-len("/stream")])
                if job is None:
                    self._send(404, {"ok": False, "error": "no such job"})
                else:
                    self._stream_job(job)
            elif p.startswith("/api/jobs/"):
                job = app.jobs.get(p[len("/api/jobs/"):])
                if job is None:
                    self._send(404, {"ok": False, "error": "no such job"})
                else:
                    self._send(200, {"ok": True, **job.to_dict()})
            else:
                self._send(404, {"ok": False, "error": "not found"})

        def _stream_job(self, job) -> None:
            """Server-Sent Events: live log chunks, then a final 'done' event."""
            self.send_response(200)
            self.send_header("Content-Type", "text/event-stream")
            self.send_header("Cache-Control", "no-cache")
            self.send_header("Connection", "close")
            self.end_headers()
            sent = 0
            try:
                while True:
                    full = job.log()
                    if len(full) > sent:
                        chunk = full[sent:]
                        sent = len(full)
                        payload = json.dumps({"log": chunk})
                        self.wfile.write(f"data: {payload}\n\n".encode())
                        self.wfile.flush()
                    if job.done.is_set() and len(job.log()) == sent:
                        payload = json.dumps(job.to_dict(include_log=False))
                        self.wfile.write(f"event: done\ndata: {payload}\n\n".encode())
                        self.wfile.flush()
                        break
                    time.sleep(0.15)
            except (BrokenPipeError, ConnectionResetError):
                pass   # client disconnected

        def do_POST(self):
            path = urlparse(self.path).path
            try:
                if path == "/api/session":
                    self._send(200, app.new_session())

                elif path == "/api/files":
                    qs = parse_qs(urlparse(self.path).query)
                    name = (qs.get("name") or [""])[0]
                    if not name:
                        self._send(400, {"ok": False, "error": "missing ?name="})
                        return
                    length = int(self.headers.get("Content-Length", 0) or 0)
                    if length > MAX_UPLOAD_BYTES:
                        self._send(413, {"ok": False, "error": f"upload "
                                         f"{length} bytes exceeds "
                                         f"{MAX_UPLOAD_BYTES}"})
                        return
                    data = self.rfile.read(length) if length else b""
                    stored = app.sandbox.put(name, data)   # may raise UnsafePath
                    self._send(200, {"ok": True, "stored": stored,
                                     "files": app.sandbox.list()})

                elif path == "/api/load":
                    body = self._read_json()
                    fname = body.get("file")
                    if not fname:
                        self._send(400, {"ok": False, "error": "missing 'file'"})
                        return
                    if app.busy():
                        self._send(409, {"ok": False, "error": "a command is "
                                         "still running in this session — wait "
                                         "for it to finish (or Stop it) first"})
                        return
                    safe = app.sandbox.confine_input(fname)  # may raise UnsafePath
                    # optional `exe` options (name parsing, autoprune/decompose…);
                    # file paths within them are confined to the sandbox.
                    opts = sanitize_options(body.get("options", ""), app.sandbox)
                    log = app.engine.load(safe, opts)
                    app.record("exe " + fname + ((" " + opts) if opts else ""))
                    app.tree_counts = parse_tree_counts(log)  # single/multi copy
                    self._send(200, {"ok": True, "log": log,
                                     "state": app.state(),
                                     "files": app.sandbox.list(),
                                     "output_files": app.output_files()})

                elif path.startswith("/api/jobs/") and path.endswith("/cancel"):
                    jid = path[len("/api/jobs/"):-len("/cancel")]
                    job = app.cancel_current(jid)
                    if job is None:
                        self._send(404, {"ok": False,
                                         "error": "no running job to cancel"})
                    else:
                        self._send(202, {"ok": True, "job_id": job.id,
                                         "status": "cancelling"})

                elif path == "/api/run":
                    body = self._read_json()
                    cmd = body.get("command")
                    if not cmd:
                        self._send(400, {"ok": False, "error": "missing 'command'"})
                        return
                    safe_cmd = sanitize_command(cmd, app.sandbox)  # may raise
                    try:
                        job = app.start_job(safe_cmd)
                    except RuntimeError:
                        self._send(409, {"ok": False, "error": "a command is "
                                         "already running in this session"})
                        return
                    app.record(cmd)   # log the user's command for reproducibility
                    self._send(202, {"ok": True, "job_id": job.id,
                                     "command": safe_cmd, "status": "running"})

                else:
                    self._send(404, {"ok": False, "error": "not found"})
            except RequestTooLarge as e:
                self._send(413, {"ok": False, "error": str(e)})
            except UnsafePath as e:
                self._send(400, {"ok": False, "error": f"unsafe path: {e}"})
            except ClannError as e:
                self._send(500, {"ok": False, "error": str(e)})

    return Handler


def make_server(host: str = "127.0.0.1", port: int = 8765,
                engine: WorkerEngine | None = None) -> ThreadingHTTPServer:
    app = _App(engine)
    httpd = ThreadingHTTPServer((host, port), make_handler(app))
    httpd.clann_app = app  # expose for tests
    return httpd


def serve(host: str = "127.0.0.1", port: int = 8765) -> None:
    httpd = make_server(host, port)
    print(f"clann-web serving on http://{host}:{port}  (session "
          f"{httpd.clann_app.session_id})")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        try:
            httpd.clann_app.engine.terminate()
        except Exception:      # noqa: BLE001
            pass
        httpd.clann_app.sandbox.destroy()
        httpd.server_close()
