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
import uuid
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs

from .engine import ClannEngine, ClannError
from .sandbox import Sandbox, UnsafePath, sanitize_command


class _App:
    """Holds the single engine + session id + sandbox (one session for now)."""

    def __init__(self, engine: ClannEngine):
        self.engine = engine
        self.sandbox = Sandbox()
        self.engine.set_workdir(self.sandbox.dir)
        self.session_id = uuid.uuid4().hex

    def new_session(self) -> dict:
        self.engine.reset()
        old = self.sandbox
        self.sandbox = Sandbox()
        self.engine.set_workdir(self.sandbox.dir)
        old.destroy()
        self.session_id = uuid.uuid4().hex
        return {"session_id": self.session_id, "state": self.engine.state(),
                "files": self.sandbox.list()}


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
            if length == 0:
                return {}
            raw = self.rfile.read(length)
            try:
                return json.loads(raw or b"{}")
            except json.JSONDecodeError:
                return {}

        def log_message(self, *args):  # keep test output quiet
            pass

        # -- routes --------------------------------------------------------
        def do_GET(self):
            if urlparse(self.path).path == "/api/session":
                self._send(200, {"session_id": app.session_id,
                                 "state": app.engine.state(),
                                 "files": app.sandbox.list()})
            else:
                self._send(404, {"ok": False, "error": "not found"})

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
                    safe = app.sandbox.confine_input(fname)  # may raise UnsafePath
                    log = app.engine.load(safe)
                    self._send(200, {"ok": True, "log": log,
                                     "state": app.engine.state()})

                elif path == "/api/run":
                    body = self._read_json()
                    cmd = body.get("command")
                    if not cmd:
                        self._send(400, {"ok": False, "error": "missing 'command'"})
                        return
                    safe_cmd = sanitize_command(cmd, app.sandbox)  # may raise
                    log = app.engine.run(safe_cmd)
                    self._send(200, {"ok": True, "log": log, "command": safe_cmd,
                                     "state": app.engine.state()})

                else:
                    self._send(404, {"ok": False, "error": "not found"})
            except UnsafePath as e:
                self._send(400, {"ok": False, "error": f"unsafe path: {e}"})
            except ClannError as e:
                self._send(500, {"ok": False, "error": str(e)})

    return Handler


def make_server(host: str = "127.0.0.1", port: int = 8765,
                engine: ClannEngine | None = None) -> ThreadingHTTPServer:
    app = _App(engine or ClannEngine())
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
        httpd.clann_app.sandbox.destroy()
        httpd.server_close()
