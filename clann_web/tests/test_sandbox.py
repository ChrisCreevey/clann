"""Step 1.2 security checks: per-session sandbox + path confinement.

Two layers are tested: the pure sanitizer (fast, no engine) and the live server
(upload rejection, output-path relocation on a real command).

Run (lib is x86_64 on Apple Silicon):
  PYTHONPATH=. arch -x86_64 /usr/bin/python3 clann_web/tests/test_sandbox.py
"""

import json
import os
import threading
import urllib.parse
import urllib.request

from clann_web.sandbox import Sandbox, UnsafePath, sanitize_command
from clann_web.server import make_server

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TUTORIAL = os.path.join(REPO, "examples", "tutorial_multicopy.ph")


def _post(base, path, payload=None):
    data = json.dumps(payload or {}).encode()
    req = urllib.request.Request(base + path, data=data,
                                 headers={"Content-Type": "application/json"},
                                 method="POST")
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            return r.status, json.loads(r.read())
    except urllib.error.HTTPError as e:
        return e.code, json.loads(e.read())


def _upload(base, name, data):
    req = urllib.request.Request(
        base + "/api/files?name=" + urllib.parse.quote(name, safe=""),
        data=data, headers={"Content-Type": "application/octet-stream"},
        method="POST")
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return r.status, json.loads(r.read())
    except urllib.error.HTTPError as e:
        return e.code, json.loads(e.read())


def _run(base, command):
    """POST /api/run (async) and poll the job to completion; returns the job."""
    import time
    _, r = _post(base, "/api/run", {"command": command})
    jid = r["job_id"]
    for _ in range(600):
        with urllib.request.urlopen(base + f"/api/jobs/{jid}") as resp:
            j = json.loads(resp.read())
        if j.get("status") in ("done", "error"):
            return j
        time.sleep(0.1)
    raise AssertionError("job did not finish")


def test_sanitizer_unit():
    sb = Sandbox()
    try:
        sb.put("trees.ph", b"();")

        # output path escaping -> relocated to a sandbox basename
        assert sanitize_command("hs savetrees=/etc/passwd", sb) == "hs savetrees=passwd"
        assert sanitize_command("hs savetrees=../../x.ph", sb) == "hs savetrees=x.ph"
        # htmlview=yes is a keyword, left alone; a path value is relocated
        assert sanitize_command("reconstruct htmlview=yes", sb) == "reconstruct htmlview=yes"
        assert sanitize_command("reconstruct htmlview=/tmp/a.html", sb) == "reconstruct htmlview=a.html"
        # sentinels untouched
        assert sanitize_command("reconstruct speciestree=memory", sb) == "reconstruct speciestree=memory"
        assert sanitize_command("hs nthreads=4 nreps=5", sb) == "hs nthreads=4 nreps=5"

        # input path must be an uploaded bare filename
        assert sanitize_command("exe trees.ph", sb) == "exe trees.ph"
        for bad in ("exe /etc/passwd", "exe ../secret.ph", "exe sub/dir.ph",
                    "usertrees file=/etc/passwd"):
            try:
                sanitize_command(bad, sb)
                raise AssertionError(f"expected UnsafePath for: {bad}")
            except UnsafePath:
                pass
    finally:
        sb.destroy()


def test_server_sandbox():
    httpd = make_server("127.0.0.1", 0)
    host, port = httpd.server_address
    base = f"http://127.0.0.1:{port}"
    sandbox_dir = httpd.clann_app.sandbox.dir
    t = threading.Thread(target=httpd.serve_forever, daemon=True)
    t.start()
    try:
        # a traversal filename on upload is rejected
        code, r = _upload(base, "../evil.ph", b"();")
        assert code == 400 and r["ok"] is False, (code, r)
        assert not os.path.exists(os.path.join(REPO, "evil.ph")), "escaped the sandbox!"

        # a normal upload lands in the sandbox
        with open(TUTORIAL, "rb") as fh:
            code, r = _upload(base, "trees.ph", fh.read())
        assert code == 200 and "trees.ph" in r["files"], (code, r)
        assert os.path.isfile(os.path.join(sandbox_dir, "trees.ph"))

        # loading a not-uploaded / traversal name is rejected
        code, r = _post(base, "/api/load", {"file": "../../etc/hosts"})
        assert code == 400 and r["ok"] is False, (code, r)

        _post(base, "/api/load", {"file": "trees.ph"})

        # an output path that tries to escape is relocated into the sandbox
        _run(base, "set criterion=dfit")
        r = _run(base, "showtrees display=no htmlview=/tmp/escape.html")
        assert r["status"] == "done", r
        assert r["command"].endswith("htmlview=escape.html"), r["command"]
        assert not os.path.exists("/tmp/escape.html"), "output escaped the sandbox!"
        assert os.path.isfile(os.path.join(sandbox_dir, "escape.html")), r
    finally:
        httpd.clann_app.sandbox.destroy()
        httpd.shutdown()
        httpd.server_close()


if __name__ == "__main__":
    test_sanitizer_unit()
    test_server_sandbox()
    print("PASS: sandbox confinement — uploads validated, inputs must be uploaded, "
          "output paths relocated into the session dir.")
