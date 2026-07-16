"""Step 4.1 security pass — the safe-defaults guarantees, as one suite.

Consolidates the properties that must hold before anyone binds beyond loopback
(PLAN_web_client.md §3.5, §3.4):

  * the server build has NO shell reachability — `libclann-server.so` carries no
    `system` reference, and a command that reaches Clann's shell choke point is
    refused at runtime instead of spawning anything;
  * every file path in a command is confined to the session sandbox — inputs must
    be uploaded (traversal rejected), outputs are relocated in;
  * oversized request bodies are rejected (413) before being read into memory;
  * the default bind address is loopback;
  * per-session temp directories are cleaned up.

Run (lib is x86_64 on Apple Silicon):
  PYTHONPATH=. arch -x86_64 /usr/bin/python3 clann_web/tests/test_security.py
"""

import glob
import inspect
import json
import os
import subprocess
import threading
import time
import urllib.error
import urllib.parse
import urllib.request

from clann_web.server import (make_server, serve,
                              MAX_UPLOAD_BYTES, MAX_JSON_BYTES)
from clann_web.engine import _default_lib_path

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TUTORIAL = os.path.join(REPO, "examples", "tutorial_multicopy.ph")


def _post(base, path, payload=None, raw=None):
    data = raw if raw is not None else json.dumps(payload or {}).encode()
    headers = {} if raw is not None else {"Content-Type": "application/json"}
    req = urllib.request.Request(base + path, data=data, headers=headers,
                                 method="POST")
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.status, json.loads(r.read())
    except urllib.error.HTTPError as e:
        try:
            return e.code, json.loads(e.read())
        except Exception:            # noqa: BLE001
            return e.code, {}


def _get(base, path):
    with urllib.request.urlopen(base + path, timeout=30) as r:
        return json.loads(r.read())


def _upload(base, name, data):
    req = urllib.request.Request(
        base + "/api/files?name=" + urllib.parse.quote(name, safe=""),
        data=data, method="POST")
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return r.status, json.loads(r.read())
    except urllib.error.HTTPError as e:
        return e.code, json.loads(e.read())


def _wait(base, jid, timeout=120):
    end = time.time() + timeout
    while time.time() < end:
        j = _get(base, f"/api/jobs/{jid}")
        if j["status"] in ("done", "error", "cancelled"):
            return j
        time.sleep(0.05)
    raise AssertionError("job did not finish")


# -- static guarantees (no server needed) ---------------------------------
def check_no_system_symbol():
    lib = _default_lib_path()
    out = subprocess.run(["nm", "-u", lib], capture_output=True, text=True)
    refs = [ln for ln in out.stdout.splitlines() if "_system" in ln
            and "system_" not in ln]
    assert not refs, f"server lib references system(): {refs}"
    print("  no system() symbol in libclann-server.so")


def check_loopback_default():
    sig = inspect.signature(serve)
    assert sig.parameters["host"].default == "127.0.0.1", sig
    ms = inspect.signature(make_server)
    assert ms.parameters["host"].default == "127.0.0.1", ms
    print("  default bind is 127.0.0.1")


# -- live server guarantees -----------------------------------------------
def check_shell_refused(base):
    """A command reaching Clann's shell choke point is refused, not run."""
    _post(base, "/api/session")
    with open(TUTORIAL, "rb") as fh:
        _upload(base, "t.ph", fh.read())
    _post(base, "/api/load", {"file": "t.ph"})
    # mrp criterion drives the parsimony step through PAUP* (clann_shell).
    _wait(base, _post(base, "/api/run", {"command": "set criterion=mrp"})[1]["job_id"])
    j = _wait(base, _post(base, "/api/run", {"command": "hs nreps=1"})[1]["job_id"])
    assert "Refused" in j.get("log", ""), (
        "expected shell refusal in log, got:\n" + j.get("log", "")[-400:])
    print("  shell/system path refused at runtime (mrp -> PAUP*)")


def check_path_confinement(base):
    _post(base, "/api/session")
    # traversal on upload
    code, _ = _upload(base, "../evil.ph", b"x")
    assert code == 400, ("upload traversal not rejected", code)
    assert not os.path.exists(os.path.join(REPO, "evil.ph"))
    # traversal on load
    code, _ = _post(base, "/api/load", {"file": "../../etc/hosts"})
    assert code == 400, ("load traversal not rejected", code)
    # load of a non-uploaded bare name is rejected (must upload first)
    code, _ = _post(base, "/api/load", {"file": "nope.ph"})
    assert code == 400, ("unknown input accepted", code)
    # output path escaping is relocated into the sandbox, original never created
    with open(TUTORIAL, "rb") as fh:
        _upload(base, "t.ph", fh.read())
    _post(base, "/api/load", {"file": "t.ph"})
    esc = "/tmp/clann_escape_probe.html"
    if os.path.exists(esc):
        os.remove(esc)
    code, r = _post(base, "/api/run", {"command": f"showtrees htmlview={esc}"})
    assert code == 202, (code, r)
    assert "htmlview=clann_escape_probe.html" in r["command"], r["command"]
    _wait(base, r["job_id"])
    assert not os.path.exists(esc), "output escaped the sandbox"
    print("  paths confined: inputs must be uploaded, outputs relocated in")


def check_size_limits(base):
    _post(base, "/api/session")
    # Oversized upload is rejected on its declared Content-Length, before the
    # (potentially huge) body is read — so we needn't actually send 64 MiB.
    req = urllib.request.Request(
        base + "/api/files?name=big.ph", method="POST", data=b"",
        headers={"Content-Length": str(MAX_UPLOAD_BYTES + 1)})
    try:
        with urllib.request.urlopen(req, timeout=30):
            rejected = False
    except urllib.error.HTTPError as e:
        rejected = (e.code == 413)
    assert rejected, "oversized upload (by Content-Length) not rejected 413"
    # oversized JSON body rejected
    code, _ = _post(base, "/api/run", {"command": "x" * (MAX_JSON_BYTES + 10)})
    assert code == 413, ("oversized JSON not rejected", code)
    print("  request-size caps enforced (413 on oversized upload + JSON)")


def _serve(base_host="127.0.0.1"):
    httpd = make_server(base_host, 0)
    port = httpd.server_address[1]
    base = f"http://{base_host}:{port}"
    t = threading.Thread(target=httpd.serve_forever, daemon=True)
    t.start()
    return httpd, base


def check_temp_cleanup():
    before = set(glob.glob("/tmp/clannweb-*"))
    httpd, base = _serve()
    sandbox_dir = httpd.clann_app.sandbox.dir
    assert os.path.isdir(sandbox_dir)
    try:
        _post(base, "/api/session")   # rotates the sandbox (old destroyed)
    finally:
        httpd.clann_app.engine.terminate()
        httpd.clann_app.sandbox.destroy()
        httpd.shutdown(); httpd.server_close()
    assert not os.path.exists(sandbox_dir), "initial sandbox not cleaned"
    after = set(glob.glob("/tmp/clannweb-*"))
    assert after <= before, f"leaked sandbox dirs: {after - before}"
    print("  per-session temp dirs cleaned up")


def run_checks():
    check_no_system_symbol()
    check_loopback_default()
    httpd, base = _serve()
    try:
        check_shell_refused(base)
        check_path_confinement(base)
        check_size_limits(base)
    finally:
        httpd.clann_app.engine.terminate()
        httpd.clann_app.sandbox.destroy()
        httpd.shutdown(); httpd.server_close()
    check_temp_cleanup()


def test_security_suite():
    run_checks()


if __name__ == "__main__":
    run_checks()
    print("PASS: no shell reachability, paths confined, request-size caps, "
          "loopback default, temp cleanup.")
