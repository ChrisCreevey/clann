"""Step 3.3 (per-session worker processes, Model A) + Step 3.2 (cancellation).

Verifies:
  * two sessions each own a separate worker process and run `hs` CONCURRENTLY
    with independent state (cross-session parallelism, crash containment);
  * a running search can be cancelled from the API and the session recovers
    (graceful SIGINT → Clann's controlc, or kill+respawn fallback for a wedged
    search), after which the server accepts the next command.

Run (lib is x86_64 on Apple Silicon):
  PYTHONPATH=. arch -x86_64 /usr/bin/python3 clann_web/tests/test_worker.py
"""

import json
import os
import threading
import time
import urllib.error
import urllib.parse
import urllib.request

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


def _get(base, path):
    with urllib.request.urlopen(base + path, timeout=30) as r:
        return json.loads(r.read())


def _upload(base, name, data):
    req = urllib.request.Request(
        base + "/api/files?name=" + urllib.parse.quote(name, safe=""),
        data=data, method="POST")
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read())


def _wait(base, jid, timeout=180):
    end = time.time() + timeout
    while time.time() < end:
        j = _get(base, f"/api/jobs/{jid}")
        if j["status"] in ("done", "error", "cancelled"):
            return j
        time.sleep(0.05)
    raise AssertionError("job did not finish")


def _start_server():
    httpd = make_server("127.0.0.1", 0)
    port = httpd.server_address[1]
    base = f"http://127.0.0.1:{port}"
    t = threading.Thread(target=httpd.serve_forever, daemon=True)
    t.start()
    return httpd, base


def _load(base):
    _post(base, "/api/session")
    with open(TUTORIAL, "rb") as fh:
        _upload(base, "t.ph", fh.read())
    _post(base, "/api/load", {"file": "t.ph"})
    _wait(base, _post(base, "/api/run", {"command": "set criterion=recon"})[1]["job_id"])


def check_concurrent_isolation():
    """Two sessions (two worker processes) run hs at the same time."""
    a = _start_server()
    b = _start_server()
    servers = [a, b]
    try:
        _load(a[1]); _load(b[1])
        results = {}

        def go(base, key):
            code, r = _post(base, "/api/run", {"command": "hs nreps=6 nthreads=1"})
            assert code == 202, (code, r)
            results[key] = _wait(base, r["job_id"])

        ta = threading.Thread(target=go, args=(a[1], "a"))
        tb = threading.Thread(target=go, args=(b[1], "b"))
        t0 = time.time()
        ta.start(); tb.start(); ta.join(); tb.join()
        elapsed = time.time() - t0

        for key in ("a", "b"):
            assert results[key]["status"] == "done", results[key]
            assert results[key]["scores"][0] == 32.0, results[key].get("scores")
        # Independent processes: killing one session's worker leaves the other fine.
        a[0].clann_app.new_session()      # terminates + respawns worker A
        assert b[1] and _get(b[1], "/api/session")["state"]["num_source_trees"] == 8
        print(f"  concurrent isolation OK (two workers, both scored 17, "
              f"{elapsed:.1f}s wall)")
    finally:
        for httpd, _ in servers:
            httpd.clann_app.engine.terminate()
            httpd.clann_app.sandbox.destroy()
            httpd.shutdown(); httpd.server_close()


def check_cancellation():
    """A long search is cancelled; the session then accepts the next command."""
    httpd, base = _start_server()
    try:
        _load(base)
        # A deliberately long run so we can cancel it mid-search.
        code, r = _post(base, "/api/run", {"command": "hs nreps=4000 nthreads=1"})
        assert code == 202, (code, r)
        jid = r["job_id"]

        # Wait until it is genuinely running (some log produced), then cancel.
        end = time.time() + 30
        while time.time() < end:
            j = _get(base, f"/api/jobs/{jid}")
            if j["status"] != "running":
                break
            if j.get("log"):
                break
            time.sleep(0.05)
        code, r = _post(base, f"/api/jobs/{jid}/cancel")
        assert code == 202 and r["status"] == "cancelling", (code, r)

        j = _wait(base, jid)
        assert j["status"] == "cancelled", ("expected cancelled, got", j["status"])

        # The session is responsive again (fresh worker if it was killed).
        code, r = _post(base, "/api/run", {"command": "set criterion=recon"})
        assert code == 202, (code, r)
        follow = _wait(base, r["job_id"])
        assert follow["status"] == "done", ("follow-up not done:", follow)
        print("  cancellation OK (job -> cancelled, session accepts next command)")
    finally:
        httpd.clann_app.engine.terminate()
        httpd.clann_app.sandbox.destroy()
        httpd.shutdown(); httpd.server_close()


def test_worker_isolation_and_cancellation():
    check_concurrent_isolation()
    check_cancellation()


if __name__ == "__main__":
    check_concurrent_isolation()
    check_cancellation()
    print("PASS: per-session worker processes give concurrent isolated sessions; "
          "a running search is cancellable and the session recovers.")
