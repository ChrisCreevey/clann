"""Step 3.1 checks: async run jobs + SSE live log streaming.

Run (lib is x86_64 on Apple Silicon):
  PYTHONPATH=. arch -x86_64 /usr/bin/python3 clann_web/tests/test_jobs.py
"""

import json
import os
import threading
import time
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
    import urllib.parse
    req = urllib.request.Request(
        base + "/api/files?name=" + urllib.parse.quote(name, safe=""),
        data=data, method="POST")
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read())


def _wait(base, jid, timeout=120):
    end = time.time() + timeout
    while time.time() < end:
        j = _get(base, f"/api/jobs/{jid}")
        if j["status"] in ("done", "error"):
            return j
        time.sleep(0.05)
    raise AssertionError("job did not finish")


def run_checks():
    httpd = make_server("127.0.0.1", 0)
    port = httpd.server_address[1]
    base = f"http://127.0.0.1:{port}"
    t = threading.Thread(target=httpd.serve_forever, daemon=True)
    t.start()
    try:
        _post(base, "/api/session")
        with open(TUTORIAL, "rb") as fh:
            _upload(base, "t.ph", fh.read())
        _post(base, "/api/load", {"file": "t.ph"})

        # /api/run returns immediately with a job id (202, status running)
        code, r = _post(base, "/api/run", {"command": "set seed=42"})
        assert code == 202 and r["job_id"] and r["status"] == "running", (code, r)
        _wait(base, r["job_id"])
        _wait(base, _post(base, "/api/run", {"command": "set criterion=recon"})[1]["job_id"])

        # stream an hs run over SSE; collect log chunks + the done event.
        # nthreads=1 and a modest nreps: single-threaded hs is stable, whereas
        # multithreaded and high-nreps-on-repeat hs have a known intermittent hang
        # (a pre-existing core bug, see PLAN Step 3.1). Still streams per-rep lines.
        code, r = _post(base, "/api/run", {"command": "hs nreps=8 nthreads=1"})
        assert code == 202, (code, r)
        jid = r["job_id"]

        chunks, done_payload = [], None
        with urllib.request.urlopen(base + f"/api/jobs/{jid}/stream", timeout=120) as s:
            event = None
            for raw in s:
                line = raw.decode().rstrip("\n")
                if line.startswith("event:"):
                    event = line.split(":", 1)[1].strip()
                elif line.startswith("data:"):
                    data = json.loads(line.split(":", 1)[1].strip())
                    if event == "done":
                        done_payload = data
                        break
                    if "log" in data:
                        chunks.append(data["log"])
                elif line == "":
                    event = None

        full = "".join(chunks)
        assert "17.000000" in full, "streamed log missing the score:\n" + full[-400:]
        assert done_payload and done_payload["status"] == "done", done_payload
        assert done_payload["scores"][0] == 17.0, done_payload.get("scores")
        assert done_payload.get("has_viewer") is True, done_payload

        # after completion the job is still queryable and carries the full log
        j = _get(base, f"/api/jobs/{jid}")
        assert j["status"] == "done" and "17.000000" in j["log"], j["status"]

        # a second run while one is active would 409 — but this one is done, so ok
        code, r = _post(base, "/api/run", {"command": "showtrees display=no"})
        assert code == 202, (code, r)
        _wait(base, r["job_id"])
    finally:
        httpd.clann_app.sandbox.destroy()
        httpd.shutdown()
        httpd.server_close()


def test_async_jobs_and_streaming():
    run_checks()


if __name__ == "__main__":
    run_checks()
    print("PASS: /api/run is async (202 + job_id); SSE streams live log chunks "
          "and a done event with structured results; job stays queryable.")
