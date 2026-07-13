"""End-to-end test for the Step 1.1 skeleton server.

Starts the real server (loading libclann-server.so) on an ephemeral loopback
port in a background thread and drives a full persistent session over HTTP:
  new session -> load tutorial -> set criterion=recon -> hs -> reconstruct
asserting the recon optimum (17) and that session state persists across calls.

Run (Apple Silicon: the lib is x86_64, so use a matching interpreter):
  arch -x86_64 /usr/bin/python3 -m pytest clann_web/tests/test_server.py -q
or without pytest:
  arch -x86_64 /usr/bin/python3 clann_web/tests/test_server.py
"""

import json
import os
import threading
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
            return json.loads(r.read())
    except urllib.error.HTTPError as e:
        # 4xx/5xx still carry a JSON body we want to inspect.
        return json.loads(e.read())


def _upload(base, name, data: bytes):
    req = urllib.request.Request(
        base + "/api/files?name=" + urllib.parse.quote(name), data=data,
        headers={"Content-Type": "application/octet-stream"}, method="POST")
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return json.loads(r.read())
    except urllib.error.HTTPError as e:
        return json.loads(e.read())


def _get(base, path):
    with urllib.request.urlopen(base + path, timeout=30) as r:
        return json.loads(r.read())


def run_checks():
    httpd = make_server("127.0.0.1", 0)  # port 0 = ephemeral
    host, port = httpd.server_address
    base = f"http://127.0.0.1:{port}"
    t = threading.Thread(target=httpd.serve_forever, daemon=True)
    t.start()
    try:
        # fresh session
        s = _post(base, "/api/session")
        assert "session_id" in s, s
        assert s["state"]["num_source_trees"] == 0, s

        # upload the tutorial into the session sandbox, then load it by name
        with open(TUTORIAL, "rb") as fh:
            up = _upload(base, "tutorial_multicopy.ph", fh.read())
        assert up["ok"] and "tutorial_multicopy.ph" in up["files"], up
        r = _post(base, "/api/load", {"file": "tutorial_multicopy.ph"})
        assert r["ok"], r
        assert r["state"]["num_source_trees"] == 8, r["state"]
        assert r["state"]["num_taxa"] == 9, r["state"]

        # persistent session: set criterion, then hs, then reconstruct
        _post(base, "/api/run", {"command": "set seed=42"})
        c = _post(base, "/api/run", {"command": "set criterion=recon"})
        assert c["state"]["criterion"] == "recon", c["state"]

        h = _post(base, "/api/run", {"command": "hs nreps=5 nthreads=1"})
        assert h["ok"], h
        assert "17.000000" in h["log"], h["log"][-500:]
        assert h["state"]["trees_in_memory"] >= 1, h["state"]
        # structured results (Step 1.3): trees + scores returned, not just log
        assert h.get("result_type") == "tree", h.get("result_type")
        assert h["trees"], h
        t0 = h["trees"][0]
        assert t0["newick"].endswith(";") and "Human" in t0["newick"], t0["newick"]
        assert t0["tree"]["children"], t0  # structured node form for the viewer
        assert h["scores"][0] == 17.0, h["scores"]

        # reconstruct uses the tree hs just left in memory (state persisted)
        rec = _post(base, "/api/run",
                    {"command": "reconstruct speciestree=memory open=no"})
        assert rec["ok"], rec
        assert "17.0000" in rec["log"], rec["log"][-500:]
        # reconciliation results carry per-tree events + dup/loss counts
        assert rec.get("result_type") == "reconciliation", rec.get("result_type")
        assert len(rec["trees"]) == 8, len(rec["trees"])
        assert all("score" in t for t in rec["trees"]), rec["trees"][0]
        # at least one duplication event somewhere in the reconciliations
        blob = json.dumps(rec["trees"])
        assert '"event": "duplication"' in blob or '"event":"duplication"' in blob, "no dup events"

        # GET session reflects the same live state
        g = _get(base, "/api/session")
        assert g["state"]["num_source_trees"] == 8, g

        # error handling: missing field
        bad = _post(base, "/api/run", {})
        assert bad.get("ok") is False, bad
    finally:
        httpd.clann_app.sandbox.destroy()
        httpd.shutdown()
        httpd.server_close()


def test_server_end_to_end():
    run_checks()


if __name__ == "__main__":
    run_checks()
    print("PASS: HTTP session load->set->hs->reconstruct; recon 17 over the wire.")
