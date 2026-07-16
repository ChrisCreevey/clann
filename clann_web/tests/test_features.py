"""Checks for the session-panel counts, exe options, and file download.

Covers the additions:
  * `exe` (load) reports single-/multi-copy tree counts in the session state;
  * the `exe` command has a typed option schema served to the form;
  * output files in the sandbox can be listed and downloaded, and download
    rejects path traversal.

Run (lib is x86_64 on Apple Silicon):
  PYTHONPATH=. arch -x86_64 /usr/bin/python3 clann_web/tests/test_features.py
"""

import json
import os
import threading
import urllib.error
import urllib.parse
import urllib.request

from clann_web.server import make_server
from clann_web.results import parse_tree_counts

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TUTORIAL = os.path.join(REPO, "examples", "tutorial_multicopy.ph")


def _post(base, path, payload=None):
    req = urllib.request.Request(base + path, data=json.dumps(payload or {}).encode(),
                                 headers={"Content-Type": "application/json"},
                                 method="POST")
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.status, json.loads(r.read())
    except urllib.error.HTTPError as e:
        return e.code, json.loads(e.read())


def _get(base, path):
    with urllib.request.urlopen(base + path, timeout=30) as r:
        return r.status, r.read()


def _upload(base, name, data):
    req = urllib.request.Request(
        base + "/api/files?name=" + urllib.parse.quote(name, safe=""),
        data=data, method="POST")
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read())


def test_parse_tree_counts_unit():
    log = ("\tNumber of single copy trees:\t3\n"
           "\tnumber of multicopy trees:\t5\n")
    assert parse_tree_counts(log) == {"num_single_copy": 3, "num_multicopy": 5}
    assert parse_tree_counts("nothing here") == {}


def run_checks():
    test_parse_tree_counts_unit()

    httpd = make_server("127.0.0.1", 0)
    base = f"http://127.0.0.1:{httpd.server_address[1]}"
    t = threading.Thread(target=httpd.serve_forever, daemon=True)
    t.start()
    try:
        _post(base, "/api/session")
        with open(TUTORIAL, "rb") as fh:
            _upload(base, "t.ph", fh.read())

        # exe has a served option schema (for the form)
        _, sch = _get(base, "/api/commands/exe/schema")
        schema = json.loads(sch)
        names = {o["name"] for o in schema["options"]}
        assert {"autoprunemono", "autodecompose", "maxnamelen"} <= names, names

        # load reports single-/multi-copy counts in the state
        code, r = _post(base, "/api/load", {"file": "t.ph"})
        assert code == 200 and r["ok"], r
        st = r["state"]
        assert st["num_single_copy"] + st["num_multicopy"] == st["num_source_trees"], st
        assert st["num_multicopy"] >= 1, st  # tutorial_multicopy has multicopy trees
        assert "t.ph" in r["output_files"]

        # download a sandbox file, and reject traversal
        code, body = _get(base, "/api/download/t.ph")
        assert code == 200 and body.startswith(open(TUTORIAL, "rb").read()[:10]), code
        try:
            _get(base, "/api/download/" + urllib.parse.quote("../../etc/hosts", safe=""))
            raise AssertionError("traversal download not rejected")
        except urllib.error.HTTPError as e:
            assert e.code in (400, 404), e.code

        print("  exe schema served; load reports single/multi counts; "
              "download works and rejects traversal")
    finally:
        httpd.clann_app.engine.terminate()
        httpd.clann_app.sandbox.destroy()
        httpd.shutdown(); httpd.server_close()


def test_features():
    run_checks()


if __name__ == "__main__":
    run_checks()
    print("PASS: single/multi-copy counts, exe options, and file download.")
