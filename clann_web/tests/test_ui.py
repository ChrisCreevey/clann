"""Step 2.1 checks: SPA shell + command palette endpoints.

Run (lib is x86_64 on Apple Silicon):
  PYTHONPATH=. arch -x86_64 /usr/bin/python3 clann_web/tests/test_ui.py
"""

import json
import threading
import urllib.request

from clann_web.server import make_server


def _get(base, path):
    with urllib.request.urlopen(base + path, timeout=30) as r:
        return r.status, r.headers.get("Content-Type", ""), r.read()


def run_checks():
    httpd = make_server("127.0.0.1", 0)
    host, port = httpd.server_address
    base = f"http://127.0.0.1:{port}"
    t = threading.Thread(target=httpd.serve_forever, daemon=True)
    t.start()
    try:
        # SPA served at /
        code, ctype, body = _get(base, "/")
        assert code == 200 and "text/html" in ctype, (code, ctype)
        html = body.decode()
        assert "<title>Clann" in html and 'id="cmd"' in html, "SPA shell missing"

        # command palette endpoint
        code, ctype, body = _get(base, "/api/commands")
        assert code == 200 and "application/json" in ctype, (code, ctype)
        names = [c["name"] for c in json.loads(body)["commands"]]
        for expected in ("hs", "nj", "showtrees", "reconstruct"):
            assert expected in names, (expected, names)
        assert all("blurb" in c for c in json.loads(body)["commands"])

        # fresh session shows an empty state
        code, ctype, body = _get(base, "/api/session")
        st = json.loads(body)["state"]
        assert st["input_file"] is None and st["num_source_trees"] == 0, st

        # viewer endpoint: placeholder before any result, real viewer after a run
        code, ctype, body = _get(base, "/api/viewer")
        assert code == 200 and "text/html" in ctype, (code, ctype)
        assert b"Run a tree command" in body, "expected placeholder"

        app = httpd.clann_app
        app.last_result_json = ('{"type":"tree","meta":{"dataset":"x"},'
                                '"trees":[{"name":"S1","tree":{"children":'
                                '[{"name":"A"},{"name":"B"}]}}]}')
        code, ctype, body = _get(base, "/api/viewer")
        html = body.decode()
        assert "const DATA =" in html and '"S1"' in html, "viewer not wrapped"
        assert html.count("/*CLANN_DATA_BEGIN*/") == 1, "template markers off"
    finally:
        httpd.clann_app.sandbox.destroy()
        httpd.shutdown()
        httpd.server_close()


def test_ui_shell():
    run_checks()


if __name__ == "__main__":
    run_checks()
    print("PASS: SPA served at /, /api/commands lists hs/nj/showtrees/reconstruct, "
          "session starts empty.")
