"""Step 2.2 checks: command schema endpoint + schema/opts_ consistency.

The consistency test guards against the form advertising an option the C command
does not actually parse (the class of bug behind `hs seed=`): every option in the
schema for a command must appear in that command's opts_<cmd>[] array in main.c.

Run (lib is x86_64 on Apple Silicon):
  PYTHONPATH=. arch -x86_64 /usr/bin/python3 clann_web/tests/test_schema.py
"""

import json
import os
import re
import threading
import urllib.request

from clann_web.commands import command_schema, _load_schema
from clann_web.server import make_server

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MAIN_C = os.path.join(REPO, "main.c")

# schema command name -> the opts_ array name in main.c
_OPTS_ARRAY = {
    "set": "opts_set", "hs": "opts_hs", "nj": "opts_nj",
    "reconstruct": "opts_reconstruct", "showtrees": "opts_showtrees",
}


def _opts_array(name: str) -> set:
    """Extract the option names (without '=') from a static opts_ array in main.c."""
    src = open(MAIN_C, encoding="utf-8").read()
    m = re.search(r"static const char \*" + re.escape(name) + r"\[\]\s*=\s*\{(.*?)\}",
                  src, re.DOTALL)
    assert m, f"{name} not found in main.c"
    body = m.group(1)
    return {tok.rstrip("=") for tok in re.findall(r'"([^"]+)"', body)}


def test_schema_matches_opts():
    schema = _load_schema()
    for cmd, arr in _OPTS_ARRAY.items():
        allowed = _opts_array(arr)
        for opt in command_schema(cmd)["options"]:
            assert opt["name"] in allowed, (
                f"schema option {cmd}.{opt['name']} is not in {arr} "
                f"(the C command does not parse it)")


def test_schema_shapes():
    # every option is well-formed; enums carry values
    for cmd, entry in _load_schema().items():
        if cmd.startswith("_"):
            continue
        for opt in entry.get("options", []):
            assert {"name", "type", "help"} <= opt.keys(), opt
            if opt["type"] == "enum":
                assert opt.get("values"), opt


def test_schema_endpoint():
    httpd = make_server("127.0.0.1", 0)
    port = httpd.server_address[1]
    base = f"http://127.0.0.1:{port}"
    t = threading.Thread(target=httpd.serve_forever, daemon=True)
    t.start()
    try:
        with urllib.request.urlopen(base + "/api/commands/hs/schema") as r:
            data = json.loads(r.read())
        names = [o["name"] for o in data["options"]]
        assert "nreps" in names and "lossmodel" in names, names
        # hs must NOT advertise seed/criterion (it does not parse them)
        assert "seed" not in names and "criterion" not in names, names
        # a command without a schema returns an empty option list
        with urllib.request.urlopen(base + "/api/commands/rfdists/schema") as r:
            assert json.loads(r.read()).get("options", []) == []
    finally:
        httpd.clann_app.sandbox.destroy()
        httpd.shutdown()
        httpd.server_close()


if __name__ == "__main__":
    test_schema_matches_opts()
    test_schema_shapes()
    test_schema_endpoint()
    print("PASS: schema options all valid for their C command; hs omits seed/"
          "criterion; endpoint serves per-command schema.")
