"""Build the interactive viewer page for the web client (Step 2.4).

Reuses the exact same self-contained viewer Clann embeds in its htmlview files
(`tools/clannview.template.html`). That template carries its data as

    const DATA = /*CLANN_DATA_BEGIN*/ { … } /*CLANN_DATA_END*/;

so we split it at those markers once and inject the JSON document Clann wrote for
the latest result (the same `{type, meta, trees:[…]}` a `resultjson=` run
produces). No duplicated rendering code — the browser viewer is identical to the
downloadable one.
"""

from __future__ import annotations

import os

_BEGIN = "/*CLANN_DATA_BEGIN*/"
_END = "/*CLANN_DATA_END*/"
_parts: tuple[str, str] | None = None


def _template_parts() -> tuple[str, str]:
    global _parts
    if _parts is None:
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        path = os.path.join(repo, "tools", "clannview.template.html")
        with open(path, encoding="utf-8") as f:
            tpl = f.read()
        b = tpl.index(_BEGIN) + len(_BEGIN)
        e = tpl.index(_END)
        _parts = (tpl[:b], tpl[e:])   # head ends after BEGIN, tail starts at END
    return _parts


def build_viewer_html(json_text: str) -> str:
    """Wrap a Clann result JSON document in the viewer template."""
    head, tail = _template_parts()
    return head + "\n" + json_text.strip() + "\n" + tail


def placeholder_html() -> str:
    return (
        "<!doctype html><meta charset='utf-8'>"
        "<body style='margin:0;font:14px -apple-system,sans-serif;color:#8a9099;"
        "display:flex;align-items:center;justify-content:center;height:100vh;"
        "background:transparent'>Run a tree command to see the interactive viewer."
        "</body>"
    )
