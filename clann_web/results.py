"""Turn a Clann command's output into structured results (Step 1.3).

Tree-producing commands (hs, nj, showtrees, reconstruct, alltrees) are run with
an injected ``resultjson=`` so Clann writes the same JSON the HTML viewer embeds
(see NOTES_html_viewer.md §3). This module reads that file and shapes an API
result: each tree keeps its structured node form (what the browser viewer
consumes) plus a derived Newick string; scores come from the JSON for
reconciliations and from the log for supertree searches.
"""

from __future__ import annotations

import json
import os
import re

# The reserved sandbox filename the server injects as resultjson=.
RESULT_JSON = "__clann_result__.json"

# Commands that produce trees we can serialise.
TREE_COMMANDS = {"hs", "hsearch", "nj", "showtrees", "reconstruct", "alltrees",
                 "consensus", "bootstrap", "boot"}


def is_tree_command(command: str) -> bool:
    toks = command.split()
    return bool(toks) and toks[0] in TREE_COMMANDS


def node_to_newick(node: dict) -> str:
    """Recursively render a viewer JSON node as Newick."""
    children = node.get("children")
    if children:
        inner = ",".join(node_to_newick(c) for c in children)
        s = "(" + inner + ")"
        name = node.get("name")
        if name:  # internal label, if any
            s += _escape(name)
    else:
        # leaf: a normal tip has a name; a reconciliation loss node may not.
        s = _escape(node.get("name") or node.get("species") or "LOST")
    length = node.get("length")
    if length is not None:
        s += f":{length}"
    return s


def _escape(name: str) -> str:
    # Newick: quote names containing whitespace or punctuation that would break parsing.
    if re.search(r"[\s(),:;']", name):
        return "'" + name.replace("'", "''") + "'"
    return name


_SINGLE_COPY = re.compile(r"single copy trees:\s*(\d+)", re.IGNORECASE)
_MULTI_COPY = re.compile(r"multicopy trees:\s*(\d+)", re.IGNORECASE)


def parse_tree_counts(log: str) -> dict:
    """Extract single-/multi-copy tree counts Clann prints on `exe` (load).

    Returns {"num_single_copy": N, "num_multicopy": M} for whichever it found,
    or {} if the log doesn't carry those lines. Mirrors tree_io.c's
    "Number of single copy trees" / "number of multicopy trees" output.
    """
    out = {}
    s = _SINGLE_COPY.search(log)
    m = _MULTI_COPY.search(log)
    if s:
        out["num_single_copy"] = int(s.group(1))
    if m:
        out["num_multicopy"] = int(m.group(1))
    return out


_SUPERTREE_SCORE = re.compile(r"Supertree\s+\d+\s+of\s+\d+.*?=\s*([-\d.]+)")


def _scores_from_log(log: str) -> list:
    return [float(m.group(1)) for m in _SUPERTREE_SCORE.finditer(log)]


def build_results(json_path: str, log: str) -> dict:
    """Return {"trees": [...], "scores": [...]} from a resultjson file + log.

    Each tree: {name, newick, tree(structured), [score, dups, losses]}.
    """
    with open(json_path) as f:
        doc = json.load(f)

    trees = doc.get("trees", [])
    log_scores = _scores_from_log(log)

    out_trees = []
    out_scores = []
    for idx, t in enumerate(trees):
        entry = {
            "name": t.get("name"),
            "newick": node_to_newick(t["tree"]) + ";",
            "tree": t["tree"],
        }
        score = None
        if "score" in t:  # reconciliation: per-tree score/dups/losses in the JSON
            score = t["score"]
            entry["score"] = t["score"]
            entry["dups"] = t.get("dups")
            entry["losses"] = t.get("losses")
        elif idx < len(log_scores):  # supertree search: score parsed from the log
            score = log_scores[idx]
            entry["score"] = score
        out_trees.append(entry)
        out_scores.append(score)

    return {"type": doc.get("type"), "meta": doc.get("meta", {}),
            "trees": out_trees, "scores": out_scores}
