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

# Default NHX filename the server injects for `reconstruct` (a real output the
# user can download and re-open), unless they pass their own nhxfile=.
RECON_NHX = "reconstructions.nhx"

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
        # internal label: an explicit name, else a support value (bootstrap/consensus)
        label = node.get("name") or node.get("support")
        if label is not None and label != "":
            s += _escape(str(label))
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


# --- Newick -> viewer node (for viewing arbitrary tree files from outputs) -----

def _parse_newick_node(s: str, i: int):
    """Recursive-descent parse of one Newick/NHX clade starting at s[i].
    Returns (node, next_index) or raises ValueError. Produces the same
    {name|children, support, length, event, species} shape the browser viewer
    consumes; NHX `[&&NHX:S=..:D=Y/N]` tags and `*LOST` leaves are interpreted
    into event/species (see _apply_nhx)."""
    node: dict = {}
    n = len(s)
    comment = None
    if i < n and s[i] == "(":
        i += 1
        children = []
        while True:
            child, i = _parse_newick_node(s, i)
            children.append(child)
            if i < n and s[i] == ",":
                i += 1
                continue
            if i < n and s[i] == ")":
                i += 1
                break
            raise ValueError("malformed clade")
        node["children"] = children
        label, i = _read_label(s, i)      # internal label = support value
        if label:
            node["support"] = label
    else:
        name, i = _read_label(s, i)        # leaf name (may be quoted)
        if name == "":
            raise ValueError("empty leaf name")
        node["name"] = name
    # NHX `[...]` comments and `:branch_length` may appear in either order.
    while i < n:
        if s[i] == "[":
            body, i = _read_comment(s, i)
            comment = body if comment is None else comment + body
            continue
        if s[i] == ":":
            i += 1
            start = i
            while i < n and s[i] not in ",():;[":
                i += 1
            try:
                node["length"] = float(s[start:i])
            except ValueError:
                pass
            continue
        break
    _apply_nhx(node, comment)
    return node, i


def _read_comment(s: str, i: int):
    """Consume a `[...]` comment starting at s[i]=='['. Returns (inner, next)."""
    j = s.find("]", i)
    if j < 0:
        raise ValueError("unterminated '['")
    return s[i + 1:j], j + 1


def _read_label(s: str, i: int):
    n = len(s)
    if i < n and s[i] == "'":              # quoted name: '...''...'
        i += 1
        out = []
        while i < n:
            if s[i] == "'":
                if i + 1 < n and s[i + 1] == "'":
                    out.append("'")
                    i += 2
                    continue
                i += 1
                break
            out.append(s[i])
            i += 1
        return "".join(out), i
    start = i
    while i < n and s[i] not in ",():;[":   # stop before an NHX comment too
        i += 1
    return s[start:i].strip(), i


def _parse_nhx_tags(comment):
    """Parse an NHX comment body (`&&NHX:S=Foo:D=Y`) into a {key: value} map."""
    tags: dict = {}
    if not comment:
        return tags
    body = comment.strip()
    if body.startswith("&&NHX"):
        body = body[len("&&NHX"):]
    for part in body.split(":"):
        if "=" in part:
            k, v = part.split("=", 1)
            tags[k] = v
    return tags


def _apply_nhx(node: dict, comment):
    """Fold NHX annotations into a parsed node (event / species), matching Clann's
    conventions: D=Y duplication, D=N speciation, S= species, `*LOST` leaf loss."""
    name = node.get("name")
    is_loss = False
    if name and name.endswith("*LOST"):
        node["name"] = name[:-len("*LOST")]
        is_loss = True
    tags = _parse_nhx_tags(comment)
    if "S" in tags:
        node["species"] = tags["S"]
    if is_loss:
        node["event"] = "loss"
    elif "children" in node:
        d = tags.get("D")
        if d == "Y":
            node["event"] = "duplication"
        elif d == "N":
            node["event"] = "speciation"


_NAME_BRACKET = re.compile(r"\[([^\]]*)\]")


def _name_from_brackets(text: str):
    """The last non-numeric [..] annotation (a tree name; numeric ones are weights)."""
    name = None
    for b in _NAME_BRACKET.findall(text):
        b = b.strip()
        if b and not _is_number(b):
            name = b
    return name


_NHX_HEADER = re.compile(r"#\s*(\S+)[^\n]*?\(score=([-+0-9.eE]+)\)")


def parse_tree_file(text: str):
    """Parse a Clann/Phylip Newick or NHX tree file into [{name, tree, ...}], or
    None if the text doesn't look like Newick trees. Handles one or more
    ;-terminated trees. Clann writes annotations AFTER the ';'
    (`<tree>;[weight][name]`), so the leading annotations of each split segment
    name the PRECEDING tree. NHX files instead carry `# name (score=..)` header
    comment lines, which are captured here and dropped before parsing."""
    if "(" not in text:
        return None
    # NHX header comments name/score the trees in order; capture then strip them.
    headers = _NHX_HEADER.findall(text)
    body = "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("#"))
    trees = []
    for part in body.split(";"):
        paren = part.find("(")
        # leading [..] before this segment's '(' annotate the previous tree
        lead = part[:paren] if paren >= 0 else part
        if trees and trees[-1]["name"].startswith("tree_"):
            nm = _name_from_brackets(lead)
            if nm:
                trees[-1]["name"] = nm
        if paren < 0:
            continue
        try:
            node, _ = _parse_newick_node(part[paren:], 0)
        except (ValueError, IndexError):
            return None
        entry = {"name": f"tree_{len(trees) + 1}", "tree": node}
        dups, losses = _count_events(node)
        if dups or losses or _has_event(node):
            entry["dups"], entry["losses"] = dups, losses
        trees.append(entry)
    # apply NHX headers (name + score) positionally when present
    for entry, (name, score) in zip(trees, headers):
        entry["name"] = name
        try:
            entry["score"] = float(score)
        except ValueError:
            pass
    return trees or None


def _count_events(node: dict):
    """(duplications, losses) among this node and its descendants.

    Losses are counted per *maximal* fully-lost subtree, not per lost leaf: a
    whole clade lost together is a single loss event. A subtree is fully lost
    when every leaf at or below it is a loss; the loss is tallied on the top
    edge of that clade (i.e. once its parent is not itself fully lost)."""
    dups, losses, _ = _count_events_rec(node)
    # Degenerate case: the entire tree is a lost clade — count it once.
    if _fully_lost(node):
        losses += 1
    return dups, losses


def _fully_lost(node: dict) -> bool:
    children = node.get("children")
    if not children:
        return node.get("event") == "loss"
    return all(_fully_lost(c) for c in children)


def _count_events_rec(node: dict):
    """Returns (dups, losses, fully_lost) for the subtree at `node`."""
    dups = 1 if node.get("event") == "duplication" else 0
    losses = 0
    children = node.get("children", ())
    if not children:
        return dups, losses, node.get("event") == "loss"
    child_fully_lost = []
    for c in children:
        d, l, fl = _count_events_rec(c)
        dups += d
        losses += l
        child_fully_lost.append(fl)
    fully_lost = all(child_fully_lost)
    if not fully_lost:
        # Each fully-lost child under this surviving node is one loss event.
        losses += sum(1 for fl in child_fully_lost if fl)
    return dups, losses, fully_lost


def _has_event(node: dict) -> bool:
    if node.get("event"):
        return True
    return any(_has_event(c) for c in node.get("children", ()))


def build_tree_view_document(trees: list, title: str) -> str:
    """A viewer result-JSON document ({type, meta, trees}) for parsed trees.
    Emits `reconciliation` (so duplication/loss glyphs render) when any tree
    carries NHX events, otherwise a plain `tree`."""
    kind = "reconciliation" if any(_has_event(t["tree"]) for t in trees) else "tree"
    meta = {"title": title}
    if kind == "reconciliation":
        meta["criterion"] = "recon"
    return json.dumps({"type": kind, "meta": meta, "trees": trees})


def _is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def looks_like_tree_file(text: str) -> bool:
    """Cheap classification: the first non-blank, non-comment content starts a
    Newick tree. Used to route an output file to the tree viewer vs text viewer."""
    for line in text.splitlines():
        t = line.strip()
        if not t or t.startswith("#"):
            continue
        return t.startswith("(") or t.lower().startswith("tree ") and "(" in t
    return False


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
