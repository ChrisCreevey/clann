"""Pure-parser checks for the output-file tree view path (no Clann worker).

Regression cover for viewing a `reconstruct` NHX file: the Python parser used by
the "view" button must interpret NHX `[&&NHX:...]` tags and `*LOST` leaves, not
leave them embedded in tip names.
"""
import json

from clann_web.results import parse_tree_file, build_tree_view_document

NHX = (
    "# tree_0  (score=2.0000)\n"
    "(((Human.a[&&NHX:S=Human:D=N],Human.b[&&NHX:S=Human:D=N])[&&NHX:D=Y:S=Human],"
    "Gorilla*LOST[&&NHX:S=Gorilla])[&&NHX:D=N],Mouse[&&NHX:S=Mouse:D=N]);\n"
)


def _first_leaf(node):
    return node if "children" not in node else _first_leaf(node["children"][0])


def test_nhx_leaf_names_have_no_annotation():
    trees = parse_tree_file(NHX)
    assert trees is not None and len(trees) == 1
    leaf = _first_leaf(trees[0]["tree"])
    assert leaf["name"] == "Human.a"
    assert leaf["species"] == "Human"
    assert "event" not in leaf  # gene-copy leaves carry no event


def test_nhx_header_name_and_score():
    t = parse_tree_file(NHX)[0]
    assert t["name"] == "tree_0"
    assert t["score"] == 2.0


def test_nhx_events_and_counts():
    t = parse_tree_file(NHX)[0]
    assert (t["dups"], t["losses"]) == (1, 1)
    # Clann's outermost paren is unannotated; the D=N speciation node is its child.
    spec = t["tree"]["children"][0]
    assert spec["event"] == "speciation"           # [&&NHX:D=N]
    dup = spec["children"][0]
    assert dup["event"] == "duplication" and dup["species"] == "Human"
    loss = spec["children"][1]
    assert loss["event"] == "loss" and loss["name"] == "Gorilla"


def test_nhx_document_is_reconciliation():
    doc = json.loads(build_tree_view_document(parse_tree_file(NHX), "reconstructions.nhx"))
    assert doc["type"] == "reconciliation"
    assert doc["meta"]["criterion"] == "recon"


def test_plain_newick_stays_tree():
    trees = parse_tree_file("((A:0.1,B:0.2)95:0.3,(C:0.4,D:0.5)80:0.6);")
    doc = json.loads(build_tree_view_document(trees, "t"))
    assert doc["type"] == "tree"
    node = trees[0]["tree"]["children"][0]
    assert node["support"] == "95" and node["length"] == 0.3
