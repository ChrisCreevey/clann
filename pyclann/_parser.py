"""Output file and stdout parsers for pyclann."""

from __future__ import annotations

import re

__all__ = [
    "parse_scored_tree_file",
    "parse_plain_tree_file",
    "parse_bootstrap_file",
    "extract_stdout_score",
]

# Matches lines like:   ((A,B),C);   [7.297618]
# or                    ((A,B),C)    [7.297618]
# (with optional semicolon before the tab/space and score in brackets)
_SCORED_LINE_RE = re.compile(
    r"^\s*(?P<newick>\([^[]+?)\s*;?\s*\[(?P<score>[0-9eE.+\-]+)\]",
)

# Matches lines like:   ((A,B),C) [0.91];   [score = 7.297618]
# (bootstrap file format: bootstrap-support in first brackets, score in second)
_BOOTSTRAP_LINE_RE = re.compile(
    r"^\s*(?P<newick>\([^[]+?)\s*\[[0-9eE.+\-]+\]\s*;?\s*\[score\s*=\s*(?P<score>[0-9eE.+\-]+)\]",
)

# Matches: "Supertree N of M score = 7.297618"
# or:      "Input tree N  score = 7.297618"
_STDOUT_SCORE_RE = re.compile(
    r"(?:Supertree\s+\d+\s+of\s+\d+|Input\s+tree\s+\d+)\s+"
    r"(?:lnL\s*=\s*|score\s*=\s*)(?P<score>[0-9eE.+\-]+)",
)


def parse_scored_tree_file(text: str) -> list[tuple[str, float]]:
    """Parse a tab-separated Newick + score file (``Heuristic_result.txt``, etc.).

    Each line has the format::

        ((A,B),C);    [7.297618]

    Returns a list of ``(newick_string, score)`` tuples in file order.
    Blank lines and lines not matching the pattern are skipped.
    """
    results: list[tuple[str, float]] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        m = _SCORED_LINE_RE.match(line)
        if m:
            newick = m.group("newick").rstrip(";").strip()
            score = float(m.group("score"))
            results.append((newick, score))
    return results


def parse_plain_tree_file(text: str) -> list[str]:
    """Parse a plain Newick file (one tree per line, no score).

    Used for ``NJ-tree.ph`` and ``consensus.ph``.
    Empty lines are skipped.  Trailing semicolons are preserved.
    """
    trees: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if line and line.startswith("("):
            trees.append(line)
    return trees


def parse_bootstrap_file(text: str) -> list[tuple[str, float]]:
    """Parse a bootstrap results file (``bootstrap.txt``).

    Each line has the format::

        ((A,B),C) [1.000000];    [score = 5.480159]

    Returns a list of ``(newick_string, score)`` tuples.
    The bootstrap-support value embedded in the Newick is preserved as-is.
    """
    results: list[tuple[str, float]] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        m = _BOOTSTRAP_LINE_RE.match(line)
        if m:
            newick = m.group("newick").rstrip(";").strip()
            score = float(m.group("score"))
            results.append((newick, score))
    return results


def extract_stdout_score(stdout: str) -> float | None:
    """Extract the first score printed to stdout by clann.

    Matches lines such as::

        Supertree 1 of 1 score = 7.297618
        Input tree 1  score = 7.297618

    Returns the score as a float, or ``None`` if no match is found.
    """
    m = _STDOUT_SCORE_RE.search(stdout)
    return float(m.group("score")) if m else None
