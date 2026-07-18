"""Command registry + option schema for the web client (Steps 2.1, 2.2).

A curated list of the analysis commands the UI offers (each with a one-line
blurb) drives the command palette; `command_schema.json` supplies the typed
option schema that drives the per-command forms. Names match Clann's REPL / the
library dispatch table and the options each command actually parses.
"""

from __future__ import annotations

import json
import os

_SCHEMA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "command_schema.json")
_schema: dict | None = None


def _load_schema() -> dict:
    global _schema
    if _schema is None:
        with open(_SCHEMA_PATH, encoding="utf-8") as f:
            _schema = json.load(f)
    return _schema


def command_schema(name: str) -> dict:
    """Return {summary, options:[…]} for a command, or {} if it has no form."""
    entry = _load_schema().get(name)
    if not entry:
        return {}
    return {"summary": entry.get("summary", ""),
            "options": entry.get("options", [])}

# Ordered by group; the web client renders each group as an <optgroup> header in
# the command dropdown. Every command here is dispatched by the library
# (clann_api.c) — the web client cannot run commands the library doesn't handle.
COMMANDS = [
    # -- Supertree reconstruction --
    {"group": "Supertree reconstruction", "name": "hs",         "blurb": "Heuristic search for the best supertree under the selected criterion"},
    {"group": "Supertree reconstruction", "name": "bootstrap",  "blurb": "Bootstrap supertree analysis under the selected criterion"},
    {"group": "Supertree reconstruction", "name": "nj",         "blurb": "Neighbour-joining supertree"},
    {"group": "Supertree reconstruction", "name": "alltrees",   "blurb": "Exhaustively search all possible supertrees (small problems)"},
    {"group": "Supertree reconstruction", "name": "usertrees",  "blurb": "Assess user-defined supertrees (from a file) to find the best scoring"},
    {"group": "Supertree reconstruction", "name": "consensus",  "blurb": "Consensus tree of all trees containing all taxa"},
    # -- Landscape analysis --
    {"group": "Landscape analysis",       "name": "recluster",  "blurb": "Cluster a landscape TSV (from hs visitedtrees=) at a given RF threshold"},
    # -- Source tree selection & modification --
    {"group": "Source tree selection & modification", "name": "savetrees",      "blurb": "Save source trees to file in Phylip format (subsets selectable by criteria)"},
    {"group": "Source tree selection & modification", "name": "showtrees",      "blurb": "Visualise selected source trees (can also save selected trees to file)"},
    {"group": "Source tree selection & modification", "name": "excludetrees",   "blurb": "Exclude source trees from analyses (by a variety of criteria)"},
    {"group": "Source tree selection & modification", "name": "includetrees",   "blurb": "Restore previously excluded source trees (by a variety of criteria)"},
    {"group": "Source tree selection & modification", "name": "deletetaxa",     "blurb": "Delete taxa from all source trees in memory (prune, preserving branch lengths)"},
    {"group": "Source tree selection & modification", "name": "restoretaxa",    "blurb": "Restore the original trees from before the last deletetaxa"},
    {"group": "Source tree selection & modification", "name": "randomisetrees", "blurb": "Randomise the source trees in memory (preserving each tree's taxa)"},
    {"group": "Source tree selection & modification", "name": "prunemonophylies",   "blurb": "Prune same-species clades to a single representative"},
    {"group": "Source tree selection & modification", "name": "decomposegenetrees", "blurb": "Decompose multi-copy gene trees into ortholog subtrees (non-destructive; writes files)"},
    # -- Miscellaneous calculations --
    {"group": "Miscellaneous calculations", "name": "rfdists",       "blurb": "Robinson–Foulds distances between all source trees"},
    {"group": "Miscellaneous calculations", "name": "generatetrees", "blurb": "Generate random supertrees & assess against the source trees in memory"},
    {"group": "Miscellaneous calculations", "name": "yaptp",         "blurb": "'Yet another PTP' — a randomisation (permutation-tail-probability) test"},
    {"group": "Miscellaneous calculations", "name": "reconstruct",   "blurb": "Gene-tree reconciliation (source trees against a species tree)"},
    # -- Experimental --
    {"group": "Experimental", "name": "sprdists",           "blurb": "Estimate SPR distances of real data vs ideal & randomised versions"},
    # -- Settings --
    {"group": "Settings", "name": "set", "blurb": "Set global options (criterion, seed, ML parameters, …)"},
]


def list_commands() -> list:
    return COMMANDS
