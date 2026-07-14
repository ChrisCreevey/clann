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

COMMANDS = [
    {"name": "hs",          "blurb": "Heuristic search for the best supertree"},
    {"name": "nj",          "blurb": "Neighbour-joining supertree (average consensus)"},
    {"name": "reconstruct", "blurb": "Reconcile gene trees onto a species tree (duplications & losses)"},
    {"name": "showtrees",   "blurb": "Display the loaded input gene trees"},
    {"name": "bootstrap",   "blurb": "Bootstrap the heuristic search"},
    {"name": "consensus",   "blurb": "Consensus of the trees in memory"},
    {"name": "alltrees",    "blurb": "Exhaustive search over all supertrees (small problems)"},
    {"name": "usertrees",   "blurb": "Score user-supplied candidate supertrees"},
    {"name": "rfdists",     "blurb": "Robinson–Foulds distances among source trees"},
    {"name": "mlscores",    "blurb": "Maximum-likelihood scores / parameter scan"},
    {"name": "set",         "blurb": "Set global options (criterion, seed, …)"},
]


def list_commands() -> list:
    return COMMANDS
