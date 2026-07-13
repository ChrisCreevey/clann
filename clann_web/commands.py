"""Command registry for the web client (Step 2.1).

A curated list of the analysis commands the UI offers, each with a one-line
blurb. Step 2.2 adds a full per-command option schema; for now this drives the
command palette (dropdown + description). Names match Clann's REPL / the
library dispatch table.
"""

from __future__ import annotations

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
