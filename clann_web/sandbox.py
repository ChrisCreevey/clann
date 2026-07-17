"""Per-session filesystem sandbox + command path confinement (Step 1.2).

Clann commands reference files freely — `exe <path>`, `speciestree=<file>`,
`savetrees=<path>`, `htmlview=<path>`, … — and the standalone tool happily reads
and writes anywhere. Behind a network port that is unacceptable. Every session
gets its own temp directory; uploads land there, the engine's working directory
is set to it, and every file path in a command is confined to it:

  * input paths (files Clann reads) must be a bare filename that has been
    uploaded into the sandbox — absolute paths, directory components, and `..`
    are rejected;
  * output paths (files Clann writes) are relocated to their basename inside the
    sandbox, so e.g. `savetrees=/etc/x` becomes `savetrees=x` in the sandbox and
    can never escape.

Non-path option values (yes/no/memory/nj/random/matrix/…) are left untouched.
"""

from __future__ import annotations

import os
import shutil
import tempfile


class UnsafePath(ValueError):
    """Raised when a command references a file it must not be allowed to."""


# Option values that are keywords, not file paths — never treated as paths.
_SENTINELS = {
    "yes", "no", "memory", "nj", "random", "all", "mindup", "none", "auto",
    "matrix", "vector", "equal", "comparisons", "splits", "taxa", "quartets",
    "4point", "ultrametric", "first", "best", "legacy", "standard",
    "score", "visits", "paper", "lust", "lnl",
}

# Options whose value is a file Clann READS (must already be in the sandbox).
_INPUT_OPTS = {
    "start", "speciestree", "file", "testsfile",
    "scorematrix", "scorematrixfile", "clanfile", "guidetree",
}

# Options whose value is a file Clann WRITES (relocated into the sandbox).
_OUTPUT_OPTS = {
    "savetrees", "htmlview", "nhxfile", "filename", "histogramfile",
    "visitedtrees", "clusteroutput", "output_file", "resultjson",
}


def safe_basename(name: str) -> str:
    """Return `name` if it is a single, safe filename component, else raise.

    Rejects absolute paths, any directory separator, and `.`/`..`.
    """
    if not name:
        raise UnsafePath("empty filename")
    if os.path.isabs(name) or "/" in name or "\\" in name:
        raise UnsafePath(f"'{name}' must be a bare filename (no directories)")
    if name in (".", "..") or name.startswith(".."):
        raise UnsafePath(f"'{name}' is not an allowed filename")
    return name


class Sandbox:
    """One session's private directory."""

    def __init__(self):
        self.dir = tempfile.mkdtemp(prefix="clannweb-")

    def put(self, name: str, data: bytes) -> str:
        """Store an uploaded file; returns the stored basename."""
        base = safe_basename(name)
        with open(os.path.join(self.dir, base), "wb") as f:
            f.write(data)
        return base

    def has(self, name: str) -> bool:
        try:
            base = safe_basename(name)
        except UnsafePath:
            return False
        return os.path.isfile(os.path.join(self.dir, base))

    def list(self) -> list[str]:
        return sorted(os.listdir(self.dir))

    def destroy(self) -> None:
        shutil.rmtree(self.dir, ignore_errors=True)

    # -- confinement helpers ------------------------------------------------
    def confine_input(self, name: str) -> str:
        """A file Clann will read: must already be uploaded into the sandbox."""
        base = safe_basename(name)  # rejects paths/traversal
        if not self.has(base):
            raise UnsafePath(
                f"input file '{base}' is not in this session — upload it first")
        return base

    def confine_output(self, name: str) -> str:
        """A file Clann will write: relocate to its basename inside the sandbox."""
        base = os.path.basename(name)
        if base in ("", ".", ".."):
            raise UnsafePath(f"invalid output filename '{name}'")
        return base


def sanitize_command(command: str, sandbox: Sandbox) -> str:
    """Confine every file path a command references to `sandbox`.

    Raises UnsafePath if an input file escapes the sandbox (or isn't uploaded).
    Output paths are rewritten to sandbox-local basenames.
    """
    tokens = command.split()
    if not tokens:
        return command

    out = [tokens[0]]
    i = 1
    # Commands whose bare first argument is an input file path:
    #   `exe <file>` / `execute <file>`   — the gene-tree file to load
    #   `usertrees <file>`                — the candidate-topologies file to score
    if tokens[0] in ("exe", "execute", "usertrees") and len(tokens) > 1 and "=" not in tokens[1]:
        out.append(sandbox.confine_input(tokens[1]))
        i = 2

    for tok in tokens[i:]:
        out.append(_confine_token(tok, sandbox))
    return " ".join(out)


def _confine_token(tok: str, sandbox: Sandbox) -> str:
    """Confine a single `key=value` token to the sandbox (or pass it through)."""
    if "=" not in tok:
        return tok
    key, val = tok.split("=", 1)
    kl = key.lower()
    if val.lower() in _SENTINELS or val == "":
        return tok
    if kl in _INPUT_OPTS:
        return f"{key}={sandbox.confine_input(val)}"
    if kl in _OUTPUT_OPTS:
        return f"{key}={sandbox.confine_output(val)}"
    return tok


def sanitize_options(options: str, sandbox: Sandbox) -> str:
    """Confine file paths in a bare option string (no leading command word).

    Used for `exe` load-time options (`clanfile=`, …) which are passed to
    clann_load_trees as parse options rather than through a full command.
    """
    return " ".join(_confine_token(t, sandbox) for t in options.split())
