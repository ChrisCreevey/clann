"""Core subprocess runner and result/error types for pyclann."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

__all__ = ["ClannResult", "ClannError", "get_clann_path", "set_clann_path", "_run_clann"]

# ---------------------------------------------------------------------------
# Binary location
# ---------------------------------------------------------------------------

_clann_path: str | None = None


def get_clann_path() -> str:
    """Return the path to the ``clann`` binary.

    Raises :class:`ClannError` if the binary cannot be found on ``PATH``.
    Override the search with :func:`set_clann_path`.
    """
    global _clann_path
    if _clann_path is not None:
        return _clann_path
    found = shutil.which("clann")
    if found is None:
        raise ClannError(
            "clann binary not found on PATH.  Install clann or call "
            "pyclann.set_clann_path('/path/to/clann') before use."
        )
    return found


def set_clann_path(path: str | os.PathLike) -> None:
    """Override the path to the ``clann`` binary for all subsequent calls."""
    global _clann_path
    _clann_path = str(path)


# ---------------------------------------------------------------------------
# Result type
# ---------------------------------------------------------------------------


@dataclass
class ClannResult:
    """Result returned by every pyclann command function.

    Attributes
    ----------
    best_tree:
        The best-scoring Newick tree string, or ``None`` if the command does
        not produce a single best tree.
    trees:
        List of all Newick tree strings returned (may be more than one when
        ``nbest > 1``).
    score:
        Score of the best tree, or ``None`` when not applicable (e.g. NJ).
    scores:
        Scores for all returned trees (parallel to *trees*).
    stdout:
        Full captured standard output from the clann process.
    stderr:
        Full captured standard error from the clann process.
    returncode:
        Process exit code (0 = success).
    command:
        The command-line string that was executed.
    """

    best_tree: str | None = None
    trees: list[str] = field(default_factory=list)
    score: float | None = None
    scores: list[float] = field(default_factory=list)
    stdout: str = ""
    stderr: str = ""
    returncode: int = 0
    command: str = ""

    def __repr__(self) -> str:  # pragma: no cover
        score_str = f"{self.score:.6f}" if self.score is not None else "None"
        tree = self.best_tree or ""
        preview = (tree[:60] + "...") if len(tree) > 60 else tree
        return (
            f"ClannResult(score={score_str}, best_tree={preview!r}, "
            f"n_trees={len(self.trees)}, returncode={self.returncode})"
        )


# ---------------------------------------------------------------------------
# Error type
# ---------------------------------------------------------------------------


class ClannError(RuntimeError):
    """Raised when a clann command fails or the binary is not configured."""


# ---------------------------------------------------------------------------
# Low-level runner
# ---------------------------------------------------------------------------


def _run_clann(
    args: Sequence[str],
    treefile: str | os.PathLike | None = None,
    output_files: dict[str, str] | None = None,
    clann_path: str | None = None,
    timeout: int | None = None,
) -> tuple[subprocess.CompletedProcess[str], dict[str, str]]:
    """Run ``clann`` with *args* and return the process result plus output file contents.

    Parameters
    ----------
    args:
        Full argument list starting with the subcommand name
        (e.g. ``["hs", "trees.ph", "nreps=10"]``).  The clann binary is
        prepended automatically.
    treefile:
        Path to the primary input tree file.  The absolute path is resolved
        before the subprocess call so that relative paths work regardless of
        the working directory used internally.
    output_files:
        Mapping of ``{key: filename_in_workdir}`` for output files to read
        after the run (e.g. ``{"result": "Heuristic_result.txt"}``).
        Values are empty strings for files that were not produced.
    clann_path:
        Override the clann binary path for this single call.
    timeout:
        Optional timeout in seconds forwarded to :func:`subprocess.run`.

    Returns
    -------
    (proc, file_contents):
        *proc* is the :class:`subprocess.CompletedProcess` instance.
        *file_contents* maps the same keys as *output_files* to their text.
    """
    binary = clann_path or get_clann_path()

    if treefile is not None:
        treefile = Path(treefile).resolve()
        if not treefile.exists():
            raise ClannError(f"Input tree file not found: {treefile}")

    with tempfile.TemporaryDirectory(prefix="pyclann_") as tmpdir:
        cmd = [binary, *args]

        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=tmpdir,
                timeout=timeout,
            )
        except FileNotFoundError as exc:
            raise ClannError(
                f"clann binary not found or not executable: {exc.filename}"
            ) from exc

        file_contents: dict[str, str] = {}
        if output_files:
            tmppath = Path(tmpdir)
            for key, fname in output_files.items():
                fpath = tmppath / fname
                file_contents[key] = fpath.read_text() if fpath.exists() else ""

    return proc, file_contents
