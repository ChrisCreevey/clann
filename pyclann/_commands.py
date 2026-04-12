"""High-level command functions for pyclann.

Each function wraps one Clann CLI command, runs the binary in a temporary
working directory, and returns a :class:`~pyclann.ClannResult`.

All keyword arguments that are not explicitly listed are forwarded to Clann as
``key=value`` option strings (e.g. ``mlbeta=2.0`` becomes the CLI token
``mlbeta=2.0``).  Boolean values are serialised as ``yes``/``no``.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

from ._core import ClannError, ClannResult, _run_clann
from ._parser import (
    extract_stdout_score,
    parse_bootstrap_file,
    parse_plain_tree_file,
    parse_scored_tree_file,
)

__all__ = ["hs", "nj", "consensus", "alltrees", "usertrees", "bootstrap", "run"]

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _opts_to_args(opts: dict[str, Any]) -> list[str]:
    """Convert a dict of keyword options to Clann ``key=value`` tokens.

    * ``None`` values are skipped.
    * ``True`` / ``False`` are converted to ``yes`` / ``no``.
    * All other values are stringified with ``str()``.
    """
    tokens: list[str] = []
    for k, v in opts.items():
        if v is None:
            continue
        if isinstance(v, bool):
            tokens.append(f"{k}={'yes' if v else 'no'}")
        else:
            tokens.append(f"{k}={v}")
    return tokens


def _best_from_scored(
    pairs: list[tuple[str, float]],
) -> tuple[str | None, float | None]:
    """Return the (best_tree, best_score) from a list of (newick, score) pairs.

    Clann minimises the score for DFIT/RF and maximises it for ML (which it
    prints as a negative lnL).  The first entry in the results file is always
    the best tree, so we simply return *pairs[0]*.
    """
    if not pairs:
        return None, None
    return pairs[0][0], pairs[0][1]


# ---------------------------------------------------------------------------
# Public command functions
# ---------------------------------------------------------------------------


def hs(
    treefile: str | os.PathLike,
    *,
    criterion: str | None = None,
    nreps: int | None = None,
    nthreads: int | None = None,
    nbest: int | None = None,
    savetrees: str | None = None,
    start: str | None = None,
    seed: int | None = None,
    clann_path: str | None = None,
    timeout: int | None = None,
    **kwargs: Any,
) -> ClannResult:
    """Run a heuristic supertree search (``clann hs``).

    Parameters
    ----------
    treefile:
        Path to the input source-tree file (Newick / Nexus).
    criterion:
        Optimality criterion: ``"dfit"`` (default), ``"ml"``, ``"rf"``,
        ``"sfit"``, ``"qfit"``, ``"avcon"``.
    nreps:
        Number of heuristic search replicates.
    nthreads:
        Number of OpenMP threads (default: all available CPUs).
    nbest:
        Number of best topologies to retain.
    savetrees:
        Override the default output filename (``Heuristic_result.txt``).
    start:
        Starting tree strategy: ``"nj"`` (default NJ tree) or ``"memory"``.
    seed:
        Random-number seed for reproducibility.
    clann_path:
        Override the path to the ``clann`` binary for this call.
    timeout:
        Optional timeout in seconds for the subprocess.
    **kwargs:
        Any additional ``key=value`` options passed directly to ``clann hs``.

    Returns
    -------
    ClannResult:
        ``.best_tree`` — best Newick string; ``.score`` — best score;
        ``.trees`` / ``.scores`` — all retained topologies.
    """
    treefile = Path(treefile).resolve()
    outfile = savetrees or "Heuristic_result.txt"

    named_opts: dict[str, Any] = {
        "criterion": criterion,
        "nreps": nreps,
        "nthreads": nthreads,
        "nbest": nbest,
        "savetrees": savetrees,
        "start": start,
        "seed": seed,
    }
    named_opts.update(kwargs)
    opt_args = _opts_to_args(named_opts)

    args = ["hs", str(treefile), *opt_args]
    proc, files = _run_clann(
        args,
        treefile=treefile,
        output_files={"result": outfile},
        clann_path=clann_path,
        timeout=timeout,
    )

    pairs = parse_scored_tree_file(files.get("result", ""))
    best_tree, best_score = _best_from_scored(pairs)
    if best_score is None:
        best_score = extract_stdout_score(proc.stdout)

    return ClannResult(
        best_tree=best_tree,
        trees=[t for t, _ in pairs],
        score=best_score,
        scores=[s for _, s in pairs],
        stdout=proc.stdout,
        stderr=proc.stderr,
        returncode=proc.returncode,
        command=" ".join(args),
    )


def nj(
    treefile: str | os.PathLike,
    *,
    missing: str | None = None,
    savetrees: str | None = None,
    clann_path: str | None = None,
    timeout: int | None = None,
    **kwargs: Any,
) -> ClannResult:
    """Construct a neighbour-joining supertree (``clann nj``).

    Parameters
    ----------
    treefile:
        Path to the input source-tree file.
    missing:
        Method for estimating missing data: ``"4point"`` (default) or
        ``"ultrametric"``.
    savetrees:
        Override the default output filename (``NJ-tree.ph``).
    clann_path:
        Override the path to the ``clann`` binary.
    timeout:
        Optional subprocess timeout in seconds.
    **kwargs:
        Any additional ``key=value`` options.

    Returns
    -------
    ClannResult:
        ``.best_tree`` — the NJ tree Newick string (with branch lengths).
    """
    treefile = Path(treefile).resolve()
    outfile = savetrees or "NJ-tree.ph"

    named_opts: dict[str, Any] = {
        "missing": missing,
        "savetrees": savetrees,
    }
    named_opts.update(kwargs)
    opt_args = _opts_to_args(named_opts)

    args = ["nj", str(treefile), *opt_args]
    proc, files = _run_clann(
        args,
        treefile=treefile,
        output_files={"result": outfile},
        clann_path=clann_path,
        timeout=timeout,
    )

    trees = parse_plain_tree_file(files.get("result", ""))
    best_tree = trees[0] if trees else None

    return ClannResult(
        best_tree=best_tree,
        trees=trees,
        score=None,
        scores=[],
        stdout=proc.stdout,
        stderr=proc.stderr,
        returncode=proc.returncode,
        command=" ".join(args),
    )


def consensus(
    treefile: str | os.PathLike,
    *,
    percentage: float | None = None,
    savetrees: str | None = None,
    clann_path: str | None = None,
    timeout: int | None = None,
    **kwargs: Any,
) -> ClannResult:
    """Compute a consensus tree (``clann consensus``).

    Parameters
    ----------
    treefile:
        Path to the input source-tree file.
    percentage:
        Minimum clade support threshold (0.5–1.0; default 0.5 = majority rule).
    savetrees:
        Override the default output filename (``consensus.ph``).
    clann_path:
        Override the path to the ``clann`` binary.
    timeout:
        Optional subprocess timeout in seconds.
    **kwargs:
        Any additional ``key=value`` options.

    Returns
    -------
    ClannResult:
        ``.best_tree`` — the consensus Newick string (support values embedded
        as branch labels).
    """
    treefile = Path(treefile).resolve()
    outfile = savetrees or "consensus.ph"

    named_opts: dict[str, Any] = {
        "percentage": percentage,
        "savetrees": savetrees,
    }
    named_opts.update(kwargs)
    opt_args = _opts_to_args(named_opts)

    args = ["consensus", str(treefile), *opt_args]
    proc, files = _run_clann(
        args,
        treefile=treefile,
        output_files={"result": outfile},
        clann_path=clann_path,
        timeout=timeout,
    )

    trees = parse_plain_tree_file(files.get("result", ""))
    best_tree = trees[0] if trees else None

    return ClannResult(
        best_tree=best_tree,
        trees=trees,
        score=None,
        scores=[],
        stdout=proc.stdout,
        stderr=proc.stderr,
        returncode=proc.returncode,
        command=" ".join(args),
    )


def alltrees(
    treefile: str | os.PathLike,
    *,
    criterion: str | None = None,
    savetrees: str | None = None,
    seed: int | None = None,
    clann_path: str | None = None,
    timeout: int | None = None,
    **kwargs: Any,
) -> ClannResult:
    """Exhaustively score all possible supertree topologies (``clann alltrees``).

    Only feasible for small taxon sets (≤ 8 taxa).

    Parameters
    ----------
    treefile:
        Path to the input source-tree file.
    criterion:
        Optimality criterion (same choices as :func:`hs`).
    savetrees:
        Override the default output filename (``top_alltrees.txt``).
    seed:
        Random-number seed.
    clann_path:
        Override the path to the ``clann`` binary.
    timeout:
        Optional subprocess timeout in seconds.
    **kwargs:
        Any additional ``key=value`` options.

    Returns
    -------
    ClannResult:
        ``.best_tree`` — best Newick; ``.trees`` / ``.scores`` — all evaluated
        topologies sorted best-first.
    """
    treefile = Path(treefile).resolve()
    # alltrees saves to top_alltrees.txt by default
    outfile = savetrees or "top_alltrees.txt"

    named_opts: dict[str, Any] = {
        "criterion": criterion,
        "savetrees": savetrees,
        "seed": seed,
    }
    named_opts.update(kwargs)
    opt_args = _opts_to_args(named_opts)

    args = ["alltrees", str(treefile), *opt_args]
    proc, files = _run_clann(
        args,
        treefile=treefile,
        output_files={"result": outfile},
        clann_path=clann_path,
        timeout=timeout,
    )

    pairs = parse_scored_tree_file(files.get("result", ""))
    best_tree, best_score = _best_from_scored(pairs)
    if best_score is None:
        best_score = extract_stdout_score(proc.stdout)

    return ClannResult(
        best_tree=best_tree,
        trees=[t for t, _ in pairs],
        score=best_score,
        scores=[s for _, s in pairs],
        stdout=proc.stdout,
        stderr=proc.stderr,
        returncode=proc.returncode,
        command=" ".join(args),
    )


def usertrees(
    treefile: str | os.PathLike,
    candidates: str | os.PathLike,
    *,
    criterion: str | None = None,
    tests: bool | None = None,
    nboot: int | None = None,
    savetrees: str | None = None,
    clann_path: str | None = None,
    timeout: int | None = None,
    **kwargs: Any,
) -> ClannResult:
    """Score user-supplied candidate topologies against source trees (``clann usertrees``).

    Parameters
    ----------
    treefile:
        Path to the source-tree file.
    candidates:
        Path to the file of candidate supertree topologies to score.
    criterion:
        Optimality criterion (same choices as :func:`hs`).
    tests:
        Whether to perform ML topology tests (Winning-Sites, KH, SH).
        Only valid with ``criterion="ml"``.
    nboot:
        Number of bootstrap replicates for topology tests.
    savetrees:
        Override the default output filename (``Usertrees_result.txt``).
    clann_path:
        Override the path to the ``clann`` binary.
    timeout:
        Optional subprocess timeout in seconds.
    **kwargs:
        Any additional ``key=value`` options.

    Returns
    -------
    ClannResult:
        ``.best_tree`` — best-scoring candidate Newick;
        ``.trees`` / ``.scores`` — all scored candidates.
    """
    treefile = Path(treefile).resolve()
    candidates = Path(candidates).resolve()
    if not candidates.exists():
        raise ClannError(f"Candidates file not found: {candidates}")

    outfile = savetrees or "Usertrees_result.txt"

    named_opts: dict[str, Any] = {
        "criterion": criterion,
        "tests": tests,
        "nboot": nboot,
        "savetrees": savetrees,
    }
    named_opts.update(kwargs)
    opt_args = _opts_to_args(named_opts)

    # CLI syntax: clann usertrees source.ph candidates.ph [options]
    args = ["usertrees", str(treefile), str(candidates), *opt_args]
    proc, files = _run_clann(
        args,
        treefile=treefile,
        output_files={"result": outfile},
        clann_path=clann_path,
        timeout=timeout,
    )

    pairs = parse_scored_tree_file(files.get("result", ""))
    best_tree, best_score = _best_from_scored(pairs)
    if best_score is None:
        best_score = extract_stdout_score(proc.stdout)

    return ClannResult(
        best_tree=best_tree,
        trees=[t for t, _ in pairs],
        score=best_score,
        scores=[s for _, s in pairs],
        stdout=proc.stdout,
        stderr=proc.stderr,
        returncode=proc.returncode,
        command=" ".join(args),
    )


def bootstrap(
    treefile: str | os.PathLike,
    *,
    criterion: str | None = None,
    nreps: int | None = None,
    nthreads: int | None = None,
    seed: int | None = None,
    savetrees: str | None = None,
    clann_path: str | None = None,
    timeout: int | None = None,
    **kwargs: Any,
) -> ClannResult:
    """Run a bootstrap analysis (``clann boot``).

    Each bootstrap replicate runs a heuristic search on a resampled set of
    source trees; a majority-rule consensus is then computed from the
    bootstrap trees.

    Parameters
    ----------
    treefile:
        Path to the input source-tree file.
    criterion:
        Optimality criterion (same choices as :func:`hs`).
    nreps:
        Number of bootstrap replicates.
    nthreads:
        Number of OpenMP threads.
    seed:
        Random-number seed for reproducibility.
    savetrees:
        Override the default bootstrap-trees filename (``bootstrap.txt``).
    clann_path:
        Override the path to the ``clann`` binary.
    timeout:
        Optional subprocess timeout in seconds.
    **kwargs:
        Any additional ``key=value`` options.

    Returns
    -------
    ClannResult:
        ``.best_tree`` — consensus Newick from ``consensus.ph``;
        ``.trees`` / ``.scores`` — individual bootstrap tree(s) and scores
        from ``bootstrap.txt``.
    """
    treefile = Path(treefile).resolve()
    boot_outfile = savetrees or "bootstrap.txt"

    named_opts: dict[str, Any] = {
        "criterion": criterion,
        "nreps": nreps,
        "nthreads": nthreads,
        "seed": seed,
        "savetrees": savetrees,
    }
    named_opts.update(kwargs)
    opt_args = _opts_to_args(named_opts)

    args = ["boot", str(treefile), *opt_args]
    proc, files = _run_clann(
        args,
        treefile=treefile,
        output_files={"bootstrap": boot_outfile, "consensus": "consensus.ph"},
        clann_path=clann_path,
        timeout=timeout,
    )

    boot_pairs = parse_bootstrap_file(files.get("bootstrap", ""))
    cons_trees = parse_plain_tree_file(files.get("consensus", ""))
    best_tree = cons_trees[0] if cons_trees else (boot_pairs[0][0] if boot_pairs else None)

    return ClannResult(
        best_tree=best_tree,
        trees=[t for t, _ in boot_pairs],
        score=None,
        scores=[s for _, s in boot_pairs],
        stdout=proc.stdout,
        stderr=proc.stderr,
        returncode=proc.returncode,
        command=" ".join(args),
    )


def run(
    command: str,
    treefile: str | os.PathLike | None = None,
    *,
    clann_path: str | None = None,
    timeout: int | None = None,
    **kwargs: Any,
) -> ClannResult:
    """Run an arbitrary Clann command string.

    This low-level function allows access to any Clann command not covered by
    the typed helpers above.  No output-file parsing is performed; results are
    available via :attr:`ClannResult.stdout`.

    Parameters
    ----------
    command:
        The Clann subcommand, e.g. ``"rfdists"``, ``"sprdists"``,
        ``"showtrees"``.
    treefile:
        Optional path to the input source-tree file.
    clann_path:
        Override the path to the ``clann`` binary.
    timeout:
        Optional subprocess timeout in seconds.
    **kwargs:
        Additional ``key=value`` options appended to the command.

    Returns
    -------
    ClannResult:
        Only ``.stdout``, ``.stderr``, and ``.returncode`` are populated;
        ``.best_tree`` and ``.score`` are ``None``.
    """
    opt_args = _opts_to_args(kwargs)
    args: list[str] = [command]
    if treefile is not None:
        args.append(str(Path(treefile).resolve()))
    args.extend(opt_args)

    proc, _ = _run_clann(
        args,
        treefile=treefile,
        output_files=None,
        clann_path=clann_path,
        timeout=timeout,
    )

    return ClannResult(
        best_tree=None,
        trees=[],
        score=None,
        scores=[],
        stdout=proc.stdout,
        stderr=proc.stderr,
        returncode=proc.returncode,
        command=" ".join(args),
    )
