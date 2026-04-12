"""pyclann — Python interface to the Clann phylogenomics tool.

This package provides a subprocess-based Python API for scripting with
Clann.  Every high-level function (``hs``, ``nj``, ``consensus``, …) runs
the ``clann`` binary in a temporary working directory, captures its output,
parses the result files, and returns a :class:`ClannResult` object.

Typical usage::

    import pyclann

    result = pyclann.hs("trees.ph", criterion="ml", nreps=10)
    print(result.best_tree)
    print(result.score)

    result = pyclann.consensus("trees.ph")
    print(result.best_tree)

The ``clann`` binary must be on the system ``PATH``, or its location can be
set explicitly::

    pyclann.set_clann_path("/path/to/clann")

For a description of all command options see the Clann ``USER_MANUAL.md``.
"""

from ._core import ClannError, ClannResult, get_clann_path, set_clann_path
from ._commands import alltrees, bootstrap, consensus, hs, nj, run, usertrees

__all__ = [
    "ClannResult",
    "ClannError",
    "get_clann_path",
    "set_clann_path",
    "hs",
    "nj",
    "consensus",
    "alltrees",
    "usertrees",
    "bootstrap",
    "run",
]
