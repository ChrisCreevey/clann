"""Setuptools shim.

Project metadata lives in ``pyproject.toml``; this file exists so that older
pip/setuptools can install the package (including editable installs), and to
**bundle the hardened server library** into the package when it has been built.

``libclann-server.so`` (a ``.dylib`` on macOS) is produced at the repository root
by ``make libclann-server.so``. Copy it into ``clann_web/`` at build time so it
ships as package-data and a non-editable ``pip install .`` can find it. If the
library hasn't been built yet, the install still succeeds — ``clann_web`` will
report a clear 'build it first' error at runtime, or honour ``$CLANN_LIB``.
"""

import os
import shutil

from setuptools import setup

_HERE = os.path.dirname(os.path.abspath(__file__))
for _name in ("libclann-server.so", "libclann-server.dylib"):
    _src = os.path.join(_HERE, _name)
    if os.path.exists(_src):
        shutil.copy2(_src, os.path.join(_HERE, "clann_web", _name))

setup()
