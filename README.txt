Clann is a program for the construction of supertrees and analysis of
phylogenomic content.


Copyright (c) 2003 - 2026 Chris Creevey
SPDX-License-Identifier: GPL-2.0-or-later

Clann is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

A copy of the GNU General Public License is included in COPYING.txt.


Installation (quick start)
--------------------------

From source:
  ./configure
  make
  sudo make install        # optional — copies clann to /usr/local/bin

See the INSTALL file for full instructions, including prerequisites (readline,
GCC/OpenMP), macOS notes, and troubleshooting.

Via Conda (Bioconda):
  conda install -c bioconda -c conda-forge clann

Via Homebrew (macOS/Linux):
  brew tap ChrisCreevey/clann
  brew install clann

See PACKAGING.md for detailed Conda and Homebrew instructions.


Dependencies
------------

Clann detects dependencies automatically at configure time:

  * GNU readline (libreadline-dev / readline-devel) — recommended for the
    interactive prompt; Clann falls back to plain input if not found.

  * GCC with OpenMP — required for multi-threaded parallel search.
    On macOS, configure auto-selects a Homebrew GCC when Apple Clang is
    detected.  Install with: brew install gcc

Install on Debian/Ubuntu:  sudo apt install build-essential libreadline-dev
Install on Red Hat/Fedora:  sudo yum install gcc readline-devel
Install on macOS:           brew install gcc readline


Documentation
-------------

  README.md     — Overview, usage examples, and command reference
  USER_MANUAL.md — Full command and option reference
  TUTORIAL.md   — Worked tutorial for new features
  PACKAGING.md  — Conda and Homebrew packaging/submission guide

For more details see the Clann home page:
  https://github.com/ChrisCreevey/clann


Contact
-------

Prof. Christopher Creevey
School of Biological Sciences
Queen's University Belfast
Northern Ireland, UK
email: chris.creevey@gmail.com

