#!/bin/bash
# Conda build script for clann (Linux and macOS)
set -euo pipefail

# On macOS, enable OpenMP via llvm-openmp if the compiler supports it.
# conda-forge's llvm-openmp package installs into $PREFIX, so we point
# the compiler there.  On Linux, libgomp is already on the search path.
if [[ "$(uname)" == "Darwin" ]]; then
    export CFLAGS="${CFLAGS} -I${PREFIX}/include"
    export LDFLAGS="${LDFLAGS} -L${PREFIX}/lib"
fi

./configure \
    --prefix="${PREFIX}" \
    --disable-silent-rules

make -j"${CPU_COUNT:-1}"
make install
