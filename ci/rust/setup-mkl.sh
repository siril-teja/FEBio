#!/usr/bin/env bash
# Install MKL 2024.2 from PyPI into $1 and lay it out the way FEBio's CMake expects
# (it only searches for lib<name>.a). The shared objects are linked under .a names, which
# GNU ld accepts. MKL 2025+ removed mkl_dcsrtrsv (used by NumCore/ILU0), hence the pin.
# Prints the cmake arguments to stdout.
set -euo pipefail
DEST=$(realpath -m "${1:-mkl}")
mkdir -p "$DEST"
python3 -m venv "$DEST/venv"
"$DEST/venv/bin/pip" install -q "mkl-devel==2024.2.*" >&2
LIB="$DEST/venv/lib"
mkdir -p "$DEST/lib"
for f in "$LIB"/libmkl_*.so.2 "$LIB"/libiomp5.so; do
  ln -sf "$f" "$DEST/lib/$(basename "$f" .2)"
done
for n in mkl_core mkl_intel_lp64 mkl_intel_thread; do
  ln -sf "$LIB/lib$n.so.2" "$DEST/lib/lib$n.a"
done
ln -sfn "$DEST/venv/include" "$DEST/include"
echo "-DUSE_MKL=ON -DMKLROOT=$DEST -DMKL_INC=$DEST/include -DMKL_LIB_DIR=$DEST/lib -DMKL_OMP_LIB=$DEST/lib/libiomp5.so"
echo "MKL_RUNTIME_DIR=$LIB" >&2
