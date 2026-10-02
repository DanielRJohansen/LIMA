#!/bin/bash
# Builds the Linux release binary and tarball in WSL, then smoke tests it on the GPU.
# Called by distribution/release.bat:  build.sh REPO_DIR COMMIT VERSION OUT_DIR
#
# Layout of the tarball (relocatable: lima finds ../share/LIMA next to its bin dir):
#   lima-<version>-linux-x86_64/
#     bin/lima
#     share/LIMA/resources/
#     LICENSE.txt, README.txt
set -euo pipefail

REPO_DIR=$1
COMMIT=$2
VERSION=$3
OUT_DIR=$4

source "$(dirname "$0")/config.sh"

# Build in the WSL filesystem, building on /mnt/c is many times slower
WORK_DIR="$HOME/lima-release"
rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"
echo "### Checking out $COMMIT"
git clone --quiet --no-checkout "$REPO_DIR" "$WORK_DIR/src"
git -C "$WORK_DIR/src" checkout --quiet "$COMMIT"

echo "### Building for CUDA architectures $CUDA_ARCHITECTURES"
cmake -S "$WORK_DIR/src" -B "$WORK_DIR/build" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=gcc-14 -DCMAKE_CXX_COMPILER=g++-14 \
    -DCMAKE_CUDA_COMPILER="$CUDA_DIR/bin/nvcc" -DCMAKE_CUDA_HOST_COMPILER=g++-14 \
    -DCMAKE_CUDA_ARCHITECTURES="$CUDA_ARCHITECTURES"
cmake --build "$WORK_DIR/build" --target lima

BUILT_VERSION=$("$WORK_DIR/build/code/LIMA/lima" --version)
[ "$BUILT_VERSION" = "lima $VERSION" ] || { echo "Built '$BUILT_VERSION', expected 'lima $VERSION'"; exit 1; }

echo "### Packaging"
NAME="lima-$VERSION-linux-x86_64"
STAGE="$WORK_DIR/$NAME"
mkdir -p "$STAGE/bin" "$STAGE/share/LIMA"
install -m755 "$WORK_DIR/build/code/LIMA/lima" "$STAGE/bin/lima"
strip "$STAGE/bin/lima"
git -C "$WORK_DIR/src" archive "$COMMIT" resources | tar -x -C "$STAGE/share/LIMA"
cp "$WORK_DIR/src/LICENSE.txt" "$WORK_DIR/src/THIRD_PARTY_NOTICES.txt" "$STAGE/"
# The citations for the bundled force fields and lipids must ship with every copy
git -C "$WORK_DIR/src" show "$COMMIT:README.md" | sed -n '/^## LIMA would not be possible/,$p' > "$STAGE/CITATIONS.md"
[ -s "$STAGE/CITATIONS.md" ] || { echo "The citations section was not found in README.md"; exit 1; }
cat > "$STAGE/README.txt" <<EOF
LIMA $VERSION for Linux (x86_64)

Requirements
  An NVIDIA GPU of the RTX 40-series or newer (or H100, B200 and similar), with a
  recent NVIDIA driver (CUDA 13 capable). Runtime libraries: glfw, OpenGL, oneTBB.
  Ubuntu/Debian and Arch users can install the .deb or PKGBUILD from the release
  instead, which handle these dependencies.

Getting started
  Extract anywhere. bin/lima finds share/LIMA next to it, so keep the layout intact.
  Run 'bin/lima --help' for the available commands.

License
  Free for small companies, noncommercial use and academia, and free to evaluate
  for 31 days, see LICENSE.txt. Other use requires a commercial license:
  daniel@lima-dynamics.com

https://github.com/$GITHUB_REPO
EOF
mkdir -p "$OUT_DIR"
tar -czf "$OUT_DIR/$NAME.tar.gz" -C "$WORK_DIR" "$NAME"
echo "Packaged $OUT_DIR/$NAME.tar.gz ($(du -h "$OUT_DIR/$NAME.tar.gz" | cut -f1))"

echo "### Smoke testing on the GPU"
# Runs the packaged binary from outside the source tree, so it must find its resources through the package layout.
# There is no display in WSL on Windows 10, so only headless commands
TEST_DIR="$WORK_DIR/smoketest"
mkdir -p "$TEST_DIR"
FIXTURES="$WORK_DIR/src/tests/clitests/fixtures"
cp "$FIXTURES"/{met_box4.gro,met_box4.top,metsol_clash.gro,metsol.gro,metsol.top,Protein_chain_A.itp,SOL.itp} "$TEST_DIR/"
LIMA="$STAGE/bin/lima"
(
    cd "$TEST_DIR"
    set -x
    "$LIMA" --help > /dev/null
    "$LIMA" makebox --box-size 5
    "$LIMA" solvate -c met_box4.gro -t met_box4.top
    "$LIMA" em -c metsol_clash.gro -t metsol.top --conf-out em.gro
    printf 'n_steps=500\ndt=2\ndata_logging_interval=50\n' > params.txt
    "$LIMA" mdrun -c metsol.gro -t metsol.top -s params.txt --conf-out md.gro --trajectory md.trr
)
for file in conf.gro met_box4_solvated.gro em.gro md.gro md.trr; do
    [ -s "$TEST_DIR/$file" ] || { echo "Smoke test did not produce $file"; exit 1; }
done
echo "Smoke tests passed"
