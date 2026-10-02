#!/bin/bash
# Builds the .deb from the tarball made by build.sh, and test-installs it in a clean Ubuntu container.
# Called by distribution/release.bat:  package-deb.sh VERSION OUT_DIR
set -euo pipefail

VERSION=$1
OUT_DIR=$2

source "$(dirname "$0")/config.sh"

WORK_DIR="$HOME/lima-release/deb"
NAME="lima-$VERSION-linux-x86_64"
DEB="lima_${VERSION}_amd64.deb"
rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"
tar -xzf "$OUT_DIR/$NAME.tar.gz" -C "$WORK_DIR"

ROOT="$WORK_DIR/root"
mkdir -p "$ROOT/DEBIAN" "$ROOT/usr/bin" "$ROOT/usr/share/doc/lima"
cp "$WORK_DIR/$NAME/bin/lima" "$ROOT/usr/bin/"
cp -r "$WORK_DIR/$NAME/share/LIMA" "$ROOT/usr/share/"
cp "$WORK_DIR/$NAME/LICENSE.txt" "$ROOT/usr/share/doc/lima/copyright"
cp "$WORK_DIR/$NAME/CITATIONS.md" "$WORK_DIR/$NAME/THIRD_PARTY_NOTICES.txt" "$ROOT/usr/share/doc/lima/"

# Dependencies come from the libraries the binary actually links, so they cannot go stale.
# dpkg-shlibdeps needs a debian/control to run in
mkdir -p "$WORK_DIR/shlibs/debian"
printf 'Source: lima\n\nPackage: lima\nArchitecture: amd64\n' > "$WORK_DIR/shlibs/debian/control"
DEPENDS=$(cd "$WORK_DIR/shlibs" && dpkg-shlibdeps -O "$ROOT/usr/bin/lima" 2>/dev/null | sed -n 's/^shlibs:Depends=//p')
[ -n "$DEPENDS" ] || { echo "dpkg-shlibdeps found no dependencies"; exit 1; }

cat > "$ROOT/DEBIAN/control" <<EOF
Package: lima
Version: $VERSION
Architecture: amd64
Maintainer: $MAINTAINER
Depends: $DEPENDS
Section: science
Priority: optional
Homepage: https://github.com/$GITHUB_REPO
Installed-Size: $(du -sk "$ROOT/usr" | cut -f1)
Description: $DESCRIPTION
 LIMA is a GPU accelerated molecular dynamics engine, with tools to build,
 solvate, minimize, simulate and render molecular systems.
 .
 Requires an NVIDIA GPU of the RTX 40-series or newer with a recent NVIDIA driver.
EOF

dpkg-deb --build --root-owner-group "$ROOT" "$OUT_DIR/$DEB" > /dev/null
echo "Packaged $OUT_DIR/$DEB ($(du -h "$OUT_DIR/$DEB" | cut -f1)), depends: $DEPENDS"

echo "### Test-installing in a clean Ubuntu container"
# No GPU in the container, so this checks installation, dependencies and resource lookup, not simulation
# (togmx needs the force field from the installed resources, and runs on the CPU)
FIXTURES="$HOME/lima-release/src/tests/clitests/fixtures"
docker run --rm -v "$OUT_DIR:/release:ro" -v "$FIXTURES:/fixtures:ro" ubuntu:24.04 bash -c "
    set -e
    export DEBIAN_FRONTEND=noninteractive
    apt-get update -qq
    apt-get install -y -qq /release/$DEB > /dev/null
    ! ldd /usr/bin/lima | grep 'not found'
    test \"\$(lima --version)\" = 'lima $VERSION'
    cd /tmp && cp /fixtures/6lzm.pdb . && lima togmx -f 6lzm.pdb && test -s conf.gro
    echo 'The .deb installs and runs'
"
