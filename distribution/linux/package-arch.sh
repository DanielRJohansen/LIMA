#!/bin/bash
# Writes the PKGBUILD for the release tarball, and test-builds and installs it in a clean Arch container.
# Called by distribution/release.bat:  package-arch.sh VERSION OUT_DIR
#
# The PKGBUILD downloads the tarball from the GitHub release, so it only works once the release is published.
# The test passes the local tarball to makepkg instead, which verifies it against the same checksum.
set -euo pipefail

VERSION=$1
OUT_DIR=$2

source "$(dirname "$0")/config.sh"

NAME="lima-$VERSION-linux-x86_64"
TARBALL="$OUT_DIR/$NAME.tar.gz"
WORK_DIR="$HOME/lima-release/arch"
rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"
tar -xzf "$TARBALL" -C "$WORK_DIR" "$NAME/bin/lima"

# Arch package for every library the binary links directly. Unknown libraries fail, rather than ship a broken package
declare -A PACKAGE_OF=(
    [libc.so.6]=glibc [libm.so.6]=glibc [libdl.so.2]=glibc [libpthread.so.0]=glibc [librt.so.1]=glibc
    [ld-linux-x86-64.so.2]=glibc [libgcc_s.so.1]=gcc-libs [libstdc++.so.6]=gcc-libs
    [libGL.so.1]=libglvnd [libGLX.so.0]=libglvnd [libOpenGL.so.0]=libglvnd [libEGL.so.1]=libglvnd
    [libGLU.so.1]=glu [libglfw.so.3]=glfw [libtbb.so.12]=onetbb [libX11.so.6]=libx11
)
DEPENDS=("nvidia-utils")	# Provides libcuda, which the CUDA runtime loads at run time
for library in $(readelf -d "$WORK_DIR/$NAME/bin/lima" | sed -n 's/.*Shared library: \[\(.*\)\]/\1/p'); do
    package=${PACKAGE_OF[$library]:-}
    [ -n "$package" ] || { echo "No Arch package known for $library, add it to package-arch.sh"; exit 1; }
    [[ " ${DEPENDS[*]} " == *" $package "* ]] || DEPENDS+=("$package")
done

SHA256=$(sha256sum "$TARBALL" | cut -d' ' -f1)
cat > "$OUT_DIR/PKGBUILD" <<EOF
# Maintainer: $MAINTAINER

pkgname=lima
pkgver=$VERSION
pkgrel=1
pkgdesc="$DESCRIPTION"
arch=('x86_64')
url="https://github.com/$GITHUB_REPO"
license=('PolyForm-Small-Business-1.0.0' 'PolyForm-Noncommercial-1.0.0' 'LicenseRef-PolyForm-Free-Trial-1.0.0')
depends=($(printf "'%s' " "${DEPENDS[@]}"))
options=('!strip' '!debug')
source=("https://github.com/$GITHUB_REPO/releases/download/v\$pkgver/lima-\$pkgver-linux-x86_64.tar.gz")
sha256sums=('$SHA256')

package() {
    cd "\$srcdir/lima-\$pkgver-linux-x86_64"
    install -Dm755 bin/lima "\$pkgdir/usr/bin/lima"
    install -d "\$pkgdir/usr/share/LIMA"
    cp -r share/LIMA/resources "\$pkgdir/usr/share/LIMA/"
    install -Dm644 LICENSE.txt "\$pkgdir/usr/share/licenses/\$pkgname/LICENSE.txt"
    install -Dm644 CITATIONS.md "\$pkgdir/usr/share/doc/\$pkgname/CITATIONS.md"
    install -Dm644 THIRD_PARTY_NOTICES.txt "\$pkgdir/usr/share/licenses/\$pkgname/THIRD_PARTY_NOTICES.txt"
}
EOF
echo "Wrote $OUT_DIR/PKGBUILD, depends: ${DEPENDS[*]}"

echo "### Test-building and installing in a clean Arch container"
# makepkg refuses to run as root, so it runs as a throwaway user. No GPU in the container, so this checks
# the package, dependencies and resource lookup (togmx reads the force field and runs on the CPU)
FIXTURES="$HOME/lima-release/src/tests/clitests/fixtures"
docker run --rm -v "$OUT_DIR:/release:ro" -v "$FIXTURES:/fixtures:ro" archlinux:base-devel bash -c "
    set -e
    pacman -Syu --noconfirm --needed > /dev/null
    useradd -m builder
    echo 'builder ALL=(ALL) NOPASSWD: ALL' > /etc/sudoers.d/builder
    cp /release/PKGBUILD /release/$NAME.tar.gz /home/builder/
    chown builder /home/builder/*
    su builder -c 'cd ~ && makepkg --syncdeps --install --noconfirm --needed' > /dev/null
    ! ldd /usr/bin/lima | grep 'not found'
    test \"\$(lima --version)\" = 'lima $VERSION'
    cd /tmp && cp /fixtures/6lzm.pdb . && lima togmx -f 6lzm.pdb && test -s conf.gro
    echo 'The PKGBUILD builds, installs and runs'
"
