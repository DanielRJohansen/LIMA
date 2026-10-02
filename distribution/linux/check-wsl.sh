#!/bin/bash
# Verifies that setup-wsl.sh has been run. Called by distribution/release.bat:  check-wsl.sh [--upload]
source "$(dirname "$0")/config.sh"

problems=()
[ -x "$CUDA_DIR/bin/nvcc" ] || problems+=("CUDA is missing from $CUDA_DIR")
command -v g++-14 > /dev/null || problems+=("g++-14 is missing")
command -v ninja > /dev/null || problems+=("ninja is missing")
command -v dpkg-shlibdeps > /dev/null || problems+=("dpkg-dev is missing")
nvidia-smi > /dev/null 2>&1 || problems+=("the GPU is not visible in WSL")
docker version > /dev/null 2>&1 || problems+=("Docker Engine is not running in this distro")
if [ "${1:-}" = "--upload" ]; then
    gh auth status > /dev/null 2>&1 || problems+=("gh is not logged in to GitHub")
fi

if [ ${#problems[@]} -gt 0 ]; then
    echo "WSL is not ready to build a release:"
    printf '  - %s\n' "${problems[@]}"
    echo "Run: wsl -d $WSL_DISTRO_NAME -- bash $(realpath "$(dirname "$0")")/setup-wsl.sh"
    exit 1
fi
