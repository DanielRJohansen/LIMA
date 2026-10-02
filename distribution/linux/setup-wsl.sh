#!/bin/bash
# One-time setup of the WSL distro that builds LIMA's Linux release. Run it inside WSL:
#   wsl -d Ubuntu-24.04 -- bash /mnt/c/<path to repo>/distribution/linux/setup-wsl.sh
# Safe to run again; it only installs what is missing.
set -euo pipefail

source "$(dirname "$0")/config.sh"
CUDA_PACKAGE_VERSION=$(basename "$CUDA_DIR" | sed 's/^cuda-//; s/\./-/')	# /usr/local/cuda-13.2 -> 13-2

echo "### Build tools and LIMA's dependencies"
sudo apt-get update
sudo apt-get install -y build-essential gcc-14 g++-14 ninja-build cmake git wget ca-certificates \
    libglfw3-dev libglm-dev libtbb-dev libgl-dev dpkg-dev gh

echo "### CUDA ${CUDA_PACKAGE_VERSION/-/.} toolkit for WSL"
# NVIDIA's WSL repository ships the toolkit without a driver; WSL uses the Windows driver
if [ ! -x "$CUDA_DIR/bin/nvcc" ]; then
    keyring=$(mktemp --suffix=.deb)
    wget -q -O "$keyring" https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-keyring_1.1-1_all.deb
    sudo dpkg -i "$keyring"
    rm "$keyring"
    sudo apt-get update
    sudo apt-get install -y "cuda-toolkit-$CUDA_PACKAGE_VERSION"
fi
"$CUDA_DIR/bin/nvcc" --version | tail -1

echo "### Checks"
nvidia-smi --query-gpu=name,driver_version --format=csv,noheader \
    || { echo "The GPU is not visible in WSL. Update the NVIDIA driver on Windows."; exit 1; }

echo "### Docker Engine, for the clean-system package tests"
# Runs inside this distro, so Docker Desktop on Windows is not needed. Requires systemd, the default in recent WSL
if [ "$(ps -p 1 -o comm=)" != "systemd" ]; then
    echo "systemd is not enabled in this distro. Add the following to /etc/wsl.conf, run 'wsl --shutdown' in Windows,"
    echo "and run this script again:"
    printf '  [boot]\n  systemd=true\n'
    exit 1
fi
sudo apt-get install -y docker.io
sudo systemctl enable --now docker
if ! id -nG | grep -qw docker; then
    sudo usermod -aG docker "$USER"
    echo "Added $USER to the docker group. Run 'wsl --shutdown' in Windows, then run this script again."
    exit 1
fi
# /usr/bin/docker must win over Docker Desktop's Windows docker.exe, which WSL may put on PATH
[ "$(command -v docker)" = "/usr/bin/docker" ] || { echo "'docker' resolves to $(command -v docker), expected /usr/bin/docker"; exit 1; }
docker version --format 'Docker {{.Server.Version}}' \
    || { echo "Docker Engine is installed but not reachable. Try 'sudo systemctl restart docker'"; exit 1; }

if ! gh auth status > /dev/null 2>&1; then
    echo
    echo "### Log in to GitHub, to upload releases"
    gh auth login --hostname github.com --git-protocol https --web
fi

echo
echo "WSL is ready to build LIMA releases"
