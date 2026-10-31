#!/bin/bash
# Builds LIMA from this checkout and installs it to /usr/bin/lima and /usr/share/LIMA
# Ubuntu/Debian. CUDA 13 must already be installed in /usr/local/cuda
set -euo pipefail
cd "$(dirname "$0")"

sudo apt-get install -y gcc-14 g++-14 cmake libglm-dev libglfw3-dev libtbb-dev libgl-dev \
    || echo "apt-get failed, continuing; cmake reports any dependency that is actually missing"

cmake -S . -B build/linux-release -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=gcc-14 -DCMAKE_CXX_COMPILER=g++-14 \
    -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc -DCMAKE_CUDA_HOST_COMPILER=g++-14
cmake --build build/linux-release --parallel
sudo cmake --install build/linux-release
