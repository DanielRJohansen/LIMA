# Shared settings for the Linux release scripts

GITHUB_REPO="DanielRJohansen/LIMA"
MAINTAINER="Daniel Johansen <daniel@lima-dynamics.com>"
DESCRIPTION="The Molecular Dynamics engine for the Generative Era"

# Keep CUDA in sync with setup-wsl.sh and the Windows release build
CUDA_DIR="/usr/local/cuda-13.2"
# RTX 40 (89), H100 (90), B200 (100), RTX 50 (120), plus PTX of the newest so future GPUs can JIT it.
# Keep in sync with distribution/release.bat
CUDA_ARCHITECTURES="89-real;90-real;100-real;120"
