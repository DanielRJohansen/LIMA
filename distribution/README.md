# Releasing LIMA

One command builds, tests and packages LIMA for Windows and Linux, and uploads it as a draft GitHub release:

```
distribution\release.bat            # the real thing
distribution\release.bat --dry-run  # everything except the upload, allows uncommitted changes
distribution\release.bat --upload-only  # retry only the upload, using the files from the last run
```

To release, bump `project(lima VERSION ...)` in the top-level `CMakeLists.txt`, commit, push, and run the script.
When it finishes, review the draft on GitHub and publish it. Publishing creates the `v<version>` tag.

## What it does

1. Checks that the tree is clean, the commit is pushed, and the version is not released yet.
2. Windows: builds a Release for all supported GPUs, packages `lima-<version>-windows-x64.zip`, and runs
   `limaclitest` against the packaged `lima.exe`.
3. Linux, in WSL: builds the same commit, packages `lima-<version>-linux-x86_64.tar.gz`, and smoke tests it on the GPU.
4. Builds `lima_<version>_amd64.deb` and test-installs it in a clean Ubuntu container.
5. Writes a `PKGBUILD` for the tarball and test-builds and installs it in a clean Arch container.
6. Uploads all of it as a draft release, with `release-notes.md` and GitHub's generated changelog.

Everything ends up in `distribution\out\<version>`.

## One-time setup

- Visual Studio 2022 with the C++ workload, and the CUDA toolkit set in `linux/config.sh` (currently 13.2).
- In WSL, run `wsl -d Ubuntu-24.04 -- bash /mnt/c/<repo>/distribution/linux/setup-wsl.sh`. It installs the Linux
  build tools, CUDA and Docker Engine (inside WSL, Docker Desktop is not needed), and logs `gh` in to GitHub.
  WSL must have systemd enabled, which is the default in recent versions.

## Supported GPUs

Release builds target RTX 40-series (sm_89), H100 (sm_90), B200 (sm_100) and RTX 50-series (sm_120), plus PTX
for newer GPUs. The list is in `release.bat` and `linux/config.sh`. Dev builds only target sm_89 to compile faster.
