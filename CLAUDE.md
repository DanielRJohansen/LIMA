# LIMA

## Coding rules

- Never use `std::format`. Use `Lima::Format` from `Format.h` instead; `std::format` costs ~3 s of compile time per translation unit (see the header comment). In `.cu` files avoid both, since nvcc can hit an internal compiler error on `<format>`; build strings with `std::to_string` there.

## Normal workflow

- Push back! If im asking you to do something, and there's a change i might not have thought through all consequences or you simply think it's a bad idea, push back and explain why. 
- If a given task is simple go ahead and implement it. Otherwise lay out a plan first for how you will solve and validate.
- Never commit. When a task is done, give a short commit message; the user reviews and commits. 
- If we're in the middle of a session with multiple items, also lay out the next item and how you'd like to proceed after writing the commit.

## Profiling & benchmarking

Workflow for GPU performance work: measure, change one isolated thing, measure again, validate.

- **Benchmark system**: STMV (`LIMA_data/benchmarking/stmv`). `limaprofile.exe --stmv <steps>` runs it headless and prints ms/step (exit code 1 just means outside the benchmark's allowed window).
  `--stmv-constrained <steps>` runs the same system with constraints = h-bonds (settled water), which takes different paths (the per-cell constrained integration, SETTLE/SHAKE). Measure it too when a change touches clustering, superclusters or integration, to catch regressions on constrained systems.
- **Build from the shell**: a .bat that calls `"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"`, then `cd /d build\x64-Release` and runs the VS-bundled ninja (`Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja\ninja.exe`) with targets `limaprofile limatest agenttesting`.
- **Keep a baseline exe**: before a change, copy `build/x64-Release/code/LIMA_TESTS/limaprofile.exe` to `limaprofile_<name>.exe`, so before/after can be compared back to back. To build a baseline of HEAD while having local changes: `git stash`, build, copy, `git stash pop`.
- **Wall-clock**: run baseline and new exe alternately, at least 2 rounds each (`--stmv 300`). GPU clocks/thermals drift ~10-20% between sessions, so only compare runs from the same session. 300-step runs are noisy (+-5%); for differences below ~3% use back-to-back `--stmv 3000` pairs. Changes to the neighbor-list/task builders only show up in wall-clock (rebuild cost is amortized over steps), not in a per-kernel ncu comparison.
- **Energy steps**: the `computePotE` kernel variants run every `nstenergy` (default 5) steps, so 20% of steps. Measure both variants; a change can speed up normal steps while slowing energy steps.
- **ncu** (Nsight Compute, call `C:/Program Files/NVIDIA Corporation/Nsight Compute <ver>/target/windows-desktop-win7-x64/ncu.exe` directly; `ncu.bat` breaks on `|` in kernel regexes). Use `--stmv 3`: step 0 is a logging step (potE computed), steps 1-2 are normal steps. limaprofile exits 1 under ncu (outside its time window); the report is still written.
  - Full report of all kernels: `ncu --set full -f -o build/x64-Release/ncu/<name> limaprofile_<name>.exe --stmv 3`
  - Quick per-kernel numbers: `ncu --metrics gpu__time_duration.sum,smsp__inst_executed.sum,launch__registers_per_thread,dram__bytes_read.sum,dram__bytes_write.sum -k "regex:NbNonlocal|SuperclusterIntegrate" ...`, then `ncu -i <rep> --page raw --csv`.
  - ncu locks clocks, so prefer its kernel durations over wall-clock for small differences.
  - Hotspots: `ncu -i <rep> --page source --csv --print-source sass -k regex:<kernel> --launch-skip 1 --launch-count 1` gives per-SASS-instruction exec counts, stall samples and L1/L2 requests per instruction. Note a barrier wait is attributed to the instruction *after* the `BAR.SYNC`.
  - First check what bounds the kernel (SpeedOfLight / SchedulerStats: issue slots busy vs DRAM throughput) before optimizing; latency fixes don't help an issue-bound kernel.
  - When "Memory Throughput" is high, find which unit: `lts__throughput`, `lts__xbar2lts_cycles_active`, `dram__throughput`, `l1tex__throughput`, and `lts__t_requests_op_{read,red}`. For global atomics the cost is the request count, not the bytes; the LSU merges the lanes of one instruction that hit the same 128B line.
  - To find the ceiling of an optimization cheaply, measure a throwaway profiling-only hack first (e.g. skip the atomics while keeping the computation alive), before designing the real change.
- **Validation** after every change:
  - `agenttesting.exe --engine-batch-tests` (batched vs standalone, bitwise)
  - `limatest.exe` full suite (RunAllUnitTests, takes no args): "Deterministic Simulations", force sanity tests, and the VC / energy-gradient drift targets.
  - VC/drift targets (`tests/automatedtests/vc_targets.csv`) are calibrated to exact numerics and exist to flag physics changes. Before a change, predict the outcome: no physics and no summation-order change => VC and gradients bitwise identical (any difference is a red flag); summation order only => small noise; physics/approximation change => a *significant* VC shift is a red flag. Compare exact values in `tests/automatedtests/vc_results.csv` (written by every limatest run, overwritten each time) against a run of the baseline `limatest_<name>.exe`; the console shows only VC and only the first failure. Never edit vc_targets.csv or benchmark targets, that is the user's call.
  - `agenttesting.exe --energy-minimization-tests` with `LIMA_EM_LABEL=<name>` when the EM path may be affected.
  - For changes that could destabilize, also run a longer STMV (`--stmv 3000`).
- **Determinism** is required: no float atomics. Deterministic accumulation uses 64-bit fixed-point integer atomics (see `ForceAccumulator`).
- **HIP backend** (`LIMA_GPU_BACKEND=HIP`, AMD RDNA only, so 32-wide wavefronts): CUDA is the reference and the only tuned path. HIP must compile and run correctly; it is never required to be fast, so a performance change is never blocked on AMD tuning, only on giving HIP a correct fallback.
  - No vendor code in kernel files: no `#if defined(LIMA_HIP)`, no inline `asm`/PTX, no raw `0xFFFFFFFFu` warp masks. Use `WarpMask`/`LIMA_FULL_WARP_MASK`, `LIMA_SHARED_UNINIT` (for `__shared__` arrays of types with member initializers, e.g. `Float3`), and `DeviceIntrinsics` (`platform/DeviceIntrinsics.cuh`, `platform/BlockCopy.cuh`). A new CUDA-only feature (PTX, cache hints, `cp.async`, ...) gets a `DeviceIntrinsics` wrapper with a simple, correct HIP fallback. Vendor `#if`s live only in `platform/`, `GpuRuntime.h` and `Engine.cu`; check with `grep -rn "LIMA_HIP\|asm(" code --include=*.cu --include=*.cuh`.
  - Data written by another block during the same kernel must be read with `DeviceIntrinsics::LoadCoherent`, never a plain load: AMD's per-CU L1 is not coherent, so a plain load works on NVIDIA and silently reads stale data on AMD.
  - Shim changes must leave the CUDA code unchanged: when touching `platform/` or `GpuRuntime.h`, diff `cuobjdump -sass LIMA_ENGINE.lib` before and after; it should be identical.
  - HIP is not built in the normal optimization loop (the AMD machine is a laptop). It is built and tested in batches, before merging dev_dj into main, so dev_dj may be HIP-broken in between. Its numerics differ (different fast-math flags), so the CUDA VC/drift targets need not hold exactly on HIP.
- **Always finish performance work with a full ncu report** of the final build for the user to inspect before they commit: `ncu --set full -f -o build/x64-Release/ncu/<name> limaprofile_<name>.exe --stmv 3`, plus one of the baseline exe when comparing. Give the report paths in the final summary, and don't delete them in cleanup.
