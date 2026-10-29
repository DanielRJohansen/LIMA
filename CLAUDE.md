# LIMA

## Coding rules

- Never use `std::format`. Use `Lima::Format` from `Format.h` instead; `std::format` costs ~3 s of compile time per translation unit (see the header comment). In `.cu` files avoid both, since nvcc can hit an internal compiler error on `<format>`; build strings with `std::to_string` there.

## Normal workflow

- Push back! If im asking you to do something, and there's a change i might not have thought through all consequences or you simply think it's a bad idea, push back and explain why. 
- If a given task is simple go ahead and implement it. Otherwise lay out a plan first for how you will solve and validate.
- Never commit. When a task is done, give a short commit message; the user reviews and commits. 
- If we're in the middle of a session with multiple items, also lay out the next item and how you'd like to proceed after writing the commit.

## Performance work

A performance change is accepted only on numbers: a baseline and a candidate exe, measured interleaved in the same session on all three system sizes, plus the full validation. One isolated change per measurement.

### The loop

1. **Hypothesis.** Before writing code, state what you expect to gain, on which size (small/medium/large), and why (which kernel, which bound). Profile first (ncu/nsys, below) to confirm the target is where the time goes and what bounds it.
2. **Ceiling.** If the change is non-trivial, first measure a throwaway profiling-only hack that removes the cost (e.g. skip the atomics but keep the computation alive). If the ceiling is below the noise, drop the idea.
3. **Baseline exe.** Before changing anything, build and copy `build/x64-Release/code/LIMA_TESTS/limaprofile.exe` (and `limatest.exe`) to `limaprofile_<name>.exe` / `limatest_<name>.exe`. Check `git status` first; to build a baseline of HEAD with local changes present: `git stash`, build, copy, `git stash pop`.
4. **Change one thing**, build, run the validation (below) before measuring, so time isn't spent measuring a broken build.
5. **Measure** with the suite protocol below, then decide with the acceptance rules.
6. **Record** the result, positive or negative, in the summary (and in memory for rejected ideas, with the numbers, so they are not retried).

### Systems: always three sizes

`limaprofile.exe --suite [runs=3]` runs all three and prints the fastest run's ms/step and ns/day for each (~25 s). All three come from `LIMA_data/salad_gromacs_bench` (the GROMACS benchmark set) and share one `lima.mdp`: PME, h-bond constraints, v-rescale, 2 fs, nstenergy 100, so the numbers are directly comparable with `lima gromacscompare`.

| size | system | atoms | ~ms/step | regime |
|---|---|---|---|---|
| small | villin | 5k | 0.058 | launch/latency-bound: kernel launches, host work, syncs, small-kernel tails |
| medium | adh_dodec | 191k | 0.50 | mixed |
| large | stmv | 1.07M | 3.3 | throughput-bound: DRAM, issue slots, atomics |

A change aimed at one regime often costs in another (more launches help STMV and hurt villin; bigger blocks do the reverse), which is why every change is measured on all three. `--salad <name> <steps>` runs one system once without warmup (for ncu/nsys); `rnase_cubic` (24k) is also available.

Additional system, when relevant:
- `--stmv <steps>`: the historical STMV benchmark *without* constraints or thermostat (`LIMA_data/benchmarking/stmv`, nstenergy 20). It takes the unconstrained integration path the suite never exercises; measure it when a change touches clustering, superclusters or integration. Exit code 1 only means outside the benchmark's allowed window. `--stmv-constrained` is the same with h-bonds.

### Measurement protocol

- Run baseline and candidate **interleaved** (A B A B A B), at least 3 rounds each of `--suite`, all in the same session. Never compare against numbers from an earlier session: GPU clocks/thermals drift 10-20% between sessions.
- Noise of the fastest-of-3 per suite run: ~1% on villin and adh_dodec, ~3% on stmv. For an effect smaller than that, use more rounds or longer runs (`--salad stmv 3000` pairs).
- Report a table per system: baseline and candidate (median of the per-round values, and min-max), and the delta in %.
- Neighbor-list/task-builder changes only show in wall-clock (rebuild cost is amortized over steps), not in a per-kernel ncu comparison.
- Energy steps: the `computePotE` kernel variants run every `nstenergy` steps (1% of suite steps, 20% in `--stmv`, 20% at the engine default of 5). Check both variants in ncu; a change can speed up normal steps while slowing energy steps.

### Acceptance

- Accept: a gain beyond the noise on at least one size, and no size regressing beyond its noise.
- A real regression on any size (consistently, across rounds) blocks the change, or must be called out explicitly with the tradeoff for the user to decide. Never hide it in an average.
- A gain within the noise is not a result; say "neutral", and only keep the change if it simplifies the code.

### Profiling tools

- **ncu** (Nsight Compute; call `C:/Program Files/NVIDIA Corporation/Nsight Compute <ver>/target/windows-desktop-win7-x64/ncu.exe` directly; `ncu.bat` breaks on `|` in kernel regexes). Run 3 steps: step 0 is a logging step (potE computed), steps 1-2 are normal steps. Under ncu, `--stmv` exits 1 (outside its time window); the report is still written.
  - Full report: `ncu --set full -f -o build/x64-Release/ncu/<name>_<system> limaprofile_<name>.exe --salad <system> 3`
  - Quick per-kernel numbers: `ncu --metrics gpu__time_duration.sum,smsp__inst_executed.sum,launch__registers_per_thread,dram__bytes_read.sum,dram__bytes_write.sum -k "regex:NbNonlocal|SuperclusterIntegrate" ...`, then `ncu -i <rep> --page raw --csv`.
  - ncu locks clocks and serializes kernels, so prefer its kernel durations over wall-clock for small per-kernel differences. It says nothing about launch gaps, so it is the wrong tool for villin-sized problems.
  - First check what bounds the kernel (SpeedOfLight / SchedulerStats: issue slots busy vs DRAM throughput) before optimizing; latency fixes don't help an issue-bound kernel.
  - Hotspots: `ncu -i <rep> --page source --csv --print-source sass -k regex:<kernel> --launch-skip 1 --launch-count 1` gives per-SASS-instruction exec counts, stall samples and L1/L2 requests per instruction. A barrier wait is attributed to the instruction *after* the `BAR.SYNC`.
  - When "Memory Throughput" is high, find which unit: `lts__throughput`, `lts__xbar2lts_cycles_active`, `dram__throughput`, `l1tex__throughput`, and `lts__t_requests_op_{read,red}`. For global atomics the cost is the request count, not the bytes; the LSU merges the lanes of one instruction that hit the same 128B line.
- **nsys** (Nsight Systems, `C:/Program Files/NVIDIA Corporation/Nsight Systems 2025.6.3/target-windows-x64/nsys.exe`) for the small system and anything launch-, sync- or host-bound: `nsys profile -o build/x64-Release/nsys/<name>_villin limaprofile_<name>.exe --salad villin 2000`. Look at GPU idle gaps between kernels, kernel count per step, and host-side API time.

### Validation (after every change)

- `agenttesting.exe --engine-batch-tests` (batched vs standalone, bitwise)
- `limatest.exe` full suite (RunAllUnitTests, takes no args): "Deterministic Simulations", force sanity tests, and the VC / energy-gradient drift targets.
- VC/drift targets (`tests/automatedtests/vc_targets.csv`) are calibrated to exact numerics and exist to flag physics changes. Before a change, predict the outcome: no physics and no summation-order change => VC and gradients bitwise identical (any difference is a red flag); summation order only => small noise; physics/approximation change => a *significant* VC shift is a red flag. Compare exact values in `tests/automatedtests/vc_results.csv` (written by every limatest run, overwritten each time) against a run of the baseline `limatest_<name>.exe`; the console shows only VC and only the first failure. Never edit vc_targets.csv or benchmark targets, that is the user's call.
- `agenttesting.exe --energy-minimization-tests` with `LIMA_EM_LABEL=<name>` when the EM path may be affected.
- For changes that could destabilize, also run a long STMV (`--salad stmv 3000` and `--stmv 3000`).
- **Determinism** is required: no float atomics. Deterministic accumulation uses 64-bit fixed-point integer atomics (see `ForceAccumulator`).

### Building

A .bat that calls `"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"`, then `cd /d build\x64-Release` and runs the VS-bundled ninja (`Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja\ninja.exe`) with targets `limaprofile limatest agenttesting`.

### Finishing

Every performance task ends with:
- The before/after table for all three sizes (plus `--stmv` if measured), with rounds, medians and spread.
- Full ncu reports of the final build and of the baseline for the system(s) the change targets (default stmv): `build/x64-Release/ncu/<name>_<system>`; an nsys report instead when the target is villin. Give the paths in the summary and don't delete them in cleanup.
- The validation results, and the VC comparison against the baseline matching the prediction.

## HIP backend

`LIMA_GPU_BACKEND=HIP`, AMD RDNA only, so 32-wide wavefronts. CUDA is the reference and the only tuned path. HIP must compile and run correctly; it is never required to be fast, so a performance change is never blocked on AMD tuning, only on giving HIP a correct fallback.
- No vendor code in kernel files: no `#if defined(LIMA_HIP)`, no inline `asm`/PTX, no raw `0xFFFFFFFFu` warp masks. Use `WarpMask`/`LIMA_FULL_WARP_MASK`, `LIMA_SHARED_UNINIT` (for `__shared__` arrays of types with member initializers, e.g. `Float3`), and `DeviceIntrinsics` (`platform/DeviceIntrinsics.cuh`, `platform/BlockCopy.cuh`). A new CUDA-only feature (PTX, cache hints, `cp.async`, ...) gets a `DeviceIntrinsics` wrapper with a simple, correct HIP fallback. Vendor `#if`s live only in `platform/`, `GpuRuntime.h` and `Engine.cu`; check with `grep -rn "LIMA_HIP\|asm(" code --include=*.cu --include=*.cuh`.
- Data written by another block during the same kernel must be read with `DeviceIntrinsics::LoadCoherent`, never a plain load: AMD's per-CU L1 is not coherent, so a plain load works on NVIDIA and silently reads stale data on AMD.
- Shim changes must leave the CUDA code unchanged: when touching `platform/` or `GpuRuntime.h`, diff `cuobjdump -sass LIMA_ENGINE.lib` before and after; it should be identical.
- HIP is not built in the normal optimization loop (the AMD machine is a laptop). It is built and tested in batches, before merging dev_dj into main, so dev_dj may be HIP-broken in between. Its numerics differ (different fast-math flags), so the CUDA VC/drift targets need not hold exactly on HIP.
