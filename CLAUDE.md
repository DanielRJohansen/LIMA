# LIMA

## Coding rules

- Never use `std::format`. Use `Lima::Format` from `Format.h` instead; `std::format` costs ~3 s of compile time per translation unit (see the header comment). In `.cu` files avoid both, since nvcc can hit an internal compiler error on `<format>`; build strings with `std::to_string` there.

## Profiling & benchmarking

Workflow for GPU performance work: measure, change one isolated thing, measure again, validate.

- **Benchmark system**: STMV (`LIMA_data/benchmarking/stmv`). `limaprofile.exe --stmv <steps>` runs it headless and prints ms/step (exit code 1 just means outside the benchmark's allowed window).
- **Keep a baseline exe**: before a change, copy `build/x64-Release/code/LIMA_TESTS/limaprofile.exe` to `limaprofile_<name>.exe`, so before/after can be compared back to back. To build a baseline of HEAD while having local changes: `git stash`, build, copy, `git stash pop`.
- **Wall-clock**: run baseline and new exe alternately, at least 2 rounds each (`--stmv 300`). GPU clocks/thermals drift ~10-20% between sessions, so only compare runs from the same session.
- **ncu** (Nsight Compute, call `target/windows-desktop-win7-x64/ncu.exe` directly; `ncu.bat` breaks on `|` in kernel regexes). Use `--stmv 3`: step 0 is a logging step (potE computed), steps 1-2 are normal steps.
  - Full report of all kernels: `ncu --set full -f -o build/x64-Release/ncu/<name> limaprofile_<name>.exe --stmv 3`
  - Quick per-kernel numbers: `ncu --metrics gpu__time_duration.sum,smsp__inst_executed.sum,launch__registers_per_thread,dram__bytes_read.sum,dram__bytes_write.sum -k "regex:NbNonlocal|SuperclusterIntegrate" ...`, then `ncu -i <rep> --page raw --csv`.
  - ncu locks clocks, so prefer its kernel durations over wall-clock for small differences.
  - Hotspots: `ncu -i <rep> --page source --csv --print-source sass -k regex:<kernel> --launch-skip 1 --launch-count 1` gives per-SASS-instruction exec counts and stall samples. Note a barrier wait is attributed to the instruction *after* the `BAR.SYNC`.
  - First check what bounds the kernel (SpeedOfLight / SchedulerStats: issue slots busy vs DRAM throughput) before optimizing; latency fixes don't help an issue-bound kernel.
- **Validation** after every change:
  - `agenttesting.exe --engine-batch-tests` (batched vs standalone, bitwise)
  - `limatest.exe` full suite (RunAllUnitTests): "Deterministic Simulations", force sanity tests, and the VC / energy-gradient drift targets. VC/drift targets are calibrated to exact numerics, so rounding-order changes can trip them; judge whether stability got *significantly* worse rather than requiring a pass.
  - `agenttesting.exe --energy-minimization-tests` with `LIMA_EM_LABEL=<name>` when the EM path may be affected.
  - For changes that could destabilize, also run a longer STMV (`--stmv 3000`).
- **Determinism** is required: no float atomics. Deterministic accumulation uses 64-bit fixed-point integer atomics (see `NbForceAccumulator`).
