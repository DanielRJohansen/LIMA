# LIMA TODO

Sections:
1. Priority: GROMACS parity (Salad benchmark)
2. Missing MD features
3. Correctness issues
4. Performance
5. Tooling, CLI and codebase
6. Ideas and research
7. AMD GPU support (HIP)

---

# 1. Priority: GROMACS parity (Salad benchmark)

## Villin (Salad, AMBER99SB) vs GROMACS 2024.2: remaining issues (2026-10-14)
- **Now matching:** bonded terms, LJ, 1-4 pairs (fudgeLJ/fudgeQQ), periodic impropers and electrostatics (1-4 Coulomb, Ewald exclusion correction).
  - Forces to 1e-4, electrostatic forces to 0.04% RMS.
  - Coulomb energy -4779.8 vs -4780.9 kJ/mol.
- **Left:** fixture `tests/automatedtests/ForceComparisons/VillinAmber` has zero charges. A charged fixture needs an RMS rather than a per-atom criterion, as PME differs by up to 0.4% per atom.
- Smaller issue (exclusions) found on the way: see [Correctness issues](#3-correctness-issues).

## Compare LIMA with GROMACS on every system in the Salad Benchmark .zip
Compare both outcome and speed.
- **Outcome:** trajectories diverge within ps (chaotic dynamics, different rounding), so compare ensemble statistics, not frames:
  - backbone RMSD vs the start structure, per-residue RMSF, radius of gyration
  - average potential energy and temperature, density for solvated boxes
  - water O-O RDF, H-bond counts, secondary structure (DSSP) over time
  - Run each system with a few seeds in both engines, so a difference can be told apart from run-to-run spread.
  - Use gmx rms/rmsf/gyrate/rdf/energy on both trajectories (LIMA writes .trr), so the analysis is identical.
- **Performance:** Salad lists GROMACS ns/day, but measure GROMACS on this PC too (WSL, GPU build).
  - Match settings: dt, cutoff, PME, nstlist, constraints, output frequency.
  - LIMA via `lima mdrun`, which prints the same ns/day table.
- **Blockers** observed in single 100 ps Villin/RNase runs (2026-10-13):
  - Once timestamps and boxes are correct, run the common RMSD/RMSF/Rg/RDF/H-bond/DSSP analyses on one LIMA and one GROMACS trajectory per system before adding replicate ensembles.

---

# 2. Missing MD features
GROMACS has these, LIMA does not (status checked 2026-10-13). Check each against what the Salad comparison actually needs, and do those first.

## Topology / input
- **`[ constraints ]` sections:** not read yet. constraints = h-bonds and `[ settles ]` are handled (see Constraints below), but topologies that list their constraints explicitly lose them.
- **MDP validation:** explicitly validate/reject unsupported settings instead of silently ignoring them: `rvdw`, PME grid settings, neighbor-list behavior, constraints, thermostat parameters, output cadence, seeds.
- **Regression coverage:** add a small AMBER99SB/TIP3P constrained-water fixture with tests for topology loading, constraint preservation, short-run stability, and readable `.trr` output. Validate against GROMACS observables, not frame-by-frame coordinates.

## Constraints
**Done** for constraints = h-bonds (2026-10-13):
- Bonds to hydrogens, and 3-atom water made rigid via its H-H distance.
- Groups (<=4 atoms, <=3 constraints) are kept within one pcluster; `SuperclusterIntegrateConstrainedKernel` solves them per pcluster with SHAKE/RATTLE (deterministic, one thread per pcluster).
- Constrained bonds stay in the bond lists (0 energy at the constrained length), so EM and MD can share a box.
- T4 in water: constraint error < 1e-5, energy conservation better than flexible at 2 fs.
- SETTLE for water done 2026-10-13 (constrained integrate kernel 36 -> 15 us on T4 in water).

**Left:**
- **LINCS / constraints = all-bonds:** postponed, h-bonds covers the usual CHARMM/AMBER setups. Needs coupled constraints across pclusters, ie. across superclusters, so a separate deterministic solver.
- **Position precision limits all constraints:** positions are absolute float coordinates in the box, so a 0.1 nm bond holds only to ~6e-6 relative whatever the solver. SETTLE assumes the water is exactly rigid before each step, so it re-snaps the rounding every step: T4 VC over 10k steps 1.50e-5 vs 1.33e-5 with SHAKE. Fix: store positions relative to the supercluster (or cell), as GROMACS keeps coordinates local. Would also help everything else.
- **Hydrogen mass repartitioning** for 4-5 fs steps. With constraints = h-bonds, T4 (vacuum and water) already runs 2000 steps at 4 fs without it; check energy conservation there first.
- **initial_temperature** draws velocities with components along the constraints, which the first step removes, so the start is somewhat colder than asked.
- Not yet compared with a constrained GROMACS run (bond lengths, water RDF, energy drift).

## Thermostats
- **Done:** v-rescale (Bussi), 2026-10-13. Deterministic (seeded, own deviates). T4 in water: mean 299.8 K, fluctuations 1.598 K vs 1.600 canonical.
- **Compatibility check:** Salad uses `tcoupl = v-rescale`, `tau_t = 0.1`, `ref_t = 300`. Compare LIMA's behavior with GROMACS v-rescale and make compatible if needed; it depends on the kinetic-energy convention used by the integrator.
- **Missing:** Nose-Hoover, Langevin/stochastic dynamics, per-group coupling (eg. protein vs solvent).

## Constant pressure (NPT)
- No barostat: C-rescale/Berendsen for equilibration, Parrinello-Rahman/MTTK for production.
- Needs box rescaling on the GPU, a virial (the fixed-point force accumulation would also have to sum pair virials), and the grid/PME to follow the box size.

## Restraints and pulling
- Position restraints are skipped (posre includes ignored in `Topologyfile.cpp`), so no restrained equilibration.
- Missing: distance/dihedral/flat-bottom restraints, pulling (umbrella sampling, steered MD).

## Other
- Virtual sites (TIP4P water, hydrogen vsites for 4-5 fs steps)
- Free-energy perturbation (lambda coupling)
- Center-of-mass motion removal
- Energy groups
- .edr-style energy output

---

# 3. Correctness issues

## Forcefields with different 1-4 scaling can't be mixed
- A topology has one `[ defaults ]` (gen-pairs, fudgeLJ, fudgeQQ), applied to the whole system.
- So merging eg. a CHARMM-style protein (fudge 1.0/1.0, as in STMV) with Slipids lipids (0.5/0.8333) gives one of them the wrong 1-4 interactions.
- `SimulationBuilder::InsertSystem` merges forcefields anyway (used by the large-EM showcase).
- Fix: per-moleculetype 1-4 scaling, or explicit `[ pairs ]` parameters with the scaling baked in.

## Exclusions differ from GROMACS in edge cases
- Exclusions come from the listed bonds/angles/dihedrals/impropers, not GROMACS' bond graph with nrexcl.
- Same result for normal topologies, but differs without angle lines, or for impropers spanning atoms more than 3 bonds apart.

## Pair list buffer: measure energy drift
Measure whether the default 0.05 nm buffer causes real drift before doing anything about it.
- GROMACS' buffer estimate (`verlet_buffer_tolerance`, implemented, off by default) asks for ~0.11-0.12 nm at nstlist 20 and 2 fs. That makes STMV ~11% slower, and by that (conservative) estimate 0.05 nm drifts well beyond GROMACS' 0.005 kJ/mol/ps per atom.
- Test: run eg. T4 in water for 100k+ steps with the default buffer and with `verlet_buffer_tolerance=0.005`, compare energy drift. Only pursue if 0.05 nm drifts noticeably.
- **Tried 2026-10-13 and reverted:** dynamic pruning as GROMACS does on GPUs (outer list with the accurate buffer every nstlist steps, a prune kernel every few steps recomputing each QuarterEntry's ownQuarterMask against an inner radius). Doesn't pay off in LIMA:
  - NbNonlocal costs the same at inner radius cutoff+0.013 as at cutoff+0.05 (2.25 ms on STMV), as nearly every 4x4 block has a pair within the cutoff itself. So pruning only removes the outer buffer's cost, while the prune kernel (0.75 ms, one lane per entry) and the larger search cost more: STMV 4.04 vs 3.84 ms/step at nstlist 20 / prune 4.
  - Longer nstlist doesn't help, as the outer buffer grows fast (0.29 nm at 40, beyond the cutoff at 80).
  - Only worth revisiting with a much cheaper prune kernel (eg. NB's 4 lanes per entry) and if the drift measurement shows a need.

---

# 4. Performance

## 3J3Q benchmarks got slower; limits raised to keep them green (2026-10-14)
Profile both, then put the limits back.
- **"3j3q load benchmark"** build step ~5.4 -> 7.5-8.9 s (limit 7.0 -> 10.0). Likely the Ewald exclusion-correction pair entries (~3x the pair entries) going through BondGroupFactory, plus building the exclusion lists per moleculetype in SuperTopology.
- **"ToGmx large CIF benchmark"** ~8.5 -> 10.2 s (limit 10 -> 12). Cause not investigated (ToGmx builds no box, so maybe the topology parser changes, or noise).

## Startup latency for small systems
Investigate time from command to first step, eg. T4 in a 7^3 nm^3 box.

## Kernel-level ideas
Not set-in-stone todos; concepts worth exploring. Some make the most sense together.
- **Engine on its own thread:** instead of Env having a thread dedicated to RunSimulation, Engine gets its own thread, so an mdrun only steps and updates the renderpipeline/simstatus. Engine would persist, so slightly lower setup latency. Engine must free its memory when idle, so we're not hogging 10 GB VRAM while Env just does preprocessing.
- **EmitQuarterEntriesKernel:** now ~0.77 ms per STMV rebuild (was 1.37). The cost was the bonded-pair exclusions (~1 active lane, dependent global loads): own particle ids now in shared memory, query ids loaded once per quarter, one pass over each bonded set for all 4 query particles. Rejected: int4 reads of the set (slower). What's left is mostly writing ~300 MB of QuarterEntries; shrinking that means changing NB's entry format.
- **FindNeighborsKernel:** now ~2.7 ms per STMV rebuild (was 3.9).
  - Done: own-particle vs candidate-sphere prefilter (196 -> 140 pair-tested candidates/SC, 90 become neighbors), branch-free pair test, `__launch_bounds__(128, 8)` for 64 regs.
  - Issue-bound (~68% issue active), so latency tricks don't help: 2 candidates per warp iteration and a candidate-particles vs own-sphere filter were both slower.
  - Left: cell enumeration (~74 cells/SC, ~7 of 32 lanes active) and the bonded check (~8%).
- **nstlist vs list-buffer tuning:** a bigger buffer and fewer rebuilds may be cheaper overall. A parameter/accuracy decision.
- **NB kernel flattened entries:** known ~4% SIMT gain left over from densetasks.

## Kernel graphs (prototyped, see stash "Tried cuda graph...")
Result: T4 -11%, STMV +2.5%.

**Recipe:**
- `Engine.cuh`: `std::array<cudaGraphExec_t, 2> stepGraphs{}` (indexed by logData) + `void LaunchStepGraph(cudaGraph_t, bool logData)`. Destroy the execs in ~Engine.
- `_deviceMaster`: after the `MakeSuperClusterTasksGPU` check and before `ForkStreamsFromMainStream()`, if `!emvariant && nSuperclusters > 0 && snf_select.empty()` (EM syncs mid-step, SNF unchecked): `cudaStreamBeginCapture(cudaStreams[0], cudaStreamCaptureModeThreadLocal)`. ThreadLocal so the log-drain thread's CUDA calls don't break the capture.
- The existing fork/join events pull the side streams (PME, bonds) into the capture; cuFFT exec captures fine. After the integrate launch: `cudaStreamEndCapture(cudaStreams[0], &graph)` -> LaunchStepGraph.
- LaunchStepGraph: if an exec exists try `cudaGraphExecUpdate(exec, graph, &info)`; on failure (topology changed after a task rebuild) `cudaGetLastError()`, destroy and re-instantiate. Then `cudaGraphLaunch(exec, cudaStreams[0])`, `cudaGraphDestroy(graph)`. Capturing every step means changed args (step, buffer pointers after Expand, live edits) need no tracking; ~10 us host/step, only ~4 re-instantiations per 2000 T4 steps.

**Measuring:** nsys needs `--cuda-graph-trace=node` to show graph kernels, but node tracing perturbs scheduling (STMV looked faster under nsys, was slower in wall-clock). Trust wall-clock A/B (`--t4 2000 1`, `--stmv 300`) of the same binary.

**Findings:**
- Launch overhead is not the issue (host runs 0.25-41 ms ahead).
- Graphs fix WDDM submit batching that idled T4 normal steps: 24/90 us -> 54 us steps.
- STMV loses because NB, PME spreading and bonds truly overlap and slow NB, which already fills the GPU (CUDA_DEVICE_MAX_CONNECTIONS=1: graph == no graph).
- Reordering launches (NB before PME) without graphs changed nothing.

**Next:** encode the order (NB waits for an event after PME's ChargeblockDistributeToGrid, bonds wait for NB) so STMV is neutral and T4 keeps most of the gain; or use graphs only below a system size.

---

# 5. Tooling, CLI and codebase

## Compile times
Compiling takes forever. Can we optimize a bit?

## Density profile analysis on the CLI
- `SimAnalysis::DensityProfile` / `CompareDensityProfiles` are only reachable through Workflow.
- Undecided shape: a gmx-style `lima density -c conf.gro -f traj.trr` from files, an `mdrun --density` flag, and how to compare/average groups of runs (eg. LABEL=GLOB arguments).
- The analyses behind ShowcaseMultisim (`SimAnalysis::WaterSelfDiffusion`) need the same before the showcase can be reproduced with lima commands.

## Solvate by tiling a pre-equilibrated water box
Possibly done by commit 78acee8b ("Solvate by tiling an equilibrated water box"); verify and remove.
- Like gmx solvate with spc216.gro, instead of a jittered cubic lattice (`SimulationBuilder::SolvateGrofile`).
- Faster, the water starts with a realistic structure, and needs less EM/equilibration. The lattice still shows faintly as a sawtooth in water density profiles of short runs.

## Codebase TODOs
Find and solve TODOs spread around the codebase.
- `BoundaryConditionPublic.cpp`: the switches on `BoundaryConditionSelect` use `case None:`, which is `SupernaturalForcesSelect::None`; it only works because both are 0. Should be `NoBC` (clang warns: -Wenum-compare-switch).

---

# 6. Ideas and research

## Implicit solvent (2026-10-12)
- **Upper bound on the speed gain:** STMV without water (167k of 1.07M atoms, same 22 nm box) runs 216.6 vs 43.7 ns/day, ~5x. Not 6.4x, because the PME grid and per-cell work stay the same. GB adds a Born-radii pass and a GB pair term, so expect maybe 2-3x net (unmeasured estimate). Large charged systems like STMV are the worst case for implicit solvent (ions, structural water, salt bridges).
- **ML implicit solvents (state as of 2026-10):** GB + GNN correction is the realistic route for an MD engine.
  - Riniker group (ETH): a GNN predicts per-atom corrections to GB-Neck2, transferable across small molecules in water (Chem. Sci. 2024), 39 organic solvents (JACS 2025), code at github.com/rinikerlab/GNNImplicitSolvent. Validated on small molecules/peptides, not large proteins.
  - Protein-trained GNNs (Airas/Ding/Zhang, ISSNet) match explicit ensembles on their training proteins; transfer is open.
  - Implicit-water ML potentials (TWIN, 2026 preprint) are about accuracy, not speed: classical explicit TIP3P is still up to 25x faster per step, tested only on ~1.3k-atom proteins, no long-range electrostatics.
- **Possible plan:**
  1. Plain GB-Neck2 in LIMA (Born radii from the supercluster neighbor lists, GB term in or beside NbNonlocalKernel). Measure against the 5x dry bound and check stability/accuracy vs explicit (eg. T4).
  2. Only if worthwhile, add the GNN per-atom correction (local neighborhoods, which the neighbor lists already provide). Its inference cost on large systems is unmeasured anywhere we found.

---

# 7. AMD GPU support (HIP)
Goal: run simulations on AMD GPUs via HIP, from the same source.

## Rules
- CUDA is the source of truth. Optimize for CUDA only; AMD gets what falls out.
- Zero NVIDIA impact: HIP plumbing compiles away under nvcc.
- Separate HIP build preset; NVIDIA builds need no ROCm.
- Hybrid: every project may use the CUDA runtime names (`cudaMalloc`, `cudaStream_t`, `CudaBuffer`...) through the shim; nothing outside Engine needs to be vendor free. Engine stays largely unchanged.
- An AMD break never blocks a CUDA commit.

## Design
One Engine source tree for both backends; every kernel is written once (no hipified copy, no `LIMA_ENGINE_HIP`).
- Files:
  - `LIMA_BASE/include/GpuRuntime.h`: host runtime names (malloc/memcpy, streams, events, errors, profiler, fp16 header). The HIP branch maps the CUDA names to HIP; grow it from `tests/hipspike/GpuCompat.h` + the report's runtime list.
  - `LIMA_ENGINE/src/platform/DeviceIntrinsics.cuh` + `BlockCopy.cuh`: named device helpers for what differs (see step 4).
  - `LIMA_ENGINE/src/platform/GpuCub.cuh`: `namespace gpucub` = `cub` / `hipcub`.
  - `LIMA_ENGINE/src/platform/GpuFft.h`: the 5 cuFFT names PME uses -> hipFFT.
- Rules:
  - CUDA builds compile to identical SASS: each helper expands to today's code under CUDA; verify every step with `cuobjdump -sass` before/after.
  - No vendor `#if` inside kernels; if a function needs two versions, both go in the platform header.
  - HIP-only paths only where needed to compile or be correct, never for tuning.
  - HIP checks at startup that the device's warpSize is 32 and throws otherwise; the hardcoded 32s stay.
- CMake: cache var `LIMA_GPU_BACKEND` = `CUDA` (default) | `HIP`; enable only that language. Under HIP the `.cu` files keep their names with `LANGUAGE HIP`, Engine links `hip::hipfft` + `hip::hipcub`, flags per backend (fast-math under HIP: decide later, determinism first).
- Preset `x64-Release-HIP`: own build dir, VS 2022 v143 (14.44) toolset, gfx1103. Same exe names. `x64-Release` (CUDA) untouched.
- ConvexHullEngine: left out of HIP builds (stub that throws "not available on HIP"); MD's hull features are CUDA only for now.
- Support RDNA wave32 only at first; wave64 (MI-series) would need generic warp-size handling.

## Tests
- Determinism is within-backend: bitwise/batched tests compare AMD to AMD.
- Cross-vendor results won't be bitwise equal: exact references and VC/drift targets are kept per backend.
- Perf limits are keyed by GPU; GPUs without limits run perf tests report-only.

## Porting order
1. Prep, no Engine changes:
   - Done (2026-10-15): dead CUDA code outside Engine removed; `LIMA_BASE/include/GpuRuntime.h` (CUDA branch only) replaces every `cuda_runtime.h`/`cuda_profiler_api.h` include, SASS unchanged.
   - Still direct: `cuda_fp16.h` in `Bodies.cuh` (HIP: `hip_fp16.h`), CUB/cuFFT in Engine.
   - Keep `CudaBuffer`, `RenderDataPipe`, `genericErrorCheck`, `GenericCopyToHost` where they are; they build for HIP through the shim. A rename (e.g. `GpuBuffer`) is optional and separate.
2. Done (2026-10-15): **go**. Card/OS support confirmed, hipify report in `docs/hip/hipify_report.txt`.
   - Setup: HIP SDK for Windows 7.2 (`C:\Program Files\AMD\ROCm\7.2`), Radeon driver 32.0.13022.3006. `hipInfo` lists the 780M as gfx1103, warpSize 32. Stock ROCm libraries, no community builds or env vars needed.
   - Spike (`tests/hipspike/`): SPIKE_PASS. R2C mode error 3.3e-8, round trip 5.4e-7.
     - hipcc 7.2 (clang 21) can't parse the VS 2026 STL (MSVC 14.51 made `isgreater` etc. `constexpr`, which HIP takes as host+device): build in a VS 2022 (14.44) developer shell (`vcvars64.bat`), then `hipcc -x hip --offload-arch=gfx1103 -O2 spike.cu -o spike.exe -lhipfft`. Run with `%HIP_PATH%bin` on PATH (else 0xC0000135, DLL not found). Don't pass `-L"%HIP_PATH%lib"`: hipcc splits it on the space.
     - The HIP preset (step 3) needs the same: the v143 (14.44) toolset, until a ROCm release supports 14.5x.
     - `__HIP_PLATFORM_AMD__` is only defined by the HIP headers, so detect hipcc with `__HIP__` (clang's own) before including anything.
   - hipify report: `hipify-clang --examine` needed workarounds for CUDA 13.1 (it supports <= 12.8); see the report header. No unconverted names. Compile probes for what hipify doesn't judge:
     - **Warp masks must be 64-bit** under HIP (static_assert): the 52 `0xFFFFFFFFu` masks need a named full-warp mask macro. `__ballot_sync` returns 64 bits; storing it in `uint32_t` is fine on wave32.
     - **`__ldcg`** exists only for half types: 6 uses in `EnergyMinimization.cuh:194`. Needs a load that is coherent across CUs, not a plain load.
     - **PTX:** `asm("trap;")` (2) -> `__builtin_trap()`; `ex2.approx.ftz.f32` in `DeviceAlgorithmsPrivate.cuh:68` -> `exp2f`/`__builtin_amdgcn_exp2f`.
     - **`cooperative_groups::memcpy_async`/`wait`** (`ParticleClusters.cu:499`) don't exist in HIP: plain block-strided copy + `__syncthreads`.
     - OK as is: all other `*_sync` warp ops incl. `__match_any_sync` and `__reduce_*_sync`, `__syncwarp`, `__popc`, `__ffs`, `__ldg`, `__half2float`, `__launch_bounds__`.
   - Runtime API list for `GpuRuntime.h`: the "converted references" table in the report (malloc/memcpy/memset (+Async), streams, events, `cudaStreamWaitEvent`, mapped host alloc, `cudaMallocManaged`, device props, errors, profiler). Libraries: CUB `DeviceReduce`/`DeviceScan`/`BlockReduce` -> hipCUB; cuFFT `PlanMany`/`SetStream`/`ExecR2C`/`ExecC2R`/`Destroy` -> hipFFT (proven by the spike).
3. Done (2026-10-15): build preset + shim. CUDA build: SASS identical, batch tests pass.
   - `LIMA_GPU_BACKEND` (CUDA|HIP) in the top-level CMake; `x64-Release-HIP` in CMakeSettings.json. Configure from a VS 2022 14.44 shell (`vcvars64.bat -vcvars_ver=14.44`); ROCm 7.2 also compiles on the Windows 10 desktop (no AMD GPU needed to build).
   - CMake can't mix MSVC and Clang, so HIP builds compile all C/C++ with ROCm's clang++ (MSVC STL/ABI). `hip::device`'s `-x hip` is cleared so only the .cu files compile as HIP.
   - HIP .cu files are C++20 (`CMAKE_HIP_STANDARD`): from C++23 the MSVC STL's constexpr isinf/isnan clash with HIP's device overloads. Engine uses no C++23 features so far.
   - `GpuRuntime.h` HIP branch, `platform/GpuFft.h`, `platform/GpuCubBlock.cuh` + `GpuCubDevice.cuh` (split so the CUDA SASS stays identical: the device headers add a kernel). ConvexHullEngine is a throwing stub under HIP.
   - Under HIP, everything compiles except the engine kernels: BASE, TOOLS, MD, TESTCASES, Engine's .cpp and CubWrappers.cu. Clang-only fixes in host code: `MDFiles.h` (inline ctors needing incomplete types), `SimulationBuilder.cpp` (`Int3 = 0`).
   - Noise to silence later: `hipError_t` is `[[nodiscard]]` (~60 warnings).
4. HIP compile done (2026-10-15): every target builds and links under HIP; CUDA SASS identical, batch tests pass.
   - `platform/DeviceIntrinsics.cuh`: `WarpMask`/`LIMA_FULL_WARP_MASK` (52 masks), `DeviceIntrinsics::Trap`, `Exp2Approx` (PTX ex2), `LoadCoherent` (`__ldcg`; HIP: agent-scope relaxed atomic load, so a stale per-CU L1 line can't be read), `LIMA_SHARED_UNINIT` (8 `__shared__ Float3` arrays: nvcc ignores their initializers, clang rejects them). `platform/BlockCopy.cuh`: `memcpy_async` (separate header: cooperative_groups adds device globals, which shifted PME's SASS constant bank offsets).
   - `LIMA_GPU_COMPILER` (in `GpuRuntime.h`) replaces `#ifdef __CUDACC__`. Kernel launches written `<< <` / `>> >` (VS autoformat) became `<<<` / `>>>`: clang needs the tokens.
   - HIP only: Engine's constructor throws unless the GPU's warpSize is 32.
   - Next, on the laptop (780M). Pull (`git fetch` + `git reset --hard origin/HIPport`), then from a VS 2022 14.44 developer shell configure `x64-Release-HIP` (VS configuration, or `cmake -S . -B build/x64-Release-HIP -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo -DLIMA_GPU_BACKEND=HIP -DCMAKE_CXX_COMPILER=<ROCm>/bin/clang++.exe -DCMAKE_C_COMPILER=<ROCm>/bin/clang.exe -DCMAKE_HIP_COMPILER=<ROCm>/bin/clang++.exe -DCMAKE_PREFIX_PATH=<ROCm>`) and build. Run with `%HIP_PATH%bin` on PATH. In order, recording pass/fail, crashes and timings here:
     1. `agenttesting.exe --engine-batch-tests`: batched vs standalone, bitwise within HIP. The first determinism check on AMD.
     2. `LIMA_EM_LABEL=hip780m agenttesting.exe --energy-minimization-tests`: compare results_hip780m with the CUDA baseline (not bitwise).
     3. `limatest.exe` (RunAllUnitTests), report only: exact references and VC/drift targets are CUDA's, so mismatches are expected; note which fail and by how much (per-backend references come later).
     4. If something crashes or hangs: first suspects are `LoadCoherent` (EM), `BlockCopyAsync` (clustering), and warp-mask/`__match_any_sync` semantics.
   - **First run on the 780M (2026-10-15, 3293b308): it runs.** Configure + build from scratch worked as written above (1091 warnings, 0 errors). No crashes or hangs.
     1. `--engine-batch-tests`: **pass**, bitwise within HIP (315 s).
     2. EM tests: not run. `LIMA_data/EnergyMinimizationTests` exists only on the desktop (not in the LIMA_data repo); copy it over to run them.
     3. `limatest.exe`: 27/39 pass (369 s). Pass: all force/potential sanity tests, **Deterministic Simulations**, Constraints, T4Lysozyme VC 2.851e-06 / 2.844e-06, PME vs Ewald (force 0.10%, energy 0.06%), angle/dihedral/Met/8Res VC within targets, the builder/ToGmx/topology tests.
        - Energy-gradient targets exceeded (CUDA-calibrated, expected): Pool 3.421e-08 / 3.421e-08 (equal to the CUDA value), PoolCompSol 2.31e-07 / 1.02e-07, Singlebond 4.04e-07 / 2.62e-07, ImproperDihedral 4.59e-05 / 4.47e-05, Solvent 7.91e-08 / 7.62e-08. PoolCompSol (2.3x) and Singlebond (1.5x) are the largest; worth a look when per-backend targets are set.
        - Perf: four batched T4 1.19x vs single (56.1 vs 47.1 ns/day), limit 1.25x. Needs a per-GPU limit.
        - Laptop data, not HIP: ForceComparisons, LoadT4, T4 throw "Unknown simulation parameter 'save_trajectory=false'" (in `LIMA_data/benchmarking/sim_params.txt` and `CompareWithOtherMdEngines/Forcecomparison1step/sim_params.txt`: the LIMA_data repo is behind the desktop's data). 3J3Q/STMV inputs missing (large systems are skipped on AMD anyway).
5. Warp-level kernels (neighbor search, NB) + determinism tests
6. PME (hipFFT), test suite (small systems)
7. Remaining features, rendering interop

## Decisions (2026-10-15)
- Target: Radeon 780M (gfx1103, RDNA3 iGPU, wave32) on Windows HIP SDK. Not officially supported by AMD, but community builds run; start with a spike (one kernel + a small hipFFT transform) before any refactor.
- Shim header, not a hipify-generated tree.
- AMD validated manually on the laptop, before releases and after big engine changes.
- First milestone: headless (`lima mdrun` + tests); rendering later.
- No HIP binaries shipped until it runs well.
- No HIP-only tuning for now.
- ConvexHullEngine stays CUDA-only and is left out of HIP builds (2026-10-15); merging it into Engine is a later sprint.
- Small systems only on AMD (no STMV/3J3Q); per-GPU test limits skip the large ones.
