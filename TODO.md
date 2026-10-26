# LIMA TODO

Sections:
1. Priority
2. Missing MD features
3. Correctness issues
4. Performance
5. Tooling, CLI and codebase
6. Ideas and research
7. AMD GPU support (HIP)

---

# 1. Priority: 

## CLI bugs from the release battery (2026-10-26)
Found with `D:\LIMA_ARTIFACTORY\release_battery` (`run_battery.sh`), lima 1.4.0.
- **Done** (2026-10-26): togmx output solvates (internal-forcefield moleculetypes stay includes, each include written once), include paths use `/`, `--set` values are case-insensitive. Covered by limaclitest togmx and wait.
- **Left:** `printToFile` still re-writes every moleculetype from a *local* include, so a local file with several moleculetypes (eg. a user-supplied `charmm27.ff/ions.itp`) is overwritten with only the last one.
- **Done** (2026-10-26): boxing a togmx protein. `solvate` takes any box (cells of length/floor(length) per axis, as the engine; integer boxes give bitwise the same output), `editconf -b/--box-size` and `-d/--distance` (diameter + 2d, cubic, as GROMACS) for the whole system, and togmx notes when the crystal cell is missing or smaller than the protein.
- **Left:** `insertmolecule` still throws `Illegal call to GetMoleculeType` for a togmx source, as its topology also has the water and ion moleculetypes. Decide whether it should insert every molecule of the source, or only the non-solvent ones.

## togmx and topology findings from the full release battery (2026-10-26)
Run `runs/20261026_190356` of `D:\LIMA_ARTIFACTORY
elease_battery`; `tools/bondscan.py DIR...` lists bonds > 0.25 nm in togmx outputs.
- **Done** (2026-10-26), tests "ToGmx damaged input" and "Topology GROMACS conventions":
  - Chains are split where C(i)-N(i+1) > 0.3 nm, each piece its own molecule (`Protein_chain_A_2`, ...) with the --ter-n/--ter-c termini, and the breaks are listed (6VXX: 33 breaks, 36 molecules; no bond > 0.25 nm after EM or MD).
  - MSE becomes MET. Other HETATM residues peptide-bonded into a chain are an error naming them; caps (ACE, NH2, NME, NMA, FOR) are dropped with a note.
  - Altlocs compare case-insensitively; a residue without the requested one keeps its most occupied one, with a note.
  - Missing heavy atoms point to PDBFixer/MODELLER; repeated atoms or residue records are an error.
  - Topology directives ignore case, `_` and `-`; `[ molecules ]` names ignore case; `[ dihedraltypes ]` of unsupported functs are skipped, a molecule using one is still an error; unsupported directives say so.
  - NoBC is rejected with a clear error (glycine_vacuo straddled the z=0 face, see [NoBC](#nobc-is-only-non-periodic-in-the-force-kernels)); with PBC the vacuum systems are stable.
- **Left:**
  - In-chain modified residues (SEP, TPO, PTR, HYP, MLY, CSO, PCA, D-amino acids, CRO) need templates to convert; for now they are an error. Capping groups could become the matching --ter-n/--ter-c instead of being dropped.
  - Missing side chains are not rebuilt (13 of ~70 real entries).
  - OPLS-AA: atomtypes whose `bond_type` differs from their name (`opls_001` -> `C`) are rejected, and Ryckaert-Bellemans dihedrals (funct 3) are unsupported; together they block every OPLS system in the battery.

## togmx: AMBER99SB force field
**Done** (2026-10-24): `lima togmx --ff amber99sb`. Villin's coordinates give exactly pdb2gmx's topology (atoms, types, bonds, pairs, angles, dihedrals, impropers), test "ToGmx AMBER99SB matches GROMACS".
- **Left:** only charged termini. Neutral and capped ends are separate residues in AMBER (ACE, NME, NHE), which togmx doesn't read or add yet.
- Not yet run on the other Salad systems.

## Villin (Salad, AMBER99SB) vs GROMACS 2024.2: remaining issues (2026-10-14)
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
- Run the common RMSD/RMSF/Rg/RDF/H-bond/DSSP analyses on one LIMA and one GROMACS trajectory per system before adding replicate ensembles.

---

# 2. Missing MD features
GROMACS has these, LIMA does not (status checked 2026-10-13). Check each against what the Salad comparison actually needs, and do those first.

## Topology / input
- **`[ constraints ]` sections:** not read yet. constraints = h-bonds and `[ settles ]` are handled (see Constraints below), but topologies that list their constraints explicitly lose them.
- **MDP validation:** integrator, pbc, coulombtype, rvdw, output cadence, thermostat/barostat and seeds are validated. `fourierspacing` sets the largest PME grid spacing. Still silently ignored: `pme-order`, `comm-mode`, and unknown keys in general.
- **Regression coverage:** add a small AMBER99SB/TIP3P constrained-water fixture with tests for topology loading, constraint preservation, short-run stability, and readable `.trr` output. Validate against GROMACS observables, not frame-by-frame coordinates.

## Constraints
constraints = h-bonds (SHAKE/RATTLE per pcluster) and SETTLE are done. Left:
- **LINCS / constraints = all-bonds:** postponed, h-bonds covers the usual CHARMM/AMBER setups. Needs coupled constraints across pclusters, ie. across superclusters, so a separate deterministic solver.
- **Position precision limits all constraints:** positions are absolute float coordinates in the box, so a 0.1 nm bond holds only to ~6e-6 relative whatever the solver. SETTLE assumes the water is exactly rigid before each step, so it re-snaps the rounding every step: T4 VC over 10k steps 1.50e-5 vs 1.33e-5 with SHAKE. Fix: store positions relative to the supercluster (or cell), as GROMACS keeps coordinates local. Would also help everything else.
- **Hydrogen mass repartitioning** for 4-5 fs steps. With constraints = h-bonds, T4 (vacuum and water) already runs 2000 steps at 4 fs without it; check energy conservation there first.
- **initial_temperature** draws velocities with components along the constraints, which the first step removes, so the start is somewhat colder than asked.
- Not yet compared with a constrained GROMACS run (bond lengths, water RDF, energy drift).

## Thermostats
- **Missing:** Nose-Hoover, Langevin/stochastic dynamics, per-group coupling (eg. protein vs solvent).

## Constant pressure (NPT)
- **Missing:** stochastic cell rescaling / a canonical NPT barostat, anisotropic or semi-isotropic coupling, and independently barostatted simulations in one batch.

## Restraints and pulling
- Position restraints are skipped (posre includes ignored in `Topologyfile.cpp`), so no restrained equilibration.
- Missing: distance/dihedral/flat-bottom restraints, pulling (umbrella sampling, steered MD).

## Other
- Virtual sites (TIP4P water, hydrogen vsites for 4-5 fs steps)
- Free-energy perturbation (lambda coupling)
- Center-of-mass motion removal
- Energy groups

---

# 3. Correctness issues

## Forcefields with different 1-4 scaling can't be mixed
- A topology has one `[ defaults ]` (gen-pairs, fudgeLJ, fudgeQQ), applied to the whole system.
- So merging eg. a CHARMM-style protein (fudge 1.0/1.0, as in STMV) with Slipids lipids (0.5/0.8333) gives one of them the wrong 1-4 interactions.
- `SimulationBuilder::InsertSystem` merges forcefields anyway (used by the large-EM showcase).
- Fix: per-moleculetype 1-4 scaling, or explicit `[ pairs ]` parameters with the scaling baked in.

## Topology sections that are silently skipped
`Topologyfile.cpp` maps `[ exclusions ]` and `[ nonbond_params ]` to `notimplemented`, so their entries are dropped without a message.
- `[ exclusions ]`: parse into the per-moleculetype exclusion lists that already drive the LJ skip and the Ewald exclusion correction. In CHARMM36 water files they duplicate what `[ settles ]` gives, so those are harmless today.
- `[ nonbond_params ]` (CHARMM NBFIX, eg. in `LIMA_data/benchmarking/aqp_ensemble`'s charmm36): the NB kernel only has per-particle combination rules, so per-pair overrides need a pair-type table (performance work). Until then, warn with the number of ignored pairs.

## NoBC is only non-periodic in the force kernels
With `boundarycondition = NoBC` the force kernels use `NoBoundaryCondition`, but the rest of the pipeline stays periodic: `ApplyBoundaryCondition` (`ParticleClusters.cu`) wraps positions, `GetPclusterPositions` wraps the cluster centers, `ClusteringPretransferKernel` hyperpositions, and the neighbor search (`SuperclusterTaskBuilder.cu`) wraps cells and uses periodic images. A particle crossing the box edge is teleported to the other side, and its neighbors across the edge are listed but get no force. No test uses NoBC.
- **Rejected for now** (2026-10-26): `verifyEngine` throws for NoBC, as glycine_vacuo blew up from wrapped bonds. To support it, make clustering and neighbor search BC-aware end to end (no wrap, edge cells absorb out-of-box clusters via `CellOf`'s clamp, no periodic cell images).
- Templating only `ApplyBoundaryCondition` is not enough: an unwrapped cluster would still be placed by its wrapped center.

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
  - Also check: `BoxImageBuilder.cpp` copies the whole `GroFile` into `BoxImage` (millions of atoms for 3J3Q). Move it if the caller doesn't need it afterwards.
- **"ToGmx large CIF benchmark"** ~8.5 -> 10.2 s (limit 10 -> 12). Cause not investigated (ToGmx builds no box, so maybe the topology parser changes, or noise).

## An open display costs ~15% simulation speed (2026-10-25)
STMV runs at ~2.87 ms/step headless and ~3.3 ms/step with the display open, with the old and the new (ray-cast + AO) shaders alike, so it is not shader cost. Find out why and fix it. Each GL frame slows the CUDA steps (~0.6 ms of GL work per frame, yet 30 fps cost ~20% and 50 fps ~37%), so suspect WDDM context switching between the GL and CUDA contexts. Note the display also redraws on a timer when nothing changed. To measure: `Benchmarks::Bench` hardcodes `job.mode = EnvMode::Headless`, so pass `envmode` through and run limaprofile with `EnvMode::Full`.

## Startup latency for small systems
Investigate time from command to first step, eg. T4 in a 7^3 nm^3 box.

## Kernel-level ideas
Not set-in-stone todos; concepts worth exploring. Some make the most sense together.
- **Engine on its own thread:** instead of Env having a thread dedicated to RunSimulation, Engine gets its own thread, so an mdrun only steps and updates the renderpipeline/simstatus. Engine would persist, so slightly lower setup latency. Engine must free its memory when idle, so we're not hogging 10 GB VRAM while Env just does preprocessing.
- **EmitQuarterEntriesKernel:** now ~0.77 ms per STMV rebuild (was 1.37). The cost was the bonded-pair exclusions (~1 active lane, dependent global loads): own particle ids now in shared memory, query ids loaded once per quarter, one pass over each bonded set for all 4 query particles. Rejected: int4 reads of the set (slower). What's left is mostly writing ~300 MB of QuarterEntries; shrinking that means changing NB's entry format.
- **FindNeighborsKernel:** now ~2.7 ms per STMV rebuild (was 3.9).
  - Issue-bound (~68% issue active), so latency tricks don't help: 2 candidates per warp iteration and a candidate-particles vs own-sphere filter were both slower.
  - Left: cell enumeration (~74 cells/SC, ~7 of 32 lanes active) and the bonded check (~8%).
- **nstlist vs list-buffer tuning:** a bigger buffer and fewer rebuilds may be cheaper overall. A parameter/accuracy decision.

## NB kernel: explored and rejected (2026-10-20)
Measured on STMV; don't retry without a new angle. Any idea that removes work below entry granularity must be judged by the warp-union cost (a warp runs its 8 entries over the union of their own-quarter masks), not per entry.
- **Combining j forces in shared memory before the atomics:** per block it would cut atomic requests ~3x, but needs ~40-55 KB smem per block (16 blocks/SM). Infeasible.
- **Ordering entries by neighbor instead of mask class:** -31% atomic requests, but +18-25% pair computations from union waste.
- **Quarter-AoS force accumulator layout:** ~2-3% ceiling, touches every ForceAccumulator user.
- **4x2 j-granularity:** packing two j-halves into a virtual quarter: NB +1%, atomics +17%. True per-(own quarter, j half) execution saves 10.9% per entry but only 0.6% under the warp union; sorting entries by the finer class makes it worse.
- **Dynamic pruning:** with compacted lists, pruning every ~5 steps at rc + 0.0125 nm saves ~7% pair work, but re-emitting the lists (~0.5 ms) costs about as much. In-place mask pruning without compaction: -0.6%.
- **Flattened entries:** each group loops over its own mask's quarters (the warp runs max popcount, not the union), entries sorted by popcount, own forces per group in shared memory. Pair work -4.7%, but the per-iteration overhead (runtime quarter, shared adds, rolled loop) costs more: NB +3% (1.66 -> 1.70 ms). Unrolling the loop spills at 64 registers: 3.4 ms. The ceiling is under 5% of the instructions and NB is half L2-bound, so not worth a cheaper variant either.

## Kernel graphs: done (2026-10-20)
MD steps are captured into a CUDA graph every step: T4 -18%, STMV neutral. Measured and rejected: encoding a launch order (NB waits for PME's spreading, bonds wait for NB) is slower on both, T4 0.082 -> 0.086/0.107 ms/step, STMV +0.8/+1.3%.
nsys needs `--cuda-graph-trace=node` to show graph kernels, but node tracing perturbs scheduling; trust wall-clock A/B.

---

# 5. Tooling, CLI and codebase

## Display: image and movie export
- **Image:** save the current view as a PNG, at a resolution above the window's (eg. 4x), with ambient occlusion. Render offscreen into a larger target, like `AmbientOcclusion`'s scene framebuffer, and read it back. Menu item plus a `lima render --screenshot PATH` option for scripted figures.
- **Movie:** a camera orbit, or a trajectory played back, written as frames or encoded to MP4 (ffmpeg if found on the PATH, else PNG frames).

## Compile times
Compiling takes forever. Can we optimize a bit?

## Density profile analysis on the CLI
- `SimAnalysis::DensityProfile` / `CompareDensityProfiles` are only reachable through Workflow.
- Undecided shape: a gmx-style `lima density -c conf.gro -f traj.trr` from files, an `mdrun --density` flag, and how to compare/average groups of runs (eg. LABEL=GLOB arguments).
- The analyses behind ShowcaseMultisim (`SimAnalysis::WaterSelfDiffusion`) need the same before the showcase can be reproduced with lima commands.

## Codebase TODOs
Left from the 2026-10-20 sweep of TODO comments:
- `TestBoxIsSavedCorrectlyBetweenSimulations` (`SetupTests.cpp`) only checks that the continued simulation exists. Make it compare frame 0 of the second run with the last frame of the first.
- `TestFilesAreCachedAsBinaries` (`FileTests.cpp`) is not registered, its data (`LIMA_data/filetests`) is gone, and topologies are no longer cached. Keep the .gro check on a temp copy of an automated-test .gro and register it.
- `Float3::rotateAroundOrigo` (`LimaTypes.cuh`) is marked "proven wrong"; it composes rotations about the fixed x, y, z axes, which is valid, but the pitch/yaw/roll naming is muddled. Only `ForceCorrectness.cpp` uses it: inline it there as `rodriguesRotatation` calls and remove it from `Float3`.

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
Runs on the Radeon 780M (gfx1103, RDNA3 iGPU, wave32) with Windows HIP SDK 7.2, from the same source: `GpuRuntime.h` + `LIMA_ENGINE/src/platform/`. Rules for kernel code are in CLAUDE.md.

## Building on the laptop
- Configure from a VS 2022 14.44 developer shell (`vcvars64.bat -vcvars_ver=14.44`): hipcc 7.2 can't parse the VS 2026 (14.5x) STL. VS configuration `x64-Release-HIP`, or `cmake -S . -B build/x64-Release-HIP -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo -DLIMA_GPU_BACKEND=HIP -DCMAKE_CXX_COMPILER=<ROCm>/bin/clang++.exe -DCMAKE_C_COMPILER=<ROCm>/bin/clang.exe -DCMAKE_HIP_COMPILER=<ROCm>/bin/clang++.exe -DCMAKE_PREFIX_PATH=<ROCm>`.
- Run with `%HIP_PATH%bin` on PATH (else 0xC0000135, DLL not found).
- ROCm 7.2 also compiles on the Windows 10 desktop, so compile checks need no AMD GPU.

## Open
- EM tests not run on HIP: `LIMA_data/EnergyMinimizationTests` exists only on the desktop; copy it over, then `LIMA_EM_LABEL=hip780m agenttesting.exe --energy-minimization-tests` and compare with the CUDA baseline (not bitwise).
- Per-backend VC/drift targets: on HIP (fast math) Anglebond 7.2e-07 / 2.5e-07 is over, Singlebond/Met/Solvent slightly over.
- Per-GPU perf limits: batched T4 on the 780M is 1.19x vs single, limit 1.25x.
- The laptop's LIMA_data is behind the desktop's (`save_trajectory=false` unknown in ForceComparisons/LoadT4/T4 sim_params).
- Silence the ~60 `[[nodiscard]]` `hipError_t` warnings.
- Remaining features and rendering interop (first milestone was headless).
- Known HIP-only performance issues, not addressed (no HIP-only tuning for now):
  - Stream fork/join costs ~300 us/step under Windows HIP (event waits, independent of stream count); one stream: +11%.
  - NB `__launch_bounds__(64, 16)` caps it at 64 VGPRs with 48 B/lane scratch spills; without the min-blocks hint +5%.

## Decisions (2026-10-15)
- AMD validated manually on the laptop, before releases and after big engine changes.
- No HIP binaries shipped until it runs well.
- ConvexHullEngine stays CUDA-only and is left out of HIP builds; merging it into Engine is a later sprint.
- Small systems only on AMD (no STMV/3J3Q); per-GPU test limits skip the large ones.
- Determinism is within-backend; cross-vendor results are not bitwise equal, so exact references and VC/drift targets are kept per backend.
- RDNA wave32 only; wave64 (MI-series) would need generic warp-size handling.
