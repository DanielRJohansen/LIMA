# Performance notes

What was learned optimizing the MD step, and what was tried and rejected with the reason, so ideas are not retried without a new
angle. Measured on an RTX 4090 (Windows, WDDM) with the CLAUDE.md protocol. Sizes: villin (5k atoms), adh_dodec (191k), stmv (1.07M),
the `lima.mdp` of `LIMA_data/salad_gromacs_bench` (PME, h-bond constraints, v-rescale every step, 2 fs, nstlist 20).

## Session 2026-10-29 (branch `perf`, from 307ecf83)

Result: villin -13.4%, adh_dodec -13.4%, stmv -4.9% (paired runs); `--stmv-constrained` -2.7%, `--stmv` neutral.

### Where the time goes (after the session)
- **stmv (2.99 ms/step):** NB ~1.7 ms, then the PME chain ~0.45 ms, mostly serial after NB (NB fills the GPU), then the constrained
  integrate ~0.23 ms. Rebuilds ~5 ms per 20 steps. The step is throughput bound: NB is at ~75% issue and ~75% L2 (atomic requests),
  the integrate kernel at 84% DRAM.
- **adh_dodec (0.43 ms/step):** NB ~180 us, the PME chain (~130 us alone) ends ~85 us after NB and is the critical path, then the
  integrate ~40 us. A rebuild costs ~680 us, mostly kernels (FindNeighbors 353 us, Emit 192 us).
- **villin (0.050 ms/step):** a plain step is ~37 us: the serial PME chain (~11 latency-bound kernels of 1.5-4 us) and the integrate
  (~7-8 us on 36 blocks); NB and bonded run alongside. A rebuild costs ~264 us of wall clock against ~145 us of kernels, but the
  difference is not host idle on the critical path: removing the waits and the recapture gained nothing (see the fixed capacity below).
- Rebuild cost was measured as T = plain + R / nstlist from runs at nstlist 20, 40 and 80 (temporary limaprofile override).

### Measuring
- Other GPU clients on the machine (the desktop apps, NZXT CAM, iCUE) time-slice the GPU: every ~16 ms a stmv step stretches by ~1.6 ms,
  and adh runs of the same exe differed by up to 8%. The 500-step stmv suite runs could not resolve 2%.
- What worked: paired alternating runs (villin 40000, adh 6000, stmv 1500 steps; A B then B A), median of the per-pair differences;
  ncu kernel durations for kernel-level decisions.
- `limaprofile --salad <system> <steps>` prints a checksum of the temperature series: equal between builds only if the dynamics are
  bitwise identical. The checksum depends on the step count.
- nsys with `--cuda-graph-trace=node` inflates host API times several fold; use it to find gaps, not to size them.
- Kernel durations in nsys include contention with concurrent kernels; check isolated ncu durations and the critical path per step.

### Committed changes and why they worked
- **Constrained integrate (904ae393).** One thread per pcluster solves RATTLE/SETTLE while the cell waits at barriers (~50% of stall
  samples). The solver indexed register arrays at runtime (group slots by constraint particle, loops to `nConstraints`), so it spilled to
  local memory. Selects and unrolled guarded loops removed that; `__launch_bounds__(256, 3)` brought registers 94 -> 80 for 3 blocks/SM.
  Kernel villin 12.5 -> 8.8 us, adh 65 -> 40 us, stmv 307 -> 232 us. Bitwise identical.
- **Thermostat kinetic energies from the integration (00d4ca5a, ce12d302).** The separate `ComputeKineticEnergyFromSlotsKernel` read all
  integration states (37 MB on stmv) every step on a side stream, competing with NB/PME: removing it entirely was -3.1% on stmv. The
  integrate kernels now write each particle's KE; the slots kernel only runs after a batch is (re)packed. Only for a measurement interval
  of 1: with nsttcouple 200 (unconstrained `--stmv`) the per-step writes cost +0.7%.
- **FindNeighbors groups (99e5591d).** The split search used groups of 8 cell offsets; villin has 36 offsets, so 5 of its 8 chunks
  worked. Group size nOffsets / nChunks and up to 32 chunks: 70.7 -> 42.9 us per rebuild.
- **SuperclusterBoundsKernel (c722d995)** and **ClusteringKernel (42035c92).** One warp per cell walked its superclusters serially;
  with 36 cells (villin) the GPU idled. 8 warps per cell for bounds (37.6 -> 10.4 us), 128 threads for clustering (43 -> 27 us); both
  also faster on stmv. Results don't depend on the block size.
- **No entry-count wait for small systems (27caacd8).** When the entries' upper bound (4 per neighbor, 512 neighbors per supercluster)
  fits in 64 MB it is allocated instead, so the host captures the next step while the GPU searches. The overflow is checked at the next
  synchronization. Pitfall: the deferred readback needs its own mapped buffer, or the next rebuild's host writes race the GPU's gathers.
- **PME z pass (b13a74c5).** Skipping ToDensity, the FFTs and the Greens multiply entirely was -21% villin, -21% adh, -5% stmv, so the
  FFT section matters at every size. `PmeFft.cuh`'s ZGreensKernel does z FFT, Greens multiply (from separable per-axis tables) and
  inverse z in one kernel, with cuFFT 2D plans for the planes: 8 kernels -> 6, villin -2.1%, adh -1.8%, stmv neutral. Not bitwise:
  VC shifts at rounding level (Solvent, PoolCompSol exceed their zero-margin targets).

### Tried and rejected
| Idea | Result | Why |
|---|---|---|
| Fewer barriers in the constrained integrate (solver computes its group's steps, 5 -> 2 syncs) | no gain | The cell waits on the serial solver work, not on the barriers themselves |
| `__launch_bounds__(256, 4)` on the constrained integrate | slower (stmv 232 -> 266 us) | 64 registers spill |
| Capturing PME first in the step graph, so its chain starts earlier | neutral (villin -0.6%) | The median critical path barely changed; launch order varies per step |
| Bonded forces as atomics into the force accumulator instead of primary/extra planes | stmv +1.1%, villin +0.7% | Bonded kernel 102 -> 161 us on stmv (atomics contend with NB's); integrate only -27 us. Confirms d9597d94 |
| Skip the primary bond planes for unbonded slots (flag in the bond references) | no gain | Never-written zero planes stay L2-resident, so they cost no DRAM |
| Skip the per-step `forcesMagnitudeSquared` write in the constrained integrate | no gain | 4 MB of scattered writes, absorbed by L2 |
| Emit with 4 warps (superclusters) per block | adh +2%, stmv -5% of Emit (-2 us/step) | Not occupancy bound: half issue bound, half load latency |
| FindNeighbors merge: load all group counts at once | 42.9 -> 41.8 us | The serial merge was not the bottleneck |
| FindNeighbors on adh with 2-5 chunks per supercluster | 359 -> 411-458 us | adh's search is throughput bound, not tail bound |
| Extra `cudaStreamQuery` flush after launching FindNeighbors | neutral | The rebuild's idle is elsewhere |
| Full custom 5-pass FFT (fused ToDensity + R2C x, y, z+Greens, y, C2R x) | ties cuFFT (stmv 195 vs 202 us) | Each own x/y pass ~1.5x slower than cuFFT's (shared-memory Stockham, ~60% occupancy) cancels the 3 saved passes. Pitfalls found on the way: bank conflicts with an even row stride, and a load loop that stored to the same array, serializing the loads |
| Hybrid FFT with the Greens function loaded per mode | stmv +0.9% | Separable on-the-fly Greens tables fixed it (neutral) |
| Double-buffered force accumulator (zero off the critical path) | not tried | ~1% stmv for even/odd step graphs |
| Store half-step velocities instead of velocity + previous force | not tried | Saves 24 B/slot, but states could no longer give full-step velocities to the host paths |
| Fixed supercluster capacity for small systems: nSuperclusters = cells x 16, superclusters past the count padded on the GPU, no clustering wait, step graph kept across rebuilds | villin +1.7%, rnase +2.5% (paired) | Bitwise identical and the graph was kept (21 captures in 2000 steps), and under nsys it looked faster (mean step 66 -> 57 us). But nsys inflates host API time: in real runs the host work around a rebuild is already hidden, so only the padding's cost showed: NB, Distribute and the search over 1.8x (rnase 2.2x) the superclusters, and fewer search chunks per supercluster |

### Open
- **villin's rebuild cycle, measured with GPU timestamps** (a `%globaltimer` stamp at each step's start and end in the step graph, and
  between the rebuild phases; no nsys). Every step takes 29.7 us on the GPU, and the gap between replayed graphs is ~2 us. The cycle's
  excess is GPU idle, not slower steps: in the rebuild step ~130 us of kernels and ~100 us of idle (37 us around the clustering wait,
  the rest at launch boundaries), plus 13 us before the step preceding a rebuild and 34 us before the second step after it. Those two
  sit at WDDM submission boundaries: a graph queued 700 us earlier still starts ~15 us late when it opens a new submission. Moving the
  boundaries (no flush at the end of the replay loop before a rebuild, no flush of the captured launch, a flush after the first
  replay) removes the gaps in the timeline, but the idle reappears inside the rebuild: wall clock neutral (+0.1%, 10 paired runs).
  The fixed capacity removed the wait's idle (37 -> 8 us) but its search got 26 us slower (fewer chunks per supercluster). So the
  idle is the cost of the rebuild's one host round trip; what is left to win is ~1-2% on villin, for a device-driven rebuild.
- **cuFFT LTO callbacks** (ToDensity into the R2C load), see PME.cuh: one more pass less, at the cost of JIT linking.
- **PME internals on stmv:** Interpolate ~125 us, Spread ~70 us (global tile-flush atomics; 2x2x2 chargeblocks per tile untried).
- **FindNeighbors cell loop** runs with ~7 of 32 lanes active (40 hit cells per supercluster on adh); flattening (cell, candidate)
  pairs might save 10-15% of the kernel, but small changes to this kernel have shifted stmv by 5% through scheduling before.
