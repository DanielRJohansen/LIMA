# Plan to reach PME performance with FMM

Analysis date: 2026-10-12. This plan uses the existing full Nsight Compute reports,
current source, existing validation results, and the publications linked below.
No new simulations or GPU profiling runs were performed for this analysis.

The acceptance criterion is complete-step FMM performance at or below PME at
matched electrostatic accuracy on STMV, on the same RTX 4090. The present result
is 12.617–12.711 ms/step against 3.742–3.750 ms/step. Parity requires removing about
8.9 ms/step. Kernel tuning alone is not an adequate strategy; the plan combines
precision and specialization fixes, a different M2L schedule/compression, and
shared near-field work. An analytical FMM representation is the planned alternative
if the Chebyshev representation cannot meet its allocated cost.

## Evidence and limits of the measurements

Inputs:

- `build/x64-Release/ncu/fmm-optimized-full.ncu-rep`
- `build/x64-Release/ncu/pme-fmm-comparison-full.ncu-rep`
- Their raw CSV exports in `build/`, and existing benchmark/test logs.
- P2M SASS exported from the existing report to `build/fmm-plan-p2m-sass.csv`.

The final normal-step FMM launches are IDs 129–182; ID 183 is the final state
store. Nsight replays kernels and affects concurrency. The following are isolated
kernel service times, not an end-to-end timeline or additive wall-clock savings.
The two previously recorded alternating wall-clock rounds remain the parity
baseline. Logging and neighbor rebuild costs must remain in the eventual result.

| Work | Saved profile time | Diagnosis |
|---|---:|---|
| M2L, all levels | 4.041 ms | Leaf contributes 3.099 ms; representation and execution both matter |
| Projection/expansion | 0.698 ms | Repeated per-cell matrix-vector products |
| P2M | 0.913 ms | FP64 pipeline 84.6%; expensive coordinate arithmetic and local arrays |
| L2P | 1.742 ms | FP64 pipeline 85.3%, issue activity 8.9%, dynamic local arrays |
| Coulomb P2P | 1.302 ms | FP64 pipeline 53.5%; directed pairs, separate traversal from LJ |
| LJ nonbonded | 2.004 ms | L2 throughput 97.5%; redundant geometry/loading relative to P2P |

The leaf M2L uses 142 registers/thread, about 24.9% achieved occupancy, 33.2%
tensor-pipeline activity, 67.1% L2 throughput, and only 2.2% DRAM throughput.
It records 497.8 million shared wavefronts against 246.7 million ideal wavefronts.
This is not a DRAM-bandwidth bottleneck. Bigger tiles alone are not the answer:
they trade improved reuse against register pressure and shared-memory costs.

L2P reports roughly 78.4 million local-load sectors, despite zero register-spill
read instructions. Runtime-indexed arrays can live in local memory without being
classified as register spills. P2M SASS contains repeated reciprocal/refinement
instructions for double division and local array stores. These observations make
compile-time order and coordinate arithmetic high-priority changes.

PME's corresponding mesh path (IDs 42–51) totals 0.681 ms; screened nonbonded work
is 2.348 ms. FMM must compete with the complete PME decomposition, including its
combined real-space LJ/Coulomb loop, not just with its FFTs.

## 1. Establish the accuracy contract before selecting cheaper settings

Keep the current order-6 implementation as a numerical reference and record two
distinct comparisons: same-order implementation improvements, and matched-PME-
accuracy parameter choices. Do not silently replace order 6 with a cheaper order.

The existing neutral test gives relative force RMS errors of about 8.65e-5 at
order 6, 4.17e-4 at order 5, and 9.93e-4 for PME. This suggests room to tune
accuracy, but these small-system results do not establish STMV accuracy.

For the implementation phase:

- Capture several real STMV configurations, including after integration and
  particles crossing leaf boundaries. Evaluate electrostatic forces separately
  from LJ/bonded forces, so large common forces cannot hide electrostatic error.
- Establish a converged, independent reference: independent periodic Ewald on
  smaller representative systems and converged independently implemented PME or
  periodic Ewald for sampled targets in the large configurations. Match boundary,
  background, self, exclusion, and pair-bond conventions explicitly.
- Compare RMS error, absolute errors and their high quantiles, and errors on
  nearly neutral/low-force particles with a denominator floor. Check energy
  gradients, total-force residual, cell-boundary behavior, and NVE drift.
- Preserve existing force/energy thresholds and bitwise standalone-versus-batch
  behavior. New precision recipes need their own error budget; a covariance
  eigenvalue tolerance is not directly a force-error tolerance.
- Select configuration per simulation, independently of other batch members or
  their retirement. Preserve 64-bit fixed-point atomics; no float atomics.

This follows the matched-accuracy comparison used in the
[GROMACS FMM study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7660746/). That paper
demonstrates important FMM advantages but also reports a dense solvated case
slower than single-GPU PME. It supports the methodology, not a guarantee that
any FMM implementation will win on STMV. Its spherical multipole order is not
numerically interchangeable with this code's Chebyshev order per axis.

## 2. Remove the FP64 and dynamic-array bottlenecks first

Primary locations: `Interpolation::Basis`, `ParticleToMultipole`, `Evaluate`, and
`Calculate` in `src/FMM.cu`.

1. Specialize the hot kernels on `Order`, as well as arithmetic type and whether
   particle potential energy is needed. Use fixed array sizes, compile-time
   indices, and bounded unrolling; inspect register pressure rather than blindly
   unrolling every loop. P2M needs basis values only, not derivative output.
2. Precompute box/level scale factors. Initially preserve double precision in
   coordinate normalization but replace per-particle double division with a
   precomputed reciprocal and multiply/FMA. This isolates the division savings
   before considering lower precision.
3. Evaluate the far-field polynomial with short FP32 partial sums or compensated
   FP32 sums for the existing low-order path, combining into FP64 only at selected
   boundaries. Hoist derivative scale factors out of the nested polynomial loop.
   Keep higher-order arithmetic and the periodic root operator conservative.
4. Add a force-only P2P/L2P specialization on nonlogging MD steps. `potE == nullptr`
   currently suppresses the output but the FMM still computes particle potential.
   Keep potential coefficients needed to obtain gradients; only remove genuinely
   unused final potential arithmetic. EM retains required energy behavior.
5. Evaluate derivative-polynomial preprocessing per leaf if repeated particle
   derivative work remains expensive. Compare the extra coefficient traffic
   against the saved arithmetic; do not assume it is beneficial.

First checkpoint: P2M+L2P below 0.9 ms combined, from 2.66 ms, while satisfying the
same numerical checks. This is a target for accepting the approach, not a
measured speedup. If arrays still generate local traffic, change the thread/data
mapping rather than merely raising register limits.

Compile-time order and register-resident coefficients are also emphasized by
[jz-fmm](https://arxiv.org/html/2609.09307v1). Its gravity benchmarks motivate the
implementation technique; they do not establish periodic molecular-dynamics
performance or accuracy for LIMA.

## 3. Reduce M2L work and redesign its data movement together

Primary locations: `Operators::Compress`, `TransformCoefficients`,
`MultipoleToLocalTensor`, and the downward pass in `Calculate`.

The current default has 216 Chebyshev coefficients and compressed rank 96.
Every leaf processes 189 interactions with a dense 96x96 matrix. For 32768 leaves,
that is approximately 114 billion operations for one matrix-product component,
or 342 billion for the three high/low products, counting FMA as two operations.
Increasing utilization alone leaves a substantial cost. Reduce the operation
count as well as improving the kernel.

Implement and compare two bounded alternatives:

**A. Improve the existing target-owned kernel.** Use a bank-conflict-free
tensor-core shared layout, a genuinely pipelined double buffer, and precomputed
interaction offsets. The current async copy is followed immediately by a wait,
so it does not overlap the next tile's load with the current tile's MMA. Tune
row/column shapes jointly with live accumulator count. Use architecture-compatible
SM80/SM89 CUTLASS/CuTe building blocks where useful; do not depend on Hopper-only
TMA or WGMMA. NVIDIA's
[GEMM guidance](https://docs.nvidia.com/cutlass/latest/media/docs/cpp/efficient_gemm.html)
specifically motivates conflict-free layouts, hierarchical reuse, and pipelining.

**B. Group interactions by translation vector and apply operator-specific
compression.** Factor each compressed operator as `C_delta ≈ U_delta V_delta^T`,
with a rank justified by a force-error budget. A two-factor application costs
approximately `2*r*r_delta` multiply-add terms instead of `r*r`; at r=96 it only
reduces arithmetic when the effective rank is below 48, before tile padding and
overhead. Use dense operators for offsets where factorization loses. Bucket
operators by padded rank and batch many source/target boxes into GEMMs. Batch
projection and expansion across cells as GEMMs too.

The [kiFMM M2L study](https://arxiv.org/html/2408.07436v2) describes individual
operator recompression and transfer-vector grouping. Its measurements are CPU
results; the GPU benefit here is a hypothesis to test. Do not copy its complete
buffering scheme: materializing every leaf interaction would create roughly
2.38 GB of FP32 partial outputs per level at rank 96. Use bounded chunks and a
deterministic reduction, or retain target-owned accumulation.

Keep the current three-product precision recipe initially. Test cheaper precision
only with explicit error evidence; dropping correction products is not a free
optimization. Try ranks 64/80/96 where supported by measured error, rather than
assuming all offsets and all levels need a rank rounded to the same 32-node tile.
Lowering rank from 96 to 64 cuts dense matrix arithmetic to 44%, but its force
accuracy is not yet established.

Intermediate checkpoint: M2L+projection/expansion below 1.5 ms combined. The final
parity budget requires approximately 0.9 ms for translations including hierarchy
work. If two substantially different schedules cannot approach that budget at
required accuracy, advance to the analytical representation in step 5. Do not
spend another round chasing a few percent in the same 96x96 loop.

## 4. Treat Coulomb near field and LJ as one pair-work problem

Primary locations: `ParticleToParticle` and exclusion handling in `src/FMM.cu`,
`NbNonlocalKernel` in `src/EngineKernels.cuh`, and backend dispatch in `Engine.cu`.

Current P2P and LJ consume 3.31 ms in isolated kernel service time. They load
positions, traverse pairs, calculate distances, and accumulate forces separately.
FMM's direct kernel also evaluates both directions of ordinary particle pairs.

First, optimize P2P's mixed-precision tile accumulation and force-only path.
Then prototype a FMM-owned near-field scheduler with two alternatives:

- Target-owned tiles with one accumulated force write per target: simplest
  deterministic reduction, but duplicated pair arithmetic.
- Unique unordered pair tiles, producing both particle forces once, using
  fixed-point tile reductions and integer atomics. Avoid one global atomic per
  elementary pair; measure whether reduced pair arithmetic pays for the writes.

Fuse LJ into the overlapping Coulomb near pairs, applying the existing LJ cutoff
and exclusions. Handle LJ-only residual pairs explicitly. On STMV, the box is
22 nm and FMM leaf width is 22/32 = 0.6875 nm, smaller than the 1.2 nm LJ cutoff.
Therefore the 27-cell FMM near list does not contain every LJ pair. Equating the
lists would silently omit forces. Conversely, FMM near pairs outside the LJ
cutoff still require full Coulomb evaluation. Use the union of the two pair
domains with distinct masks and image identities.

Reuse the existing LJ pair formula and bonded conventions. Preserve periodic
images and nearest-image exclusion subtraction, especially on small trees where
different image offsets wrap to the same cell. Maintain EM's separate handling
of large forces. Prove complete, disjoint ownership of near and far interactions.

Keep PME and FMM separate: backend policy can declare ownership of near-field
work, with shared low-level LJ math. FMM does not call PME or use its mesh. A
FMM-specific pair scheduler may replace the shared LJ kernel only in FMM mode.

Checkpoint: combined near-field work below 2 ms, with a final target around
1.4 ms. If the separate LJ loop remains at 2 ms, consider this an unresolved
architectural cost, not an unrelated kernel outside the optimization scope.

## 5. Analytical Laplace FMM is the deliberate alternative

If the Chebyshev branch misses its translation/particle budget, replace its
expansion representation inside FMM with normalized real solid harmonics,
including P2M/M2M/M2L/L2L/L2P and a corresponding periodic root operator.
Retain the existing backend interface and the validated Chebyshev version as a
reference during development.

A degree-p harmonic expansion has `(p+1)^2` real degrees of freedom. Degree 8
has 81 versus the current 216 volume interpolation nodes; this count alone
does not imply equal accuracy or a speedup. It does offer a representation
tailored to the harmonic Coulomb kernel, with analytic gradients and symmetry
relations. Compare a symmetry-aware direct M2L against rotation-based M2L at
the required order; lower asymptotic complexity is not automatically faster on
this GPU. Re-derive periodic background/self/exclusion behavior, not just the
free-space operator.

The [CUDA M2L paper](https://www.mpinat.mpg.de/634623/Kohnke_2021_IJHPCA.pdf)
demonstrates presorted interaction lists and symmetry-aware spherical-harmonic
operators. Use its structural ideas rather than transplanting assumptions about
GPU generation or expansion order. Mutual pair symmetry halves direct pair
arithmetic; reusing a reverse M2L matrix by transposition does not, by itself,
halve the two directed matrix products.

A surface-equivalent kiFMM with FFT-accelerated local translations is a secondary
alternative if harmonic arithmetic is unfavorable. Such local FFTs would remain
an FMM algorithm, but this would be a representation change: the current
Chebyshev nodes cannot simply be passed to a uniform-grid convolution. Benchmark
only one serious alternative after the compressed branch's checkpoint, rather
than maintaining three unfinished implementations.

## 6. Tune geometry and scheduling after the kernels change

For STMV's 1,066,628 particles, the current depth rule selects 32768 leaves for
integer occupancy targets 33 through 260. Sweeping 64/128/256 does no work at all
to the geometry. The actual mean occupancy is about 32.55, not 64. Adjacent tree
depths change leaf count by a factor of eight.

After measuring the revised near/far balance, use a model such as
`T(B) ≈ a*B + b*N*N/B + c*N` to select promising geometries. If its optimum falls
between octree depths, consider a small multi-root forest or mixed-axis levels.
These require compatible periodic/root operators and shape-specific transfers;
they are not a parameter-only tweak. Do not undertake this redesign before the
new near/far cost ratio is known.

P2P can run after sorting/range construction, independently of the FMM tree
passes. PME and FMM already overlap their electrostatics stream with other force
kernels; there is no entire missing 2 ms of guaranteed overlap to recover.
Consider overlap only when the revised kernels use complementary resources.
CUDA graphs, less sorting, and buffer-lifetime reductions are finishing work.
The existing sort is about 0.1 ms, so none is the principal route to parity.

## Final acceptance and decision budget

The following are engineering allocations needed to support parity, not forecasts
or claimed speedups. Avoid adding independent maximum savings to predict a result.

| Component | Target allocation |
|---|---:|
| All translation, projection, and hierarchy work | 0.90 ms |
| P2M and L2P | 0.55 ms |
| Combined LJ and Coulomb near field | 1.40 ms |
| Gather/sort/metadata | 0.15 ms |
| Bonded/integration/common work | 0.35 ms |
| Remaining overhead allowance | 0.35 ms |
| Total planning budget | 3.70 ms |

The kernel allocations guide which design to keep. The authoritative gate is
end-to-end time, because these kernels overlap and compete for GPU resources.

Implementation order: accuracy contract; precision/order specialization;
M2L compression/scheduling; unified near-field; analytical representation if
the translation budget is missed; final geometry/scheduling tuning. Record
isolated effects and repeat full validation when changes affect shared paths.

For acceptance, alternate PME and FMM for at least two 300-step rounds with no
competing GPU workload, same settings and logging, and then repeat the longer
3000-step workload. Require FMM no slower than PME within repeatability, and
prefer a margin rather than a single winning sample. Include initialization
time/memory separately. Run independent accuracy, NVE, determinism, batch
retirement, EM, and sanitizer checks; produce final full Nsight reports as
required by CLAUDE.md. A fixed-order gain and a matched-accuracy parameter gain
must be reported separately.

## Implementation checkpoint — 2026-10-12

The first implementation pass keeps order 6 and the leaf target of 64. It adds
compile-time interpolation specialization, host coordinate scales, separable
compensated L2P, tensor P2M and shared-basis transforms, force-only nonlogging
steps, per-offset rank-16/32 recompression with dense fallback, and geometry-based
parent-neighborhood partitioning. Five-level independent Ewald tests, same-depth
batch retirement, logging-cadence determinism and extreme charge scales were
added. Existing force and energy thresholds were not relaxed.

Fresh full Nsight reports compare the preserved previous optimized executable
with this implementation. Last nonlogging-step kernel durations on the RTX 4090:

| Work | Previous ms | Current ms |
|---|---:|---:|
| M2L, all levels | 4.031 | 3.089 |
| Projection and expansion, all levels | 0.698 | 0.542 |
| P2M | 0.911 | 0.328 |
| L2P and exclusions | 1.749 | 0.487 |
| Direct P2P | 1.301 | 1.112 |
| Shared LJ | 1.992 | 2.000 |

The serialized kernel sum fell from 11.668 to 8.554 ms (26.7%). This sum includes
final state storage and does not account for overlap or host/launch overhead;
wall-clock benchmarks remain the acceptance gate. PME's corresponding sum is
3.377 ms. Reports are `build/x64-Release/ncu/fmm-plan-pass1-full.ncu-rep`,
`fmm-plan-pass1-baseline-full.ncu-rep` and `pme-plan-pass1-full.ncu-rep` in that
directory. Detailed measurements and validation are recorded in
`build/FMM-plan-implementation-results.md`.

The translation milestone is missed: M2L plus basis transforms alone still cost
3.63 ms, above the 1.5 ms intermediate gate. Additional shared-memory padding,
splitting the leaf M2L into three partitions, and overly strict recompression
were measured and rejected. The retained compression uses an explicit additional
relative Frobenius residual budget, documented in ELECTROSTATICS.md.

The next work should therefore follow the structural branch of this plan:
prototype analytical spherical-harmonic translations and compare their full
particle/translation cost at independently matched accuracy; in parallel design
the union of LJ and Coulomb near-field neighborhoods for a shared traversal.
The present leaf neighborhood does not cover the entire LJ cutoff, so simply
adding LJ to P2P would omit interactions. Preserve exclusions, periodic image
multiplicity and deterministic accumulation. Finish the production STMV sampled
force comparison before selecting a lower order or a different leaf target.
The new sparse five-level reference tests exercise production geometry but do
not establish matched-accuracy performance on the full STMV charge distribution.

PME parity remains unachieved. No parameter change or kernel-sum comparison is
being used to claim that the end-to-end acceptance gate has passed.
The final two alternating 300-step rounds measured previous FMM at 12.679 and
12.762 ms/step, current FMM at 10.121 and 9.084, and PME at 3.609 and 3.725.
Mean FMM time improved 24.5%, but remains 2.62 times PME. A 3000-step current
FMM run completed at 10.308 ms/step. An earlier unstable timing set and its
unexplained 64.942 ms/step baseline outlier are retained in the results document.
