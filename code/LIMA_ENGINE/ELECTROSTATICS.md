# Electrostatics backends

Select the solver in `sim_params.txt`:

```ini
enable_electrostatics=true
longrange_method=FMM
fmm_order=6
fmm_leaf_particles=64
```

`longrange_method=PME` is the default, including for existing parameter files.
The C++ setting is `params.longrange_method = LongRangeElectrostaticsMethod::FMM`.
MDP input accepts `coulombtype=PME` or `coulombtype=FMM`, plus the two optional
`fmm_...` settings above. Both solvers currently require periodic boundaries and
the engine's existing orthorhombic, integer-nanometer box dimensions.

## Ownership

`Electrostatics::Create` is the only place that selects an implementation.
`Electrostatics::Controller` accepts current superclusters and writes the shared
force accumulator in MD, or the persistent-cluster force/energy buffer in EM.
The factory also selects the short-range policy; the controller reports optional
capacity status. The engine
owns scheduling, stream synchronization, and batch lifecycle.

`PME.cu/.cuh` own the grid, charge blocks, cuFFT plans, Ewald screening parameter,
and reciprocal/self terms. `FMM.cu/.cuh` own their tree, interpolation operators,
particle sorting, near-field pairs, exclusions, and periodic correction. FMM
does not include or call PME. In FMM mode the shared nonbonded kernel computes
Lennard-Jones interactions; FMM computes the complete Coulomb interaction.
`cutoff_nm` still controls LJ neighbors, but does not truncate FMM electrostatics.

## FMM calculation

This is a uniform octree, tensor Chebyshev interpolation FMM, following the
interpolation formulation described by
[Fong and Darve (2009)](https://doi.org/10.1016/j.jcp.2009.08.031).
Each timestep performs:

1. Wrap positions and stably sort particles into leaf cells on the GPU.
2. Form equivalent node charges (P2M) and translate them up the tree (M2M).
3. Apply multipole-to-local translations between children of neighboring parents
   that are not immediate neighbors (M2L).
4. Translate local potentials down the tree (L2L), then interpolate potentials
   and their analytic derivatives at the particles (L2P).
5. Sum the 27 neighboring leaf cells directly (P2P).

The tree accounts for the nearest 27 copies of the root box. At initialization,
a double-precision Ewald calculation builds a smooth root operator for the
periodic Green function minus those 27 copies. This includes conducting
(tin-foil) boundary conditions, the regular self-image potential, and a uniform
neutralizing background for a net charge. The root operator is reused; no
particle Ewald sum or FFT runs during timesteps.

Topology exclusions subtract the nearest-image Coulomb interaction of excluded
pairs. Their periodic images remain. The existing bonded kernels, including
their current pair-bond behavior, remain responsible for bonded forces. The
PME implementation retains its existing conventions; exact energy equality
between methods is not promised for net-charged or excluded-pair systems.

`fmm_order` is the number of interpolation nodes **per axis**, from 4 through 8.
It controls approximation accuracy and cost, and is not an error tolerance.
`fmm_leaf_particles` is a positive target mean occupancy, not a hard capacity.
The smallest tree has two cells per axis. Dense or uneven distributions can
require a smaller occupancy target to reduce direct work.

Simulations are grouped internally by tree depth and share immutable operators.
Each simulation uses the same depth and particle summation order when run alone
or in a batch. Retiring a member does not change a survivor's approximation.
Different methods or FMM settings are incompatible scheduler batches.

## Cost and validation

The uniform tree uses a shared compressed basis for M2L and separable one-axis
M2M/L2L transfers. On Ampere and newer GPUs, M2L uses tensor cores with scaled
high/low FP16 components and compensated accumulation. For shared rank 96 at
orders 5–6, individual operators use two factors of padded rank 16 or 32 when
complete-pivot cross approximation meets a relative Frobenius residual bound
of `10^(-(order+6)/2)`. Other offsets retain dense compressed operators. This is
an additional approximation, checked against independent force/energy references.
The factors are reused across all tree levels. Geometry selects a fixed partition
of parent neighborhoods, with deterministic reduction between partitions.

Orders 4–6 also use high/low tensor products for tiled P2M and for shared-basis
projection/expansion at larger levels. Power-of-two input normalization preserves
charge range. Older GPUs retain scalar deposition and dense translations.
Interpolation order is specialized at compile time; particle coordinate scales
are computed on the host. L2P uses separable contractions and compensated float
sums at orders 4–6; orders 7–8 retain double interpolation arithmetic. P2P shares
source particles and keeps double accumulation between short float pair tiles.
Nonlogging MD steps skip final potential-energy work; EM always computes energy.
Forces use deterministic fixed-point accumulation; batch membership and logging
cadence do not change the trajectory.

Dense operator tables are retained for the fallback and occupy about 5.5 MiB at
order 4, 63 MiB at order 6, and 352 MiB at order 8. Compressed operators, packed
tensor data, tree coefficients, and particle buffers require additional memory.
Increasing order remains expensive. Strongly anisotropic boxes or clustered
systems need their own accuracy and performance checks.

On an RTX 4090, the latest alternating 300-step STMV runs measured 9.08–10.12
ms/step for order 6 with a leaf target of 64, versus 3.61–3.73 ms/step for PME.
The previous optimized FMM measured 12.68–12.76 ms/step in those same rounds;
this pass reduced mean wall time by 24.5%. FMM remains about 2.6 times slower
than PME and has **not reached performance parity**. Translation and the
separate LJ/Coulomb near-field work remain substantial costs. These timings
exclude operator initialization and do not establish matched-accuracy STMV
performance. A 3000-step FMM run completed at 10.31 ms/step.

Build `agenttesting` and run `agenttesting --fmm-tests`. The suite uses a separate
double-precision periodic Ewald reference and checks forces, energies, increasing
order, rectangular/net-charged boxes, exclusions across periodic boundaries,
MD/EM paths, batch retirement with equal and different tree depths, parameter
roundtrips, disabled electrostatics, invalid orders, and NVE energy conservation.
It samples independent Ewald forces on five-level neutral and charged rectangular
trees, and checks extreme charge scales and energy-logging cadence. It also
checks PME against the same neutral reference. `--engine-batch-tests` exercises
the existing PME and disabled-electrostatics paths after the interface refactor.
