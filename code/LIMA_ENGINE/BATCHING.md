# Engine batches

`Engine({ simulation.get(), otherSimulation.get() })` accepts nonowning pointers
that must outlive the engine. `Environment` groups prepared jobs into batches of
up to four simulations and keeps their owners alive until the engine is destroyed.

The scheduler compares the actual simulations after preprocessing and configuration.
It starts with the oldest prepared job and selects compatible jobs from a bounded
eight-job preparation lookahead. The scheduler counts submitted jobs until their
preprocessing finishes. When compatible candidates do not fill a batch, it waits
for that count to reach zero or for the bounded preparation buffer to fill.
`mustRunAlone`, displayed, stepwise and preparation-only jobs form standalone
scheduling boundaries. A single submission never waits for future work.
Preparation and postprocessing still overlap the single GPU worker.

Results report the actual `batchId` and `batchSize`. Each member's `engineTime`
ends at its own retirement (zero for zero-step members), excluding engine setup.
Results are delivered to independent postprocessing callbacks after the entire
batch finishes. Preprocessing or postprocessing failures affect only that job;
an engine execution failure is reported to every member of its batch.

## Storage and indexing

`EngineBatchData` owns the GPU allocations. `EngineSimulationData` holds each
member's host bookkeeping, status and ranges. `SimulationDeviceData` supplies
the small amount of simulation-specific state used by kernels.

- Dense particle IDs, pcluster IDs, bond-group references, supercluster IDs,
  query offsets and result offsets refer to batch-global storage.
- Particle slots (`pcluster * 4 + lane`) and dense particle IDs are distinct.
- Particle indices inside a pcluster and bond indices inside a bond group remain
  local. Invalid IDs remain invalid when packing.
- Packing rebases device metadata and exclusion sets without changing the
  original host topology.
- Every simulation owns a disjoint range of spatial bins. Periodic neighbor
  lookup wraps local box coordinates before adding the simulation's bin offset.
- Logging uses concatenated per-simulation ring buffers, each with its own stride.

All members must start at step zero, have equal box dimensions and agree on
`SimParams` except `dt` and `n_steps`. Particle counts, topology, degrees of freedom
and PME self-energy corrections are simulation-specific. The engine and scheduler
use the same compatibility check in `BatchCompatibility.h`.

## Execution and completion

The existing pipeline streams remain shared by the batch. Force kernels and
integration launch over active work IDs rather than launching once per simulation.
The timestep is shared; integration reads each member's `dt` and thermostat scalar.

Retirement happens after the member's last integration and host measurements.
It drains pending log entries, copies final state once, and removes the member
from the active pcluster, bond-group and supercluster work lists. Its particle and
topology ranges remain allocated and stable. A subsequent neighbor-list rebuild
also drops its superclusters and NB queries. Repeated finalization and stepping a
finished batch are harmless.

PME assigns dense grid slots to active members, uses batched R2C/C2R plans, and
normalizes by one simulation's grid size. Slot mapping and FFT plans are rebuilt
when active membership changes. Equal box sizes currently imply equal PME grids.
Green's-function values are shared; charges, interpolation and corrections are not.

Clustering/task sorting, force gathering and integer PME charge accumulation keep
their existing ordering. Kinetic energies are computed in one batch launch, then
reduced separately with the existing Thrust reduction for each member. EM stopping
also examines each member separately. Further reduction/launch optimizations can
be made independently of this storage migration.

`GetRunStatus(id)`, `OffloadPclusterState(id)` and particle-control setters address
individual members; the default ID is zero. `IsFinished()` describes the whole
batch. The live editor explicitly selects `EngineRunMode::Interactive`, restricted
to one member, so its zero-step setting and changes between MD and EM continue
to work without terminating the engine.

## Validation

Run `limatest --engine-batch-tests` for direct engine regressions. These compare
coordinates, velocities, forces, energies, trajectories and temperatures exactly
against independent runs, with different particle counts, timesteps and stopping
steps. They cover PME on/off, MD/EM, cross-pcluster bonds, zero-step members,
partial log buffers, repeated finalization, frozen retired GPU coordinates and
the interactive path.

Run the same command under CUDA Compute Sanitizer's `memcheck` tool to check GPU
accesses. The existing full `limatest` suite remains the physics regression suite.
Run `limatest --environment-batch-tests` for scheduler regressions, including four
T4 simulations matching an isolated reference, different particle/step counts and
timesteps, zero-step members, incompatible parameters, standalone/preparation-only
jobs, bounded queue draining and callback failures.
