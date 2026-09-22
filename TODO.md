# Multi-simulation batching
Goal: Improve GPU utilization and aggregate simulation throughput for small systems (<20³ nm³) by processing multiple independent simulations in lockstep within the same kernel launches.
Approach: Treat multiple simulations as one batch. Rather than separate CUDA streams or launches, the Engine flatten equivalent work across simulations into shared workloads:
-  Concatenate pclusters, superclusters, bond groups, NB tasks, results, etc. across simulations. 
-  Prefer batch-global indices/offsets rather than adding a literal [batchIndex] dimension to every buffer. 
-  Launch kernels over the total work across all simulations, e.g. NbNonlocalKernel<<<totalNTasks,...>>>. 
-  Work items that require simulation-specific state carry/derive a simulationId. 
-  Keep simulations synchronized at the same pipeline/timestep stage. 
PME: Group simulations with compatible PME grid dimensions and use batched cuFFT transforms. Reciprocal-space calculations and normalization remain independent per simulation.
Expected benefits: Better SM utilization for small simulations, fewer total kernel launches, reduced launch/synchronization overhead, and potentially much better utilization of small FFTs. 
Key constraints: Preserve deterministic results independently for every simulation. No interactions or reductions may cross simulation boundaries. Handle different particle/task counts efficiently without padding everything to the largest simulation. 
Instead of initializing Engine with 1 simulation, it will be initialized with a batch of simulations. These simulations may not terminate at the same time, at which point engine should remove nb-tasks related to a terminated simulation, and should no longer compute forces or integrate that simulation.

Batching rules:
For sim's to be batched they must share the following criteria:
boxsize, emvariant, stepspernlistupdate, steps_per_temperature_measurement, bcselect, enableEs, cutoffNM, applythermostat, ref_t, save_energy and more.
In other words, the only params that'll really be different for each sims are the particles, topology, dt and nSteps. All but nSteps should be trivial to implement, whereas nSteps means some sim's needs to stop before other.

First milestone: Get engine up and running with this new paradigm. Dont focus on the batching logic in environment, start by just making batches of size1.


## Add simulation tabs to Display. Tabs automatically appear/disappear as rendercontexts are made/freed. 