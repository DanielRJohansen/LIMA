## Compiling takes forever. Can we optimize a bit?

## For Windows (and linux) installs, it'd be nice to add lima to path. This would mean an installer for windows atleast, no?


## Investigate startup latency for small systems (say t4 in 7^3nm^3 box). Defined as time from command to first step.

## The tests: Test T4 (MD performance) and LoadT4 (setup performance) are always noisy. Find a better way to isolate them from other tests.

# Performance ideas. These are not set in stone todo's, rather concepts that are interesting to explore to see if they give any speed. Some makes the most sense together
## Instead of Env having a thread dedicated to RunSimulation, Engine instead get's it's own thread, so when engine starts an mdrun it literally only steps & updates the renderpipelin/simstatus. Another benefit is that Engine would persist, meaning slightly lower latency for setup. It's import that Engine free's its memory when it has no work, such that we are hogging 10 GB VRAM, while env is just doing preprocessing.
## Look into all structs/classes used used in Kernels, and see if they could be improved. Especially with advanced stuff such as alignment, __restrict__, or simply ordering/datatypes of members
## EmitQuarterEntriesKernel: now ~0.77 ms per STMV rebuild (was 1.37). The cost was the bonded-pair exclusions (~1 active lane, dependent global loads): own particle ids now in shared memory, query ids loaded once per quarter, one pass over each bonded set for all 4 query particles. Rejected: int4 reads of the set (slower). What's left is mostly writing ~300 MB of QuarterEntries; shrinking that means changing NB's entry format.
## FindNeighborsKernel: now ~2.7 ms per STMV rebuild (was 3.9). Done: own-particle vs candidate-sphere prefilter (196 -> 140 pair-tested candidates/SC, 90 become neighbors), branch-free pair test, __launch_bounds__(128, 8) for 64 regs. The kernel is issue-bound (~68% issue active), so latency tricks don't help: 2 candidates per warp iteration and a candidate-particles vs own-sphere filter were both slower. Left: cell enumeration (~74 cells/SC, ~7 of 32 lanes active) and the bonded check (~8%).
## nstlist vs list-buffer tuning: a bigger buffer and fewer rebuilds may be cheaper overall. Parameter/accuracy call for you.
## NB kernel flattened entries: known ~4% SIMT gain left over from densetasks.
## Rebuild-step host syncs: ~1.1 ms GPU idle per nlist rebuild (~55 us/step on STMV). genericErrorCheck(stream) after FindNeighbors/Emit syncs, CubWrappers::ExclusiveScan always syncs (its shared TempStorage relies on it), and count readbacks (nSuperclusters, boundaries, nQuarterEntries, overflow) each force a round trip. Use NoSync checks, a stream-ordered scan, and fewer/merged readbacks (or capacity-sized buffers so counts stay on GPU).


# Kernel graphs (prototyped, see stash "Tried cuda graph..."): T4 -11%, STMV +2.5%. Recipe:
- Engine.cuh: `std::array<cudaGraphExec_t, 2> stepGraphs{}` (indexed by logData) + `void LaunchStepGraph(cudaGraph_t, bool logData)`. Destroy the execs in ~Engine.
- _deviceMaster: after the `MakeSuperClusterTasksGPU` check and before `ForkStreamsFromMainStream()`, if `!emvariant && nSuperclusters > 0 && snf_select.empty()` (EM syncs mid-step, SNF unchecked): `cudaStreamBeginCapture(cudaStreams[0], cudaStreamCaptureModeThreadLocal)`. ThreadLocal so the log-drain thread's CUDA calls don't break the capture.
- The existing fork/join events pull the side streams (PME, bonds) into the capture; cuFFT exec captures fine. After the integrate launch: `cudaStreamEndCapture(cudaStreams[0], &graph)` -> LaunchStepGraph.
- LaunchStepGraph: if an exec exists try `cudaGraphExecUpdate(exec, graph, &info)`; on failure (topology changed after a task rebuild) `cudaGetLastError()`, destroy and re-instantiate. Then `cudaGraphLaunch(exec, cudaStreams[0])`, `cudaGraphDestroy(graph)`. Capturing every step means changed args (step, buffer pointers after Expand, live edits) need no tracking; ~10 us host/step, only ~4 re-instantiations per 2000 T4 steps.
- Measuring: nsys needs `--cuda-graph-trace=node` to show graph kernels, but node tracing perturbs scheduling (STMV looked faster under nsys, was slower in wall-clock). Trust wall-clock A/B (`--t4 2000 1`, `--stmv 300`) of the same binary.
- Findings: launch overhead is not the issue (host runs 0.25-41 ms ahead). Graphs fix WDDM submit batching that idled T4 normal steps 24/90 us -> 54 us steps. STMV loses because NB, PME spreading and bonds truly overlap and slow NB, which already fills the GPU (CUDA_DEVICE_MAX_CONNECTIONS=1: graph == no graph). Reordering launches (NB before PME) without graphs changed nothing.
- Next: encode the order (NB waits for an event after PME's ChargeblockDistributeToGrid, bonds wait for NB) so STMV is neutral and T4 keeps most of the gain; or graphs only below a system size.