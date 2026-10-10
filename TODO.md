## Compiling takes forever. Can we optimize a bit?

## For Windows (and linux) installs, it'd be nice to add lima to path. This would mean an installer for windows atleast, no?


## Investigate startup latency for small systems (say t4 in 7^3nm^3 box). Defined as time from command to first step.

## The tests: Test T4 (MD performance) and LoadT4 (setup performance) are always noisy. Find a better way to isolate them from other tests.

# Performance ideas. These are not set in stone todo's, rather concepts that are interesting to explore to see if they give any speed. Some makes the most sense together
## Instead of Env having a thread dedicated to RunSimulation, Engine instead get's it's own thread, so when engine starts an mdrun it literally only steps & updates the renderpipelin/simstatus. Another benefit is that Engine would persist, meaning slightly lower latency for setup. It's import that Engine free's its memory when it has no work, such that we are hogging 10 GB VRAM, while env is just doing preprocessing.
## Rework the Kernel dispatch ordering, or perhaps even test a Kernel-graph dispatch. I dont know if we'd need 2 graphs, 1 for the logging step and one for the basic step?
## Look into all structs/classes used used in Kernels, and see if they could be improved. Especially with advanced stuff such as alignment, __restrict__, or simply ordering/datatypes of members
## Figure out if the GPU-busy time is near 100%, if not figure out what is causing the GPU to be idle.
## Drop the end-of-step host sync in MD, so the GPU no longer idles between steps; other streams wait on a GPU event instead. ~−130 µs/step (~4%), low risk.
## ClusteringKernel occupancy: 178 regs, 16% occupancy, 32-thread blocks, 1.5 ms per rebuild. Cut registers or widen blocks.
## EmitQuarterEntriesKernel: only 6.8/32 threads active, 1.4 ms per rebuild. Give each warp several superclusters, or a flatter loop.
## Rebuild host round trips: ~1 ms GPU idle per rebuild from count readbacks. Size buffers to capacity, or chain on the GPU.
## FindNeighborsKernel algorithm: 54% of its 16×16 particle tests find nothing, 5.5 ms per rebuild. Reuse the previous list as the candidate set, or pre-filter at particle level.
## nstlist vs list-buffer tuning: a bigger buffer and fewer rebuilds may be cheaper overall. Parameter/accuracy call for you.
## BondgroupsKernel packing: 18/32 threads active, latency-bound. Fill warps with several bondgroups per block so position loads overlap. ~−50 to 100 µs.
## NB kernel flattened entries: known ~4% SIMT gain left over from densetasks.
## Rebuild-step host syncs: ~1.1 ms GPU idle per nlist rebuild (~55 us/step on STMV). genericErrorCheck(stream) after FindNeighbors/Emit syncs, CubWrappers::ExclusiveScan always syncs (its shared TempStorage relies on it), and count readbacks (nSuperclusters, boundaries, nQuarterEntries, overflow) each force a round trip. Use NoSync checks, a stream-ordered scan, and fewer/merged readbacks (or capacity-sized buffers so counts stay on GPU).
