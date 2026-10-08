## Compiling takes forever. Can we optimize a bit?

## For Windows (and linux) installs, it'd be nice to add lima to path. This would mean an installer for windows atleast, no?


## Investigate startup latency for small systems (say t4 in 7^3nm^3 box). Defined as time from command to first step.

## Give EM the 4x4 quarter skipping that MD has. MD's NbNonlocalKernel only computes the 4x4 blocks of particle pairs (a quarter of each supercluster against a quarter of the other) that have a pair within range, roughly half of them. NbNonlocalEmKernel still computes all 256 pairs of every listed supercluster pair, which made EM ~13% slower per step once the neighbor search stopped missing pairs. It can't reuse the MD kernel directly, since EM forces can exceed NbForceAccumulator's fixed-point range, so it needs a variant that stores and sums its results in a fixed order.