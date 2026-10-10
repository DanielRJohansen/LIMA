## Compiling takes forever. Can we optimize a bit?

## For Windows (and linux) installs, it'd be nice to add lima to path. This would mean an installer for windows atleast, no?


## Investigate startup latency for small systems (say t4 in 7^3nm^3 box). Defined as time from command to first step.

## The tests: Test T4 (MD performance) and LoadT4 (setup performance) are always noisy. Find a better way to isolate them from other tests.

## PME self-energy correction is applied once per particle instead of once per system. PME::Controller::CalcEnergyCorrection returns the whole system's Ewald self-energy (-kappa/sqrt(pi) * C * sum q^2), and InterpolateForcesKernel adds it to every charged particle's potE (then halves it), so the logged potential energy is off by roughly N/2 times the self-energy. Forces are unaffected, but energies and the VC/drift values computed from them are. The per particle correction should be -kappa/sqrt(pi) * C * q_i^2, added after the halving.
