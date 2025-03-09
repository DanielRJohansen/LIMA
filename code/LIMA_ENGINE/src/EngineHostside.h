#include "KernelConstants.cuh"
#include "PhysicsUtils.cuh"







std::unique_ptr<Simulation> Engine::takeBackSim() {
	assert(sim_dev);
	sim_dev->boxState.CopyDataToHost(*simulation->box_host);
	return std::move(simulation);
}

void Engine::verifyEngine() {
	LIMA_UTILS::genericErrorCheck("Error before engine initialization.\n");

	const int nBlocks = simulation->box_host->boxparams.boxSize;
	assert(nBlocks * nBlocks * nBlocks < INT32_MAX && "Neighborlist cannot handle such large gridnode_ids");

	if constexpr (ENABLE_ES_LR) {
		if (simulation->simparams_host.enable_electrostatics && simulation->simparams_host.bc_select != PBC) {
			throw std::invalid_argument("Electrostatics only supported with PBC at the current time");
		}
	}
}





ForceEnergyInterims::ForceEnergyInterims(int nCompounds, int nTinymols, int nSolventblocks, int nBondgroups) {
	if (nCompounds > 0) {
		const size_t byteSize = sizeof(ForceEnergy) * nCompounds * MAX_COMPOUND_PARTICLES;
		cudaMalloc(&forceEnergyFarneighborShortrange, byteSize);
		cudaMalloc(&forceEnergyImmediateneighborShortrange, byteSize);
		cudaMalloc(&forceEnergyBonds, byteSize);
		cudaMalloc(&forceEnergiesPME, byteSize);
		
		cudaMemset(forceEnergyFarneighborShortrange, 0, byteSize);
		cudaMemset(forceEnergyImmediateneighborShortrange, 0, byteSize);
		cudaMemset(forceEnergyBonds, 0, byteSize);
		cudaMemset(forceEnergiesPME, 0, byteSize);		
	}

	if (nBondgroups > 0) {
		cudaMalloc(&forceEnergiesBondgroups, sizeof(ForceEnergy) * SolventBlock::MAX_SOLVENTS_IN_BLOCK * nSolventblocks);
		cudaMemset(forceEnergiesBondgroups, 0, sizeof(ForceEnergy) * SolventBlock::MAX_SOLVENTS_IN_BLOCK * nSolventblocks);
	}

	if (nTinymols > 0) {
		const size_t byteSize = sizeof(ForceEnergy) * SolventBlock::MAX_SOLVENTS_IN_BLOCK * nSolventblocks;
		cudaMalloc(&forceEnergiesCompoundinteractions, byteSize);
		cudaMalloc(&forceEnergiesTinymolinteractions, byteSize);
		cudaMalloc(&forceEnergiesTinymolBondgroups, byteSize);

		cudaMemset(forceEnergiesCompoundinteractions, 0, byteSize);
		cudaMemset(forceEnergiesTinymolinteractions, 0, byteSize);
		cudaMemset(forceEnergiesTinymolBondgroups, 0, byteSize);
	}
}

void ForceEnergyInterims::Free() const {
	if (forceEnergyFarneighborShortrange != nullptr) {
		cudaFree(forceEnergyFarneighborShortrange);
		cudaFree(forceEnergyImmediateneighborShortrange);
		cudaFree(forceEnergyBonds);
		cudaFree(forceEnergiesPME);
		cudaFree(forceEnergiesBondgroups);
	}

	if (forceEnergiesCompoundinteractions != nullptr) { // The buffers are never allocated in some sims
		cudaFree(forceEnergiesCompoundinteractions);
		cudaFree(forceEnergiesTinymolinteractions);
		cudaFree(forceEnergiesTinymolBondgroups);
	}

	LIMA_UTILS::genericErrorCheck("Error during CompoundForceEnergyInterims destruction");
}








constexpr std::array<float, 2 * DeviceConstants::BSPLINE_LUT_SIZE> PrecomputeBsplineTable()
{
    const int N = DeviceConstants::BSPLINE_LUT_SIZE;
    std::array<float, 2 * N> result{};

    for (int i = 0; i < N; i++)
    {
        const double f = static_cast<double>(i) / static_cast<double>(N-1);

        // w0 = (1 - f)^3 / 6
        const double w0 = (1 - f) * (1 - f) * (1. - f) / 6.;
        const double w1 = (4. - 6. * f * f + 3. * f * f * f) / 6.;

        // Store in array: [w0, w1]
        result[i] = static_cast<float>(w0);
        result[N + i] = static_cast<float>(w1);
    }

    return result;
}

// Precomputes ERFC-related scalars from 0 to cutoffNM
std::array<float, DeviceConstants::ERFC_LUT_SIZE> PrecomputeErfcForcescalarTable(float cutoffNM) {
	const float ewaldKappa = PhysicsUtils::CalcEwaldkappa(cutoffNM);

	std::array<float, DeviceConstants::ERFC_LUT_SIZE> result{};

	for (int i = 0; i < DeviceConstants::ERFC_LUT_SIZE; i++)	{
		const double fraction = static_cast<double>(i) / static_cast<double>(DeviceConstants::ERFC_LUT_SIZE - 1);
		const double correspondingDistance = fraction * cutoffNM;

        const double erfcTerm = erfc(correspondingDistance * ewaldKappa);
        const float scalar = erfcTerm + 2. * ewaldKappa / sqrt(PI) * correspondingDistance * exp(-ewaldKappa * ewaldKappa * (correspondingDistance* correspondingDistance));

		result[i] = scalar;
	}
	return result;
}

std::array<float, DeviceConstants::ERFC_LUT_SIZE> PrecomputeErfcPotentialscalarTable(float cutoffNM) {
	const float ewaldKappa = PhysicsUtils::CalcEwaldkappa(cutoffNM);

	std::array<float, DeviceConstants::ERFC_LUT_SIZE> result{};

	for (int i = 0; i < DeviceConstants::ERFC_LUT_SIZE; i++)	{
		const double fraction = static_cast<double>(i) / static_cast<double>(DeviceConstants::ERFC_LUT_SIZE - 1);
		const double correspondingDistance = fraction * cutoffNM;

		const float scalar = erfc(correspondingDistance * ewaldKappa);
		result[i] = scalar;
	}
	return result;
}

