#include "KernelConstants.cuh"
#include "PhysicsUtils.cuh"







void Engine::CopySimulationToHost() {
	assert(sim_dev);
	sim_dev->boxState.CopyDataToHost(*simulation->box);	
}

void Engine::verifyEngine() {
	LIMA_UTILS::genericErrorCheck("Error before engine initialization.\n");

	Int3 dim = simulation->box->boxparams.boxSize;
	assert(dim.x < 1024 && dim.y < 1024 && dim.z < 1024 && "Neighborlist cannot handle such large gridnode_ids");

	if constexpr (ENABLE_ES_LR) {
		if (simulation->simParams.enable_electrostatics && simulation->simParams.bc_select != PBC) {
			throw std::invalid_argument("Electrostatics only supported with PBC at the current time");
		}
	}
}





ForceEnergyInterims::ForceEnergyInterims(int nBondgroups, int nParticles, int nPclusters) {
	if (nPclusters > 0) {
		const size_t byteSize = sizeof(ForceEnergy) * nPclusters * PersistentCluster::maxParticles;
		cudaMalloc(&bonded, byteSize);
		cudaMalloc(&snf, byteSize);
		cudaMalloc(&pme, byteSize);

		cudaMemset(bonded, 0, byteSize);
		cudaMemset(snf, 0, byteSize);
		cudaMemset(pme, 0, byteSize);
	}

	if (nBondgroups > 0) {
		cudaMalloc(&forceEnergiesBondgroups, sizeof(ForceEnergy) * BondGroup::maxParticles * nBondgroups);
		cudaMemset(forceEnergiesBondgroups, 0, sizeof(ForceEnergy) * BondGroup::maxParticles * nBondgroups);
	}

	if (nParticles > 0) {
		//cudaMalloc(&nbNonlocal, sizeof(ForceEnergy) * nParticles);
		//cudaMalloc(&bonded, sizeof(ForceEnergy) * nParticles);
		//cudaMalloc(&forceEnergySNF, sizeof(ForceEnergy) * nParticles);

		//cudaMemset(nbNonlocal, 0, sizeof(ForceEnergy) * nParticles);
		//cudaMemset(bonded, 0, sizeof(ForceEnergy) * nParticles);
		//cudaMemset(forceEnergySNF, 0, sizeof(ForceEnergy) * nParticles);
	}
}

void ForceEnergyInterims::Free() const {
	if (forceEnergiesBondgroups != nullptr) {
		cudaFree(forceEnergiesBondgroups);
	}

	if (bonded != nullptr) {
		cudaFree(bonded);
		cudaFree(snf);
		cudaFree(pme);
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

