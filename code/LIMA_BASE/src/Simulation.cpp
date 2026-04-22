#include "Simulation.cuh"

#include <functional>
#include <type_traits> // For std::is_integral, std::is_floating_point, and static_assert
#include <filesystem>
#include "Filehandling.h"
#include "MDFiles.h"
#include <fstream>

#include <concepts>









Box::Box(Float3 boxSizeNM) {
	if ((NodeIndex(boxSizeNM.ToInt3()).toFloat3() - boxSizeNM).len() > 0.0000001) {
		throw std::invalid_argument("Boxsize was not an integer");
	}

	boxparams.boxSize = boxSizeNM.ToInt3();
	//solventblockgrid_circularqueue = SolventBlocksCircularQueue::createQueue(boxparams.boxSize);
}





Simulation::Simulation(const SimParams& params) :
	simParams{ params }
{
	box = std::make_unique<Box>();
}

Simulation::Simulation(const SimParams& params, std::unique_ptr<Box> _box) :
	simParams{ params }
{
	box = std::move(_box);
}

void Simulation::PrepareDataBuffers() {
	// Allocate buffers. We need to allocate for atleast 1 step, otherwise the bootstrapping mechanism will fail.
	const auto n_steps = std::max(simParams.n_steps, uint64_t{ 1 });
	// Standard Data Buffers 
	{
		// Permanent Outputs for energy & trajectory analysis
		const int nPclusters = box->persistentClusters.size();
		const int particlesUpperbound = nPclusters * PersistentCluster::maxParticles;
		const size_t n_datapoints = particlesUpperbound * n_steps / simParams.data_logging_interval;
		const auto datasize_str = std::to_string((float)((2. * sizeof(float) * n_datapoints + sizeof(Float3) * n_datapoints) * 1e-6));
		
		//m_logger->print("Malloc " + datasize_str + " MB on host for data buffers\n");


		potE_buffer = std::make_unique<ParticleDataBuffer<float>>(particlesUpperbound, n_steps, simParams.data_logging_interval, nPclusters);
		vel_buffer = std::make_unique<ParticleDataBuffer<float>>(particlesUpperbound, n_steps, simParams.data_logging_interval, nPclusters);
		forceBuffer = std::make_unique<ParticleDataBuffer<Float3>>(particlesUpperbound, n_steps, simParams.data_logging_interval, nPclusters);
		traj_buffer = std::make_unique<ParticleDataBuffer<Float3>>(particlesUpperbound, n_steps, simParams.data_logging_interval, nPclusters);
		
		temperature_buffer.reserve(n_steps / simParams.steps_per_temperature_measurement + 1);
	}

	// Trainingdata buffers
	{
#ifdef GENERATETRAINDATA
		uint64_t n_loggingdata_host = 10 * n_steps;
		uint64_t n_traindata_host = n_steps * N_DATAGAN_VALUES * MAX_COMPOUND_PARTICLES * simulation.boxparams_host.n_compounds;
		auto datasize_str = std::to_string((float)(sizeof(Float3) * n_traindata_host + sizeof(float) * n_loggingdata_host) * 1e-9);
		m_logger->print("Reserving " + datasize_str + "GB host mem for logging and training data\n");

		simulation.loggingdata.resize(n_loggingdata_host);
		simulation.trainingdata.resize(n_traindata_host);
#endif
	}
}
