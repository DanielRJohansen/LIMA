#include "Simulation.cuh"

#include <functional>
#include <type_traits> // For std::is_integral, std::is_floating_point, and static_assert
#include <filesystem>
#include "Filehandling.h"
#include "MDFiles.h"
#include <fstream>

#include <concepts>









Box::Box(Float3 boxSizeNM) {
	if ((NodeIndex(boxSizeNM.ToInt3()).toFloat3() - boxSizeNM).len() > 0.0000001
		|| boxSizeNM.x != boxSizeNM.y || boxSizeNM.y != boxSizeNM.z) {
		throw std::invalid_argument("Boxsize was not an integer, or was not cubic");
	}

	boxparams.boxSize = static_cast<int>(boxSizeNM.x);
	solventblockgrid_circularqueue = SolventBlocksCircularQueue::createQueue(static_cast<int>(boxSizeNM.x));
}





Simulation::Simulation(const SimParams& params) :
	simparams_host{ params }
{
	box_host = std::make_unique<Box>();
}

Simulation::Simulation(const SimParams& params, std::unique_ptr<Box> box) :
	simparams_host{ params }
{
	box_host = std::move(box);
}

void Simulation::PrepareDataBuffers() {
	// Allocate buffers. We need to allocate for atleast 1 step, otherwise the bootstrapping mechanism will fail.
	const auto n_steps = std::max(simparams_host.n_steps, uint64_t{ 1 });
	// Standard Data Buffers 
	{
		// Permanent Outputs for energy & trajectory analysis
		const int particlesUpperbound = box_host->boxparams.total_particles_upperbound;
		const size_t n_datapoints = particlesUpperbound * n_steps / simparams_host.data_logging_interval;
		const auto datasize_str = std::to_string((float)((2. * sizeof(float) * n_datapoints + sizeof(Float3) * n_datapoints) * 1e-6));
		//m_logger->print("Malloc " + datasize_str + " MB on host for data buffers\n");


		potE_buffer = std::make_unique<ParticleDataBuffer<float>>(particlesUpperbound, box_host->boxparams.n_compounds, n_steps, simparams_host.data_logging_interval);
		vel_buffer = std::make_unique<ParticleDataBuffer<float>>(particlesUpperbound, box_host->boxparams.n_compounds, n_steps, simparams_host.data_logging_interval);
		forceBuffer = std::make_unique<ParticleDataBuffer<Float3>>(particlesUpperbound, box_host->boxparams.n_compounds, n_steps, simparams_host.data_logging_interval);
		traj_buffer = std::make_unique<ParticleDataBuffer<Float3>>(particlesUpperbound, box_host->boxparams.n_compounds, n_steps, simparams_host.data_logging_interval);
		
		temperature_buffer.reserve(n_steps / simparams_host.steps_per_temperature_measurement + 1);
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
