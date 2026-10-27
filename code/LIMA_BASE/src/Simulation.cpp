#include "Simulation.cuh"
#include "Format.h"

#include <functional>
#include <type_traits> // For std::is_integral, std::is_floating_point, and static_assert
#include <filesystem>
#include "Filehandling.h"
#include "MDFiles.h"
#include <fstream>

#include <concepts>









Box::Box(Float3 boxSizeNM) {
	if (!(boxSizeNM.x > 0.f && boxSizeNM.y > 0.f && boxSizeNM.z > 0.f))
		throw std::invalid_argument("Box size must be positive");
	boxparams.SetBoxLength(boxSizeNM);
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

struct TrajectoryStream {
	TrajectoryStream(const std::filesystem::path& path, std::vector<int> atomOfSlot, int nAtoms, int interval, float dt)
		: writer{ path, nAtoms }, atomOfSlot(std::move(atomOfSlot)), frame(nAtoms), interval(interval),
		// As Trajectory, so a frame's time is the same float as when the whole trajectory was kept and written at the end
		frameSpacingPs(static_cast<float>(dt * interval * NANO_TO_PICO)), dtPs(static_cast<float>(dt * NANO_TO_PICO)) {}

	MDFiles::TrajectoryWriter writer;
	std::vector<int> atomOfSlot;	// By logging slot: the atom's index in the input coordinates, or -1 for a padding lane
	std::vector<Float3> frame;		// [nm] In input atom order
	int interval = 0;				// [steps]
	float frameSpacingPs = 0.f;
	float dtPs = 0.f;
	int64_t lastStep = -1;
};

Simulation::~Simulation() = default;

void Simulation::StartTrajectoryStream(const std::filesystem::path& path) {
	const int loggingInterval = simParams.nstenergy;
	const int interval = simParams.nstxout_compressed;
	if (loggingInterval <= 0)
		throw std::invalid_argument("A trajectory needs nstenergy > 0, as its frames come from the logged data");
	if (interval <= 0 || interval % loggingInterval != 0)
		throw std::invalid_argument(Lima::Format("nstxout_compressed ({}) must be a positive multiple of nstenergy ({})",
			interval, loggingInterval));

	std::vector<int> atomOfSlot;
	atomOfSlot.reserve(box->persistentClustersMetadata.size() * PersistentCluster::maxParticles);
	int nAtoms = 0;
	for (const auto& metadata : box->persistentClustersMetadata)
		for (const int atom : metadata.particleIdsGlobal) {
			atomOfSlot.push_back(atom);
			nAtoms = std::max(nAtoms, atom + 1);
		}
	trajectoryStream = std::make_unique<TrajectoryStream>(path, std::move(atomOfSlot), nAtoms, interval, simParams.dt.value);
}

bool Simulation::TrajectoryFrameAt(int64_t step) const {
	return trajectoryStream && step % trajectoryStream->interval == 0;
}

void Simulation::WriteTrajectoryFrame(int64_t step, Float3 boxSize, const Float3* slots) {
	TrajectoryStream& stream = *trajectoryStream;
	if (step <= stream.lastStep)
		return;	// The final frame, when it is also a regular one
	for (size_t slot = 0; slot < stream.atomOfSlot.size(); ++slot)
		if (const int atom = stream.atomOfSlot[slot]; atom >= 0)
			stream.frame[atom] = slots[slot];
	const float timePs = step % stream.interval == 0 ? stream.frameSpacingPs * static_cast<int>(step / stream.interval) : step * stream.dtPs;
	stream.writer.WriteFrame(step, timePs, boxSize, stream.frame.data());
	stream.lastStep = step;
}

void Simulation::FinishTrajectoryStream() {
	trajectoryStream.reset();
}

void Simulation::PrepareDataBuffers() {
	// Allocate buffers. We need to allocate for atleast 1 step, otherwise the bootstrapping mechanism will fail.
	const auto n_steps = std::max(simParams.nsteps.value, uint64_t{ 1 });
	// Standard Data Buffers
	{
		// Permanent Outputs for energy & trajectory analysis
		const int nPclusters = box->persistentClusters.size();
		const int particlesUpperbound = nPclusters * PersistentCluster::maxParticles;
		const int interval = simParams.nstenergy;
		// As ParticleDataBuffer sizes itself: one entry per logged step, and at least one
		const uint64_t nEntries = interval > 0 ? std::max(n_steps / interval + (n_steps % interval != 0), uint64_t{ 1 }) : 1;
		const uint64_t bytesPerEntry = static_cast<uint64_t>(particlesUpperbound) * ((logged.trajectory ? sizeof(Float3) : 0)
			+ (logged.energy ? 2 * sizeof(float) : 0) + (logged.forces ? sizeof(Float3) : 0));
		const uint64_t bytes = bytesPerEntry * nEntries;
		// Fail before allocating, rather than swapping for hours or being killed partway into the run
		if (const uint64_t available = FileUtils::AvailablePhysicalMemory(); available > 0 && bytes > available) {
			std::string selected;
			for (const auto& [on, name] : { std::pair{ logged.trajectory, "trajectory" }, { logged.energy, "energies" }, { logged.forces, "forces" } })
				if (on) selected += (selected.empty() ? "" : ", ") + std::string{ name };
			throw std::runtime_error(Lima::Format("Logging this run ({}) needs {:.1f} GB of RAM, {} frames of {:.1f} MB, but {:.1f} GB is available. "
				"Raise nstenergy (now {}), run fewer steps, or leave out the outputs that need it (eg. --trajectory, save_energy)",
				selected, bytes / 1e9, nEntries, bytesPerEntry / 1e6, available / 1e9, interval));
		}

		if (logged.energy) {
			potE_buffer = std::make_unique<ParticleDataBuffer<float>>(particlesUpperbound, n_steps, interval, nPclusters);
			vel_buffer = std::make_unique<ParticleDataBuffer<float>>(particlesUpperbound, n_steps, interval, nPclusters);
		}
		if (logged.forces)
			forceBuffer = std::make_unique<ParticleDataBuffer<Float3>>(particlesUpperbound, n_steps, interval, nPclusters);
		if (logged.trajectory)
			traj_buffer = std::make_unique<ParticleDataBuffer<Float3>>(particlesUpperbound, n_steps, interval, nPclusters);

		temperature_buffer.reserve(n_steps / simParams.nsttcouple + 1);
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
