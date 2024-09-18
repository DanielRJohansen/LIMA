#include "Simulation.cuh"

#include <map>
#include <functional>
#include <type_traits> // For std::is_integral, std::is_floating_point, and static_assert
#include <filesystem>
#include "Filehandling.h"
#include <fstream>

using std::string;
namespace fs = std::filesystem;

template <typename T>
constexpr T convertStringvalueToValue(std::vector<std::pair<string, T>> pairs, const string& key_str, const string& val_str) {
	for (auto& pair : pairs) {
		if (pair.first == val_str) {
			return pair.second;
		}
	}

	throw std::runtime_error("Illegal key-value pair in sim_params.txt: " + key_str + " " + val_str);
}

// Helper function for overwriting templates
template <typename T>
constexpr void overwriteParamNonNumbers(std::map<std::string, std::string>& dict, const std::string& key, T& val, std::function<T(const string&)> transform) {
	if (dict.count(key)) {
		val = transform(dict[key]);
	}
}

template <typename T>
constexpr void overloadParamNumber(std::map<std::string, std::string>& dict, T& param,
	std::string key, std::function<T(const T&)> transform = [](const T& v) {return v; })
{
	static_assert(std::is_integral_v<T> || std::is_floating_point_v<T>, "T must be an integral or floating-point type.");

	if (dict.count(key)) {
		if constexpr (std::is_integral_v<T>) {
			// For integral types (e.g., int)
			param = static_cast<T>(std::stoi(dict[key]));
		}
		else if constexpr (std::is_floating_point_v<T>) {
			// For floating-point types (e.g., float, double)
			param = static_cast<T>(transform(std::stof(dict[key])));
		}
	}
}

SimParams::SimParams(const fs::path& path) {
	auto dict = Filehandler::parseINIFile(path.string());
	overloadParamNumber<float>(dict, dt, "dt", [](const float& val) {return val * FEMTO_TO_LIMA; });
	overloadParamNumber(dict, n_steps, "n_steps");
	overloadParamNumber(dict, data_logging_interval, "data_logging_interval");

	overwriteParamNonNumbers<bool>(dict, "em", em_variant, 
		[](const string& value) {return convertStringvalueToValue<bool>({ {"true", true }, {"false", false}}, "em", value); }
	);
	overwriteParamNonNumbers<bool>(dict, "enable_electrostatics", enable_electrostatics,
				[](const string& value) {return convertStringvalueToValue<bool>({ {"true", true }, {"false", false}}, "enable_electrostatics", value); }
	);


	overwriteParamNonNumbers<BoundaryConditionSelect>(dict, "boundarycondition", bc_select,
		[](const string& value) {return convertStringvalueToValue<BoundaryConditionSelect>({ {"PBC", PBC }, {"NoBC", NoBC}}, "boundarycondition", value); }
	);
}


void SimParams::dumpToFile(const fs::path& filename) {
	std::ofstream file(filename);
	if (!file.is_open()) {
		// Handle the error, e.g., by throwing an exception or logging an error message
		throw std::runtime_error("Unable to open file: " + filename.string());
	}

	file << "n_steps=" << n_steps << "\n";
	file << "dt=" << static_cast<int>(std::round(dt * LIMA_TO_FEMTO)) << " # [femto seconds]\n";
	file << "em=" << (em_variant ? "true" : "false") << " # Is an energy-minimization sim\n";
	file << "boundarycondition=" << (bc_select == PBC ? "PBC" : "No Boundary Condition") << " # (PBC, NoBC)\n";
	//file << "Supernatural Forces: " << (snf_select == HorizontalSqueeze ? "Horizontal Squeeze" : "None") << "\n";
	file << "data_logging_interval=" << data_logging_interval << " # [steps]\n";
	file << "enable_electrostatics=" << (enable_electrostatics ? "true" : "false") << "\n";
	file.close();
}










Box::Box(int boxSizeNM) {
	boxparams.boxSize = boxSizeNM;
	compoundcoordsCircularQueue = CompoundcoordsCircularQueue::CreateQueue();
	solventblockgrid_circularqueue = SolventBlocksCircularQueue::createQueue(boxSizeNM);
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

		temperature_buffer.reserve(n_steps / STEPS_PER_THERMOSTAT + 1);
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

//void InputSimParams::overloadParams(std::map<std::string, double>& dict) {
//	overloadParam(dict, &dt, "dt", FEMTO_TO_LIMA);	// convert [fs] to [ls]
//	overloadParam(dict, &n_steps, "n_steps");
//}



//SimParams::SimParams(const InputSimParams& ip) : 
//	n_steps(ip.n_steps), dt(ip.dt),em_variant(ip.em_variant), bc_select(ip.boundarycondition)
//{}

DatabuffersDevice::DatabuffersDevice(int total_particles_upperbound, int n_compounds, int loggingInterval) :
	total_particles_upperbound{ total_particles_upperbound }
{
	// Permanent Outputs for energy & trajectory analysis
	{
		const size_t n_datapoints = total_particles_upperbound * nStepsInBuffer;
		const size_t bytesize_mb = (sizeof(float) * n_datapoints + sizeof(Float3) * n_datapoints) / 1'000'000;
		assert(n_datapoints && "Tried creating traj or potE buffers with 0 datapoints");
		assert(bytesize_mb < 6'000 && "Tried reserving >6GB data on device");

		
		cudaMallocManaged(&potE_buffer, sizeof(*potE_buffer) * n_datapoints);
		cudaMallocManaged(&traj_buffer, sizeof(*traj_buffer) * n_datapoints);
		cudaMallocManaged(&vel_buffer, sizeof(*vel_buffer) * n_datapoints);
		cudaMallocManaged(&forceBuffer, sizeof(*forceBuffer) * n_datapoints);

		//cudaMemset // TODO: switch to memset
		std::vector<float>potE_zero(n_datapoints, 0);
		std::vector<Float3>traj_zero(n_datapoints, Float3{});

		cudaMemcpy(potE_buffer, potE_zero.data(), sizeof(float) * n_datapoints, cudaMemcpyHostToDevice);
		cudaMemcpy(traj_buffer, traj_zero.data(), sizeof(Float3) * n_datapoints, cudaMemcpyHostToDevice);
	}
	

	// TRAINING DATA and TEMPRARY OUTPUTS
	{
#ifdef GENERATETRAINDATA
		size_t n_outdata = 10 * STEPS_PER_LOGTRANSFER;
		size_t n_traindata = std::max(static_cast<size_t>(N_DATAGAN_VALUES) * MAX_COMPOUND_PARTICLES * n_compounds * STEPS_PER_TRAINDATATRANSFER, size_t{ 1 });	// This is a hack; can't allocate 0 bytes.
		assert(n_outdata && "Tried to create outdata buffer with 0 datapoints");
		assert(n_traindata && "Tried to create traindata buffer with 0 datapoints");
		
		size_t bytesize_mb = (sizeof(float) * n_outdata + sizeof(Float3) * n_traindata) / 1'000'000;
		assert(bytesize_mb < 6'000 && "Tried reserving >6GB data on device");

		cudaMallocManaged(&outdata, sizeof(float) * n_outdata);	// 10 data streams for 10k steps. 1 step at a time.
		cudaMallocManaged(&data_GAN, sizeof(Float3) * n_traindata);
#endif
	}
	
}

void DatabuffersDevice::freeMembers() {
	cudaFree(potE_buffer);
	cudaFree(traj_buffer);
	cudaFree(vel_buffer);
	cudaFree(forceBuffer);

#ifdef GENERATETRAINDATA
	cudaFree(outdata);
	cudaFree(data_GAN);
#endif
}
