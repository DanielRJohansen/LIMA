#include "Simulation.cuh"

#include <functional>
#include <type_traits> // For std::is_integral, std::is_floating_point, and static_assert
#include <filesystem>
#include "Filehandling.h"
#include "MDFiles.h"
#include <fstream>

#include <concepts>

using std::string;
namespace fs = std::filesystem;

using Dictionary = std::unordered_map<string, string>;

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
constexpr void overwriteParamNonNumbers(std::unordered_map<std::string, std::string>& dict, const std::string& key, T& val, std::function<T(const string&)> transform) {
	if (dict.count(key)) {
		val = transform(dict[key]);
	}
}

void Readb(const Dictionary& dict, bool& value, const std::string& key_name) {
	if (!dict.contains(key_name))
		return;

	std::string value_str = dict.at(key_name);
	if (value_str != "true" && value_str != "false") {
		throw std::runtime_error("Illegal key-value pair in sim_params.txt: " + key_name + " " + value_str);
	}
	
	value = (value_str == "true");
}
template <std::integral T>
constexpr void Readi(std::unordered_map<std::string, std::string>& dict, T& param,
	const std::string& key, std::function<T(const T&)> transform = [](const T& v) { return v; })
{
	if (dict.count(key)) {
		param = static_cast<T>(transform(std::stoll(dict[key])));
	}
}
template <std::floating_point T, typename Transform = std::function<T(const T&)>>
constexpr void Readf(std::unordered_map<std::string, std::string>& dict, T& param,
	const std::string& key, Transform transform = [](const T& v) { return v; })
{
	if (dict.count(key)) {
		param = static_cast<T>(transform(std::stod(dict[key])));
	}
}

SimParams::SimParams(const fs::path& path) {
	const bool forceKeysAndValuesLowercase = true;
	auto dict = FileUtils::parseINIFile(path.string(), forceKeysAndValuesLowercase);

	Readi(dict, n_steps, "n_steps");
	Readf(dict, dt, "dt", [](auto val) {return val * FEMTO_TO_LIMA; });
	Readb(dict, em_variant, "em");
	Readf(dict, em_force_tolerance, "em_force_tolerance");

	overwriteParamNonNumbers<BoundaryConditionSelect>(dict, "boundarycondition", bc_select,
		[](const string& value) {return convertStringvalueToValue<BoundaryConditionSelect>({ {"pbc", PBC }, {"nobc", NoBC} }, "boundarycondition", value); }
	);
	Readb(dict, enable_electrostatics, "enable_electrostatics");
	Readf(dict, cutoff_nm, "cutoff_nm");
	// Skip SNF for now	

	Readi(dict, data_logging_interval, "data_logging_interval");	
	Readb(dict, save_trajectory, "save_trajectory");
	Readb(dict, save_energy, "save_energy");
	// Skip colormethod for now

	Readi(dict, steps_per_temperature_measurement, "steps_per_temperature_measurement");
	Readb(dict, apply_thermostat, "apply_thermostat");
}


void SimParams::dumpToFile(const fs::path& filename) {
	std::ofstream file(filename);
	if (!file.is_open()) {
		// Handle the error, e.g., by throwing an exception or logging an error message
		throw std::runtime_error("Unable to open file: " + filename.string());
	}

	file << "// Main params" << "\n";
	file << "n_steps=" << n_steps << "\n";
	file << "dt=" << static_cast<int>(std::round(dt * LIMA_TO_FEMTO)) << " # [fs]\n";
	file << "em=" << (em_variant ? "true" : "false") << " # Is an energy-minimization sim\n";
	file << "em_force_tolerance=" << em_force_tolerance << " # [kJ/mol/nm]\n"; // TODO add units

	file << "// Physics params" << "\n";
	file << "boundarycondition=" << (bc_select == PBC ? "PBC" : "No Boundary Condition") << " # (PBC, NoBC)\n";
	file << "enable_electrostatics=" << (enable_electrostatics ? "true" : "false") << "\n";
	file << "cutoff_nm=" << cutoff_nm << " # [nm]\n";
	// Skip SNF for now

	file << "// Output params" << "\n";
	file << "data_logging_interval=" << data_logging_interval << " # [steps]\n";
	file << "save_trajectory=" << (save_trajectory ? "true" : "false") << "\n";
	file << "save_energy=" << (save_energy ? "true" : "false") << " # (Save kinetic and potential energy to file)" << "\n";
	// Skip colormethod for now
	
	file << "// Temperature params" << "\n";
	file << "steps_per_temperature_measurement=" << steps_per_temperature_measurement << " # [steps]\n";
	file << "apply_thermostat=" << (apply_thermostat ? "true" : "false") << " # will speed up or slow down particles to achieve the desired temperature\n";
	file << "# desired_temperature - not currently available forced to be 300 [k]\n";

	file.close();
}










Box::Box(int boxSizeNM) {
	boxparams.boxSize = boxSizeNM;
	compoundcoordsCircularQueue = CompoundcoordsCircularQueue_Host::CreateQueue();
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

std::unique_ptr<MDFiles::TrrFile> Simulation::ToTracjectoryFile() const{	
	auto trrFile = std::make_unique<MDFiles::TrrFile>(Float3{ box_host->boxparams.boxSize });

	assert(simparams_host.data_logging_interval == traj_buffer->GetLoggingInterval());
	const int nLoggedSteps = step / simparams_host.data_logging_interval;
	trrFile->positions.reserve(nLoggedSteps);

	for (int64_t step = 0; step < step; step += simparams_host.data_logging_interval) {

		std::vector<Float3> row(box_host->boxparams.total_particles);
	    int index = 0; 

	    // Todo: this can be optimized with some insert magic, but i do not have the brain capacity. Ask gpt?
	    for (int compound_id = 0; compound_id < box_host->boxparams.n_compounds; compound_id++) {
	        for (int pid = 0; pid < box_host->compounds[compound_id].n_particles; pid++) {
				row[index++] = traj_buffer->GetMostRecentCompoundparticleDatapoint(compound_id, pid, step);
	        }
	    }

	    for (int solvent_index = 0; solvent_index < box_host->boxparams.n_solvents; solvent_index++) {
	        row[index++] = traj_buffer->GetMostRecentSolventparticleDatapointAtIndex(solvent_index, step);
	    }

		trrFile->positions.emplace_back(std::move(row));
	}

	return std::move(trrFile);
}
