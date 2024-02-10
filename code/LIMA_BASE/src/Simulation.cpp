#include <assert.h>
#include <map>
#include <functional>
#include <type_traits> // For std::is_integral, std::is_floating_point, and static_assert

#include "Simulation.cuh"
#include "Forcefield.cuh"

using std::string;

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

SimParams::SimParams(const std::string& path) {
	auto dict = Filehandler::parseINIFile(path);
	overloadParamNumber<float>(dict, dt, "dt", [](const float& val) {return val * FEMTO_TO_LIMA; });
	overloadParamNumber(dict, n_steps, "n_steps");
	overloadParamNumber(dict, box_size, "boxlen");

	overwriteParamNonNumbers<bool>(dict, "em", em_variant, 
		[](const string& value) {return convertStringvalueToValue<bool>({ {"true", true }, {"false", false}}, "em", value); }
	);

	overwriteParamNonNumbers<BoundaryConditionSelect>(dict, "boundarycondition", bc_select,
		[](const string& value) {return convertStringvalueToValue<BoundaryConditionSelect>({ {"PBC", PBC }, {"NoBC", NoBC}}, "boundarycondition", value); }
	);
}


void SimParams::dumpToFile(const std::string& filename) {
	std::ofstream file(filename);
	if (!file.is_open()) {
		// Handle the error, e.g., by throwing an exception or logging an error message
		throw std::runtime_error("Unable to open file: " + filename);
	}

	file << "n_steps=" << n_steps << "\n";
	file << "dt=" << static_cast<int>(std::round(dt * LIMA_TO_FEMTO)) << " # [femto seconds]\n";
	file << "em=" << (em_variant ? "true" : "false") << " # Is an energy-minimization sim\n";
	file << "boundarycondition=" << (bc_select == PBC ? "PBC" : "No Boundary Condition") << " # (PBC, NoBC)\n";
	//file << "Supernatural Forces: " << (snf_select == HorizontalSqueeze ? "Horizontal Squeeze" : "None") << "\n";
	file << "boxlen=" << box_size << " # [nm]\n";

	file.close();
}



Box::~Box() {
	if (owns_members) { deleteMembers(); }
}

void Box::moveToDevice() {
	/*const size_t bytes_total = sizeof(Compound) * boxparams.n_compounds
		+ sizeof(CompoundState) * MAX_COMPOUNDS * 3u
		+ sizeof(NeighborList) * (MAX_SOLVENTS + MAX_COMPOUNDS);*/
	//printf("BOX: moving %.2f MB to device\n", (float)bytes_total * 1e-6);

	compounds = genericMoveToDevice(compounds, MAX_COMPOUNDS);
	bridge_bundle = genericMoveToDevice(bridge_bundle, 1);

	solvents = genericMoveToDevice(solvents, boxparams.n_solvents);

	coordarray_circular_queue = genericMoveToDevice(coordarray_circular_queue, MAX_COMPOUNDS * STEPS_PER_LOGTRANSFER);
	solventblockgrid_circularqueue = solventblockgrid_circularqueue->moveToDevice();


	bonded_particles_lut_manager = genericMoveToDevice(bonded_particles_lut_manager, 1);

	forcefield = genericMoveToDevice(forcefield, 1);

	cudaDeviceSynchronize();
	is_on_device = true;
}

void Box::deleteMembers() {
	if (is_on_device) {
		cudaFree(compounds);
		cudaFree(coordarray_circular_queue);
		cudaFree(solventblockgrid_circularqueue);

		cudaFree(forcefield);

		cudaFree(bridge_bundle);
		cudaFree(bonded_particles_lut_manager);

		if (boxparams.n_solvents > 0) {
			cudaFree(solvents);
		}		
	}
	else {
		delete[] compounds;
		delete[] coordarray_circular_queue;
		delete solventblockgrid_circularqueue;


		delete[] bridge_bundle;
		delete[] bonded_particles_lut_manager;

		// TMP, forcefield should maybe come with other members?
		if (forcefield) { delete forcefield; }

		if (boxparams.n_solvents > 0) {
			delete[] solvents;
		}
	}
	owns_members = false;	
	cudaDeviceSynchronize();
	auto cuda_status = cudaGetLastError();
	if (cuda_status != cudaSuccess) {
		std::cout << "\nCuda error code: " << cuda_status << " - " << cudaGetErrorString(cuda_status) << std::endl;
		exit(1);
	}
}

std::unique_ptr<Box> SimUtils::copyToHost(Box* box_dev) {
	auto box = std::make_unique<Box>();

	// First copy the box meta info, which we need to some of the dataptr transfers
	cudaMemcpy(box.get(), box_dev, sizeof(Box), cudaMemcpyDeviceToHost);

	genericCopyToHost(&box->compounds, MAX_COMPOUNDS);
	genericCopyToHost(&box->coordarray_circular_queue, MAX_COMPOUNDS * STEPS_PER_LOGTRANSFER);

	genericCopyToHost(&box->solvents, box->boxparams.n_solvents);
	box->solventblockgrid_circularqueue = box->solventblockgrid_circularqueue->copyToHost();

	
	genericCopyToHost(&box->forcefield, 1);

	genericCopyToHost(&box->bridge_bundle, 1);
	genericCopyToHost(&box->bonded_particles_lut_manager, 1);
	


	box->owns_members = true;
	box->is_on_device = false;
	//printf("Box copied to host\n");
	return box;
}


Simulation::Simulation(const SimParams& params, const std::string& molecule_path, EnvMode envmode) :
	simparams_host{ params }
{
	box_host = std::make_unique<Box>();
	//forcefield = std::make_unique<Forcefield>(envmode == Headless ? SILENT : V1, molecule_path);
}

void Simulation::copyBoxVariables() {
	boxparams_host = box_host->boxparams;
	//n_compounds = box_host->boxparams.n_compounds;
	//n_bridges = box_host->bridge_bundle->n_bridges;


	//n_solvents = box_host->boxparams.n_solvents;
	//blocks_per_solventkernel = (int)ceil((float)n_solvents / (float)THREADS_PER_SOLVENTBLOCK);

	compounds_host.resize(boxparams_host.n_compounds);
	for (int i = 0; i < boxparams_host.n_compounds; i++)
		compounds_host[i] = box_host->compounds[i];

	// Need this variable both on host and device
	//total_particles_upperbound = box_host->boxparams.n_compounds * MAX_COMPOUND_PARTICLES + SolventBlockGrid::blocks_total * MAX_SOLVENTS_IN_BLOCK;
	//box_host->boxparams.total_particles_upperbound = total_particles_upperbound;
}

//void InputSimParams::overloadParams(std::map<std::string, double>& dict) {
//	overloadParam(dict, &dt, "dt", FEMTO_TO_LIMA);	// convert [fs] to [ls]
//	overloadParam(dict, &n_steps, "n_steps");
//}



//SimParams::SimParams(const InputSimParams& ip) : 
//	n_steps(ip.n_steps), dt(ip.dt),em_variant(ip.em_variant), bc_select(ip.boundarycondition)
//{}

DatabuffersDevice::DatabuffersDevice(size_t total_particles_upperbound, int n_compounds) {
	// Permanent Outputs for energy & trajectory analysis
	{
		const size_t n_datapoints = total_particles_upperbound * STEPS_PER_LOGTRANSFER / LOG_EVERY_N_STEPS;
		const size_t bytesize_mb = (sizeof(float) * n_datapoints + sizeof(Float3) * n_datapoints) / 1'000'000;
		assert(n_datapoints && "Tried creating traj or potE buffers with 0 datapoints");
		assert(bytesize_mb < 6'000 && "Tried reserving >6GB data on device");

		
		cudaMallocManaged(&potE_buffer, sizeof(*potE_buffer) * n_datapoints);
		cudaMallocManaged(&traj_buffer, sizeof(*traj_buffer) * n_datapoints);
		cudaMallocManaged(&vel_buffer, sizeof(*vel_buffer) * n_datapoints);

		std::vector<float>potE_zero(n_datapoints, 0);
		std::vector<Float3>traj_zero(n_datapoints, Float3{});

		cudaMemcpy(potE_buffer, potE_zero.data(), sizeof(float) * n_datapoints, cudaMemcpyHostToDevice);
		cudaMemcpy(traj_buffer, traj_zero.data(), sizeof(Float3) * n_datapoints, cudaMemcpyHostToDevice);
	}
	

#ifdef USEDEBUGF3
	uint64_t bytes_for_debugf3 = sizeof(Float3) * DEBUGDATAF3_NVARS * simulation->total_particles_upperbound * simulation->n_steps;
	cudaMallocManaged(&simulation->box->debugdataf3, bytes_for_debugf3);
#endif

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

#ifdef GENERATETRAINDATA
	cudaFree(outdata);
	cudaFree(data_GAN);
#endif
}
