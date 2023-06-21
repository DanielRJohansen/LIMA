#include "LIMA_MD/include/Analyzer.cuh"
#include "LIMA_BASE/include/Printer.h"
#include "LIMA_BASE/include/Utilities.h"

#include <algorithm>

using namespace LIMA_Print;

template<typename T>
void __device__ distributedSummation(T* arrayptr, int array_len) {				// Places the result at pos 0 of input_array
	T temp;			// This is a lazy soluation, but maybe it is also fast? Definitely simple..
	for (int i = 1; i < array_len; i *= 2) {	// Distributed averaging							// Make a generic and SAFER function for this, PLEASE OK??
		if ((threadIdx.x + i) < array_len) {
			temp = arrayptr[threadIdx.x] + arrayptr[threadIdx.x + i];// *0.5f;	// easier to just divide by sum of solvents at host
		}
		__syncthreads();
		arrayptr[threadIdx.x] = temp;
		__syncthreads();
	}
}

void __global__ monitorCompoundEnergyKernel(Box* box, const SimParams* simparams, Float3* traj_buffer, float* potE_buffer, Float3* data_out) {		// everything here breaks if not all compounds are identical in particle count and particle mass!!!!!!!
	__shared__ Float3 energy[MAX_COMPOUND_PARTICLES];
	__shared__ Compound compound;


	const uint64_t step = blockIdx.x + uint64_t{ 1 };
	const int compound_index = blockIdx.y;
	energy[threadIdx.x] = Float3(0.f);


	if (threadIdx.x == 0) {
		data_out[compound_index + (step - 1) * box->boxparams.n_compounds] = Float3{};
		compound = box->compounds[compound_index];
	}
	__syncthreads();

	if (threadIdx.x >= compound.n_particles) {
		return;
	}
	__syncthreads();

	const uint8_t& atom_type = compound.atom_types[threadIdx.x];
	const float mass = box->forcefield->particle_parameters[atom_type].mass;

	const uint32_t compound_offset = compound_index * MAX_COMPOUND_PARTICLES;
	const int step_offset = step * box->boxparams.total_particles_upperbound;
	const float potE = potE_buffer[threadIdx.x + compound_offset + step_offset];

	const Float3 pos_tsub1 = traj_buffer[threadIdx.x + compound_offset + (step - uint64_t{ 1 }) * box->boxparams.total_particles_upperbound];
	const Float3 pos_tadd1 = traj_buffer[threadIdx.x + compound_offset + (step + uint64_t{ 1 }) * box->boxparams.total_particles_upperbound];
	const float n_steps = 2.f;

	const float kinE = EngineUtils::calcKineticEnergy(&pos_tadd1, &pos_tsub1, mass, simparams->constparams.dt * n_steps / NANO_TO_LIMA);	// convert dt from [ls] to [ns]
	const float totalE = potE + kinE;

	energy[threadIdx.x] = Float3(potE, kinE, totalE);
	__syncthreads();

	distributedSummation(energy, MAX_COMPOUND_PARTICLES);
	__syncthreads();

	if (threadIdx.x == 0) {
		data_out[compound_index + (step - 1) * box->boxparams.n_compounds] = energy[0];
	}
}





void __global__ monitorSolventEnergyKernel(Box* box, const SimParams* simparams, Float3* traj_buffer, float* potE_buffer, Float3* data_out) {
	__shared__ Float3 energy[THREADS_PER_SOLVENTBLOCK];



	int solvent_index = threadIdx.x + blockIdx.y * THREADS_PER_SOLVENTBLOCK;
	int step = blockIdx.x + 1;
	int compounds_offset = box->boxparams.n_compounds * MAX_COMPOUND_PARTICLES;


	energy[threadIdx.x] = Float3(0.f);
	if (threadIdx.x == 0) {
		data_out[(step - 1) * gridDim.y + blockIdx.y] = energy[0];
	}
	if (solvent_index >= box->boxparams.n_solvents) { return; }


	
	Float3 pos_tsub1 = traj_buffer[compounds_offset + solvent_index + (step - 1) * box->boxparams.total_particles_upperbound];
	Float3 pos_tadd1 = traj_buffer[compounds_offset + solvent_index + (step + 1) * box->boxparams.total_particles_upperbound];

	float mass = SOLVENT_MASS;
	float kinE = EngineUtils::calcKineticEnergy(&pos_tadd1, &pos_tsub1, mass, simparams->constparams.dt*2.f / NANO_TO_LIMA);	// convert [ls] to [ns]

	float potE = potE_buffer[compounds_offset + solvent_index + step * box->boxparams.total_particles_upperbound];

	if (potE != 0.f) {
		//printf("step %04d solvate %04d pot %f, compound_offset %d, step_offset  %d\n", step, solvent_index, potE, compounds_offset, step*box->total_particles_upperbound);
	}

	float totalE = potE + kinE;

	energy[threadIdx.x] = Float3(potE, kinE, totalE);
	__syncthreads();
	distributedSummation(energy, THREADS_PER_SOLVENTBLOCK);
	if (threadIdx.x == 0) {
		data_out[(step - 1) * gridDim.y + blockIdx.y] = energy[0];
	}
}



Analyzer::AnalyzedPackage Analyzer::analyzeEnergy(Simulation* simulation) {	// Calculates the avg J/mol // calculate energies separately for compounds and solvents. weigh averages based on amount of each
	LIMA_UTILS::genericErrorCheck("Cuda error before analyzeEnergy\n");

	const auto n_steps = simulation->getStep();
	if (simulation->getStep() < 3) { return Analyzer::AnalyzedPackage(); }

	std::vector<Float3> average_energy;
	average_energy.resize(n_steps - 2);	// Ignore first and last step

	// We need to split up the analyser into steps, as we cannot store all positions traj on device at once.
	uint64_t max_steps_per_kernel = 1000;
	uint64_t particles_per_step = simulation->boxparams_host.total_particles_upperbound;
	uint64_t max_values_per_kernel = (max_steps_per_kernel + 2) * particles_per_step;							// Pad steps with 2 for vel calculation
	const std::string bytesize = std::to_string((sizeof(Float3) + sizeof(double)) * (max_values_per_kernel) * 1e-6);
	m_logger->print("Analyzer malloc " + bytesize + " MB on device\n");
	cudaMalloc(&traj_buffer_device, sizeof(Float3) * max_values_per_kernel);
	cudaMalloc(&potE_buffer_device, sizeof(float) * max_values_per_kernel);

	for (int i = 0; i < ceil((double)n_steps / (double)max_steps_per_kernel); i++) {
		uint64_t step_offset = i * max_steps_per_kernel;												// offset one since we can't analyse step 1
		uint64_t steps_in_kernel = std::min(max_steps_per_kernel, n_steps - step_offset);

		// Create a array of len 1002, where index 0 and 1001 are padded values
		moveAndPadData(simulation, steps_in_kernel, step_offset);

		std::vector<Float3> average_solvent_energy = analyzeSolvateEnergy(simulation, steps_in_kernel);
		std::vector<Float3> average_compound_energy = analyzeCompoundEnergy(simulation, steps_in_kernel);

		for (uint64_t ii = 0; ii < steps_in_kernel; ii++) {
			int64_t step = step_offset + ii - 1;	// -1 because index 0 is unused
			if (step == -1 || step >= n_steps-2u) { continue; }	// Dont save first step, as the kinE is slightly wrong
			average_energy[step] = (average_solvent_energy[ii] + average_compound_energy[ii]);
		}
	}

	cudaFree(traj_buffer_device);
	cudaFree(potE_buffer_device);

	m_logger->finishSection("Finished analyzing energies");
	return AnalyzedPackage(average_energy, simulation->temperature_buffer);
}

void Analyzer::moveAndPadData(Simulation* simulation, uint64_t steps_in_kernel, uint64_t step_offset) {
	const uint64_t particles_per_step = simulation->boxparams_host.total_particles_upperbound;

	// First move the middle bulk to device
	//cudaMemcpy(&traj_buffer_device[1 * particles_per_step], &simulation->traj_buffer[step_offset * particles_per_step], sizeof(Float3) * steps_in_kernel * particles_per_step, cudaMemcpyHostToDevice);
	cudaMemcpy(&traj_buffer_device[1 * particles_per_step], &simulation->traj_buffer->data()[step_offset * particles_per_step], sizeof(Float3) * steps_in_kernel * particles_per_step, cudaMemcpyHostToDevice);
	cudaMemcpy(&potE_buffer_device[1 * particles_per_step], &simulation->potE_buffer[step_offset * particles_per_step], sizeof(float) * steps_in_kernel * particles_per_step, cudaMemcpyHostToDevice);
	LIMA_UTILS::genericErrorCheck("Cuda error during analyzer transfer2\n");

	// Then pad the front. If step 0, we pad with zero. If step n we pad with n-1
	uint64_t paddingSrcIndex = step_offset == 0 ? 0 : step_offset - 1;
	//cudaMemcpy(&traj_buffer_device[0], &simulation->traj_buffer[paddingSrcIndex * particles_per_step], sizeof(Float3) * particles_per_step, cudaMemcpyHostToDevice);
	cudaMemcpy(&traj_buffer_device[0], &simulation->traj_buffer->data()[paddingSrcIndex * particles_per_step], sizeof(Float3) * particles_per_step, cudaMemcpyHostToDevice);
	cudaMemcpy(&potE_buffer_device[0], &simulation->potE_buffer[paddingSrcIndex * particles_per_step], sizeof(float) * particles_per_step, cudaMemcpyHostToDevice);
	LIMA_UTILS::genericErrorCheck("Cuda error during analyzer transfer1\n");

	// Then pad the end. if step
	const uint64_t end_index = step_offset + steps_in_kernel;
	paddingSrcIndex = end_index == simulation->getStep() ? end_index - 1 : end_index;
	//cudaMemcpy(&traj_buffer_device[(steps_in_kernel + 1) * particles_per_step], &simulation->traj_buffer[paddingSrcIndex * particles_per_step], sizeof(Float3) * particles_per_step, cudaMemcpyHostToDevice);
	cudaMemcpy(&traj_buffer_device[(steps_in_kernel + 1) * particles_per_step], &simulation->traj_buffer->data()[paddingSrcIndex * particles_per_step], sizeof(Float3) * particles_per_step, cudaMemcpyHostToDevice);
	cudaMemcpy(&potE_buffer_device[(steps_in_kernel + 1) * particles_per_step], &simulation->potE_buffer[paddingSrcIndex * particles_per_step], sizeof(float) * particles_per_step, cudaMemcpyHostToDevice);

	LIMA_UTILS::genericErrorCheck("Cuda error during analyzer transfer\n");
}

std::vector<Float3> Analyzer::analyzeSolvateEnergy(Simulation* simulation, uint64_t n_steps) {
	// Start by creating array of energies of value 0
	std::vector<Float3> average_solvent_energy(n_steps);

	int blocks_per_solventkernel = (int)ceil((float)simulation->boxparams_host.n_solvents / (float)THREADS_PER_SOLVENTBLOCK);

	// If any solvents are present, fill above array
	if (simulation->boxparams_host.n_solvents > 0) {

		std::vector<Float3> average_solvent_energy_blocked(n_steps * blocks_per_solventkernel);
		Float3* data_out;
		cudaMalloc(&data_out, sizeof(Float3) * blocks_per_solventkernel * n_steps);

		dim3 block_dim(n_steps, blocks_per_solventkernel, 1);
		monitorSolventEnergyKernel << < block_dim, THREADS_PER_SOLVENTBLOCK >> > (simulation->sim_dev->box, simulation->sim_dev->params, traj_buffer_device, potE_buffer_device, data_out);
		LIMA_UTILS::genericErrorCheck("Cuda error during analyzeSolvateEnergy\n");

		cudaMemcpy(average_solvent_energy_blocked.data(), data_out, sizeof(Float3) * blocks_per_solventkernel * n_steps, cudaMemcpyDeviceToHost);
		cudaDeviceSynchronize();
		cudaFree(data_out);

		for (uint64_t step = 0; step < n_steps; step++) {
			average_solvent_energy[step] = Float3(0.f);
			for (int block = 0; block < blocks_per_solventkernel; block++) {
				average_solvent_energy[step] += average_solvent_energy_blocked[block + step * blocks_per_solventkernel];
			}
			//average_solvent_energy[step] *= (1.f / simulation->boxparams_host.n_solvents);
		}

	}

	return average_solvent_energy;
}


std::vector<Float3> Analyzer::analyzeCompoundEnergy(Simulation* simulation, uint64_t steps_in_kernel) {
	const uint64_t n_datapoints = simulation->boxparams_host.n_compounds * steps_in_kernel;

	std::vector<Float3> total_compound_energy(steps_in_kernel);

	if (simulation->extraparams.total_compound_particles > 0) {
		std::vector<Float3> host_data(n_datapoints);

		Float3* data_out;
		cudaMalloc(&data_out, sizeof(Float3) * n_datapoints);

		dim3 block_dim(static_cast<uint32_t>(steps_in_kernel), simulation->sim_dev->box->boxparams.n_compounds, 1);
		monitorCompoundEnergyKernel << < block_dim, MAX_COMPOUND_PARTICLES >> > (simulation->sim_dev->box, simulation->sim_dev->params, traj_buffer_device, potE_buffer_device, data_out);
		cudaDeviceSynchronize();
		LIMA_UTILS::genericErrorCheck("Cuda error during analyzeCompoundEnergy\n");

		cudaMemcpy(host_data.data(), data_out, sizeof(Float3) * n_datapoints, cudaMemcpyDeviceToHost);
		cudaFree(data_out);


		for (uint64_t step = 0; step < steps_in_kernel; step++) {
			for (uint64_t i = 0; i < simulation->sim_dev->box->boxparams.n_compounds; i++) {
				total_compound_energy[step] += host_data[i + step * simulation->sim_dev->box->boxparams.n_compounds];
			}
			//total_compound_energy[step] *= (1.f / (simulation->total_compound_particles));
		}

	}

	return total_compound_energy;
}

float getMin(const std::vector<float>& vec) {
	return *std::min_element(vec.begin(), vec.end());
}

float getMax(const std::vector<float>& vec) {
	return *std::max_element(vec.begin(), vec.end());
}

float getMean(const std::vector<float>& vec)
{
	double sum = 0.;
	for (auto elem : vec) { sum += static_cast<double>(elem); }	
	return static_cast<float>(sum / static_cast<double>(vec.size()));
}

float getStdDev(const std::vector<float>& vec) {
	if (vec.size() == 0) { return 0.f; }

	const double mean = getMean(vec);

	double variance = 0;
	for (auto elem : vec) { variance += (elem - mean) * (elem - mean); }

	const double deviation = variance / static_cast<double>(vec.size());
	return static_cast<float>(std::abs(std::sqrt(deviation)));
}

float Analyzer::getVarianceCoefficient(const std::vector<float>& vec) {
	if (vec.empty()) { return 0.f; } 
	const float stddev = getStdDev(vec);
	const float mean = getMean(vec);

	if (stddev == 0.f && mean == 0.f) { return 0.f; }
	return  stddev / mean;
}

void printRow(string title, std::vector<float>& vec) {
	if (vec.empty()) { return; }
	LIMA_Printer::printTableRow(
		title, { 
			getMin(vec), 
			getMax(vec), 
			getStdDev(vec),
			(vec.back() - vec.front()) / vec.front() });
}

void Analyzer::printEnergy(AnalyzedPackage* package) {
	LIMA_Printer::printTableRow({ "", "min", "max", "Std. deviation", "Change 0->n"});
	printRow("potE", package->pot_energy);
	printRow("kinE", package->kin_energy);
	printRow("totalE", package->total_energy);
}