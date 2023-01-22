#include "Analyzer.cuh"
#include "Printer.h"

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

void __global__ monitorCompoundEnergyKernel(Box* box, Float3* traj_buffer, float* potE_buffer, Float3* data_out) {		// everything here breaks if not all compounds are identical in particle count and particle mass!!!!!!!
	__shared__ Float3 energy[MAX_COMPOUND_PARTICLES];
	__shared__ Compound compound;


	uint64_t step = blockIdx.x + uint64_t{ 1 };
	const int compound_index = blockIdx.y;
	energy[threadIdx.x] = Float3(0.f);


	if (threadIdx.x == 0) {
		//printf("index: %d\n", compound_index + (step - 1) * N_MONITORBLOCKS_PER_STEP);
		data_out[compound_index + (step - 1) * box->n_compounds] = Float3{};
		//mass = box->compounds[compound_index].particles[0]
		compound = box->compounds[compound_index];
	}
	__syncthreads();

	if (threadIdx.x >= compound.n_particles) {
		return;
	}
	__syncthreads();

	const uint8_t& atom_type = compound.atom_types[threadIdx.x];
	const float mass = box->forcefield_device_box->particle_parameters[atom_type].mass;

	const uint32_t compound_offset = compound_index * MAX_COMPOUND_PARTICLES;
	const int step_offset = step * box->total_particles_upperbound;
	const float potE = potE_buffer[threadIdx.x + compound_offset + step_offset];

	const Float3 pos_tsub1 = traj_buffer[threadIdx.x + compound_offset + (step - uint64_t{ 1 }) * box->total_particles_upperbound];
	const Float3 pos_tadd1 = traj_buffer[threadIdx.x + compound_offset + (step + uint64_t{ 1 }) * box->total_particles_upperbound];
	const float n_steps = 2.f;

	float kinE = EngineUtils::calcKineticEnergy(&pos_tadd1, &pos_tsub1, mass, box->dt * n_steps);
	float totalE = potE + kinE;

	energy[threadIdx.x] = Float3(potE, kinE, totalE);
	__syncthreads();

	distributedSummation(energy, MAX_COMPOUND_PARTICLES);
	__syncthreads();

	if (threadIdx.x == 0) {
		data_out[compound_index + (step - 1) * box->n_compounds] = energy[0];
	}
}





void __global__ monitorSolventEnergyKernel(Box* box, Float3* traj_buffer, float* potE_buffer, Float3* data_out) {
	__shared__ Float3 energy[THREADS_PER_SOLVENTBLOCK];



	int solvent_index = threadIdx.x + blockIdx.y * THREADS_PER_SOLVENTBLOCK;
	int step = blockIdx.x + 1;
	int compounds_offset = box->n_compounds * MAX_COMPOUND_PARTICLES;


	energy[threadIdx.x] = Float3(0.f);
	if (threadIdx.x == 0) {
		data_out[(step - 1) * gridDim.y + blockIdx.y] = energy[0];
	}
	if (solvent_index >= box->n_solvents) { return; }


	
	Float3 pos_tsub1 = traj_buffer[compounds_offset + solvent_index + (step - 1) * box->total_particles_upperbound];
	Float3 pos_tadd1 = traj_buffer[compounds_offset + solvent_index + (step + 1) * box->total_particles_upperbound];

	float mass = 12.011000 * 1e-3f;
	float kinE = EngineUtils::calcKineticEnergy(&pos_tadd1, &pos_tsub1, mass, box->dt*2.f);

	float potE = potE_buffer[compounds_offset + solvent_index + step * box->total_particles_upperbound];

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
	EngineUtils::genericErrorCheck("Cuda error before analyzeEnergy\n");

	const auto n_steps = simulation->getStep();
	if (simulation->getStep() < 1) { return Analyzer::AnalyzedPackage(); }

	//printf("Analyzer malloc %.4f GB on host\n", sizeof(Float3) * analysable_steps * 1e-9);
	//Float3* average_energy = new Float3[analysable_steps];
	std::vector<Float3> average_energy;
	average_energy.resize(n_steps - 2);	// Ignore first and last step

	// We need to split up the analyser into steps, as we cannot store all positions traj on device at once.
	uint64_t max_steps_per_kernel = 1000;
	uint64_t particles_per_step = simulation->total_particles_upperbound;
	uint64_t max_values_per_kernel = (max_steps_per_kernel + 2) * particles_per_step;							// Pad steps with 2 for vel calculation
	printf("Analyzer malloc %.2f MB on device\n", (sizeof(Float3) + sizeof(double)) * (max_values_per_kernel) * 1e-6);
	cudaMalloc(&traj_buffer_device, sizeof(Float3) * max_values_per_kernel);
	cudaMalloc(&potE_buffer_device, sizeof(float) * max_values_per_kernel);

	for (int i = 0; i < ceil((double)n_steps / (double)max_steps_per_kernel); i++) {
		uint64_t step_offset = i * max_steps_per_kernel;												// offset one since we can't analyse step 1
		uint64_t steps_in_kernel = std::min(max_steps_per_kernel, n_steps - step_offset);

		// Create a array of len 1002, where index 0 and 1001 are padded values
		moveAndPadData(simulation, steps_in_kernel, step_offset);

		Float3* average_solvent_energy = analyzeSolvateEnergy(simulation, steps_in_kernel);
		std::vector<Float3> average_compound_energy = analyzeCompoundEnergy(simulation, steps_in_kernel);

		for (uint64_t ii = 0; ii < steps_in_kernel; ii++) {
			int64_t step = step_offset + ii - 1;	// -1 because index 0 is unused
			if (step == -1 || step >= n_steps-2u) { continue; }	// Dont save first step, as the kinE is slightly wrong
			average_energy[step] = (average_solvent_energy[ii] + average_compound_energy[ii]);
		}
		delete[] average_solvent_energy;
	}

	cudaFree(traj_buffer_device);
	cudaFree(potE_buffer_device);

	printH2("Finished analyzing energies", false, true);
	return AnalyzedPackage(average_energy, simulation->temperature_buffer, simulation->n_temp_values);
}

void Analyzer::moveAndPadData(Simulation* simulation, uint64_t steps_in_kernel, uint64_t step_offset) {
	const uint64_t particles_per_step = simulation->total_particles_upperbound;

	// First move the middle bulk to device
	cudaMemcpy(&traj_buffer_device[1 * particles_per_step], &simulation->traj_buffer[step_offset * particles_per_step], sizeof(Float3) * steps_in_kernel * particles_per_step, cudaMemcpyHostToDevice);
	cudaMemcpy(&potE_buffer_device[1 * particles_per_step], &simulation->potE_buffer[step_offset * particles_per_step], sizeof(float) * steps_in_kernel * particles_per_step, cudaMemcpyHostToDevice);
	EngineUtils::genericErrorCheck("Cuda error during analyzer transfer2\n");

	// Then pad the front. If step 0, we pad with zero. If step n we pad with n-1
	uint64_t paddingSrcIndex = step_offset == 0 ? 0 : step_offset - 1;
	cudaMemcpy(&traj_buffer_device[0], &simulation->traj_buffer[paddingSrcIndex * particles_per_step], sizeof(Float3) * particles_per_step, cudaMemcpyHostToDevice);
	cudaMemcpy(&potE_buffer_device[0], &simulation->potE_buffer[paddingSrcIndex * particles_per_step], sizeof(float) * particles_per_step, cudaMemcpyHostToDevice);
	EngineUtils::genericErrorCheck("Cuda error during analyzer transfer1\n");

	// Then pad the end. if step
	const uint64_t end_index = step_offset + steps_in_kernel;
	paddingSrcIndex = end_index == simulation->getStep() ? end_index - 1 : end_index;
	cudaMemcpy(&traj_buffer_device[(steps_in_kernel + 1) * particles_per_step], &simulation->traj_buffer[paddingSrcIndex * particles_per_step], sizeof(Float3) * particles_per_step, cudaMemcpyHostToDevice);
	cudaMemcpy(&potE_buffer_device[(steps_in_kernel + 1) * particles_per_step], &simulation->potE_buffer[paddingSrcIndex * particles_per_step], sizeof(float) * particles_per_step, cudaMemcpyHostToDevice);

	cudaDeviceSynchronize();
	EngineUtils::genericErrorCheck("Cuda error during analyzer transfer\n");
}

// TODO: Fix this fucntion like the compound one
Float3* Analyzer::analyzeSolvateEnergy(Simulation* simulation, uint64_t n_steps) {
	// Start by creating array of energies of value 0
	Float3* average_solvent_energy = new Float3[n_steps];
	for (int i = 0; i < n_steps; i++)
		average_solvent_energy[i] = Float3(0.f);

	// If any solvents are present, fill above array
	if (simulation->n_solvents > 0) {

		Float3* average_solvent_energy_blocked = new Float3[n_steps * simulation->blocks_per_solventkernel];
		Float3* data_out;
		cudaMalloc(&data_out, sizeof(Float3) * simulation->blocks_per_solventkernel * n_steps);

		dim3 block_dim(n_steps, simulation->blocks_per_solventkernel, 1);
		monitorSolventEnergyKernel << < block_dim, THREADS_PER_SOLVENTBLOCK >> > (simulation->box, traj_buffer_device, potE_buffer_device, data_out);
		cudaDeviceSynchronize();
		EngineUtils::genericErrorCheck("Cuda error during analyzeSolvateEnergy\n");

		cudaMemcpy(average_solvent_energy_blocked, data_out, sizeof(Float3) * simulation->blocks_per_solventkernel * n_steps, cudaMemcpyDeviceToHost);
		cudaDeviceSynchronize();

		for (uint64_t step = 0; step < n_steps; step++) {
			average_solvent_energy[step] = Float3(0.f);
			for (int block = 0; block < simulation->blocks_per_solventkernel; block++) {
				average_solvent_energy[step] += average_solvent_energy_blocked[block + step * simulation->blocks_per_solventkernel];
			}
			average_solvent_energy[step] *= (1.f / simulation->n_solvents);
		}

		cudaFree(data_out);
		delete[] average_solvent_energy_blocked;
	}

	return average_solvent_energy;
}


std::vector<Float3> Analyzer::analyzeCompoundEnergy(Simulation* simulation, uint64_t steps_in_kernel) {
	uint64_t n_datapoints = simulation->n_compounds * steps_in_kernel;

	//Float3* average_compound_energy = new Float3[n_steps];
	std::vector<Float3> total_compound_energy(steps_in_kernel);
	//total_compound_energy.resize(steps_in_kernel);

	for (int i = 0; i < steps_in_kernel; i++)
		total_compound_energy[i] = Float3(0.f);


	if (simulation->total_compound_particles > 0) {
		Float3* host_data = new Float3[n_datapoints];

		Float3* data_out;
		cudaMalloc(&data_out, sizeof(Float3) * n_datapoints);

		dim3 block_dim(steps_in_kernel, simulation->box->n_compounds, 1);
		monitorCompoundEnergyKernel << < block_dim, MAX_COMPOUND_PARTICLES >> > (simulation->box, traj_buffer_device, potE_buffer_device, data_out);
		cudaDeviceSynchronize();
		EngineUtils::genericErrorCheck("Cuda error during analyzeCompoundEnergy\n");

		cudaMemcpy(host_data, data_out, sizeof(Float3) * n_datapoints, cudaMemcpyDeviceToHost);
		cudaDeviceSynchronize();


		for (uint64_t step = 0; step < steps_in_kernel; step++) {
			for (uint64_t i = 0; i < simulation->box->n_compounds; i++) {
				total_compound_energy[step] += host_data[i + step * simulation->box->n_compounds];
			}
			//total_compound_energy[step] *= (1.f / (simulation->total_compound_particles));
		}

		cudaFree(data_out);
		delete[] host_data;
	}

	return total_compound_energy;
}

float getMin(std::vector<float>& vec) {
	return *std::min_element(vec.begin(), vec.end());
}

float getMax(std::vector<float>& vec) {
	return *std::max_element(vec.begin(), vec.end());
}

float getMean(vector<float> vec)
{
	double sum = 0.;
	for (auto elem : vec) { sum += static_cast<double>(elem); }	
	return static_cast<float>(sum / static_cast<double>(vec.size()));
}

float getStdDev(std::vector<float>& vec) {
	if (vec.size() == 0) { return 0.f; }

	double mean = getMean(vec);

	double err_sum = 0;
	for (auto elem : vec) { err_sum += (elem - mean) * (elem - mean); }

	double mean_err_sq = err_sum / static_cast<double>(vec.size());
	return static_cast<float>(std::abs(std::sqrt(mean_err_sq)));
}

float Analyzer::getStdDevNorm(std::vector<float>& vec) {

	return vec.front() != 0.f ? getStdDev(vec) / vec.front() : 0.f;
}

void printRow(string title, vector<float>& vec) {
	LIMA_Printer::printTableRow(
		title, { 
			getMin(vec), 
			getMax(vec), 
			getStdDev(vec) / vec.front(),
			(vec.back() - vec.front()) / vec.front() });
}

void Analyzer::printEnergy(AnalyzedPackage* package) {
	LIMA_Printer::printTableRow({ "", "min", "max", "Std. deviation", "Change 0->n"});
	printRow("potE", package->pot_energy);
	printRow("kinE", package->kin_energy);
	printRow("totalE", package->total_energy);
}