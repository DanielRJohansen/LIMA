#include "Analyzer.cuh"
#include "Printer.h"
#include "Forcefield.cuh"
#include "EngineUtils.cuh"

#include <algorithm>

using namespace LIMA_Print;

const int THREADS_PER_SOLVENTBLOCK_ANALYZER = 128;


template<typename T>
void __device__ distributedSummation(T* arrayptr, int array_len) {				// Places the result at pos 0 of input_array
	T temp;			// This is a lazy soluation, but maybe it is also fast? Definitely simple..
	for (int i = 1; i < array_len; i *= 2) {	// Distributed averaging							// Make a generic and SAFER function for this, PLEASE OK??
		if ((threadIdx.x + i) < array_len) {
			temp = arrayptr[threadIdx.x] + arrayptr[threadIdx.x + i];
		}
		__syncthreads();
		arrayptr[threadIdx.x] = temp;
		__syncthreads();
	}
}

void __global__ monitorCompoundEnergyKernel(Compound* compounds, const ForceField_NB* const forcefield, const BoxParams boxparams, float* potE_buffer, float* vel_buffer, Float3* data_out) {		// everything here breaks if not all compounds are identical in particle count and particle mass!!!!!!!
	__shared__ Float3 energy[MAX_COMPOUND_PARTICLES];
	__shared__ Compound compound;


	const int64_t step = blockIdx.x;	// Step relative to current batch
	const int64_t compound_index = blockIdx.y;
	const int64_t particle_index = threadIdx.x;
	energy[particle_index] = Float3(0.f);


	if (particle_index == 0) {
		data_out[compound_index + (step) * boxparams.n_compounds] = Float3{};
		compound = compounds[compound_index];
	}
	__syncthreads();

	if (particle_index >= compound.n_particles) {
		return;
	}
	__syncthreads();

	const uint8_t atom_type = compound.atom_types[particle_index];
	const float mass = forcefield->particle_parameters[atom_type].mass;

	const int64_t compound_offset = compound_index * MAX_COMPOUND_PARTICLES;
	const int64_t step_offset = step * boxparams.total_particles_upperbound;
	const float potE = potE_buffer[particle_index + compound_offset + step_offset];

	const float speed = vel_buffer[particle_index + compound_offset + step_offset];
	const float kinE = EngineUtils::calcKineticEnergy(speed, mass);	// remove direction from vel

	const float totalE = potE + kinE;

	energy[particle_index] = Float3(potE, kinE, totalE);
	__syncthreads();

	distributedSummation(energy, MAX_COMPOUND_PARTICLES);
	__syncthreads();

	if (particle_index == 0) {
		data_out[compound_index + (step) * boxparams.n_compounds] = energy[0];
	}
}





void __global__ monitorSolventEnergyKernel(const BoxParams boxparams, float* potE_buffer, float* vel_buffer, Float3* data_out) {
	__shared__ Float3 energy[THREADS_PER_SOLVENTBLOCK_ANALYZER];


	const int solvent_index = threadIdx.x + blockIdx.y * THREADS_PER_SOLVENTBLOCK_ANALYZER;
	const int step = blockIdx.x;
	const int compounds_offset = boxparams.n_compounds * MAX_COMPOUND_PARTICLES;
	const int step_offset = step * boxparams.total_particles_upperbound;

	energy[threadIdx.x] = Float3(0.f);
	if (threadIdx.x == 0) {
		data_out[(step) * gridDim.y + blockIdx.y] = energy[0];
	}
	if (solvent_index >= boxparams.n_solvents) { return; }


	const float velocity = vel_buffer[step_offset + compounds_offset + solvent_index];
	const float kinE = EngineUtils::calcKineticEnergy(velocity, SOLVENT_MASS);	// remove direction from vel
	float potE = potE_buffer[compounds_offset + solvent_index + step * boxparams.total_particles_upperbound];

	const float totalE = potE + kinE;

	energy[threadIdx.x] = Float3(potE, kinE, totalE);
	__syncthreads();
	distributedSummation(energy, THREADS_PER_SOLVENTBLOCK_ANALYZER);
	if (threadIdx.x == 0) {
		data_out[(step) * gridDim.y + blockIdx.y] = energy[0];
	}
}



Analyzer::AnalyzedPackage Analyzer::analyzeEnergy(Simulation* simulation) {	// Calculates the avg J/mol // calculate energies separately for compounds and solvents. weigh averages based on amount of each
	LIMA_UTILS::genericErrorCheck("Cuda error before analyzeEnergy\n");

	//const int64_t n_steps = simulation->getStep();
	const int64_t n_entryindices = LIMALOGSYSTEM::getMostRecentDataentryIndex(simulation->getStep());

	if (n_entryindices < 2) { return Analyzer::AnalyzedPackage(); }


	// First set up some stuff needed on device, that is currently on host
	cudaMalloc(&forcefield_device, sizeof(ForceField_NB));
	cudaMemcpy(forcefield_device, &simulation->forcefield->getNBForcefieldRef(), sizeof(ForceField_NB), cudaMemcpyHostToDevice);

	if (simulation->boxparams_host.n_compounds > 0) {
		cudaMalloc(&compounds_device, sizeof(Compound) * simulation->compounds_host.size());
		cudaMemcpy(compounds_device, simulation->compounds_host.data(), sizeof(Compound) * simulation->compounds_host.size(), cudaMemcpyHostToDevice);
	}
	


	std::vector<Float3> average_energy;
	average_energy.resize(n_entryindices - 2);	// Ignore first and last step	// TODO: Rework this, no longer necessary as we use VVS

	// We need to split up the analyser into steps, as we cannot store all positions traj on device at once.
	int64_t max_steps_per_kernel = 100;
	int64_t particles_per_step = simulation->boxparams_host.total_particles_upperbound;
	int64_t max_values_per_kernel = max_steps_per_kernel * particles_per_step;							// Pad steps with 2 for vel calculation

	const std::string bytesize = std::to_string((sizeof(Float3) + sizeof(double)) * (max_values_per_kernel) * 1e-6);
	m_logger->print("Analyzer malloc " + bytesize + " MB on device\n");
	cudaMalloc(&potE_buffer_device, sizeof(float) * max_values_per_kernel);
	cudaMalloc(&vel_buffer_device, sizeof(float) * max_values_per_kernel);

	for (int64_t i = 0; i < ceil((double)n_entryindices / (double)max_steps_per_kernel); i++) {
		const int64_t step_offset = i * max_steps_per_kernel;												// offset one since we can't analyse step 1
		const int64_t steps_in_kernel = std::min(max_steps_per_kernel, n_entryindices - step_offset);

		cudaMemcpy(potE_buffer_device, &simulation->potE_buffer->data()[step_offset * particles_per_step], sizeof(float) * steps_in_kernel * particles_per_step, cudaMemcpyHostToDevice);
		cudaMemcpy(vel_buffer_device, &simulation->vel_buffer->data()[step_offset * particles_per_step], sizeof(float) * steps_in_kernel * particles_per_step, cudaMemcpyHostToDevice);
		LIMA_UTILS::genericErrorCheck("Cuda error during analyzer transfer2\n");

		std::vector<Float3> average_solvent_energy = analyzeSolvateEnergy(simulation, steps_in_kernel);
		std::vector<Float3> average_compound_energy = analyzeCompoundEnergy(simulation, steps_in_kernel);

		for (int64_t ii = 0; ii < steps_in_kernel; ii++) {
			int64_t step = step_offset + ii - 1;	// -1 because index 0 is unused
			if (step == -1 || step >= n_entryindices -2u) { continue; }	// Dont save first step, as the kinE is slightly wrong
			average_energy[step] = (average_solvent_energy[ii] + average_compound_energy[ii]);
		}
	}

	cudaFree(potE_buffer_device);
	cudaFree(vel_buffer_device);
	cudaFree(forcefield_device);
	if (simulation->boxparams_host.n_compounds > 0) {
		cudaFree(compounds_device);
	}

	m_logger->finishSection("Finished analyzing energies");
	return AnalyzedPackage(average_energy, simulation->temperature_buffer);
}

std::vector<Float3> Analyzer::analyzeSolvateEnergy(Simulation* simulation, uint64_t n_steps) {
	// Start by creating array of energies of value 0
	std::vector<Float3> average_solvent_energy(n_steps);

	int blocks_per_solventkernel = (int)ceil((float)simulation->boxparams_host.n_solvents / (float)THREADS_PER_SOLVENTBLOCK_ANALYZER);

	// If any solvents are present, fill above array
	if (simulation->boxparams_host.n_solvents > 0) {

		std::vector<Float3> average_solvent_energy_blocked(n_steps * blocks_per_solventkernel);
		Float3* data_out;
		cudaMalloc(&data_out, sizeof(Float3) * blocks_per_solventkernel * n_steps);

		dim3 block_dim(n_steps, blocks_per_solventkernel, 1);
		monitorSolventEnergyKernel << < block_dim, THREADS_PER_SOLVENTBLOCK_ANALYZER >> > (simulation->boxparams_host, potE_buffer_device, vel_buffer_device, data_out);	// TODO: FIx
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

		dim3 block_dim(static_cast<uint32_t>(steps_in_kernel), simulation->boxparams_host.n_compounds, 1);
		monitorCompoundEnergyKernel << < block_dim, MAX_COMPOUND_PARTICLES >> > (compounds_device, forcefield_device, simulation->boxparams_host, potE_buffer_device, vel_buffer_device, data_out);
		cudaDeviceSynchronize();
		LIMA_UTILS::genericErrorCheck("Cuda error during analyzeCompoundEnergy\n");

		cudaMemcpy(host_data.data(), data_out, sizeof(Float3) * n_datapoints, cudaMemcpyDeviceToHost);
		cudaFree(data_out);


		for (uint64_t step = 0; step < steps_in_kernel; step++) {
			for (uint64_t i = 0; i < simulation->boxparams_host.n_compounds; i++) {
				total_compound_energy[step] += host_data[i + step * simulation->boxparams_host.n_compounds];
			}
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

float getVarianceCoefficient(const std::vector<float>& vec) {
	if (vec.empty()) { return 0.f; } 
	const float stddev = getStdDev(vec);
	const float mean = getMean(vec);

	if (stddev == 0.f && mean == 0.f) { return 0.f; }
	return  stddev / std::abs(mean);
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





float calculateSlopeLinearRegression(const std::vector<float>& y_values, const float mean) {
	size_t n = y_values.size();
	float sum_x = 0;
	float sum_y = 0;
	float sum_xy = 0;
	float sum_xx = 0;

	for (size_t i = 0; i < n; ++i) {
		sum_x += i;
		sum_y += y_values[i];
		sum_xy += i * y_values[i];
		sum_xx += i * i;
	}

	const float slope = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x);
	const float slope_coefficient = slope / mean;
	return slope_coefficient;
}

Analyzer::AnalyzedPackage::AnalyzedPackage(std::vector<Float3>& avg_energy, std::vector<float> temperature) {
	energy_data = avg_energy;
	//auto e_cnt = energy_data.size();

	temperature_data = temperature;
	//memcpy(temperature_data.data(), t_ptr, t_cnt);

	auto e_cnt = energy_data.size();
	pot_energy.resize(e_cnt);
	kin_energy.resize(e_cnt);
	total_energy.resize(e_cnt);
	for (int i = 0; i < e_cnt; i++) {
		pot_energy[i] = energy_data[i].x;
		kin_energy[i] = energy_data[i].y;
		total_energy[i] = energy_data[i].z;
	}

	mean_energy = getMean(total_energy);

	energy_gradient = calculateSlopeLinearRegression(total_energy, mean_energy);
	variance_coefficient = getVarianceCoefficient(total_energy);
}


























void Analyzer::findAndDumpPiecewiseEnergies(const Simulation& sim, const std::string& workdir) {
	std::vector<float> energies;
	
	for (auto entryindex = 0; entryindex < LIMALOGSYSTEM::getMostRecentDataentryIndex(sim.getStep()-1); entryindex++) {

		for (int compound_id = 0; compound_id < sim.boxparams_host.n_compounds; compound_id++) {
			for (int particle_id = 0; particle_id < MAX_COMPOUND_PARTICLES; particle_id++) {
				
				const float potE = sim.potE_buffer->getCompoundparticleDatapointAtIndex(compound_id, particle_id, entryindex);

				const uint8_t& atom_type = sim.compounds_host[compound_id].atom_types[particle_id];
				const float mass = sim.forcefield->getNBForcefield().particle_parameters[atom_type].mass;
				const float vel = sim.vel_buffer->getCompoundparticleDatapointAtIndex(compound_id, particle_id, entryindex);
				const float kinE = EngineUtils::calcKineticEnergy(vel, mass);
				
				energies.emplace_back(potE);
				energies.emplace_back(kinE);
			}
		}

		for (int solvent_id = 0; solvent_id < sim.boxparams_host.n_solvents; solvent_id++) {

			const float potE = sim.potE_buffer->getSolventparticleDatapointAtIndex(solvent_id, entryindex);

			const float mass = sim.forcefield->getNBForcefield().particle_parameters[ATOMTYPE_SOLVENT].mass;
			const float vel = sim.vel_buffer->getSolventparticleDatapointAtIndex(solvent_id, entryindex);
			const float kinE = EngineUtils::calcKineticEnergy(vel, mass);			

			energies.emplace_back(potE);
			energies.emplace_back(kinE);
		}
	}

	Filehandler::dumpToFile(energies.data(), energies.size(), workdir + "/PiecewiseEnergy.bin");
}

