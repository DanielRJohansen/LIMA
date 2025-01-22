#include "Analyzer.cuh"
#include "DeviceAlgorithms.cuh"
#include "PhysicsUtils.cuh"


#include "Constants.h"
#include "Printer.h"
#include "Statistics.h"


#include <algorithm>
#include <numeric>

using namespace LIMA_Print;

const int THREADS_PER_SOLVENTBLOCK_ANALYZER = 128;




void __global__ monitorCompoundEnergyKernel(Compound* compounds, const ForceField_NB* const forcefield, const BoxParams boxparams, float* potE_buffer, float* vel_buffer, Float3* data_out) {		// everything here breaks if not all compounds are identical in particle count and particle mass!!!!!!!
	__shared__ Float3 energy[MAX_COMPOUND_PARTICLES];
	__shared__ Compound compound;


	const int64_t step = blockIdx.x;	// Step relative to current batch
	const int64_t compound_index = blockIdx.y;
	const int64_t particle_index = threadIdx.x;
	energy[particle_index] = Float3(0.f);


	if (particle_index == 0) {
		data_out[compound_index + (step) * boxparams.n_compounds] = Float3{};
		compound = compounds[compound_index]; // TODO: this is silly, just use the global one directly
	}
	__syncthreads();

	if (particle_index >= compound.n_particles) {
		return;
	}
	__syncthreads();

	const float mass = compound.atomMasses[particle_index];

	const int64_t compound_offset = compound_index * MAX_COMPOUND_PARTICLES;
	const int64_t step_offset = step * boxparams.total_particles_upperbound;
	const float potE = potE_buffer[particle_index + compound_offset + step_offset];

	const float speed = vel_buffer[particle_index + compound_offset + step_offset];
	const float kinE = PhysicsUtils::calcKineticEnergy(speed, mass);	// remove direction from vel

	const float totalE = potE + kinE;

	energy[particle_index] = Float3(potE, kinE, totalE);
	__syncthreads();

	LAL::distributedSummation(energy, MAX_COMPOUND_PARTICLES);
	__syncthreads();

	if (particle_index == 0) {
		data_out[compound_index + (step) * boxparams.n_compounds] = energy[0];
	}
}





void __global__ monitorSolventEnergyKernel(const BoxParams boxparams, float* potE_buffer, float* vel_buffer, Float3* data_out, const TinyMolState* const tinyMols, const ForcefieldTinymol* const forcefield) {
	__shared__ Float3 energy[THREADS_PER_SOLVENTBLOCK_ANALYZER];


	const int solvent_index = threadIdx.x + blockIdx.y * THREADS_PER_SOLVENTBLOCK_ANALYZER;
	const int64_t step = blockIdx.x;
	const int compounds_offset = boxparams.n_compounds * MAX_COMPOUND_PARTICLES;
	const int64_t step_offset = step * boxparams.total_particles_upperbound;

	energy[threadIdx.x] = Float3(0.f);
	if (threadIdx.x == 0) {
		data_out[(step) * gridDim.y + blockIdx.y] = energy[0];
	}
	if (solvent_index >= boxparams.n_solvents) { return; }

	const float mass = forcefield->types[tinyMols[solvent_index].tinymolTypeIndex].mass;
	const float velocity = vel_buffer[step_offset + compounds_offset + solvent_index];
	const float kinE = PhysicsUtils::calcKineticEnergy(velocity, mass);	// remove direction from vel
	float potE = potE_buffer[compounds_offset + solvent_index + step * boxparams.total_particles_upperbound];

	const float totalE = potE + kinE;

	energy[threadIdx.x] = Float3(potE, kinE, totalE);
	__syncthreads();
	LAL::distributedSummation(energy, THREADS_PER_SOLVENTBLOCK_ANALYZER);
	if (threadIdx.x == 0) {
		data_out[(step) * gridDim.y + blockIdx.y] = energy[0];
	}
}


const int NUM_BINS = 16;
const int BIN_BASE = 10;
__device__ int getBinIndex(float value) {
	int binIndex = (value > 0) ? log10f(value) / log10f(BIN_BASE) : -(log10f(-value) / log10f(BIN_BASE));
	binIndex += NUM_BINS / 2;  // Center the bins around zero
	return std::min(std::max(binIndex, 0), NUM_BINS - 1);  // Clamp to valid range
}

__global__ void potEHistogramKernel(Compound* compounds, int total_particles_upperbound, float* potE_buffer, int* histogramData, int64_t step) {
	__shared__ int shared_histogram[NUM_BINS];

	const int compound_index = blockIdx.x;   // Unique index for each compound
	const int particle_index = threadIdx.x;  // Unique index for each particle within a compound

	// Initialize shared histogram to zero
	if (particle_index < NUM_BINS) {
		shared_histogram[particle_index] = 0;
	}
	__syncthreads();
	if (particle_index >= compounds[compound_index].n_particles) {
		return;
	}

	// Calculate the offsets
	const int64_t compound_offset = compound_index * MAX_COMPOUND_PARTICLES;
	const float potE = potE_buffer[particle_index + compound_offset];
	// Determine the bin index for the current potential energy
	int binIndex = getBinIndex(potE);

	// Atomically increment the appropriate bin in the shared histogram
	atomicAdd(&shared_histogram[binIndex], 1);
	__syncthreads();

	// First thread writes the shared histogram to the global histogram
	if (particle_index == 0) {
		for (int i = 0; i < NUM_BINS; ++i) {
			if (shared_histogram[i] > 0) {
				atomicAdd(&histogramData[i], shared_histogram[i]);
			}
		}
	}
}











std::vector<Float3> analyzeSolvateEnergy(Simulation* simulation, uint64_t n_steps, float* potE_buffer_device, float* vel_buffer_device, const ForcefieldTinymol* const forcefield_device, const TinyMolState* const tinyMols) {
	// Start by creating array of energies of value 0
	std::vector<Float3> average_solvent_energy(n_steps);

	int blocks_per_solventkernel = (int)ceil((float)simulation->box_host->boxparams.n_solvents / (float)THREADS_PER_SOLVENTBLOCK_ANALYZER);

	// If any solvents are present, fill above array
	if (simulation->box_host->boxparams.n_solvents > 0) {

		std::vector<Float3> average_solvent_energy_blocked(n_steps * blocks_per_solventkernel);
		Float3* data_out;
		cudaMalloc(&data_out, sizeof(Float3) * blocks_per_solventkernel * n_steps);

		dim3 block_dim(n_steps, blocks_per_solventkernel, 1);
		monitorSolventEnergyKernel << < block_dim, THREADS_PER_SOLVENTBLOCK_ANALYZER >> > (simulation->box_host->boxparams, potE_buffer_device, vel_buffer_device, data_out, tinyMols, forcefield_device);
		LIMA_UTILS::genericErrorCheck("Cuda error during analyzeSolvateEnergy\n");

		cudaMemcpy(average_solvent_energy_blocked.data(), data_out, sizeof(Float3) * blocks_per_solventkernel * n_steps, cudaMemcpyDeviceToHost);
		cudaDeviceSynchronize();
		cudaFree(data_out);

		for (uint64_t step = 0; step < n_steps; step++) {
			average_solvent_energy[step] = Float3(0.f);
			for (int block = 0; block < blocks_per_solventkernel; block++) {
				average_solvent_energy[step] += average_solvent_energy_blocked[block + step * blocks_per_solventkernel];
			}
			//average_solvent_energy[step] *= (1.f / simulation->box_host->boxparams.n_solvents);
		}

	}

	return average_solvent_energy;
}

std::vector<Float3> analyzeCompoundEnergy(Simulation* simulation, uint64_t steps_in_kernel, float* potE_buffer_device, float* vel_buffer_device, Compound* compounds_device, ForceField_NB* forcefield_device) {
	const uint64_t n_datapoints = simulation->box_host->boxparams.n_compounds * steps_in_kernel;

	std::vector<Float3> total_compound_energy(steps_in_kernel);

	if (simulation->box_host->boxparams.total_compound_particles > 0) {
		std::vector<Float3> host_data(n_datapoints);

		Float3* data_out;
		cudaMalloc(&data_out, sizeof(Float3) * n_datapoints);

		dim3 block_dim(static_cast<uint32_t>(steps_in_kernel), simulation->box_host->boxparams.n_compounds, 1);
		monitorCompoundEnergyKernel << < block_dim, MAX_COMPOUND_PARTICLES >> > (compounds_device, forcefield_device, simulation->box_host->boxparams, potE_buffer_device, vel_buffer_device, data_out);
		cudaDeviceSynchronize();
		LIMA_UTILS::genericErrorCheck("Cuda error during analyzeCompoundEnergy\n");

		cudaMemcpy(host_data.data(), data_out, sizeof(Float3) * n_datapoints, cudaMemcpyDeviceToHost);
		cudaFree(data_out);


		for (uint64_t step = 0; step < steps_in_kernel; step++) {
			for (uint64_t i = 0; i < simulation->box_host->boxparams.n_compounds; i++) {
				total_compound_energy[step] += host_data[i + step * simulation->box_host->boxparams.n_compounds];
			}
		}

	}

	return total_compound_energy;
}

SimAnalysis::AnalyzedPackage SimAnalysis::analyzeEnergy(Simulation* simulation) {	// Calculates the avg J/mol // calculate energies separately for compounds and solvents. weigh averages based on amount of each
	LIMA_UTILS::genericErrorCheck("Cuda error before analyzeEnergy\n");

	const int64_t n_entryindices = LIMALOGSYSTEM::getMostRecentDataentryIndex(simulation->getStep(), simulation->simparams_host.data_logging_interval);

	if (n_entryindices < 2) { return AnalyzedPackage(); }




	// First set up some stuff needed on device, that is currently on host
	float* potE_buffer_device = nullptr;
	float* vel_buffer_device = nullptr;
	Compound* compounds_device = nullptr;	
	ForceField_NB* forcefield_device = nullptr;
	cudaMalloc(&forcefield_device, sizeof(ForceField_NB));
	cudaMemcpy(forcefield_device, &simulation->forcefield, sizeof(ForceField_NB), cudaMemcpyHostToDevice);

	ForcefieldTinymol* tinymolForcefield_device = GenericCopyToDevice(&simulation->forcefieldTinymol, 1);
	TinyMolState* tinymols = GenericCopyToDevice(simulation->box_host->tinyMols);

	if (simulation->box_host->boxparams.n_compounds > 0) {
		cudaMalloc(&compounds_device, sizeof(Compound) * simulation->box_host->boxparams.n_compounds);
		cudaMemcpy(compounds_device, simulation->box_host->compounds.data(), sizeof(Compound) * simulation->box_host->boxparams.n_compounds, cudaMemcpyHostToDevice);
	}



	std::vector<Float3> average_energy;
	average_energy.resize(n_entryindices - 2);	// Ignore first and last step	// TODO: Rework this, no longer necessary as we use VVS

	// We need to split up the analyser into steps, as we cannot store all positions traj on device at once.
	int64_t max_steps_per_kernel = 100;
	int64_t particles_per_step = simulation->box_host->boxparams.total_particles_upperbound;
	int64_t max_values_per_kernel = max_steps_per_kernel * particles_per_step;							// Pad steps with 2 for vel calculation

	const std::string bytesize = std::to_string((sizeof(Float3) + sizeof(double)) * (max_values_per_kernel) * 1e-6);
	//m_logger->print("Analyzer malloc " + bytesize + " MB on device\n");
	cudaMalloc(&potE_buffer_device, sizeof(float) * max_values_per_kernel);
	cudaMalloc(&vel_buffer_device, sizeof(float) * max_values_per_kernel);

	for (int64_t i = 0; i < ceil((double)n_entryindices / (double)max_steps_per_kernel); i++) {
		const int64_t step_offset = i * max_steps_per_kernel;												// offset one since we can't analyse step 1
		const int64_t steps_in_kernel = std::min(max_steps_per_kernel, n_entryindices - step_offset);

		cudaMemcpy(potE_buffer_device, &simulation->potE_buffer->data()[step_offset * particles_per_step], sizeof(float) * steps_in_kernel * particles_per_step, cudaMemcpyHostToDevice);
		cudaMemcpy(vel_buffer_device, &simulation->vel_buffer->data()[step_offset * particles_per_step], sizeof(float) * steps_in_kernel * particles_per_step, cudaMemcpyHostToDevice);
		LIMA_UTILS::genericErrorCheck("Cuda error during analyzer transfer2\n");

		std::vector<Float3> average_solvent_energy = analyzeSolvateEnergy(simulation, steps_in_kernel, potE_buffer_device, vel_buffer_device, tinymolForcefield_device, tinymols);
		std::vector<Float3> average_compound_energy = analyzeCompoundEnergy(simulation, steps_in_kernel, potE_buffer_device, vel_buffer_device, compounds_device, forcefield_device);

		for (int64_t ii = 0; ii < steps_in_kernel; ii++) {
			int64_t step = step_offset + ii - 1;	// -1 because index 0 is unused
			if (step == -1 || step >= n_entryindices - 2u) { continue; }	// Dont save first step, as the kinE is slightly wrong
			average_energy[step] = (average_solvent_energy[ii] + average_compound_energy[ii]);
		}
	}

	cudaFree(potE_buffer_device);
	cudaFree(vel_buffer_device);
	cudaFree(forcefield_device);
	if (simulation->box_host->boxparams.n_compounds > 0) {
		cudaFree(compounds_device);
	}

	cudaFree(tinymolForcefield_device);
	cudaFree(tinymols);

	//m_logger->finishSection("Finished analyzing energies");
	return AnalyzedPackage(average_energy, simulation->temperature_buffer);
}

float getMin(const std::vector<float>& vec) {
	return *std::min_element(vec.begin(), vec.end());
}

float getMax(const std::vector<float>& vec) {
	return *std::max_element(vec.begin(), vec.end());
}

float getVarianceCoefficient(const std::vector<float>& vec) {
	if (vec.empty()) { return 0.f; } 
	const float stddev = Statistics::StdDev(vec);
	const float mean = Statistics::Mean(vec);

	if (stddev == 0.f && mean == 0.f) { return 0.f; }
	return  stddev / std::abs(mean);
}

void printRow(string title, const std::vector<float>& vec) {
	if (vec.empty()) { return; }
	LIMA_Printer::printTableRow(
		title, { 
			getMin(vec), 
			getMax(vec), 
			Statistics::StdDev(vec),
			(vec.back() - vec.front()) / vec.front() });
}

void SimAnalysis::AnalyzedPackage::Print() const {
	LIMA_Printer::printTableRow({ "", "min", "max", "Std. deviation", "Change 0->n"});
	printRow("potE", pot_energy);
	printRow("kinE", kin_energy);
	printRow("totalE", total_energy);
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

SimAnalysis::AnalyzedPackage::AnalyzedPackage(std::vector<Float3>& avg_energy, std::vector<float> temperature) {
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

	mean_energy = Statistics::Mean(total_energy);

	energy_gradient = calculateSlopeLinearRegression(total_energy, mean_energy);
	variance_coefficient = getVarianceCoefficient(total_energy);
}






















//
//
//
//
//void Analyzer::findAndDumpPiecewiseEnergies(const Simulation& sim, const std::string& workdir) {
//	std::vector<float> energies;
//	
//	for (auto entryindex = 0; entryindex < LIMALOGSYSTEM::getMostRecentDataentryIndex(sim.getStep()-1, sim.simparams_host.data_logging_interval); entryindex++) {
//
//		for (int compound_id = 0; compound_id < sim.box_host->boxparams.n_compounds; compound_id++) {
//			for (int particle_id = 0; particle_id < MAX_COMPOUND_PARTICLES; particle_id++) {
//				
//				const float potE = sim.potE_buffer->getCompoundparticleDatapointAtIndex(compound_id, particle_id, entryindex);
//
//				const uint8_t& atom_type = sim.box_host->compounds[compound_id].atom_types[particle_id];
//				const float mass = sim.forcefield.particle_parameters[atom_type].mass;
//				const float vel = sim.vel_buffer->getCompoundparticleDatapointAtIndex(compound_id, particle_id, entryindex);
//				const float kinE = PhysicsUtils::calcKineticEnergy(vel, mass);
//				
//				energies.emplace_back(potE);
//				energies.emplace_back(kinE);
//			}
//		}
//
//		for (int solvent_id = 0; solvent_id < sim.box_host->boxparams.n_solvents; solvent_id++) {
//
//			const float potE = sim.potE_buffer->getSolventparticleDatapointAtIndex(solvent_id, entryindex);
//
//			const float mass = sim.forcefield.particle_parameters[ATOMTYPE_SOLVENT].mass;
//			const float vel = sim.vel_buffer->getSolventparticleDatapointAtIndex(solvent_id, entryindex);
//			const float kinE = PhysicsUtils::calcKineticEnergy(vel, mass);
//
//			energies.emplace_back(potE);
//			energies.emplace_back(kinE);
//		}
//	}
//
//	FileUtils::dumpToFile(energies.data(), energies.size(), workdir + "/PiecewiseEnergy.bin");
//}




std::vector<int64_t> MakeBinLabels() {
	std::vector<int64_t> bins;

	int64_t current_bin = 10;

	while (bins.size() < NUM_BINS / 2) {
		bins.push_back(current_bin);
		current_bin = (current_bin == 0) ? 10 : current_bin * 10;
	}

	std::vector<double> negative_bins;
	current_bin = -10;
	while (negative_bins.size() < NUM_BINS / 2) {
		negative_bins.push_back(current_bin);
		current_bin *= 10;
	}

	std::reverse(negative_bins.begin(), negative_bins.end());
	bins.insert(bins.begin(), negative_bins.begin(), negative_bins.end());

	return bins;
}

void SimAnalysis::PlotPotentialEnergyDistribution(const Simulation& simulation, const std::filesystem::path& dir, const std::vector<int>& stepsToPlot) {
	int* histogramDataDevice;
	cudaMalloc(&histogramDataDevice, NUM_BINS * sizeof(int));
		
	float* energyBufferDevice;	
	cudaMalloc(&energyBufferDevice, sizeof(float) * simulation.box_host->boxparams.total_particles_upperbound);

	Compound* compoundsDevice;
	cudaMalloc(&compoundsDevice, sizeof(Compound) * simulation.box_host->boxparams.n_compounds);
	cudaMemcpy(compoundsDevice, simulation.box_host->compounds.data(), sizeof(Compound) * simulation.box_host->boxparams.n_compounds, cudaMemcpyHostToDevice);

	std::ofstream out_file(dir / "histogram_data.bin", std::ios::binary);
	int nPlots = stepsToPlot.size();
	out_file.write(reinterpret_cast<char*>(&nPlots), sizeof(int));
	for (int64_t step : stepsToPlot) {
		cudaMemcpy(energyBufferDevice, simulation.potE_buffer->GetBufferAtStep(step), sizeof(float) * simulation.box_host->boxparams.total_particles_upperbound, cudaMemcpyHostToDevice);
		cudaMemset(histogramDataDevice, 0, NUM_BINS * sizeof(int));

		cudaDeviceSynchronize();
		potEHistogramKernel << <simulation.box_host->boxparams.n_compounds, MAX_COMPOUND_PARTICLES >> > (compoundsDevice, simulation.box_host->boxparams.total_particles_upperbound, energyBufferDevice, histogramDataDevice, step);
		cudaDeviceSynchronize();

		std::vector<int> histogramDataHost;
		GenericCopyToHost(histogramDataDevice, histogramDataHost, NUM_BINS);

		std::vector<int64_t> bins = MakeBinLabels();
		
		out_file.write(reinterpret_cast<char*>(bins.data()), bins.size() * sizeof(int64_t));
		out_file.write(reinterpret_cast<char*>(histogramDataHost.data()), histogramDataHost.size() * sizeof(int));
	}
	out_file.close();

	cudaFree(energyBufferDevice);
	cudaFree(compoundsDevice);
}







int SimAnalysis::CountOscillations(std::vector<float>& data) {
	if (data.size() < 2) return 0;

	// Calculate mean value
	const float mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();

	int count = 0;
	bool aboveMean = data[0] > mean;

	for (size_t i = 1; i < data.size(); ++i) {
		bool currentAboveMean = data[i] > mean;
		if (!aboveMean && currentAboveMean) {
			++count;
		}
		aboveMean = currentAboveMean;
	}

	return count;
}





