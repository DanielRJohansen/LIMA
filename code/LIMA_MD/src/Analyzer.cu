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



// everything here breaks if not all compounds are identical in particle count and particle mass!!!!!!!
// blockdim = (16,4,1)
void __global__ MonitorPclusterEnergy(const PersistentClusterMeta* const pcMeta, float* potE_buffer, float* vel_buffer, double3* data_out /*TODO: Its silly to store total, store just the float2 instead!*/, int nPclusters) {
	const int pcId = blockIdx.x * blockDim.x + threadIdx.x;
	const int64_t step = blockIdx.y;	// Step relative to current batch	
	const int pid = threadIdx.y;

	//if (pcId == 0 && pid == 0)


	if (pcId >= nPclusters)
		return;
	
	const int pidGlobal = pcMeta[pcId].particleIdsGlobal[pid];
	if (pidGlobal == -1)
		return;


	
	const float mass = pcMeta[pcId].mass[pid];

	const int64_t step_offset = step * nPclusters * PersistentCluster::maxParticles;
	const int64_t bufferIndex = pid + pcId * PersistentCluster::maxParticles + step_offset;

	const float potE = potE_buffer[bufferIndex];
	const float speed = vel_buffer[bufferIndex];
	const float kinE = PhysicsUtils::calcKineticEnergy(speed, mass);	// remove direction from vel

	const float totalE = potE + kinE;

	/*energy[particle_index] = Float3(potE, kinE, totalE);
	__syncthreads();

	LAL::distributedSummation(energy, MAX_COMPOUND_PARTICLES);
	__syncthreads();*/

	data_out[bufferIndex] = double3{ potE, kinE, totalE };
}










SimAnalysis::AnalyzedPackage SimAnalysis::analyzeEnergy(Simulation* simulation) {	// Calculates the avg J/mol // calculate energies separately for compounds and solvents. weigh averages based on amount of each
	LIMA_UTILS::genericErrorCheck("Cuda error before analyzeEnergy\n");

	const std::vector<PersistentClusterMeta>& pcMetaHost = simulation->box->persistentClustersMetadata;
	if (pcMetaHost.empty()) {
		return SimAnalysis::AnalyzedPackage{};
	}
	const int64_t n_entryindices = LIMALOGSYSTEM::getMostRecentDataentryIndex(simulation->getStep(), simulation->simParams.data_logging_interval);
	if (n_entryindices < 2) { return AnalyzedPackage(); }

	const int nParticlesUpperbound = pcMetaHost.size() * PersistentCluster::maxParticles;
	int64_t max_steps_per_kernel = 100;

	// First set up some stuff needed on device, that is currently on host
	PersistentClusterMeta* pcMetaDevice = GenericCopyToDevice(pcMetaHost);
	float* potE_buffer_device = nullptr;
	float* vel_buffer_device = nullptr;
	double3* energiesDev = nullptr;
	std::vector<double3> energiesHost(nParticlesUpperbound * max_steps_per_kernel, double3{});

	cudaMalloc(&potE_buffer_device, sizeof(float) * max_steps_per_kernel * nParticlesUpperbound);
	cudaMalloc(&vel_buffer_device, sizeof(float) * max_steps_per_kernel * nParticlesUpperbound);
	cudaMalloc(&energiesDev, sizeof(double3) * max_steps_per_kernel * nParticlesUpperbound);



	std::vector<Float3> average_energy(n_entryindices);


	// We need to split up the analyser into steps, as we cannot store all positions traj on device at once.
	for (int64_t i = 0; i < ceil((double)n_entryindices / (double)max_steps_per_kernel); i++) {
		const int64_t step_offset = i * max_steps_per_kernel;												// offset one since we can't analyse step 1
		const int64_t steps_in_kernel = std::min(max_steps_per_kernel, n_entryindices - step_offset);		

		cudaMemcpy(potE_buffer_device, &simulation->potE_buffer->data()[step_offset * nParticlesUpperbound], sizeof(float) * steps_in_kernel * nParticlesUpperbound, cudaMemcpyHostToDevice);
		cudaMemcpy(vel_buffer_device, &simulation->vel_buffer->data()[step_offset * nParticlesUpperbound], sizeof(float) * steps_in_kernel * nParticlesUpperbound, cudaMemcpyHostToDevice);
		cudaMemset(energiesDev, 0, sizeof(double3) * max_steps_per_kernel * nParticlesUpperbound);
		LIMA_UTILS::genericErrorCheck("Cuda error during analyzer transfer2\n");


		int pClustersPerCudablock = 16;
		int nCudablocks = (pcMetaHost.size() + pClustersPerCudablock-1) / pClustersPerCudablock;

		dim3 gridDim(
			nCudablocks,
			static_cast<uint32_t>(steps_in_kernel), 
			1
		);
		dim3 block_dim(
			pClustersPerCudablock,
			4, 
			1
		);

		MonitorPclusterEnergy<<< gridDim, block_dim>>> (pcMetaDevice, potE_buffer_device, vel_buffer_device, energiesDev, pcMetaHost.size());
		LIMA_UTILS::genericErrorCheck("Cuda error during MonitorPclusterEnergy\n");

		cudaMemcpy(energiesHost.data(), energiesDev, sizeof(double3) * steps_in_kernel * nParticlesUpperbound, cudaMemcpyDeviceToHost);


		for (uint64_t stepRelative = 0; stepRelative < steps_in_kernel; stepRelative++) {
			const int absStep = stepRelative + step_offset;

			double3 sum{};
			int cnt = 0;
			for (int pcid = 0; pcid < pcMetaHost.size(); pcid++) {
				for (int pid = 0; pid < PersistentCluster::maxParticles; pid++) {
					// TODO: Use precise particle count and std::accumulate here. Or reduce on gpu...
					if (pcMetaHost[pcid].particleIdsGlobal[pid] == -1) { continue; }

					const int64_t bufferIndex = pid + pcid * PersistentCluster::maxParticles + stepRelative * nParticlesUpperbound;
					sum.x += energiesHost[bufferIndex].x;
					sum.y += energiesHost[bufferIndex].y;
					sum.z += energiesHost[bufferIndex].z;
					cnt++;
				}
			}


			average_energy[absStep] = Float3(sum.x / static_cast<double>(cnt), sum.y / static_cast<double>(cnt), sum.z / static_cast<double>(cnt));
		}

		

		//std::vector<Float3> average_solvent_energy = analyzeSolvateEnergy(simulation, steps_in_kernel, potE_buffer_device, vel_buffer_device, tinymolForcefield_device, tinymols);
		//std::vector<Float3> average_compound_energy = analyzeCompoundEnergy(simulation, steps_in_kernel, potE_buffer_device, vel_buffer_device, compounds_device, forcefield_device);

		//for (int64_t ii = 0; ii < steps_in_kernel; ii++) {
		//	int64_t step = step_offset + ii - 1;	// -1 because index 0 is unused
		//	if (step == -1 || step >= n_entryindices - 2u) { continue; }	// Dont save first step, as the kinE is slightly wrong
		//	average_energy[step] = (average_solvent_energy[ii] + average_compound_energy[ii]);
		//}
	}

	cudaFree(potE_buffer_device);
	cudaFree(vel_buffer_device);
	cudaFree(energiesDev);
	cudaFree(pcMetaDevice);

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
//	for (auto entryindex = 0; entryindex < LIMALOGSYSTEM::getMostRecentDataentryIndex(sim.getStep()-1, sim.simParams.data_logging_interval); entryindex++) {
//
//		for (int compound_id = 0; compound_id < sim.box->boxparams.n_compounds; compound_id++) {
//			for (int particle_id = 0; particle_id < MAX_COMPOUND_PARTICLES; particle_id++) {
//				
//				const float potE = sim.potE_buffer->getCompoundparticleDatapointAtIndex(compound_id, particle_id, entryindex);
//
//				const uint8_t& atom_type = sim.box->compounds[compound_id].atom_types[particle_id];
//				const float mass = sim.forcefield.particle_parameters[atom_type].mass;
//				const float vel = sim.vel_buffer->getCompoundparticleDatapointAtIndex(compound_id, particle_id, entryindex);
//				const float kinE = PhysicsUtils::calcKineticEnergy(vel, mass);
//				
//				energies.emplace_back(potE);
//				energies.emplace_back(kinE);
//			}
//		}
//
//		for (int solvent_id = 0; solvent_id < sim.box->boxparams.n_solvents; solvent_id++) {
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

	//int64_t current_bin = 10;

	//while (bins.size() < NUM_BINS / 2) {
	//	bins.push_back(current_bin);
	//	current_bin = (current_bin == 0) ? 10 : current_bin * 10;
	//}

	//std::vector<double> negative_bins;
	//current_bin = -10;
	//while (negative_bins.size() < NUM_BINS / 2) {
	//	negative_bins.push_back(current_bin);
	//	current_bin *= 10;
	//}

	//std::reverse(negative_bins.begin(), negative_bins.end());
	//bins.insert(bins.begin(), negative_bins.begin(), negative_bins.end());

	return bins;
}

void SimAnalysis::PlotPotentialEnergyDistribution(const Simulation& simulation, const std::filesystem::path& dir, const std::vector<int>& stepsToPlot) {
	//int* histogramDataDevice;
	//cudaMalloc(&histogramDataDevice, NUM_BINS * sizeof(int));
	//	
	//float* energyBufferDevice;	
	//cudaMalloc(&energyBufferDevice, sizeof(float) * simulation.box->boxparams.total_particles_upperbound);

	//Compound* compoundsDevice;
	//cudaMalloc(&compoundsDevice, sizeof(Compound) * simulation.box->boxparams.n_compounds);
	//cudaMemcpy(compoundsDevice, simulation.box->compounds.data(), sizeof(Compound) * simulation.box->boxparams.n_compounds, cudaMemcpyHostToDevice);

	//std::ofstream out_file(dir / "histogram_data.bin", std::ios::binary);
	//int nPlots = stepsToPlot.size();
	//out_file.write(reinterpret_cast<char*>(&nPlots), sizeof(int));
	//for (int64_t step : stepsToPlot) {
	//	cudaMemcpy(energyBufferDevice, simulation.potE_buffer->GetBufferAtStep(step), sizeof(float) * simulation.box->boxparams.total_particles_upperbound, cudaMemcpyHostToDevice);
	//	cudaMemset(histogramDataDevice, 0, NUM_BINS * sizeof(int));

	//	cudaDeviceSynchronize();
	//	potEHistogramKernel << <simulation.box->boxparams.n_compounds, MAX_COMPOUND_PARTICLES >> > (compoundsDevice, simulation.box->boxparams.total_particles_upperbound, energyBufferDevice, histogramDataDevice, step);
	//	cudaDeviceSynchronize();

	//	std::vector<int> histogramDataHost;
	//	GenericCopyToHost(histogramDataDevice, histogramDataHost, NUM_BINS);

	//	std::vector<int64_t> bins = MakeBinLabels();
	//	
	//	out_file.write(reinterpret_cast<char*>(bins.data()), bins.size() * sizeof(int64_t));
	//	out_file.write(reinterpret_cast<char*>(histogramDataHost.data()), histogramDataHost.size() * sizeof(int));
	//}
	//out_file.close();

	//cudaFree(energyBufferDevice);
	//cudaFree(compoundsDevice);
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





