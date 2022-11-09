#include "Analyzer.cuh"






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

/*
double __device__ calcKineticEnergy(Float3* pos1, Float3* pos2, double mass, double dt) {
	EngineUtils::applyHyperpos(pos1, pos2);


	double vel = (*pos1 - *pos2).len() * 0.5f / dt;
	double kinE = 0.5 * mass * vel * vel;
	return kinE;
}
*/


void __global__ monitorCompoundEnergyKernel(Box* box, Float3* traj_buffer, float* potE_buffer, Float3* data_out) {		// everything here breaks if not all compounds are identical in particle count and particle mass!!!!!!!
	__shared__ Float3 energy[MAX_COMPOUND_PARTICLES];
	__shared__ Compound compound;


	unsigned long long int step = blockIdx.x + 1;
	const int compound_index = blockIdx.y;
	energy[threadIdx.x] = Float3(0.f);








	if (threadIdx.x == 0) {
		//printf("index: %d\n", compound_index + (step - 1) * N_MONITORBLOCKS_PER_STEP);
		data_out[compound_index + (step - 1) * box->n_compounds] = Float3(0, 0, 0);
		//mass = box->compounds[compound_index].particles[0]
		compound = box->compounds[compound_index];
	}

	__syncthreads();

	if (threadIdx.x >= compound.n_particles) {
		return;
	}
	__syncthreads();


	float potE = potE_buffer[threadIdx.x + compound_index * MAX_COMPOUND_PARTICLES + step * box->total_particles_upperbound];

	Float3 pos_tsub1 = traj_buffer[threadIdx.x + compound_index * MAX_COMPOUND_PARTICLES + (step - 1) * box->total_particles_upperbound];
	Float3 pos_tadd1 = traj_buffer[threadIdx.x + compound_index * MAX_COMPOUND_PARTICLES + (step + 1) * box->total_particles_upperbound];
	//EngineUtils::applyHyperpos(&pos_tadd1, &pos_tsub1);
	//testspace::testerfn(4);
	//double kinE = calcKineticEnergy(&pos_tadd1, &pos_tsub1, compound.particles[threadIdx.x].mass, box->dt);
	float kinE = EngineUtils::calcKineticEnergy(&pos_tadd1, &pos_tsub1, SOLVENT_MASS, box->dt);



	double totalE = potE + kinE;


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
	int step = blockIdx.x + 1;																	// In cuda uint64_t isn't necessarily == llu
	int compounds_offset = box->n_compounds * MAX_COMPOUND_PARTICLES;


	energy[threadIdx.x] = Float3(0.f);
	if (threadIdx.x == 0) {
		data_out[(step - 1) * gridDim.y + blockIdx.y] = energy[0];
	}
	if (solvent_index >= box->n_solvents) { return; }


	
	Float3 pos_tsub1 = traj_buffer[compounds_offset + solvent_index + (step - 1) * box->total_particles_upperbound];
	Float3 pos_tadd1 = traj_buffer[compounds_offset + solvent_index + (step + 1) * box->total_particles_upperbound];

	float kinE = EngineUtils::calcKineticEnergy(&pos_tadd1, &pos_tsub1, SOLVENT_MASS, box->dt);

	float potE = potE_buffer[compounds_offset + solvent_index + step * box->total_particles_upperbound];

	if (potE != 0.f) {
		//printf("step %04d solvate %04d pot %f, compound_offset %d, step_offset  %d\n", step, solvent_index, potE, compounds_offset, step*box->total_particles_upperbound);
	}

	double totalE = potE + kinE;

	energy[threadIdx.x] = Float3(potE, kinE, totalE);
	__syncthreads();
	distributedSummation(energy, THREADS_PER_SOLVENTBLOCK);
	if (threadIdx.x == 0) {
		data_out[(step - 1) * gridDim.y + blockIdx.y] = energy[0];
	}
}



Analyzer::AnalyzedPackage Analyzer::analyzeEnergy(Simulation* simulation) {	// Calculates the avg J/mol // calculate energies separately for compounds and solvents. weigh averages based on amount of each
	int64_t analysable_steps = simulation->getStep() - 3;
	EngineUtils::genericErrorCheck("Cuda error before analyzeEnergy\n");
	if (analysable_steps < 1) {
		printf("FATAL ERROR, no steps to analyze");
		exit(1);
	}
	
	printf("Analyzer malloc %.4f GB on host\n", sizeof(Float3) * analysable_steps * 1e-9);
	Float3* average_energy = new Float3[analysable_steps];


	//uint64_t n_values = simulation->total_particles_upperbound * simulation->getStep();
	//printf("Analyzer malloc %f MB on device\n", (sizeof(Float3) + sizeof(double)) * n_values * 1e-6);

	/*
	//cudaMalloc(&traj_buffer_device, sizeof(Float3) * n_values);
	//cudaMalloc(&potE_buffer_device, sizeof(double) * n_values);
	//cudaMemcpy(traj_buffer_device, simulation->traj_buffer, sizeof(Float3) * n_values, cudaMemcpyHostToDevice);
	//cudaMemcpy(potE_buffer_device, simulation->potE_buffer, sizeof(double) * n_values, cudaMemcpyHostToDevice);


	Float3* average_solvent_energy = analyzeSolvateEnergy(simulation, analysable_steps);
	Float3* average_compound_energy = analyzeCompoundEnergy(simulation, analysable_steps);	//avg energy PER PARTICLE in compound

	for (uint64_t i = 0; i < analysable_steps; i++) {
		average_energy[i] = (average_solvent_energy[i] * simulation->box->n_solvents * 1 + average_compound_energy[i] * simulation->total_compound_particles) * (1.f/ (simulation->total_particles));
		//average_energy[i].print('E');
	}
	*/

	// We need to split up the analyser into steps, as we cannot store all positions traj on device at once.
	uint64_t max_steps_per_kernel = 1000;
	uint64_t particles_per_step = simulation->total_particles_upperbound;
	uint64_t max_values_per_kernel = (max_steps_per_kernel +2) * particles_per_step;							// Pad steps with 2 for vel calculation
	printf("Analyzer malloc %.2f MB on device\n", (sizeof(Float3) + sizeof(double)) * (max_values_per_kernel) * 1e-6);
	cudaMalloc(&traj_buffer_device, sizeof(Float3) * max_values_per_kernel);
	cudaMalloc(&potE_buffer_device, sizeof(float) * max_values_per_kernel);

	for (int i = 0; i < ceil((double)analysable_steps / (double)max_steps_per_kernel); i++) {
		uint64_t step_offset = i * max_steps_per_kernel + 1;												// offset one since we can't analyse step 1
		uint64_t steps_in_kernel = min(max_steps_per_kernel, analysable_steps - step_offset);
		uint64_t values_in_kernel = (steps_in_kernel + 2) * particles_per_step;

		cudaMemcpy(traj_buffer_device, &simulation->traj_buffer[(step_offset-1) * particles_per_step], sizeof(Float3) * values_in_kernel, cudaMemcpyHostToDevice);		// -1 to additional datap for vel calculation. +2 for -1 and -1 for vel calculation
		cudaMemcpy(potE_buffer_device, &simulation->potE_buffer[(step_offset-1) * particles_per_step], sizeof(float) * values_in_kernel, cudaMemcpyHostToDevice);



		cudaDeviceSynchronize();
		EngineUtils::genericErrorCheck("Cuda error during analyzer transfer\n");

		Float3* average_solvent_energy = analyzeSolvateEnergy(simulation, steps_in_kernel);
		Float3* average_compound_energy = analyzeCompoundEnergy(simulation, steps_in_kernel);	//avg energy PER PARTICLE in compound
		for (uint64_t ii = 0; ii < steps_in_kernel; ii++) {
			uint64_t step = step_offset - 1 + ii;	// -1 because index 0 is unused
			average_energy[step] = (average_solvent_energy[ii] * simulation->box->n_solvents + average_compound_energy[ii] * simulation->total_compound_particles) * (1.f / (simulation->total_particles));
			//printf("avg %f\n", average_energy[step].x);
			//average_energy[i].print('E');
		}
		delete[] average_solvent_energy, average_compound_energy;
	}


	cudaFree(traj_buffer_device);
	cudaFree(potE_buffer_device);


	printf("\n########## Finished analyzing energies ##########\n\n");
	return AnalyzedPackage(average_energy, analysable_steps, simulation->temperature_buffer, simulation->n_temp_values);;
}

Float3* Analyzer::analyzeSolvateEnergy(Simulation* simulation, uint64_t n_steps) {
	dim3 block_dim(n_steps, simulation->blocks_per_solventkernel, 1);
	Float3* average_solvent_energy = new Float3[n_steps];
	Float3* average_solvent_energy_blocked = new Float3[n_steps * simulation->blocks_per_solventkernel];
	Float3* data_out;
	cudaMalloc(&data_out, sizeof(Float3) * simulation->blocks_per_solventkernel * n_steps);

	for (int i = 0; i < n_steps; i++)
		average_solvent_energy[i] = Float3(0.f);

	
	monitorSolventEnergyKernel << < block_dim, THREADS_PER_SOLVENTBLOCK >> > (simulation->box, traj_buffer_device, potE_buffer_device, data_out);
	cudaDeviceSynchronize();
	cudaMemcpy(average_solvent_energy_blocked, data_out, sizeof(Float3) * simulation->blocks_per_solventkernel * n_steps, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();

	for (uint64_t step = 0; step < n_steps; step++) {
		average_solvent_energy[step] = Float3(0.f);
		for (int block = 0; block < simulation->blocks_per_solventkernel; block++) {
			average_solvent_energy[step] += average_solvent_energy_blocked[block + step * simulation->blocks_per_solventkernel];
		}
		average_solvent_energy[step] *= (1.f / simulation->n_solvents);
		//if (average_solvent_energy[step].x != 0.f)
		//printf("avg sol pot %f\n", average_solvent_energy[step].x);
	}


	cudaFree(data_out);
	delete[] average_solvent_energy_blocked;
	EngineUtils::genericErrorCheck("Cuda error during analyzeSolvateEnergy\n");
	return average_solvent_energy;
}


Float3* Analyzer::analyzeCompoundEnergy(Simulation* simulation, uint64_t n_steps) {
	uint64_t n_datapoints = simulation->n_compounds * n_steps;
	dim3 block_dim(n_steps, simulation->box->n_compounds, 1);
	Float3* average_compound_energy = new Float3[n_steps];

	for (int i = 0; i < n_steps; i++)
		average_compound_energy[i] = Float3(0.f);


	Float3* host_data = new Float3[n_datapoints];
	Float3* data_out;

	cudaMalloc(&data_out, sizeof(Float3) * n_datapoints * 2);
	

	



	monitorCompoundEnergyKernel << < block_dim, MAX_COMPOUND_PARTICLES >> > (simulation->box, traj_buffer_device, potE_buffer_device, data_out);
	cudaDeviceSynchronize();
	cudaMemcpy(host_data, data_out, sizeof(Float3) * n_datapoints, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();


	for (uint64_t step = 0; step < n_steps; step++) {
		for (uint64_t i = 0; i < simulation->box->n_compounds; i++) {
			//if (host_data[i + step * 256].x > 10000)
				//printf("Block: %d energy: %f\n", i, host_data[i + step * 256].x);
			average_compound_energy[step] += host_data[i + step * simulation->box->n_compounds];
		}
		average_compound_energy[step] *= (1.f / (simulation->total_compound_particles));
	}


	cudaFree(data_out);
	delete[] host_data;
	EngineUtils::genericErrorCheck("Cuda error during analyzeCompoundEnergy\n");

	return average_compound_energy;
}






