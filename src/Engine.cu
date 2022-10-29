#include "Engine.cuh"

#include "EngineUtils.cu"


__constant__ ForceField forcefield_device;
//__constant__ ForceField forcefield_nb_device;


Engine::Engine() {}
Engine::Engine(Simulation* simulation, ForceField forcefield_host) {
	EngineUtils::genericErrorCheck("Error before engine initialization.\n");
	this->simulation = simulation;
	//nlist_data_collection = new NListDataCollection(simulation);
	nlist_manager = new NListManager(simulation);



	int Ckernel_shared_mem = sizeof(Compound) + sizeof(CompoundState) + sizeof(NeighborList) + sizeof(Float3) * NEIGHBORLIST_MAX_SOLVENTS + sizeof(uint8_t) * NEIGHBORLIST_MAX_SOLVENTS;	
	int Skernel_shared_mem = sizeof(Float3) * MAX_COMPOUND_PARTICLES + sizeof(uint8_t) * MAX_COMPOUND_PARTICLES + sizeof(Solvent) * THREADS_PER_SOLVENTBLOCK;
	printf("Compoundkernel shared mem. size: %d B\n", Ckernel_shared_mem);
	printf("Solventkernel shared mem. size: %d B\n", Skernel_shared_mem);

	EngineUtils::genericErrorCheck("Error before moving forcefield to device\n");



	//ForceField forcefield_host = FFM.getForcefield();
	//ForceField forcefield_host = FFM.getNBForcefield();
	this->forcefield_host = forcefield_host;
	cudaMemcpyToSymbol(forcefield_device, &forcefield_host, sizeof(ForceField), 0, cudaMemcpyHostToDevice);	// So there should not be a & before the device __constant__

	printf("Forcefield size: %d bytes\n", sizeof(ForceField));

	//ForceField forcefield_nb_host = FFM.getNBForcefield();
	//cudaMemcpyToSymbol(forcefield_nb_device, &forcefield_nb_host, sizeof(ForceField), 0, cudaMemcpyHostToDevice);
	//cudaDeviceSynchronize();
	EngineUtils::genericErrorCheck("Error while moving forcefield to device\n");

	handleNLISTS(simulation, false, true);				// Fix neighborlists before running

	printf("Engine ready\n\n\n");
}





void Engine::deviceMaster() {
	EngineUtils::genericErrorCheck("Error before step!");
	step();
	EngineUtils::genericErrorCheck("Error after step!");
}

void Engine::hostMaster() {						// This is and MUST ALWAYS be called after the deviceMaster, and AFTER intStep()!
	auto t0 = std::chrono::high_resolution_clock::now();
	if ((simulation->getStep() % STEPS_PER_LOGTRANSFER) == 0) {
		offloadLoggingData();
		//offloadPositionData();

		if ((simulation->getStep() % STEPS_PER_THERMOSTAT) == 0 && ENABLE_BOXTEMP) {
			handleBoxtemp();
		}
	}
	if ((simulation->getStep() % STEPS_PER_TRAINDATATRANSFER) == 0) {
		offloadTrainData();
	}

	
	handleNLISTS(simulation);
	

	if ((simulation->getStep() % STEPS_PER_THERMOSTAT) == 1) {	// So this runs 1 step AFTER handleBoxtemp
		simulation->box->thermostat_scalar = 1.f;
	}

	auto t1 = std::chrono::high_resolution_clock::now();

	int cpu_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
	timings = timings + Int3(0,0,cpu_duration);
}









//--------------------------------------------------------------------------	CPU workload --------------------------------------------------------------//


void Engine::handleNLISTS(Simulation* simulation, bool async, bool force_update) {

	if (((nlist_manager->stepsSinceUpdate(simulation->getStep() ) >= STEPS_PER_NLIST_UPDATE) || force_update) && !updatenlists_mutexlock) {
		updatenlists_mutexlock = 1;
		
		nlist_manager->offloadPositionDataNLIST(simulation);
		// Lots of waiting time here...
		cudaDeviceSynchronize();

		nlist_manager->updateNeighborLists(simulation, &updatenlists_mutexlock, force_update, async, &timings.z);
	}

	if (nlist_manager->updated_neighborlists_ready) {
		nlist_manager->pushNlistsToDevice(simulation);
	}
}


void Engine::offloadLoggingData() {
	uint64_t step_offset = (simulation->getStep() - STEPS_PER_LOGTRANSFER) ;	// Tongue in cheek here, i think this is correct...

	cudaMemcpy(&simulation->potE_buffer[step_offset * simulation->total_particles_upperbound], simulation->box->potE_buffer, sizeof(float) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER, cudaMemcpyDeviceToHost);
	
	cudaMemcpy(&simulation->traj_buffer[step_offset * simulation->total_particles_upperbound], simulation->box->traj_buffer, sizeof(Float3) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER, cudaMemcpyDeviceToHost);

	cudaMemcpy(&simulation->logging_data[step_offset * 10], simulation->box->outdata, sizeof(float) * 10 * STEPS_PER_LOGTRANSFER, cudaMemcpyDeviceToHost);
}

void Engine::offloadPositionData() {
	//uint64_t step_offset = (simulation->getStep() - STEPS_PER_LOGTRANSFER) * simulation->total_particles_upperbound;	// Tongue in cheek here, i think this is correct...
	//cudaMemcpy(&simulation->traj_buffer[step_offset], simulation->box->traj_buffer, sizeof(Float3) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER, cudaMemcpyDeviceToHost);
}

void Engine::offloadTrainData() {
	uint64_t values_per_step = N_DATAGAN_VALUES * MAX_COMPOUND_PARTICLES * simulation->n_compounds;
	uint64_t step_offset = (simulation->getStep() - STEPS_PER_TRAINDATATRANSFER) * values_per_step;	// fix max_compound to the actual count save LOTS of space!. Might need a file in simout that specifies cnt for loading in other programs...
	cudaMemcpy(&simulation->traindata_buffer[step_offset], simulation->box->data_GAN, sizeof(Float3) * values_per_step * STEPS_PER_TRAINDATATRANSFER, cudaMemcpyDeviceToHost);
	EngineUtils::genericErrorCheck("Cuda error during traindata offloading\n");
}



																																			// THIS fn requires mallocmanaged!!   // HARD DISABLED HERE
void Engine::handleBoxtemp() {
	const float target_temp = 310.f;				// [k]
	Float3 temp_package = EngineUtils::getBoxTemperature(simulation, forcefield_host);
	float temp = temp_package.x;
	float biggest_contribution = temp_package.y;

	simulation->temperature_buffer[simulation->n_temp_values++] = temp;

	// So we avoid dividing by 0
	float temp_safe = temp == 0.f ? 1 : temp;
	float temp_scalar = target_temp / temp;
		// I just added this to not change any temperatures too rapidly
	temp_scalar = min(temp_scalar, 1.01f);
	temp_scalar = max(temp_scalar, 0.99f);

	if (PRINT_TEMP || temp > 500.f || temp < 100.f) { printf("\n %d Temperature: %.1f Biggest contrib: %.0f avg kinE %.0f\n", (simulation->getStep() - 1) / STEPS_PER_THERMOSTAT, temp, biggest_contribution, temp_package.z); }
		
	if (temp > target_temp/4.f && temp < target_temp*4.f || true) {
		if (APPLY_THERMOSTAT && simulation->getStep() > 10) {
			
			simulation->box->thermostat_scalar = temp_scalar;

			if (temp_scalar != temp_scalar ){//} || abs(temp_scalar) == "inf") {
				printf("Scalar: %f\n", simulation->box->thermostat_scalar);
				exit(0);
			}			
		}		
	}
	else {
		printf("Critical temperature encountered (%0.02f [k])\n", temp);
		simulation->box->critical_error_encountered = true;
	}
}



//--------------------------------------------------------------------------	SIMULATION BEGINS HERE --------------------------------------------------------------//
void Engine::step() {
	auto t0 = std::chrono::high_resolution_clock::now();
	cudaDeviceSynchronize();

	if (simulation->box->bridge_bundle->n_bridges > 0) {																		// TODO: Illegal access to device mem!!
		compoundBridgeKernel << < simulation->box->bridge_bundle->n_bridges, MAX_PARTICLES_IN_BRIDGE >> > (simulation->box);	// Must come before compoundKernel()		// DANGER
	}
		
	cudaDeviceSynchronize();
	if (simulation->n_compounds > 0) {
		compoundKernel << < simulation->box->n_compounds, THREADS_PER_COMPOUNDBLOCK >> > (simulation->box);
	}
	cudaDeviceSynchronize();

#ifdef ENABLE_SOLVENTS
	if (simulation->n_solvents > 0) { 
		solventForceKernel << < simulation->blocks_per_solventkernel, THREADS_PER_SOLVENTBLOCK >> > (simulation->box); 
	}
#endif
	cudaDeviceSynchronize();

	auto t1 = std::chrono::high_resolution_clock::now();


	
	CompoundState* temp = simulation->box->compound_state_array;
	simulation->box->compound_state_array = simulation->box->compound_state_array_next;
	simulation->box->compound_state_array_next = temp;
	
	
	Solvent* temp_s = simulation->box->solvents;
	simulation->box->solvents = simulation->box->solvents_next;
	simulation->box->solvents_next = temp_s;
	
	
	//cudaMemcpy(simulation->box->compound_state_array, simulation->box->compound_state_array_next, sizeof(CompoundState) * MAX_COMPOUNDS, cudaMemcpyDeviceToDevice);	// Update all positions, after all forces have been calculated
	
	
	
	//cudaMemcpy(simulation->box->solvents, simulation->box->solvents_next, sizeof(Solvent) * MAX_SOLVENTS, cudaMemcpyDeviceToDevice);
	cudaDeviceSynchronize();
	EngineUtils::genericErrorCheck("Error during step or state_transfer\n");		// Temp, we want to do host stuff while waiting for async GPU operations...	// SLOW
	auto t2 = std::chrono::high_resolution_clock::now();


	simulation->incStep();


	int force_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
	int copy_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	timings = timings + Int3(force_duration, copy_duration, 0);
}





// ------------------------------------------------------------------------------------------- DEVICE FUNCTIONS -------------------------------------------------------------------------------------------//



__device__ void applyPBC(Float3* current_position) {	// Only changes position if position is outside of box;
	for (int dim = 0; dim < 3; dim++) {
		*current_position->placeAt(dim) += BOX_LEN * (current_position->at(dim) < 0.f);
		*current_position->placeAt(dim) -= BOX_LEN * (current_position->at(dim) > BOX_LEN);
	}
}



/*
constexpr double sigma = 0.3923f;										// nm, point at which d/dE or force is 0.
constexpr double epsilon = 0.5986 * 1'000.f;							// J/mol
*/

__device__ inline float calcSigma(uint8_t atomtype1, uint8_t atomtype2) {
	return (forcefield_device.particle_parameters[atomtype1].sigma + forcefield_device.particle_parameters[atomtype2].sigma) * 0.5;
}
__device__ inline float calcEpsilon(uint8_t atomtype1, uint8_t atomtype2) {
	return sqrtf(forcefield_device.particle_parameters[atomtype1].epsilon * forcefield_device.particle_parameters[atomtype2].epsilon);
}

__device__ Float3 calcLJForce(Float3* pos0, Float3* pos1, float* data_ptr, float* potE, float sigma, float epsilon, int type1=-1, int type2=-1) {
	// Calculates LJ force on p0	(attractive to p1. Negative values = repulsion )//
	// input positions in cartesian coordinates [nm]
	// sigma [nm]
	// epsilon [J/mol]->[(kg*nm^2)/(ns^2*mol)]
	//	Returns force in J/mol*M		?????????????!?!?//


	//return Float3(sigma);	// Test for forcing LJ calc, but always small force, even with bonded particles!


	// Directly from book
	float dist_sq = (*pos1 - *pos0).lenSquared();
	float s = sigma * sigma / dist_sq;								// [nm^2]/[nm^2] -> unitless
	s = s * s * s;
	float force_scalar = 24.f * epsilon * s / dist_sq * (1.f - 2.f * s);	// Attractive. Negative, when repulsive		[(kg*nm^2)/(nm^2*ns^2*mol)] ->----------------------	[(kg)/(ns^2*mol)]	

#ifdef LIMA_VERBOSE
	Float3 ddd = (*pos1 - *pos0) * force_scalar;
	if (ddd.x != ddd.x) {
		printf("\nError is here\n");
		printf("Scalar %f dist_sq %f\n", force_scalar, dist_sq);
		pos0->print('0');
		pos0->print('1');
	}
	//if (threadIdx.x == 32 && blockIdx.x == 28 && force < -600000.f && force > -700000) {
	//if (abs(force) > 20e+8 && type1 == 200 && type2 == 1) {
	if (abs(force_scalar) > 20e+9 || *potE != *potE ) {
		//printf("\nDist %f   D2 %f    force %.1f [MN]     sigma: %.3f \tepsilon %.0f  t1 %d t2 %d\n", (float)sqrt(dist_sq), (float) dist_sq, (float)force_scalar *1e-6, sigma, epsilon, type1, type2);
		//pos0->print('1');
		//pos1->print('2');
		//printf("Block %d Thread %d\n", blockIdx.x, threadIdx.x);
		//printf("Thread %d Block %d self %f %f %f other %f %f %f\n", threadIdx.x, blockIdx.x, pos0->x, pos0->y, pos0->z, pos1->x, pos1->y, pos1->z);
		//printf("Force %f %f %f\n", force_unit_vector.x * force, force_unit_vector.y * force, force_unit_vector.z * force);
	}

	{
		//Float3 force_vec = force_unit_vector * force;
		/*if (force_vec.x != force_vec.x) {
			//printf("Force: %f\n", force);
			//force_unit_vector.print('u');
			printf("Thread %d Block %d self %f %f %f other %f %f %f\n", threadIdx.x, blockIdx.x, pos0->x, pos0->y, pos0->z, pos1->x, pos1->y, pos1->z);

		}*/
	}
#endif
	return (*pos1 - *pos0) * force_scalar;										// GN/mol [(kg*nm)/(ns^2*mol)]
}

//constexpr double kb = 17.5 * 1e+6;		//	J/(mol*nm^2)	/ kg/(ns^2 * mol)
__device__ void calcPairbondForces(Float3* pos_a, Float3* pos_b, PairBond* bondtype, Float3* results, float* potE) {
	// Calculates bond force on both particles					//
	// Calculates forces as J/mol*M								//

	Float3 difference = *pos_a - *pos_b;						//	[nm]
	double error = difference.len() - bondtype->b0;			//	[nm]

	*potE += 0.5 * bondtype->kb * (error * error);					// [J/mol]

	difference = difference.norm();								// dif_unit_vec, but shares register with dif
	double force_scalar = -bondtype->kb * error;				//	[J/(mol*nm^2)]*nm =>	kg*nm^2*ns^-2/(mol*nm^2)*nm = kg*nm/(mol*ns^2)				 [N/mol] directionless, yes?

	results[0] = difference * force_scalar;				// [GN]
	results[1] = difference * force_scalar * -1;		// [GN]

#ifdef LIMA_VERBOSE
	if (force_scalar > 2e+7) {
		printf("thread %d  error %f ref %f force %f\n", threadIdx.x, error, bondtype->b0, force_scalar);
		//pos_a->print('a');
		//pos_b->print('b');
	}
#endif
		
}


//constexpr double ktheta = 65 * 1e+3;	// J/mol
__device__ void calcAnglebondForces(Float3* pos_left, Float3* pos_middle, Float3* pos_right, AngleBond* angletype, Float3* results, float* potE) {
	Float3 v1 = *pos_left - *pos_middle;
	Float3 v2 = *pos_right - *pos_middle;
	Float3 normal1 = v1.cross(v2);
	Float3 normal2 = v2.cross(v1);

	Float3 inward_force_direction1 = (v1.cross(normal2)).norm();
	Float3 inward_force_direction2 = (v2.cross(normal1)).norm();

	double angle = Float3::getAngle(v1, v2);
	double error = angle - angletype->theta_0;// *reference_angle;

	//error = (error / (2.f * PI)) * 360.f;

	//*potE += 0.5 * ktheta * error * error * 0.5;
	//double force_scalar = ktheta * (error);
	*potE += 0.5 * angletype->k_theta * error * error * 0.5;
	double force_scalar = angletype->k_theta * (error);

	results[0] = inward_force_direction1 * force_scalar;
	results[2] = inward_force_direction2 * force_scalar;
	results[1] = (results[0] + results[2]) * -1;


	//printf("\nangle %f error %f force %f t0 %f kt %f\n", angle, error, force_scalar, angletype->theta_0, angletype->k_theta);

}
__device__ void calcDihedralbondForces(Float3* pos_left, Float3* pos_lm, Float3* pos_rm, Float3* pos_right, DihedralBond* dihedral, Float3* results, float* potE) {
	Float3 normal1 = (*pos_left-*pos_lm).cross((*pos_rm-*pos_lm)).norm();		// Should this not be normalized????
	Float3 normal2 = (*pos_lm-*pos_rm).cross((*pos_right-*pos_rm)).norm();			// Is inward or outward? Refactor!!!
	// Both vectors point "left". If pos_left+n1 is closer to pos_right than pos_left-n1, then n1 is pointing inwards, and n2 outwards. Also vice-versa.
	float torsion = Float3::getAngle(normal1, normal2);




	if (((*pos_left + normal1) - *pos_right).len() < ((*pos_left - normal1) - *pos_right).len())	// Wait, can't i just check if torsion > pi | torsion < 0`???????????????????????
		normal2 *= -1;
	else
		normal1 *= -1;
	// Now both normals are pointing inwards

	//printf("Torsion %f\n", normal1.len());
	//float error = (torsion - dihedral->phi_0);	// *360.f / 2.f * PI;	// Check whether k_phi is in radians or degrees

	//error = (error / (2.f * PI)) * 360.f;

	
	//*potE += 0.5 * dihedral->k_phi * error * error * 0.5;
	//double force_scalar = dihedral->k_phi * (error);

	//double force_scalar = sinf(torsion - dihedral->phi_0);		// Should be -sinf? (cos)' = -sin??!?!

	float force_scalar = dihedral->k_phi * sinf(dihedral->n * torsion - dihedral->phi_0);
	*potE = dihedral->k_phi * (1 - cosf(dihedral->n * torsion - dihedral->phi_0));

	//printf("Torsion %f ref %f force_scalar %f\n", torsion, dihedral->phi_0, force_scalar);
	//force_scalar *= dihedral->k_phi;
	if (abs(force_scalar) > 183300) {
		pos_left->print('L');
		pos_lm->print('l');
		pos_rm->print('r');
		pos_right->print('R');
		//printf("torsion %f      ref %f     error %f     force: %f\n", torsion, dihedral->phi_0, error, force_scalar);
		printf("torsion %f      ref %f     error %f     force: %f\n", torsion, dihedral->phi_0, 0.f, force_scalar);
	}
		



	//results[0] = normal1 * force_scalar * -1.f;
	//results[3] = normal2 * force_scalar;
	results[0] = normal1 * force_scalar;
	results[3] = normal2 * force_scalar;
	// Not sure about the final two forces, for now we'll jsut split the sum of opposite forces between them.
	results[1] = (results[0] + results[3]) * -1.f * 0.5;
	results[2] = (results[0] + results[3]) * -1.f * 0.5;


}




//__device__ Float3 computeIntermolecularLJForces(Float3* self_pos, int n_particles, Float3* positions, float* data_ptr, float* potE_sum, uint8_t atomtype_self, uint8_t* atomtypes_others, LJ_Ignores* lj_ignore_list, uint32_t* global_ids) {	// Assumes all positions are 
__device__ Float3 computerIntermolecularLJForces(Float3* self_pos, uint8_t atomtype_self, float* potE_sum, uint32_t global_id_self, float* data_ptr,
	Compound* neighbor_compound, Float3* neighbor_positions, int neighborcompound_id, BondedParticlesLUT& bonded_particles_lut) {	
	Float3 force(0.f);

	for (int neighborparticle_id = 0; neighborparticle_id < neighbor_compound->n_particles; neighborparticle_id++) {

		// If thread's assoc. particle is bonded to the particle in neighborcompound, continue
		if (*bonded_particles_lut.get(threadIdx.x, neighborparticle_id)) { continue; }

		const int neighborparticle_atomtype = neighbor_compound->atom_types[neighborparticle_id];	//// TEMPORARY, this is waaaay to many global mem accesses

		force += calcLJForce(self_pos, &neighbor_positions[neighborparticle_id], data_ptr, potE_sum,
			calcSigma(atomtype_self, neighborparticle_atomtype), calcEpsilon(atomtype_self, neighborparticle_atomtype),
			global_id_self, neighbor_compound->particle_global_ids[neighborparticle_id]
		);
	}	
	return force;// *24.f * 1e-9;
}

__device__ Float3 computeIntramolecularLJForces(Compound* compound, CompoundState* compound_state, float* potE_sum, float* data_ptr, BondedParticlesLUT& bonded_particles_lut) {
	Float3 force(0.f);
	for (int i = 0; i < compound->n_particles; i++) {

		if (i != threadIdx.x && !(* bonded_particles_lut.get(threadIdx.x, i))) {

			force += calcLJForce(&compound_state->positions[threadIdx.x], &compound_state->positions[i], data_ptr, potE_sum,
				calcSigma(compound->atom_types[threadIdx.x], compound->atom_types[i]), calcEpsilon(compound->atom_types[threadIdx.x], compound->atom_types[i]),
				compound->particle_global_ids[threadIdx.x], compound->particle_global_ids[i]
			);// *24.f * 1e-9;
		}
	}
	return force;
}





__device__ Float3 computeSolventToSolventLJForces(Float3* self_pos, NeighborList* nlist, Solvent* solvents, float* data_ptr, float* potE_sum) {	// Specific to solvent kernel
	Float3 force(0.f);
	for (int i = 0; i < nlist->n_solvent_neighbors; i++) {
		Solvent neighbor = solvents[nlist->neighborsolvent_ids[i]];			
		EngineUtils::applyHyperpos(&neighbor.pos, self_pos);
		force += calcLJForce(self_pos, &neighbor.pos, data_ptr, potE_sum, 
			forcefield_device.particle_parameters[ATOMTYPE_SOL].sigma, 
			forcefield_device.particle_parameters[ATOMTYPE_SOL].epsilon);
		/*if (abs(force.len()) > 1000000) {
			(*self_pos - neighbor.pos).print('d');
			printf("F %f %f %f n %d  solvent_id %d neighbor_id %d\n", force.x, force.y, force.z, nlist->n_solvent_neighbors, threadIdx.x + blockIdx.x*blockDim.x, nlist->neighborsolvent_ids[i]);
		}*/
	}
	return force;// *24.f * 1e-9;
}
__device__ Float3 computeSolventToCompoundLJForces(Float3* self_pos, int n_particles, Float3* positions, float* data_ptr, float* potE_sum, uint8_t atomtype_self) {	// Specific to solvent kernel
	Float3 force(0, 0, 0);
	for (int i = 0; i < n_particles; i++) {
		force += calcLJForce(self_pos, &positions[i], data_ptr, potE_sum,
			calcSigma(atomtype_self, ATOMTYPE_SOL), calcEpsilon(atomtype_self, ATOMTYPE_SOL)
			//(forcefield_device.particle_parameters[atomtype_self].sigma + forcefield_device.particle_parameters[ATOMTYPE_SOL].sigma) * 0.5f,
			//sqrtf(forcefield_device.particle_parameters[atomtype_self].epsilon * forcefield_device.particle_parameters[ATOMTYPE_SOL].epsilon)
		);
	}
	return force;// *24.f * 1e-9;
}
__device__ Float3 computeCompoundToSolventLJForces(Float3* self_pos, int n_particles, Float3* positions, float* data_ptr, float* potE_sum, uint8_t atomtype_self, uint8_t* atomtypes_others) {	// Assumes all positions are 
	Float3 force(0.f);
	for (int i = 0; i < n_particles; i++) {
		force += calcLJForce(self_pos, &positions[i], data_ptr, potE_sum,
			(forcefield_device.particle_parameters[atomtype_self].sigma + forcefield_device.particle_parameters[atomtypes_others[i]].sigma) * 0.5f,
			(forcefield_device.particle_parameters[atomtype_self].epsilon + forcefield_device.particle_parameters[atomtypes_others[i]].epsilon) * 0.5,
			atomtype_self, atomtypes_others[i]
		);
	}
	return force;// *24.f * 1e-9;
}

/*
__device__ Float3 computeCompoundToSolventLJForces(Float3* self_pos, int n_particles, Float3* positions, double* data_ptr, double* potE_sum) {	// Specific to solvent kernel
	Float3 force(0, 0, 0);
	for (int i = 0; i < n_particles; i++) {
		Float3 hyperpos = positions[i];			// copy, DONT ref it as all threads will cann applyHyperpos
		EngineUtils::applyHyperpos(self_pos, &hyperpos);
		force += calcLJForce(self_pos, &hyperpos, data_ptr, potE_sum);	
	}
	return force;
}
*/
template <typename T>	// Can either be Compound or CompoundBridgeCompact
__device__ Float3 computePairbondForces(T* entity, Float3* positions, Float3* utility_buffer, float* potE) {	// only works if n threads >= n bonds
//__device__ Float3 computePairbondForces(CompoundBridgeCompact* entity, Float3* positions, Float3* utility_buffer, float* potE) {	// only works if n threads >= n bonds
	utility_buffer[threadIdx.x] = Float3(0.f);
	for (int bond_offset = 0; (bond_offset * blockDim.x) < entity->n_singlebonds; bond_offset++) {
		PairBond* pb = nullptr;
		Float3 forces[2];
		int bond_index = threadIdx.x + bond_offset * blockDim.x;

		if (bond_index < entity->n_singlebonds) {
			pb = &entity->singlebonds[bond_index];

			calcPairbondForces(
				&positions[pb->atom_indexes[0]],
				&positions[pb->atom_indexes[1]],
				//&pb->reference_dist,
				pb,
				forces, potE
			);

			//forces[0] = Float3(0.f); forces[1] = Float3(0.f);
#ifndef INWD
			if (pb->invertLJ) {
				float temp = 0.f;
				//float sigma = (forcefield_device.particle_parameters[entity->atom_types[pb->atom_indexes[0]]].sigma + forcefield_device.particle_parameters[entity->atom_types[pb->atom_indexes[1]]].sigma) * 0.5;
				
				float epsilon = sqrtf(forcefield_device.particle_parameters[entity->atom_types[pb->atom_indexes[0]]].epsilon * forcefield_device.particle_parameters[entity->atom_types[pb->atom_indexes[1]]].epsilon);
				//Float3 anti_lj_force = calcLJForce(&positions[pb->atom_indexes[0]], &positions[pb->atom_indexes[1]], &temp, potE, sigma, epsilon, -1,-1);
				Float3 anti_lj_force = calcLJForce(&positions[pb->atom_indexes[0]], &positions[pb->atom_indexes[1]], &temp, potE, 
					calcSigma(entity->atom_types[pb->atom_indexes[0]], entity->atom_types[pb->atom_indexes[1]]), 
					calcEpsilon(entity->atom_types[pb->atom_indexes[0]], entity->atom_types[pb->atom_indexes[1]]),
					-1, -1);
				forces[0] -= anti_lj_force; forces[1] += anti_lj_force;
				if (blockIdx.x == 0 && (pb->atom_indexes[0] == 1 || pb->atom_indexes[1] == 1)) {
					//printf("ANTI LJ %d %d\n", pb->atom_indexes[0], pb->atom_indexes[1]);
					//anti_lj_force.print('A');
				}

				//if (entity->particle_refs[pb->atom_indexes[0]].global_id == 649 || entity->particle_refs[pb->atom_indexes[1]].global_id == 649) {
					//anti_lj_force.print('B');
					//printf("B dist %.08f\n", (positions[pb->atom_indexes[0]] - positions[pb->atom_indexes[1]]).len());
				//}
					
				
			}
#endif
		}


		for (int i = 0; i < blockDim.x; i++) {
			if (threadIdx.x == i && pb != nullptr) {
				utility_buffer[pb->atom_indexes[0]] += forces[0];
				utility_buffer[pb->atom_indexes[1]] += forces[1];
			}
			__syncthreads();
		}
	}
	return utility_buffer[threadIdx.x];
}

template <typename T>	// Can either be Compound or CompoundBridgeCompact
__device__ Float3 computeAnglebondForces(T* entity, Float3* positions, Float3* utility_buffer, float* potE) {
	utility_buffer[threadIdx.x] = Float3(0, 0, 0);
	for (int bond_offset = 0; (bond_offset * blockDim.x) < entity->n_anglebonds; bond_offset++) {
		AngleBond* ab = nullptr;
		Float3 forces[3];
		int bond_index = threadIdx.x + bond_offset * blockDim.x;

		if (bond_index < entity->n_anglebonds) {
			ab = &entity->anglebonds[bond_index];


			calcAnglebondForces(
				&positions[ab->atom_indexes[0]],
				&positions[ab->atom_indexes[1]],
				&positions[ab->atom_indexes[2]],
				//&ab->reference_angle,
				ab,
				forces, potE
			);

			//forces[0] = Float3(0.f); forces[1] = Float3(0.f); forces[2] = Float3(0.f);
#ifndef INWD
			if (ab->invertLJ) {
				float temp = 0.f;
				//float sigma = (forcefield_device.particle_parameters[entity->atom_types[ab->atom_indexes[0]]].sigma + forcefield_device.particle_parameters[entity->atom_types[ab->atom_indexes[2]]].sigma) * 0.5;
				//float epsilon = sqrtf(forcefield_device.particle_parameters[entity->atom_types[ab->atom_indexes[0]]].epsilon * forcefield_device.particle_parameters[entity->atom_types[ab->atom_indexes[2]]].epsilon);
				//Float3 anti_lj_force = calcLJForce(&positions[ab->atom_indexes[0]], &positions[ab->atom_indexes[2]], &temp, potE, sigma, epsilon, -1, -1);
				Float3 anti_lj_force = calcLJForce(&positions[ab->atom_indexes[0]], &positions[ab->atom_indexes[2]], &temp, potE, 
					calcSigma(entity->atom_types[ab->atom_indexes[0]], entity->atom_types[ab->atom_indexes[2]]), 
					calcEpsilon(entity->atom_types[ab->atom_indexes[0]], entity->atom_types[ab->atom_indexes[2]]),
					-1, -1);
				forces[0] -= anti_lj_force; forces[2] += anti_lj_force;
				if (blockIdx.x == 0 && (ab->atom_indexes[0] == 1 || ab->atom_indexes[2] == 1)) {
					//printf("ANTI LJ %d %d %d\n", ab->atom_indexes[0], ab->atom_indexes[1], ab->atom_indexes[2]);
					//anti_lj_force.print('A');
				}
			}
#endif
		}


		for (int i = 0; i < blockDim.x; i++) {
			if (threadIdx.x == i && ab != nullptr) {
				utility_buffer[ab->atom_indexes[0]] += forces[0];
				utility_buffer[ab->atom_indexes[1]] += forces[1];
				utility_buffer[ab->atom_indexes[2]] += forces[2];
			}
			__syncthreads();
		}
	}

	return utility_buffer[threadIdx.x];
}

template <typename T>	// Can either be Compound or CompoundBridgeCompact
__device__ Float3 computeDihedralForces(T* entity, Float3* positions, Float3* utility_buffer, float* potE) {
	utility_buffer[threadIdx.x] = Float3(0, 0, 0);
	for (int bond_offset = 0; (bond_offset * blockDim.x) < entity->n_dihedrals; bond_offset++) {
		DihedralBond* db = nullptr;
		Float3 forces[4];
		int bond_index = threadIdx.x + bond_offset * blockDim.x;

		if (bond_index < entity->n_dihedrals) {
			db = &entity->dihedrals[bond_index];
			//printf("Firing %d of %d\n", bond_index, entity);
			calcDihedralbondForces(
				&positions[db->atom_indexes[0]],
				&positions[db->atom_indexes[1]],
				&positions[db->atom_indexes[2]],
				&positions[db->atom_indexes[3]],
				db,
				forces, 
				potE
			);

#ifndef INWD
			//forces[0] = Float3(0.f); forces[1] = Float3(0.f); forces[2] = Float3(0.f); forces[3] = Float3(0.f);
			if (db->invertLJ) {
				float temp = 0.f;
				//float sigma = (forcefield_device.particle_parameters[entity->atom_types[db->atom_indexes[0]]].sigma + forcefield_device.particle_parameters[entity->atom_types[db->atom_indexes[3]]].sigma) * 0.5;
				//float epsilon = sqrtf(forcefield_device.particle_parameters[entity->atom_types[db->atom_indexes[0]]].epsilon * forcefield_device.particle_parameters[entity->atom_types[db->atom_indexes[3]]].epsilon);
				//Float3 anti_lj_force = calcLJForce(&positions[db->atom_indexes[0]], &positions[db->atom_indexes[3]], &temp, potE, sigma, epsilon, db->atom_indexes[0], db->atom_indexes[3]);
				Float3 anti_lj_force = calcLJForce(&positions[db->atom_indexes[0]], &positions[db->atom_indexes[3]], &temp, potE, 
					calcSigma(entity->atom_types[db->atom_indexes[0]], entity->atom_types[db->atom_indexes[3]]), 
					calcEpsilon(entity->atom_types[db->atom_indexes[0]], entity->atom_types[db->atom_indexes[3]]),
					db->atom_indexes[0], db->atom_indexes[3]);

				//Float3 anti_lj_force = calcLJForce(&positions[db->atom_indexes[0]], &positions[db->atom_indexes[3]], &temp, potE, sigma, epsilon, );
				forces[0] -= anti_lj_force; forces[3] += anti_lj_force;

			}
#endif

		}

		for (int i = 0; i < blockDim.x; i++) {
			if (threadIdx.x == i && db != nullptr) {
				utility_buffer[db->atom_indexes[0]] += forces[0];
				utility_buffer[db->atom_indexes[1]] += forces[1];
				utility_buffer[db->atom_indexes[2]] += forces[2];
				utility_buffer[db->atom_indexes[3]] += forces[3];
			}
			__syncthreads();
		}
	}
	//if (utility_buffer[threadIdx.x].x != utility_buffer[threadIdx.x].x)
		//utility_buffer[threadIdx.x].print('d');
	return utility_buffer[threadIdx.x];
}

__device__ void integratePosition(Float3* pos, Float3* pos_tsub1, Float3* force, const double mass, const double dt, float* thermostat_scalar, int p_index, bool issolvent =false) {
	// Force is in ??Newton, [kg * nm /(mol*ns^2)] //

	float prev_vel = (*pos - *pos_tsub1).len();
		
	Float3 temp = *pos;
	EngineUtils::applyHyperpos(pos, pos_tsub1);
	//printf("dt %.08f\n", dt);
	//*pos = *pos * 2 - *pos_tsub1 + *force * (1.f / mass) * dt * dt;	// no *0.5?	// [nm] - [nm] + [kg*m*s^-2]/[kg] * [ns]^2
	//double magnitude_equalizer = 1e+4;		// Due to force being very large and dt being very small
	//double masstime_scaling = ((dt*magnitude_equalizer)*(dt * magnitude_equalizer)) / mass;

	//*pos = *pos * 2 - *pos_tsub1 + *force * masstime_scaling;		// [nm] - [nm] + [kg/mol*m*/s^2]/[kg/mol] * [s]^2 * (1e-9)^2	=> [nm]-[nm]+[]
	//pos = *pos * 2 - *pos_tsub1 + (*force * (1./(magnitude_equalizer * magnitude_equalizer))) * masstime_scaling;		// [nm] - [nm] + [kg/mol*m * /s^2]/[kg/mol] * [s]^2 * (1e-9)^2	=> [nm]-[nm]+[]
	*pos = *pos * 2 - *pos_tsub1 + *force * (dt / mass) * dt;		// [nm] - [nm] + [kg/mol*m*/s ^ 2] / [kg / mol] * [s] ^ 2 * (1e-9) ^ 2 = > [nm] - [nm] + []
	*pos_tsub1 = temp;

	
	return;
	Float3 delta_pos = *pos - *pos_tsub1;
	*pos = *pos_tsub1 + delta_pos * *thermostat_scalar;

	if (delta_pos.len() > 0.05)
		printf("\nSol: %d b %d t %d.      Distance/step: %f prev: %f    Force: %f\n", issolvent, blockIdx.x, threadIdx.x, delta_pos.len(), prev_vel, force->len());
#ifdef LIMA_VERBOSE
	if ((*pos-*pos_tsub1).len() > 0.1) {
		printf("\nP_index %d Thread %d blockId %d\tForce %f mass  %f \Dist %f\n", p_index, threadIdx.x, blockIdx.x, force->len(), mass, (*pos - *pos_tsub1).len());
		//printf("\nP_index %d Thread %d blockId %d\tForce %f mass  %f \tFrom %f %f %f\tTo %f %f %f\n", p_index, threadIdx.x, blockIdx.x, force->len(), mass, pos_tsub1->x, pos_tsub1->y, pos_tsub1->z, pos->x, pos->y, pos->z);		
	}
#endif
}

__device__ void integratePositionRampUp(Float3* pos, Float3* pos_tsub1, Float3* force, const double mass, const double dt, float rampup_scalar) {
	// Force is in ??Newton, [kg * nm /(mol*ns^2)] //

	//float prev_vel = (*pos - *pos_tsub1).len();

	Float3 temp = *pos;
	EngineUtils::applyHyperpos(pos, pos_tsub1);
	//printf("dt %.08f\n", dt);
	//*pos = *pos * 2 - *pos_tsub1 + *force * (1.f / mass) * dt * dt;	// no *0.5?	// [nm] - [nm] + [kg*m*s^-2]/[kg] * [ns]^2
	//double magnitude_equalizer = 1e+4;		// Due to force being very large and dt being very small
	//double masstime_scaling = ((dt * magnitude_equalizer) * (dt * magnitude_equalizer)) / mass;

	//Float3 vel_scaled = (*pos - *pos_tsub1) * (1. / rampup_scalar);
	//Float3 acc_scaled = (*force * (1. / (magnitude_equalizer * magnitude_equalizer))) * masstime_scaling * (1./rampup_scalar);

	//*pos = *pos + vel_scaled + acc_scaled;
	//*pos = *pos * 2 - *pos_tsub1 + *force * masstime_scaling;		// [nm] - [nm] + [kg/mol*m*/s^2]/[kg/mol] * [s]^2 * (1e-9)^2	=> [nm]-[nm]+[]
	//*pos = *pos * 2 - *pos_tsub1 + (*force * (1. / (magnitude_equalizer * magnitude_equalizer))) * masstime_scaling;		// [nm] - [nm] + [kg/mol*m*/s^2]/[kg/mol] * [s]^2 * (1e-9)^2	=> [nm]-[nm]+[]


	*pos = *pos * 2. - *pos_tsub1 + *force * (dt / mass) * dt * 0.f;		// [nm] - [nm] + [kg/mol*m*/s ^ 2] / [kg / mol] * [s] ^ 2 * (1e-9) ^ 2 = > [nm] - [nm] + []
	*pos_tsub1 = temp;

	return;

	float vel_scalar = log2f(force->len()/1000.);
	vel_scalar = max(vel_scalar, 1.f);
//	rampup_scalar = min(rampup_scalar, 2.f);
	Float3 delta_pos = *pos - *pos_tsub1;
	*pos = *pos_tsub1 + delta_pos * (1. / vel_scalar);
}




// ------------------------------------------------------------------------------------------- KERNELS -------------------------------------------------------------------------------------------//




#define compound_index blockIdx.x
__global__ void compoundKernel(Box* box) {
	__shared__ Compound compound;
	__shared__ CompoundState compound_state;
	__shared__ NeighborList neighborlist;
	__shared__ BondedParticlesLUT bonded_particles_lut;
	__shared__ Float3 utility_buffer[THREADS_PER_COMPOUNDBLOCK];


	return;

	if (threadIdx.x == 0) {
		compound.loadMeta(&box->compounds[blockIdx.x]);
		compound_state.setMeta(compound.n_particles);
		neighborlist.loadMeta(&box->compound_neighborlists[blockIdx.x]);

	}


	__syncthreads();
	compound.loadData(&box->compounds[blockIdx.x]);
	compound_state.loadData(&box->compound_state_array[blockIdx.x]);
	neighborlist.loadData(&box->compound_neighborlists[blockIdx.x]);
	__syncthreads();


	float potE_sum = 0;
	float data_ptr[4];
	for (int i = 0; i < 4; i++)
		data_ptr[i] = 0;
	data_ptr[2] = 9999;
	data_ptr[3] = box->step + 1;
	

	Float3 force = compound.forces[threadIdx.x];
	
	//force = Float3(0.f);
	Float3 force_LJ_sol(0.f);
	


	// ------------------------------------------------------------ Intramolecular Operations ------------------------------------------------------------ //
	{
		bonded_particles_lut.load(*box->bonded_particles_lut_manager->get(compound_index, compound_index));
		EngineUtils::applyHyperpos(&compound_state.positions[0], &compound_state.positions[threadIdx.x]);
		__syncthreads();

		force += computePairbondForces(&compound, compound_state.positions, utility_buffer, &potE_sum);
		force += computeAnglebondForces(&compound, compound_state.positions, utility_buffer, &potE_sum);
		force += computeDihedralForces(&compound, compound_state.positions, utility_buffer, &potE_sum);
		force += computeIntramolecularLJForces(&compound, &compound_state, &potE_sum, data_ptr, bonded_particles_lut);
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------- //
	//force = Float3(0.f);
	// --------------------------------------------------------------- Intermolecular forces --------------------------------------------------------------- //	
	for (int i = 0; i < neighborlist.n_compound_neighbors; i++) {
		int neighborcompound_id = neighborlist.neighborcompound_ids[i];
		if (neighborcompound_id == blockIdx.x) {
			printf("i: %d  Block %d neighbors nc %d of %d\n",i, blockIdx.x, neighborcompound_id, neighborlist.n_compound_neighbors);
		}

		int neighborcompound_particles = box->compound_state_array[neighborcompound_id].n_particles;

		if (threadIdx.x < neighborcompound_particles) {
			utility_buffer[threadIdx.x] = box->compound_state_array[neighborcompound_id].positions[threadIdx.x];
			//utility_buffer_small[threadIdx.x] = box->compounds[neighborcompound_id].atom_types[threadIdx.x];				/// HEEEEY WHY ARE WE NOT USING THIS?!??!?!?!?!
			EngineUtils::applyHyperpos(&compound_state.positions[0], &utility_buffer[threadIdx.x]);
		}

		BondedParticlesLUT* compoundpair_lut = box->bonded_particles_lut_manager->get(compound_index, neighborcompound_id);
		bonded_particles_lut.load(*compoundpair_lut);
		__syncthreads();

		if (threadIdx.x < compound.n_particles) {
			force += computerIntermolecularLJForces(&compound_state.positions[threadIdx.x], compound.atom_types[threadIdx.x], &potE_sum, compound.particle_global_ids[threadIdx.x], data_ptr,
			//anti_inter += computerIntermolecularLJForces(&compound_state.positions[threadIdx.x], compound.atom_types[threadIdx.x], &compound.lj_ignore_list[threadIdx.x], &potE_sum, compound.particle_global_ids[threadIdx.x], data_ptr,
				&box->compounds[neighborcompound_id], utility_buffer, neighborcompound_id, bonded_particles_lut);
		}
		__syncthreads();				// CRITICAL			
	}
	//force += anti_inter;
	//force = Float3(0.f);
	// ------------------------------------------------------------------------------------------------------------------------------------------------------ //




	// --------------------------------------------------------------- Solvation forces --------------------------------------------------------------- //
#ifdef ENABLE_SOLVENTS
	*/
	for (int offset = 0; offset * blockDim.x < neighborlist.n_solvent_neighbors; offset += blockDim.x) {
		int solvent_nlist_index = offset + threadIdx.x; // index in neighborlist

		if (solvent_nlist_index < neighborlist.n_solvent_neighbors) {
			utility_buffer[threadIdx.x] = box->solvents[neighborlist.neighborsolvent_ids[solvent_nlist_index]].pos;
			EngineUtils::applyHyperpos(&compound_state.positions[0], &utility_buffer[threadIdx.x]);
		}
		__syncthreads();

		if (threadIdx.x < compound.n_particles) {
			force_LJ_sol += computeSolventToCompoundLJForces(&compound_state.positions[threadIdx.x], blockDim.x, utility_buffer, data_ptr, &potE_sum, compound.atom_types[threadIdx.x]);			
		}
	}
	force += force_LJ_sol;
#endif
	// ------------------------------------------------------------------------------------------------------------------------------------------------ //


	//box->potE_buffer[threadIdx.x + blockIdx.x * MAX_COMPOUND_PARTICLES * N_DATAGAN_VALUES + (box->step % STEPS_PER_TRAINDATATRANSFER) * gridDim.x * MAX_COMPOUND_PARTICLES * N_DATAGAN_VALUES] = force.len();
	box->traj_buffer[threadIdx.x + blockIdx.x * MAX_COMPOUND_PARTICLES * N_DATAGAN_VALUES + (box->step % STEPS_PER_TRAINDATATRANSFER) * gridDim.x * MAX_COMPOUND_PARTICLES * N_DATAGAN_VALUES] = Float3(0.);// compound_state.positions[threadIdx.x];	// N DATAGAN???!?!


	// ------------------------------------------------------------ Integration ------------------------------------------------------------ //
	if (threadIdx.x < compound.n_particles) {
		//Float3 force_(force.x, force.y, force.z);
		//integratePosition(&compound_state.positions[threadIdx.x], &compound.prev_positions[threadIdx.x], &force_, forcefield_device.particle_parameters[compound.atom_types[threadIdx.x]].mass, box->dt, &box->thermostat_scalar, threadIdx.x, false);		
		if (box->step >= RAMPUP_STEPS || INTEGRATION_RAMPUP_ENABLED) {
			integratePosition(&compound_state.positions[threadIdx.x], &compound.prev_positions[threadIdx.x], &force, forcefield_device.particle_parameters[compound.atom_types[threadIdx.x]].mass, box->dt, &box->thermostat_scalar, threadIdx.x, false);
		}
		else {
			integratePositionRampUp(&compound_state.positions[threadIdx.x], &compound.prev_positions[threadIdx.x], &force, forcefield_device.particle_parameters[compound.atom_types[threadIdx.x]].mass, box->dt, RAMPUP_STEPS-box->step);
		}
		
		box->compounds[blockIdx.x].prev_positions[threadIdx.x] = compound.prev_positions[threadIdx.x];
	}
	__syncthreads();
	// ------------------------------------------------------------------------------------------------------------------------------------- //





	// ------------------------------------ PERIODIC BOUNDARY CONDITION ------------------------------------------------------------------------------------------------- // 
	if (threadIdx.x == 0) {
		applyPBC(&compound_state.positions[threadIdx.x]);
	}
	EngineUtils::applyHyperpos(&compound_state.positions[0], &compound_state.positions[threadIdx.x]);	// So all particles follows p0
	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------ //
	





	

	
	// ------------------------------------ DATA LOG ------------------------------- //	
	{
		if (threadIdx.x < compound.n_particles || 1) {																							// TODO: Remove || 1
			int step_offset = (box->step % STEPS_PER_LOGTRANSFER) * box->total_particles_upperbound;
			int compound_offset = blockIdx.x * MAX_COMPOUND_PARTICLES;
			//box->potE_buffer[threadIdx.x + compound_offset + step_offset] = potE_sum;			// TODO: Should be += since bridge sets value first

		}		
		__syncthreads();


		if (blockIdx.x == 0 && threadIdx.x == 0) {			
			Float3 pos_prev_temp = compound.prev_positions[threadIdx.x];			
			EngineUtils::applyHyperpos(&compound_state.positions[threadIdx.x], &pos_prev_temp);
			
			int step_offset = (box->step % STEPS_PER_LOGTRANSFER) * 10;

			box->outdata[0 + step_offset] = (compound_state.positions[threadIdx.x] - pos_prev_temp).len() / box->dt;
			box->outdata[1 + step_offset] = EngineUtils::calcKineticEnergy(&compound_state.positions[threadIdx.x], &pos_prev_temp, forcefield_device.particle_parameters[compound.atom_types[threadIdx.x]].mass, box->dt);
			box->outdata[2 + step_offset] = potE_sum;																											// This does NOT take bridge potE into account!!!!!!!
			box->outdata[3 + step_offset] = force.len();

			//box->outdata[5 + box->step * 10] = data_ptr[2];// closest particle
			//box->outdata[6 + box->step * 10] = data_ptr[1];// force.len();
		}


		int n_compounds_total = gridDim.x;
		int step_offset = (box->step % STEPS_PER_TRAINDATATRANSFER) * n_compounds_total * MAX_COMPOUND_PARTICLES * N_DATAGAN_VALUES;
		int compound_offset = blockIdx.x * MAX_COMPOUND_PARTICLES * N_DATAGAN_VALUES;
		int particle_offset = threadIdx.x * N_DATAGAN_VALUES;
		box->data_GAN[0 + particle_offset + compound_offset + step_offset] = compound_state.positions[threadIdx.x];
		box->data_GAN[1 + particle_offset + compound_offset + step_offset] = force_LJ_sol;
		box->data_GAN[2 + particle_offset + compound_offset + step_offset] = force;

		if (threadIdx.x >= compound.n_particles)
			box->data_GAN[0 + particle_offset + compound_offset + step_offset] = Float3(-1.f);
//		box->data_GAN[1 + threadIdx.x * 6 + step_offset] = force_bond + force_angle;
//		box->data_GAN[2 + threadIdx.x * 6 + step_offset] = force_LJ_com;
//		box->data_GAN[3 + threadIdx.x * 6 + step_offset] = force_LJ_sol;
		//box->outdata[0 + threadIdx.x * 10]
	}
	// ----------------------------------------------------------------------------- //

	if (force.len() > 2e+10) {
		printf("\n\nCritical force %.0f           block %d thread %d\n\n\n", force.len(), blockIdx.x, threadIdx.x);
		box->critical_error_encountered = true;
	}
		
	
	box->compound_state_array_next[blockIdx.x].loadData(&compound_state);
}
#undef compound_id










__global__ void solventForceKernel(Box* box) {
#define solvent_index (blockIdx.x * blockDim.x + threadIdx.x)
#define thread_active (solvent_index < box->n_solvents)
	//__shared__ Float3 solvent_positions[THREADS_PER_SOLVENTBLOCK];
	__shared__ Float3 utility_buffer[MAX_COMPOUND_PARTICLES];
	__shared__ uint8_t utility_buffer_small[MAX_COMPOUND_PARTICLES];
	__shared__ Solvent solvents[THREADS_PER_SOLVENTBLOCK];

	float potE_sum = 0;
	float data_ptr[4];	// Pot, force, closest particle, ?
	for (int i = 0; i < 4; i++)
		data_ptr[i] = 0;
	data_ptr[2] = 9999.f;
	Float3 force(0,0,0);

	Float3 fff;


	//Solvent* solvent = &solvents[threadIdx.x];
	Solvent solvent;
	//Float3 solvent_pos;
	if (thread_active) {
		solvent = box->solvents[solvent_index];
	}	

	// --------------------------------------------------------------- Molecule Interactions --------------------------------------------------------------- //
	for (int i = 0; i < box->n_compounds; i++) {
		//continue;		 // DANGER
		int n_compound_particles = box->compound_state_array[i].n_particles;
		// First all threads help loading the molecule
		if (threadIdx.x < n_compound_particles) {
			utility_buffer[threadIdx.x] = box->compound_state_array[i].positions[threadIdx.x];
			utility_buffer_small[threadIdx.x] = box->compounds[i].atom_types[threadIdx.x];
		}
		__syncthreads();
			

		if (thread_active) {
			EngineUtils::applyHyperpos(&utility_buffer[0], &solvent.pos);									// Move own particle in relation to compound-key-position
			force += computeCompoundToSolventLJForces(&solvent.pos, n_compound_particles, utility_buffer, data_ptr, &potE_sum, ATOMTYPE_SOL, utility_buffer_small);
		}
		__syncthreads();


	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //



	// --------------------------------------------------------------- Solvent Interactions --------------------------------------------------------------- //
	if (thread_active) {
		// DANGER
		force += computeSolventToSolventLJForces(&solvent.pos, &box->solvent_neighborlists[solvent_index], box->solvents, data_ptr, &potE_sum);
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //


	if (solvent_index < box->n_solvents)
		box->traj_buffer[box->n_compounds * MAX_COMPOUND_PARTICLES + solvent_index + (box->step % STEPS_PER_LOGTRANSFER) * box->total_particles_upperbound] = Float3(0.);
	//__syncthreads();
	// ------------------------------------ DATA LOG ------------------------------- //
	if (thread_active) {
		int compounds_offset = box->n_compounds * MAX_COMPOUND_PARTICLES;
		int step_offset = (box->step % STEPS_PER_LOGTRANSFER) * box->total_particles_upperbound;

		//box->potE_buffer[compounds_offset + solvent_index + step_offset] = (double) potE_sum;	//  data_ptr[0];
		//box->potE_buffer[compounds_offset + solvent_index + step_offset] = 0;
		//printf("pot: %f\n", box->potE_buffer[compounds_offset + solvent_index + (box->step) * box->total_particles]);
		box->potE_buffer[compounds_offset + solvent_index + step_offset] = force.len();
		box->traj_buffer[compounds_offset + solvent_index + step_offset] = solvent.pos;
		/*if ((solvent.pos.x > 7 || solvent.pos.y > 7 || solvent.pos.z > 7) && box->step == 0)
			solvent.pos.print('y');*/
		box->traj_buffer[compounds_offset + solvent_index + step_offset] = solvent.pos;
	}

	// ----------------------------------------------------------------------------- //

	if (thread_active) {
//		printf("%f\n", (solvent.pos - box->solvents[solvent_index].pos).len());

		int p_index = MAX_COMPOUND_PARTICLES + solvent_index;
		if (box->step >= RAMPUP_STEPS || !INTEGRATION_RAMPUP_ENABLED) {
			integratePosition(&solvent.pos, &solvent.pos_tsub1, &force, forcefield_device.particle_parameters[ATOMTYPE_SOL].mass, box->dt, &box->thermostat_scalar, p_index, true);
		}
		else {
			integratePositionRampUp(&solvent.pos, &solvent.pos_tsub1, &force, forcefield_device.particle_parameters[ATOMTYPE_SOL].mass, box->dt, RAMPUP_STEPS - box->step);

		}
		float len = (solvent.pos - solvent.pos_tsub1).len();


		applyPBC(&solvent.pos);	


		
	}



	if (thread_active) {
		if (solvent.pos.x != solvent.pos.x) {
			solvent.pos.print('s');
			box->critical_error_encountered = true;
		}
		box->solvents_next[solvent_index] = solvent;
	}




#undef solvent_index
#undef thread_active
}






__global__ void compoundBridgeKernel(Box* box) {
//#define compound_id blockIdx.x
#define particle_id_bridge threadIdx.x
//#define particle_active particle_id_bridge < bridge.n_particles
	__shared__ CompoundBridgeCompact bridge;
	__shared__ Float3 positions[MAX_PARTICLES_IN_BRIDGE];
	__shared__ Float3 utility_buffer[MAX_PARTICLES_IN_BRIDGE];							// waaaaay too biggg
	__shared__ uint8_t utility_buffer_small[MAX_PARTICLES_IN_BRIDGE];

	if (threadIdx.x == 0) {
		bridge.loadMeta(&box->bridge_bundle->compound_bridges[blockIdx.x]);
	}
	__syncthreads();
	bridge.loadData(&box->bridge_bundle->compound_bridges[blockIdx.x]);
	positions[particle_id_bridge] = Float3(0.f);
	if (particle_id_bridge < bridge.n_particles) {
		ParticleRefCompact* p_ref = &bridge.particle_refs[particle_id_bridge];
		positions[particle_id_bridge] = box->compound_state_array[p_ref->compound_id].positions[p_ref->local_id_compound];
	}
	__syncthreads();



	float potE_sum = 0;

	Float3 force(0.f);

	// ------------------------------------------------------------ Intercompund Operations ------------------------------------------------------------ //
	{											// So for the very first step, these ´should all be 0, but they are not??										TODO: Look into this at some point!!!! Also, can cant really apply hyperpos here without breaking stuff, sindce
																																								// The bridge spans such a big volume! :/
		EngineUtils::applyHyperpos(&positions[0], &positions[particle_id_bridge]);
		__syncthreads();	
		force += computePairbondForces(&bridge, positions, utility_buffer, &potE_sum);
		force += computeAnglebondForces(&bridge, positions, utility_buffer, &potE_sum);
		force += computeDihedralForces(&bridge, positions, utility_buffer, &potE_sum);
	}

	__syncthreads();
	// --------------------------------------------------------------------------------------------------------------------------------------------------- //

	if (particle_id_bridge < bridge.n_particles) {
		ParticleRefCompact* p_ref = &bridge.particle_refs[particle_id_bridge];
		box->compounds[p_ref->compound_id].forces[p_ref->local_id_compound] = force;



		int compound_offset = p_ref->compound_id * MAX_COMPOUND_PARTICLES;
		int step_offset = (box->step % STEPS_PER_LOGTRANSFER) * box->total_particles_upperbound;

		box->potE_buffer[p_ref->local_id_compound + compound_offset + step_offset] = potE_sum;
		//force.print(char((int)threadIdx.x+48));
	}


}



















/*	VELOCITY VERLET STORMER
__device__ void integratePosition(CompactParticle* particle, Float3* particle_pos, Float3* particle_force, double* dt) {
	*particle_pos = *particle_pos + (particle->vel + *particle_force * (0.5 / particle->mass) * *dt) * *dt;
}
__device__ void integrateVelocity(CompactParticle* particle, Float3* particle_force, double* dt) {
	particle->vel = particle->vel + (*particle_force + particle->force_prev) * (0.5 / particle->mass) * *dt;
}
*/