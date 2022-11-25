#include "Engine.cuh"


// Pre-calculate a solvent-X paired forcefield, to save ALOT of calc in kernels
__constant__ ForceField_NB forcefield_device;




void Engine::setDeviceConstantMemory() {
	cudaMemcpyToSymbol(forcefield_device, &forcefield_host, sizeof(ForceField_NB), 0, cudaMemcpyHostToDevice);	// So there should not be a & before the device __constant__
	cudaDeviceSynchronize();
	EngineUtils::genericErrorCheck("Error while moving forcefield to device\n");
}




// ------------------------------------------------------------------------------------------- DEVICE FUNCTIONS -------------------------------------------------------------------------------------------//

__device__ inline float calcSigma(uint8_t atomtype1, uint8_t atomtype2) {
	return (forcefield_device.particle_parameters[atomtype1].sigma + forcefield_device.particle_parameters[atomtype2].sigma) * 0.5;
}
__device__ inline float calcEpsilon(uint8_t atomtype1, uint8_t atomtype2) {
	return sqrtf(forcefield_device.particle_parameters[atomtype1].epsilon * forcefield_device.particle_parameters[atomtype2].epsilon);
}

__device__ Float3 calcLJForce(Float3* pos0, Float3* pos1, float* data_ptr, float* potE, float sigma, float epsilon, int type1 = -1, int type2 = -1) {
	// Calculates LJ force on p0	(attractive to p1. Negative values = repulsion )//
	// input positions in cartesian coordinates [nm]
	// sigma [nm]
	// epsilon [J/mol]->[(kg*nm^2)/(ns^2*mol)]
	// Returns force in J/mol*M		?????????????!?!?//

	// Directly from book
	float dist_sq = (*pos1 - *pos0).lenSquared();
	float s = sigma * sigma / dist_sq;								// [nm^2]/[nm^2] -> unitless
	s = s * s * s;
	float force_scalar = 24.f * epsilon * s / dist_sq * (1.f - 2.f * s);	// Attractive. Negative, when repulsive		[(kg*nm^2)/(nm^2*ns^2*mol)] ->----------------------	[(kg)/(ns^2*mol)]	

	*potE += 4. * epsilon * s * (s - 1.f);

	//printf("eps %f sig %f s6 %f dist %f dist2 %f\n", epsilon, sigma, s, (*pos1 - *pos0).len(), dist_sq);
#ifdef LIMA_VERBOSE
	Float3 ddd = (*pos1 - *pos0) * force_scalar;
	if (ddd.x != ddd.x) {
		printf("\nError is here\n");
		printf("Scalar %f dist_sq %f\n", force_scalar, dist_sq);
		pos0->print('0');
		pos0->print('1');
	}
#endif
	return (*pos1 - *pos0) * force_scalar;										// GN/mol [(kg*nm)/(ns^2*mol)]
}


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
	Float3 normal1 = (*pos_left - *pos_lm).cross((*pos_rm - *pos_lm)).norm();		// Should this not be normalized????
	Float3 normal2 = (*pos_lm - *pos_rm).cross((*pos_right - *pos_rm)).norm();			// Is inward or outward? Refactor!!!
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


__device__ Float3 computerIntercompoundLJForces(Float3* self_pos, uint8_t atomtype_self, float* potE_sum, uint32_t global_id_self, float* data_ptr,
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

__device__ Float3 computeIntracompoundLJForces(Compound* compound, CompoundState* compound_state, float* potE_sum, float* data_ptr, BondedParticlesLUT* bonded_particles_lut) {
	Float3 force(0.f);
	if (threadIdx.x >= compound->n_particles) { return force; }
	for (int i = 0; i < compound->n_particles; i++) {
		if (i != threadIdx.x && !(*bonded_particles_lut->get(threadIdx.x, i))) {
			//*potE_sum = 16;
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
	//if (blockIdx.x + threadIdx.x == 0) printf("n neighbor solvents %d\n", nlist->n_solvent_neighbors);
	for (int i = 0; i < nlist->n_solvent_neighbors; i++) {
		Solvent neighbor = solvents[nlist->neighborsolvent_ids[i]];
		EngineUtils::applyHyperpos(&neighbor.pos, self_pos);
		force += calcLJForce(self_pos, &neighbor.pos, data_ptr, potE_sum,
			forcefield_device.particle_parameters[ATOMTYPE_SOL].sigma,
			forcefield_device.particle_parameters[ATOMTYPE_SOL].epsilon);
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

template <typename T>	// Can either be Compound or CompoundBridgeCompact
__device__ Float3 computePairbondForces(T* entity, Float3* positions, Float3* utility_buffer, float* potE) {	// only works if n threads >= n bonds
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
	utility_buffer[threadIdx.x] = Float3(0.f);

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

	return utility_buffer[threadIdx.x];
}

__device__ void integratePosition(Float3* pos, Float3* pos_tsub1, Float3* force, const float mass, const double dt, float* thermostat_scalar, int p_index, bool print = false) {
	// Force is in ??Newton, [kg * nm /(mol*ns^2)] //	



	// Always do hyperpos first!
	EngineUtils::applyHyperpos(pos, pos_tsub1);

	Float3 temp = *pos;
	float prev_vel = (*pos - *pos_tsub1).len();

	// [nm] - [nm] + [kg/mol*m*/s ^ 2] / [kg / mol] * [s] ^ 2 * (1e-9) ^ 2 = > [nm] - [nm] + []
	Float3 dx = *pos - *pos_tsub1;
	Float3 ddx = force->mul_highres(dt) * (dt / static_cast<double>(mass));
	Float3 x_ = (dx).add_highres(ddx);

	if (print) {
	//	printf("Thread %d    prev_vel %f    Force scalar %f    Acc %.10f    Change %.10f\n", threadIdx.x, prev_vel, force->len(), acc, (*pos - *pos_tsub1).len() - prev_vel);
		printf("x  %.8f  dx %.8f   ddx %.8f    x_ %.8f   dif %.8f\n", pos->len(), dx.len(), ddx.len(), x_.len(), ((*pos + x_) - *pos).len());
	}
	*pos += x_;
	*pos_tsub1 = temp;

	double d = force->len();
	double acc = d* dt / mass * dt;


	
	Float3 delta_pos = *pos - *pos_tsub1;
	*pos = *pos_tsub1 + delta_pos * *thermostat_scalar;

	
	{
		float mass2 = 2.f * mass;
		float mass_inv = 1.f / mass;
		float dtdivm = dt / static_cast<double>(mass);
		float fdt = force->mul_highres(dt).len();
		//printf("dt/mass %.8f   fdt %.8f ddx %.8f    x_ %f\n", dtdivm, fdt, ddx.len(), x_.len());
		//printf("mass %f dt^2/mass %.8f, dt %.8f ddx %.8f   \n", mass, dt / mass * dt, dt, ddx);
	}
	

	//if (delta_pos.len() > 0.05) {
	//	printf("\nSol: %d b %d t %d.      Distance/step: %f prev: %f    Force: %f\n", issolvent, blockIdx.x, threadIdx.x, delta_pos.len(), prev_vel, force->len());
	//}

#ifdef LIMA_VERBOSE
	if ((*pos - *pos_tsub1).len() > 0.1) {
		printf("\nP_index %d Thread %d blockId %d\tForce %f mass  %f \Dist %f\n", p_index, threadIdx.x, blockIdx.x, force->len(), mass, (*pos - *pos_tsub1).len());
		//printf("\nP_index %d Thread %d blockId %d\tForce %f mass  %f \tFrom %f %f %f\tTo %f %f %f\n", p_index, threadIdx.x, blockIdx.x, force->len(), mass, pos_tsub1->x, pos_tsub1->y, pos_tsub1->z, pos->x, pos->y, pos->z);		
	}
#endif
}

__device__ void integratePositionRampUp(Float3* pos, Float3* pos_tsub1, Float3* force, const double mass, const double dt, int step) {
	EngineUtils::applyHyperpos(pos, pos_tsub1);

	Float3 temp = *pos;
	float prev_vel = (*pos - *pos_tsub1).len();

	*pos = *pos * 2.f - *pos_tsub1 + *force * (dt / mass) * dt;		// [nm] - [nm] + [kg/mol*m*/s ^ 2] / [kg / mol] * [s] ^ 2 * (1e-9) ^ 2 = > [nm] - [nm] + []
	*pos_tsub1 = temp;

	Float3 delta_pos = *pos - *pos_tsub1;
	float dist = delta_pos.len();
	//float vel_scalar = std::min(1.f, prev_vel / new_vel);

	float rampup_progress = static_cast<float>(step + 1) / static_cast<float>(RAMPUP_STEPS);
	float vel_scalar = MAX_RAMPUP_DIST / dist;				// Scalar so dist of vector will be MAX_RAMPUP_DIST
	vel_scalar = vel_scalar <= 1.f ? vel_scalar : 1.f;		// Constrain scalar so we never accelerate any vector
	vel_scalar += (1.f - vel_scalar) * rampup_progress;		// Move the scalar closer to 1., as we get closer to finishing rampup
	*pos = *pos_tsub1 + delta_pos * vel_scalar;

	{
	// Ensure that particles cannot increase their momentum during rampup
	//Float3 delta_pos = *pos - *pos_tsub1;
	//float new_vel = delta_pos.len();
	////float vel_scalar = std::min(1.f, prev_vel / new_vel);
	//float vel_scalar = prev_vel / new_vel;
	//vel_scalar = vel_scalar <= 1.01f ? vel_scalar : 1.01f;
	//*pos = *pos_tsub1 + delta_pos * vel_scalar;
	}
}


__device__ inline void LogCompoundData(Compound& compound, Box* box, CompoundState& compound_state, float* potE_sum, Float3& force, Float3& force_LJ_sol) {
	const int compound_index = blockIdx.x;

	{
		const int steps_since_transfer = (box->step % STEPS_PER_LOGTRANSFER);
		const int step_offset = steps_since_transfer * box->total_particles_upperbound;
		const int compound_offset = compound_index * MAX_COMPOUND_PARTICLES;
		const int index = step_offset + compound_offset + threadIdx.x;
		// Log the previous pos, as this comes after integration, and we do want that juicy pos(t0) ;)
		box->traj_buffer[index] = compound.prev_positions[threadIdx.x];
		box->potE_buffer[index] = *potE_sum;
	}
	//if (blockIdx.x == 0 && threadIdx.x == 0) {
	//	Float3 pos_prev_temp = compound.prev_positions[threadIdx.x];
	//	EngineUtils::applyHyperpos(&compound_state.positions[threadIdx.x], &pos_prev_temp);
	//	const int step_offset = (box->step % STEPS_PER_LOGTRANSFER) * 10;
	//	box->outdata[0 + step_offset] = (compound_state.positions[threadIdx.x] - pos_prev_temp).len() / box->dt;
	//	box->outdata[1 + step_offset] = EngineUtils::calcKineticEnergy(&compound_state.positions[threadIdx.x], &pos_prev_temp, forcefield_device.particle_parameters[compound.atom_types[threadIdx.x]].mass, box->dt);
	//	box->outdata[2 + step_offset] = potE_sum;																											// This does NOT take bridge potE into account!!!!!!!
	//	box->outdata[3 + step_offset] = force.len();
	//}
	//const int step_offset = (box->step % STEPS_PER_TRAINDATATRANSFER) * n_compounds_total * MAX_COMPOUND_PARTICLES * N_DATAGAN_VALUES;
	//int compound_offset = blockIdx.x * MAX_COMPOUND_PARTICLES * N_DATAGAN_VALUES;
	//int particle_offset = threadIdx.x * N_DATAGAN_VALUES;
	//box->data_GAN[0 + particle_offset + compound_offset + step_offset] = compound_state.positions[threadIdx.x];
	//box->data_GAN[1 + particle_offset + compound_offset + step_offset] = force_LJ_sol;
	//box->data_GAN[2 + particle_offset + compound_offset + step_offset] = force;
	//if (threadIdx.x >= compound.n_particles)
	//	box->data_GAN[0 + particle_offset + compound_offset + step_offset] = Float3(-1.f);
	//		box->data_GAN[1 + threadIdx.x * 6 + step_offset] = force_bond + force_angle;
	//		box->data_GAN[2 + threadIdx.x * 6 + step_offset] = force_LJ_com;
	//		box->data_GAN[3 + threadIdx.x * 6 + step_offset] = force_LJ_sol;
			//box->outdata[0 + threadIdx.x * 10]
}

__device__ inline void LogSolventData(bool thread_active, Box* box, int solvent_index, Float3& force, Solvent& solvent) {
	if (thread_active) {
		int compounds_offset = box->n_compounds * MAX_COMPOUND_PARTICLES;
		int step_offset = (box->step % STEPS_PER_LOGTRANSFER) * box->total_particles_upperbound;

		box->potE_buffer[compounds_offset + solvent_index + step_offset] = force.len();
		box->traj_buffer[compounds_offset + solvent_index + step_offset] = solvent.pos;
	}
}


// ------------------------------------------------------------------------------------------- KERNELS -------------------------------------------------------------------------------------------//




#define compound_index blockIdx.x
__global__ void compoundKernel(Box* box) {
	__shared__ Compound compound;
	__shared__ CompoundState compound_state;
	__shared__ NeighborList neighborlist;
	__shared__ BondedParticlesLUT bonded_particles_lut;
	__shared__ Float3 utility_buffer[THREADS_PER_COMPOUNDBLOCK];


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
	Float3 force_LJ_sol(0.f);


	// ------------------------------------------------------------ Intramolecular Operations ------------------------------------------------------------ //
	{
		bonded_particles_lut.load(*box->bonded_particles_lut_manager->get(compound_index, compound_index));
		EngineUtils::applyHyperpos(&compound_state.positions[0], &compound_state.positions[threadIdx.x]);
		__syncthreads();

		force += computePairbondForces(&compound, compound_state.positions, utility_buffer, &potE_sum);
		force += computeAnglebondForces(&compound, compound_state.positions, utility_buffer, &potE_sum);
		force += computeDihedralForces(&compound, compound_state.positions, utility_buffer, &potE_sum);
		force += computeIntracompoundLJForces(&compound, &compound_state, &potE_sum, data_ptr, &bonded_particles_lut);
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------- //

	// --------------------------------------------------------------- Intermolecular forces --------------------------------------------------------------- //	
	for (int i = 0; i < neighborlist.n_compound_neighbors; i++) {
		int neighborcompound_id = neighborlist.neighborcompound_ids[i];
		int neighborcompound_particles = box->compound_state_array[neighborcompound_id].n_particles;

		if (threadIdx.x < neighborcompound_particles) {
			utility_buffer[threadIdx.x] = box->compound_state_array[neighborcompound_id].positions[threadIdx.x];
			EngineUtils::applyHyperpos(&compound_state.positions[0], &utility_buffer[threadIdx.x]);
		}

		BondedParticlesLUT* compoundpair_lut = box->bonded_particles_lut_manager->get(compound_index, neighborcompound_id);
		bonded_particles_lut.load(*compoundpair_lut);
		__syncthreads();

		if (threadIdx.x < compound.n_particles) {
			force += computerIntercompoundLJForces(&compound_state.positions[threadIdx.x], compound.atom_types[threadIdx.x], &potE_sum, compound.particle_global_ids[threadIdx.x], data_ptr,
				&box->compounds[neighborcompound_id], utility_buffer, neighborcompound_id, bonded_particles_lut);
		}
		__syncthreads();
	}
	// ------------------------------------------------------------------------------------------------------------------------------------------------------ //



	// --------------------------------------------------------------- Solvation forces --------------------------------------------------------------- //
#ifdef ENABLE_SOLVENTS
	for (int i = 0; i * blockDim.x < neighborlist.n_solvent_neighbors; i++) {
		int solvent_nlist_index = i * blockDim.x + threadIdx.x; // index in neighborlist

		if (solvent_nlist_index < neighborlist.n_solvent_neighbors) {
			utility_buffer[threadIdx.x] = box->solvents[neighborlist.neighborsolvent_ids[solvent_nlist_index]].pos;
			EngineUtils::applyHyperpos(&compound_state.positions[0], &utility_buffer[threadIdx.x]);
		}
		__syncthreads();

		if (threadIdx.x < compound.n_particles) {
			force_LJ_sol += computeSolventToCompoundLJForces(&compound_state.positions[threadIdx.x], blockDim.x, utility_buffer, data_ptr, &potE_sum, compound.atom_types[threadIdx.x]);
		}
		__syncthreads();
	}
	force += force_LJ_sol;
#endif
	// ------------------------------------------------------------------------------------------------------------------------------------------------ //



	// ------------------------------------------------------------ Integration ------------------------------------------------------------ //
	if (threadIdx.x < compound.n_particles) {
		if (box->step >= RAMPUP_STEPS) {
			integratePosition(&compound_state.positions[threadIdx.x], &compound.prev_positions[threadIdx.x], &force, forcefield_device.particle_parameters[compound.atom_types[threadIdx.x]].mass, box->dt, &box->thermostat_scalar, threadIdx.x, box->step > 810);
		}
		else {
			integratePositionRampUp(&compound_state.positions[threadIdx.x], &compound.prev_positions[threadIdx.x], &force, forcefield_device.particle_parameters[compound.atom_types[threadIdx.x]].mass, box->dt, box->step);
		}
		box->compounds[blockIdx.x].prev_positions[threadIdx.x] = compound.prev_positions[threadIdx.x];
	}
	__syncthreads();
	// ------------------------------------------------------------------------------------------------------------------------------------- //





	// ------------------------------------------------------ PERIODIC BOUNDARY CONDITION ------------------------------------------------------------------------------- // 
	if (threadIdx.x == 0) {
		EngineUtils::applyPBC(&compound_state.positions[threadIdx.x]);
	}
	__syncthreads();
	EngineUtils::applyHyperpos(&compound_state.positions[0], &compound_state.positions[threadIdx.x]);	// So all particles follows p0
	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------ //


	LogCompoundData(compound, box, compound_state, &potE_sum, force, force_LJ_sol);

	if (force.len() > 2e+10) {
		printf("\n\nCritical force %.0f           block %d thread %d\n\n\n", force.len(), blockIdx.x, threadIdx.x);
		box->critical_error_encountered = true;
	}


	box->compound_state_array_next[blockIdx.x].loadData(&compound_state);
}
#undef compound_id









#define solvent_index (blockIdx.x * blockDim.x + threadIdx.x)
#define thread_active (solvent_index < box->n_solvents)
__global__ void solventForceKernel(Box* box) {
	__shared__ Float3 utility_buffer[MAX_COMPOUND_PARTICLES];
	__shared__ uint8_t utility_buffer_small[MAX_COMPOUND_PARTICLES];

	float potE_sum = 0;
	float data_ptr[4];	// Pot, force, closest particle, ?
	for (int i = 0; i < 4; i++)
		data_ptr[i] = 0;
	data_ptr[2] = 9999.f;


	Float3 force(0.f);
	Solvent solvent;

	if (thread_active) {
		solvent = box->solvents[solvent_index];
	}

	// --------------------------------------------------------------- Molecule Interactions --------------------------------------------------------------- //
	for (int i = 0; i < box->n_compounds; i++) {
		int n_compound_particles = box->compound_state_array[i].n_particles;
		// First all threads help loading the molecule
		if (threadIdx.x < n_compound_particles) {
			utility_buffer[threadIdx.x] = box->compound_state_array[i].positions[threadIdx.x];
			utility_buffer_small[threadIdx.x] = box->compounds[i].atom_types[threadIdx.x];
		}
		__syncthreads();

		//  We can optimize here by loading and calculate the paired sigma and eps, jsut remember to loop threads, if there are many aomttypes.


		if (thread_active) {
			EngineUtils::applyHyperpos(&utility_buffer[0], &solvent.pos);									// Move own particle in relation to compound-key-position
			force += computeCompoundToSolventLJForces(&solvent.pos, n_compound_particles, utility_buffer, data_ptr, &potE_sum, ATOMTYPE_SOL, utility_buffer_small);
		}
		__syncthreads();
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //



	// --------------------------------------------------------------- Solvent Interactions --------------------------------------------------------------- //
	if (thread_active) {
		force += computeSolventToSolventLJForces(&solvent.pos, &box->solvent_neighborlists[solvent_index], box->solvents, data_ptr, &potE_sum);
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //


	//if (solvent_index < box->n_solvents)
	//	box->traj_buffer[box->n_compounds * MAX_COMPOUND_PARTICLES + solvent_index + (box->step % STEPS_PER_LOGTRANSFER) * box->total_particles_upperbound] = Float3(0.);

	LogSolventData(thread_active, box, solvent_index, force, solvent);


	if (thread_active) {
		int p_index = MAX_COMPOUND_PARTICLES + solvent_index;
		if (box->step >= RAMPUP_STEPS) {
			integratePosition(&solvent.pos, &solvent.pos_tsub1, &force, forcefield_device.particle_parameters[ATOMTYPE_SOL].mass, box->dt, &box->thermostat_scalar, p_index, true);
		}
		else {
			integratePositionRampUp(&solvent.pos, &solvent.pos_tsub1, &force, forcefield_device.particle_parameters[ATOMTYPE_SOL].mass, box->dt, box->step);
		}

		EngineUtils::applyPBC(&solvent.pos);
	}


	if (thread_active) {
		if (EngineUtils::calcHyperDist(&solvent.pos, &solvent.pos_tsub1) > 0.05) { 
			printf("Solvent was too fast: %f\n", (solvent.pos - solvent.pos_tsub1).len());
			box->critical_error_encountered = true; 
		}
	}


	if (thread_active) {
		if (solvent.pos.x != solvent.pos.x) {
			solvent.pos.print('s');
			box->critical_error_encountered = true;
		}
		box->solvents_next[solvent_index] = solvent;
	}
}
#undef solvent_index
#undef thread_active




//#define compound_id blockIdx.x
#define particle_id_bridge threadIdx.x
__global__ void compoundBridgeKernel(Box* box) {
	__shared__ CompoundBridgeCompact bridge;
	__shared__ Float3 positions[MAX_PARTICLES_IN_BRIDGE];
	__shared__ Float3 utility_buffer[MAX_PARTICLES_IN_BRIDGE];							// waaaaay too biggg

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
	{											// So for the very first step, these ´should all be 0, but they are not??										TODO: Look into this at some point!!!! 
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

		const int compound_offset = p_ref->compound_id * MAX_COMPOUND_PARTICLES;
		const int step_offset = (box->step % STEPS_PER_LOGTRANSFER) * box->total_particles_upperbound;

		box->potE_buffer[p_ref->local_id_compound + compound_offset + step_offset] = potE_sum;
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