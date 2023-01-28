#include "Engine.cuh"

#pragma warning ( push )
#pragma warning ( disable: E0020 )

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

__device__ Float3 calcLJForce(Float3* pos0, Float3* pos1, float* data_ptr, float* potE, const float sigma, const float epsilon, int type1 = -1, int type2 = -1) {
	// Calculates LJ force on p0	(attractive to p1. Negative values = repulsion )//
	// input positions in cartesian coordinates [nm]
	// sigma [nm]
	// epsilon [J/mol]->[(kg*nm^2)/(ns^2*mol)]
	// Returns force in J/mol*M		?????????????!?!?//

	// Directly from book
	float dist_sq = (*pos1 - *pos0).lenSquared();
	float s = (sigma * sigma) / dist_sq;								// [nm^2]/[nm^2] -> unitless	// OPTIM: Only calculate sigma_squared, since we never use just sigma
	s = s * s * s;
	float force_scalar = 24.f * epsilon * s / dist_sq * (1.f - 2.f * s) ;// *FEMTO_TO_LIMA* FEMTO_TO_LIMA;	// Attractive. Negative, when repulsive		[(kg*nm^2)/(nm^2*ns^2*mol)] ->----------------------	[(kg)/(ns^2*mol)]	

	*potE += 4. * epsilon * s * (s - 1.f) * 0.5;

#ifdef LIMA_VERBOSE
		//if (threadIdx.x == 0 && blockIdx.x == 0) {
		//	float dist = (*pos0 - *pos1).len() * NORMALIZER;
		//	printf("\ndist %f force %f pot %f, sigma %f, s %f distsq %f eps %.10f\n", dist, force_scalar * (*pos1 - *pos0).len(), *potE, sigma * NORMALIZER, s, dist_sq, epsilon);
		//}
		pos0->print('0');
		pos0->print('1');
	}
#endif

	return (*pos1 - *pos0) * force_scalar;										// GN/mol [(kg*nm)/(ns^2*mol)]
}


__device__ void calcPairbondForces(Float3* pos_a, Float3* pos_b, PairBond* bondtype, Float3* results, float* potE) {
	// Calculates bond force on both particles					//
	// Calculates forces as J/mol*M								//
	// kb [J/(mol*lm^2)]
	Float3 difference = *pos_a - *pos_b;						//	[lm]
	float error = difference.len() - bondtype->b0;				//	[lm]

	*potE += 0.5f * bondtype->kb * (error * error);				// [J/mol]
	
	float force_scalar = -bondtype->kb * error;				//	[J/(mol*lm)] = [kg/(mol*s^2)]

	difference = difference.norm();								// dif_unit_vec, but shares variable with dif
	results[0] = difference * force_scalar;						// [kg * lm / (mol*ls^2)] = [lN]
	results[1] = difference * force_scalar * -1;				// [kg * lm / (mol*ls^2)] = [lN]

	if (0.5f * bondtype->kb * (error * error) > 100000) {
		//printf("error %f    pote %f    force %f\n", error, 0.5f * bondtype->kb * (error * error), force_scalar);
	}
}


__device__ void calcAnglebondForces(Float3* pos_left, Float3* pos_middle, Float3* pos_right, AngleBond* angletype, Float3* results, float* potE) {
	const Float3 v1 = (*pos_left - *pos_middle).norm();
	const Float3 v2 = (*pos_right - *pos_middle).norm();
	const Float3 normal = v1.cross(v2).norm();	// Poiting towards y, when right is pointing toward x

	const Float3 inward_force_direction1 = (v1.cross(normal * -1.f)).norm();
	const Float3 inward_force_direction2 = (v2.cross(normal)).norm();

	const float angle = Float3::getAngle(v1, v2);					// [rad]
	const float error = angle - angletype->theta_0;				// [rad]

	// Simple implementation
	//*potE += angletype->k_theta * error * error * 0.5f;		// Energy [J/mol]0
	//float torque = angletype->k_theta * (error);				// Torque [J/(mol*rad)]

	// Correct implementation
	*potE += -angletype->k_theta * (cosf(error) - 1.f);		// Energy [J/mol]
	const float torque = angletype->k_theta * sinf(error);	// Torque [J/(mol*rad)]

	results[0] = inward_force_direction1 * (torque / (*pos_left - *pos_middle).len());
	results[2] = inward_force_direction2 * (torque / (*pos_right - *pos_middle).len());
	results[1] = (results[0] + results[2]) * -1.f;


	if (threadIdx.x == 2 && blockIdx.x == 0 || 1) {
		//printf("\nangle %f error %f force %f t0 %f kt %f\n", angle, error, force_scalar, angletype->theta_0, angletype->k_theta);
	}	
}
__device__ void calcDihedralbondForces(Float3* pos_left, Float3* pos_lm, Float3* pos_rm, Float3* pos_right, DihedralBond* dihedral, Float3* results, float* potE) {
	Float3 normal1 = (*pos_left - *pos_lm).cross((*pos_rm - *pos_lm)).norm();		
	Float3 normal2 = (*pos_lm - *pos_rm).cross((*pos_right - *pos_rm)).norm();	
	// Both vectors point "left" (looking from lm to rm). 

	float torsion = Float3::getAngle(normal2, normal1);
	//const float torsion = normal1.getAngleSigned(normal2);

	const bool angle_is_negative = (normal2.dot(*pos_left - *pos_lm)) > 0.f;
	torsion = angle_is_negative ? torsion * -1.f : torsion;

	normal2 *= -1;
	// Now  normal2 is flipped meaning both vectors point inward when 0 < torsion < 3.14, and outwards otherwise

	*potE += dihedral->k_phi * (1.f + cosf(dihedral->n * torsion - dihedral->phi_0));

	float torque = -dihedral->k_phi * (dihedral->n * sinf(dihedral->n * torsion - dihedral->phi_0));

	if (0) {
		normal1.print('1');
		normal2.print('2');
		pos_left->print('L');
		pos_lm->print('l');
		pos_rm->print('r');
		pos_right->print('R');
		printf("angle neg %d\n", angle_is_negative);
		//printf("torsion %f      ref %f     error %f     force: %f\n", torsion, dihedral->phi_0, error, force_scalar);
		float pot = dihedral->k_phi * (1 + cosf(dihedral->n * torsion - dihedral->phi_0));
		printf("torsion %f     torque: %f    pot %f     phi_0 %f k_phi %f\n", torsion, torque, pot, dihedral->phi_0, dihedral->k_phi);
	}



	results[0] = normal1 * (torque / (*pos_left-*pos_lm).len());
	results[3] = normal2 * (torque / (*pos_right-*pos_rm).len());
	// Not sure about the final two forces, for now we'll jsut split the sum of opposite forces between them.
	//results[1] = (results[0] + results[3]) * -1.f * 0.5;
	//results[2] = (results[0] + results[3]) * -1.f * 0.5;
	results[1] = (results[0]) * -1.f;
	results[2] = (results[3]) * -1.f;
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

//__device__ Float3 computeSolventToSolventLJForces(Float3* self_pos, NeighborList* nlist, Solvent* solvents, float* data_ptr, float* potE_sum) {	// Specific to solvent kernel
//	Float3 force(0.f);
//	//if (blockIdx.x + threadIdx.x == 0) printf("n neighbor solvents %d\n", nlist->n_solvent_neighbors);
//	for (int i = 0; i < nlist->n_solvent_neighbors; i++) {
//		Solvent neighbor = solvents[nlist->neighborsolvent_ids[i]];
//		EngineUtils::applyHyperpos(&neighbor.pos, self_pos);
//		force += calcLJForce(self_pos, &neighbor.pos, data_ptr, potE_sum,
//			forcefield_device.particle_parameters[ATOMTYPE_SOL].sigma,
//			forcefield_device.particle_parameters[ATOMTYPE_SOL].epsilon);
//	}
//	return force;// *24.f * 1e-9;
//}
__device__ Float3 computeSolventToCompoundLJForces(Float3* self_pos, int n_particles, Float3* positions, float* data_ptr, float* potE_sum, uint8_t atomtype_self) {	// Specific to solvent kernel
	Float3 force{};
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
	utility_buffer[threadIdx.x] = Float3(0);
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

__device__ void integratePosition(Coord& coord, Coord& coord_tsub1, Float3* force, const float mass, const double dt, const float thermostat_scalar) {
	//EngineUtils::applyHyperpos(pos, pos_tsub1);

	//Float3 temp = *pos;
	//float prev_vel = (*pos - *pos_tsub1).len();
	Coord prev_vel = coord - coord_tsub1;

	// [nm] - [nm] + [kg/mol*m*/s ^ 2] / [kg / mol] * [s] ^ 2 * (1e-9) ^ 2 = > [nm] - [nm] + []
	//int32_t x = coord.x;
	//int32_t dx = coord.x - coord_tsub1.x;
	//int32_t ddx = static_cast<int32_t>(force->x * dt * dt / static_cast<double>(mass));
	//*pos += (*pos - *pos_tsub1) + *force * dt * (dt / static_cast<double>(mass));

	coord += (coord - coord_tsub1) * thermostat_scalar + *force * dt * dt / mass;
	//pos->x = x + dx + ddx;
	//coord.x = x + dx + ddx;

	//Double3 pos_d{ *pos };

	if (threadIdx.x == 0 && blockIdx.x == 0) {
		//force->print('f');
		//printf("dt %f mass %f\n", dt, mass);
		//(*force*dt*dt / mass).print('F');
		//uint32_t diff = coord.x - x - dx;
		//printf("x  %d  dx %d  force %.10f ddx %d    x_ %d   dif %d\n", x, dx, force->x, ddx, dx + ddx, diff);
	}
}

__device__ void integratePosition(Float3* pos, Float3* pos_tsub1, Float3* force, const float mass, const float dt, float* thermostat_scalar, int p_index, bool print = false) {
	// Force is in ??Newton, [kg * nm /(mol*ns^2)] //	

	// Always do hyperpos first!				TODO: Remove, as this is now done by moving the origo, right?
	//EngineUtils::applyHyperpos(pos, pos_tsub1);

	Float3 temp = *pos;
	float prev_vel = (*pos - *pos_tsub1).len();

	// [nm] - [nm] + [kg/mol*m*/s ^ 2] / [kg / mol] * [s] ^ 2 * (1e-9) ^ 2 = > [nm] - [nm] + []
	/*double x = pos->x;
	double dx = pos->x - pos_tsub1->x;
	double ddx = force->x * dt * dt / static_cast<double>(mass);*/

	// Without thermostat
	//*pos += (*pos - *pos_tsub1) + *force * dt * (dt / static_cast<double>(mass));
	*pos += (*pos - *pos_tsub1) * thermostat_scalar + *force * dt * (dt / static_cast<double>(mass));
	//pos->x = x + dx + ddx;


	//Double3 pos_d{ *pos };

	if (threadIdx.x + blockIdx.x == 0) {
		//printf("step %d f %f\n", 0, force->x);
		//float actual_x = x * NORMALIZER * LIMA_SCALE;
		//double diff = (double)pos->x - x - dx;
		//printf("x  %.10f  dx %.10f  force %.10f ddx %.10f    x_ %.10f   dif %.10f\n", x, dx, force->x, ddx, dx + ddx, diff);
	}
	////*pos += x_;
	////*pos = Float3{ pos_new.x, pos_new.y, pos_new.z };
	//if (print && threadIdx.x == 1) {
	//	//printf("Thread %d    prev_vel %f    Force scalar %f    Acc %.10f    Change %.10f\n", threadIdx.x, prev_vel, force->len(), acc, (*pos - *pos_tsub1).len() - prev_vel);
	//	//x.print();
	
	//		//printf("x  %.8f  dx %.8f   ddx %.8f    x_ %.8f   dif %.8f\n", x_rel.len(), dx.len(), ddx.len(), x_.len(), ((x_new_abs) - *pos).len());
	//	//pos->print('f');
	//	//pos_new.print('d');
	//}

	*pos_tsub1 = temp;

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


__device__ inline void LogCompoundData(Compound& compound, Box* box, CompoundCoords& compound_coords, float* potE_sum, Float3& force, Float3& force_LJ_sol) {
	const uint32_t index = EngineUtils::getLoggingIndexOfParticle(box->step, box->total_particles_upperbound, blockIdx.x, threadIdx.x);
	box->traj_buffer[index] = LIMAPOSITIONSYSTEM::getGlobalPosition(compound_coords);
	box->potE_buffer[index] = *potE_sum;
}

//__device__ inline void LogSolventData(bool thread_active, Box* box, int solvent_index, Float3& force, Solvent& solvent) {
//	if (thread_active) {
//		//int compounds_offset = box->n_compounds * MAX_COMPOUND_PARTICLES;
//		//int step_offset = (box->step % STEPS_PER_LOGTRANSFER) * box->total_particles_upperbound;
//			//box->potE_buffer[compounds_offset + solvent_index + step_offset] = force.len();
//		//box->traj_buffer[compounds_offset + solvent_index + step_offset] = solvent.pos;
//		const uint32_t index = EngineUtils::getLoggingIndexOfParticle(box->step, box->total_particles_upperbound, box->n_compounds, solvent_index);
//		box->potE_buffer[index] = force.len();
//		box->traj_buffer[index] = solvent.pos;
//	}
//}

__device__ void getCompoundHyperpositionsAsFloat3(const Coord& origo_self, const CompoundCoords* querycompound, Float3* output_buffer, Coord* utility_coord, int step) { 
	if (threadIdx.x == 0) {
		Coord querycompound_hyperorigo = LIMAPOSITIONSYSTEM::getHyperOrigo(origo_self, querycompound->origo);

		// calc Relative Position Shift from the origo-shift
		*utility_coord = LIMAPOSITIONSYSTEM::getRelShift(origo_self, querycompound_hyperorigo);
	}
	__syncthreads();

	//Coord prev_rel_pos = box->compound_coord_array_prev[blockIdx.x].rel_positions[threadIdx.x] - rel_pos_shift;
	Coord queryparticle_coord = querycompound->rel_positions[threadIdx.x] - *utility_coord;
	output_buffer[threadIdx.x] = queryparticle_coord.toFloat3();
}


// ------------------------------------------------------------------------------------------- KERNELS -------------------------------------------------------------------------------------------//




#define compound_index blockIdx.x
__global__ void compoundKernel(Box* box) {
	__shared__ Compound compound;				// Mostly bond information
	__shared__ CompoundState compound_state;	// Relative position in [lm]
	__shared__ CompoundCoords compound_coords;	// Global positions in [lm]
	__shared__ NeighborList neighborlist;		
	__shared__ BondedParticlesLUT bonded_particles_lut;
	__shared__ Float3 utility_buffer[THREADS_PER_COMPOUNDBLOCK];
	__shared__ Coord utility_coord;
	__shared__ Coord rel_pos_shift;	// Maybe not needed, jsut use the utility one above?

	if (threadIdx.x == 0) {
		compound.loadMeta(&box->compounds[blockIdx.x]);
		compound_state.setMeta(compound.n_particles);
		neighborlist.loadMeta(&box->compound_neighborlists[blockIdx.x]);
	}
	__syncthreads();




	compound.loadData(&box->compounds[blockIdx.x]);
	auto* coordarray_ptr = CoordArrayQueueHelpers::getCoordarrayPtr(box->coordarray_circular_queue, box->step, blockIdx.x);
	compound_coords.loadData(*coordarray_ptr);
	neighborlist.loadData(&box->compound_neighborlists[blockIdx.x]);
	__syncthreads();

	LIMAPOSITIONSYSTEM::getRelativePositions(compound_coords, compound_state);
	__syncthreads();


	float potE_sum = 0;
	float data_ptr[4]{};
	for (int i = 0; i < 4; i++)
		data_ptr[i] = 0;
	data_ptr[2] = 9999;
	data_ptr[3] = box->step + 1;


	Float3 force = compound.forces[threadIdx.x];
	Float3 force_LJ_sol(0.f);


	if (blockIdx.x == 1 && threadIdx.x == 0) {
		compound_state.positions[0].print('S');
		force.print('F');
	}
	// ------------------------------------------------------------ Intramolecular Operations ------------------------------------------------------------ //
	{
		bonded_particles_lut.load(*box->bonded_particles_lut_manager->get(compound_index, compound_index));
		__syncthreads();

		force += computePairbondForces(&compound, compound_state.positions, utility_buffer, &potE_sum);
		force += computeAnglebondForces(&compound, compound_state.positions, utility_buffer, &potE_sum);
		force += computeDihedralForces(&compound, compound_state.positions, utility_buffer, &potE_sum);
		force += computeIntracompoundLJForces(&compound, &compound_state, &potE_sum, data_ptr, &bonded_particles_lut);
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------- //


	// --------------------------------------------------------------- Intermolecular forces --------------------------------------------------------------- //	
	for (int i = 0; i < neighborlist.n_compound_neighbors; i++) {
		const uint16_t neighborcompound_id = neighborlist.neighborcompound_ids[i];
		


		const auto coords_ptr = CoordArrayQueueHelpers::getCoordarrayPtr(box->coordarray_circular_queue, box->step, neighborcompound_id);
		getCompoundHyperpositionsAsFloat3(compound_coords.origo, coords_ptr, utility_buffer, &utility_coord, box->step);

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
	//for (int i = 0; i * blockDim.x < neighborlist.n_solvent_neighbors; i++) {
	//	int solvent_nlist_index = i * blockDim.x + threadIdx.x; // index in neighborlist

	//	if (solvent_nlist_index < neighborlist.n_solvent_neighbors) {
	//		utility_buffer[threadIdx.x] = box->solvents[neighborlist.neighborsolvent_ids[solvent_nlist_index]].pos;
	//		EngineUtils::applyHyperpos(&compound_state.positions[0], &utility_buffer[threadIdx.x]);
	//	}
	//	__syncthreads();

	//	if (threadIdx.x < compound.n_particles) {
	//		force_LJ_sol += computeSolventToCompoundLJForces(&compound_state.positions[threadIdx.x], blockDim.x, utility_buffer, data_ptr, &potE_sum, compound.atom_types[threadIdx.x]);
	//	}
	//	__syncthreads();
	//}
	//force += force_LJ_sol;
#endif
	// ------------------------------------------------------------------------------------------------------------------------------------------------ //


	// ------------------------------------------------------------ Integration ------------------------------------------------------------ //
	{
		// For the very first step, we need to fetch the prev position from the last index of the circular queue! 
		const int actual_stepindex_of_prev = box->step == 0 ? STEPS_PER_LOGTRANSFER - 1 : box->step - 1;
		const auto* coordarray_prev_ptr = CoordArrayQueueHelpers::getCoordarrayPtr(box->coordarray_circular_queue, actual_stepindex_of_prev, blockIdx.x);
		if (threadIdx.x == 0) {
			Coord prev_hyper_origo = LIMAPOSITIONSYSTEM::getHyperOrigo(compound_coords.origo, coordarray_prev_ptr->origo);
			rel_pos_shift = LIMAPOSITIONSYSTEM::getRelShift(compound_coords.origo, prev_hyper_origo);
		}
		__syncthreads();

		Coord prev_rel_pos = coordarray_prev_ptr->rel_positions[threadIdx.x] - rel_pos_shift;


		if (threadIdx.x < compound.n_particles) {
			integratePosition(compound_coords.rel_positions[threadIdx.x], prev_rel_pos, &force, forcefield_device.particle_parameters[compound.atom_types[threadIdx.x]].mass, box->dt, box->thermostat_scalar);
		}
		__syncthreads();
	}
	// ------------------------------------------------------------------------------------------------------------------------------------- //




	// ------------------------------------------------------ PERIODIC BOUNDARY CONDITION ------------------------------------------------------------------------------- // 
	LIMAPOSITIONSYSTEM::applyPBC(compound_coords);
	__syncthreads();

	
	{
		__shared__ Coord shift_lm;
		if (threadIdx.x == 0) {
			shift_lm = LIMAPOSITIONSYSTEM::shiftOrigo(compound_coords);
		}
		__syncthreads();

		//if (threadIdx.x == 0 && blockIdx.x == 1) {
		//	shift_lm.print('S');
		//}

		LIMAPOSITIONSYSTEM::shiftRelPos(compound_coords, shift_lm);
		__syncthreads();
	}
	//if (threadIdx.x == 0 && blockIdx.x == 0) LIMAPOSITIONSYSTEM::getGlobalPositionNM(compound_coords).print('a');
	//EngineUtils::applyHyperpos(&compound_state.positions[0], &compound_state.positions[threadIdx.x]);	// So all particles follows p0
	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------ //
	if (threadIdx.x == 0 && blockIdx.x == 0) {
		//compound_coords.origo.print('a');
		//compound_coords.rel_positions[0].print('A');
	}

	LogCompoundData(compound, box, compound_coords, &potE_sum, force, force_LJ_sol);

	if (force.len() > 2e+10) {
		printf("\n\nCritical force %.0f           block %d thread %d\n\n\n", force.len(), blockIdx.x, threadIdx.x);
		box->critical_error_encountered = true;
	}

	auto* coordarray_next_ptr = CoordArrayQueueHelpers::getCoordarrayPtr(box->coordarray_circular_queue, box->step + 1, blockIdx.x);
	coordarray_next_ptr->loadData(compound_coords);
}
#undef compound_id









#define solvent_index (blockIdx.x * blockDim.x + threadIdx.x)
#define thread_active (solvent_index < box->n_solvents)
__global__ void solventForceKernel(Box* box) {
	__shared__ Float3 utility_buffer[MAX_COMPOUND_PARTICLES];	
	__shared__ uint8_t utility_buffer_small[MAX_COMPOUND_PARTICLES];
	__shared__ CompoundCoords compound_coords;	// For neighboring compounds
	float potE_sum = 0;
	float data_ptr[4];	// Pot, force, closest particle, ?
	for (int i = 0; i < 4; i++)
		data_ptr[i] = 0;
	data_ptr[2] = 9999.f;


	//Float3 force(0.f);
	//Solvent solvent;

	//if (thread_active) {
	//	solvent = box->solvents[solvent_index];
	//}

	//// --------------------------------------------------------------- Molecule Interactions --------------------------------------------------------------- //
	//for (int i = 0; i < box->n_compounds; i++) {
	//	int n_compound_particles = box->compounds[i].n_particles;

	//	// First all threads help loading the molecule
	//	if (threadIdx.x < n_compound_particles) {
	//		// First load particles of neighboring compound
	//		auto* coordarray_ptr = CoordArrayQueueHelpers::getCoordarrayPtr(box->coordarray_circular_queue, box->step, blockIdx.x);
	//		compound_coords.loadData(*coordarray_ptr);
	//		LIMAPOSITIONSYSTEM::getRelativePositions(compound_coords.rel_positions, utility_buffer);

	//		// Then load atomtypes of neighboring compound
	//		utility_buffer_small[threadIdx.x] = box->compounds[i].atom_types[threadIdx.x];
	//	}
	//	__syncthreads();

	//	//  We can optimize here by loading and calculate the paired sigma and eps, jsut remember to loop threads, if there are many aomttypes.


	//	if (thread_active) {
	//		EngineUtils::applyHyperpos(&utility_buffer[0], &solvent.pos);									// Move own particle in relation to compound-key-position
	//		force += computeCompoundToSolventLJForces(&solvent.pos, n_compound_particles, utility_buffer, data_ptr, &potE_sum, ATOMTYPE_SOL, utility_buffer_small);
	//	}
	//	__syncthreads();
	//}
	//// ----------------------------------------------------------------------------------------------------------------------------------------------------- //



	//// --------------------------------------------------------------- Solvent Interactions --------------------------------------------------------------- //
	//if (thread_active) {
	//	force += computeSolventToSolventLJForces(&solvent.pos, &box->solvent_neighborlists[solvent_index], box->solvents, data_ptr, &potE_sum);
	//}
	//// ----------------------------------------------------------------------------------------------------------------------------------------------------- //


	////if (solvent_index < box->n_solvents)
	////	box->traj_buffer[box->n_compounds * MAX_COMPOUND_PARTICLES + solvent_index + (box->step % STEPS_PER_LOGTRANSFER) * box->total_particles_upperbound] = Float3(0.);

	//LogSolventData(thread_active, box, solvent_index, force, solvent);


	//if (thread_active) {
	//	int p_index = MAX_COMPOUND_PARTICLES + solvent_index;
	//	if (box->step >= RAMPUP_STEPS) {
	//		integratePosition(&solvent.pos, &solvent.pos_tsub1, &force, forcefield_device.particle_parameters[ATOMTYPE_SOL].mass, box->dt, &box->thermostat_scalar, p_index, true);
	//	}
	//	else {
	//		integratePositionRampUp(&solvent.pos, &solvent.pos_tsub1, &force, forcefield_device.particle_parameters[ATOMTYPE_SOL].mass, box->dt, box->step);
	//	}

	//	EngineUtils::applyPBC(&solvent.pos);
	//}


	//if (thread_active) {
	//	if (EngineUtils::calcHyperDist(&solvent.pos, &solvent.pos_tsub1) > 0.05) { 
	//		printf("Solvent was too fast: %f\n", (solvent.pos - solvent.pos_tsub1).len());
	//		box->critical_error_encountered = true; 
	//	}
	//}


	//if (thread_active) {
	//	if (solvent.pos.x != solvent.pos.x) {
	//		solvent.pos.print('s');
	//		box->critical_error_encountered = true;
	//	}
	//	box->solvents_next[solvent_index] = solvent;
	//}
}
#undef solvent_index
#undef thread_active













#define particle_id_bridge threadIdx.x
__global__ void compoundBridgeKernel(Box* box) {
	__shared__ CompoundBridgeCompact bridge;
	__shared__ Float3 positions[MAX_PARTICLES_IN_BRIDGE];
	__shared__ Float3 utility_buffer[MAX_PARTICLES_IN_BRIDGE];							// waaaaay too biggg
	__shared__ Coord particle_coords[MAX_PARTICLES_IN_BRIDGE];
	__shared__ Coord utility_coord;

	if (threadIdx.x == 0) {
		bridge.loadMeta(&box->bridge_bundle->compound_bridges[blockIdx.x]);

		// Calculate necessary shift in relative positions for right, so right share the origo with left.
		utility_coord = LIMAPOSITIONSYSTEM::getRelativeShiftBetweenCoordarrays(box->coordarray_circular_queue, box->step, bridge.compound_id_left, bridge.compound_id_right);
	}
	__syncthreads();

	bridge.loadData(&box->bridge_bundle->compound_bridges[blockIdx.x]);
	__syncthreads();

	if (particle_id_bridge < bridge.n_particles) {
		ParticleRefCompact& p_ref = bridge.particle_refs[particle_id_bridge];
		particle_coords[particle_id_bridge] = CoordArrayQueueHelpers::getCoordarrayPtr(box->coordarray_circular_queue, box->step, p_ref.compound_id)->origo;

		// If we are on the right side, we need to shift 
		if (p_ref.compound_id == bridge.compound_id_right) { particle_coords[particle_id_bridge] += utility_coord; }
	}
	__syncthreads();

	// Load the now shifted relative coords into float3 positions for force calcs.
	LIMAPOSITIONSYSTEM::getRelativePositions(particle_coords, positions);

	float potE_sum = 0;
	Float3 force{};

	// ------------------------------------------------------------ Intercompund Operations ------------------------------------------------------------ //
	{											// So for the very first step, these ´should all be 0, but they are not??										TODO: Look into this at some point!!!! 
		force += computePairbondForces(&bridge, positions, utility_buffer, &potE_sum);
		force += computeAnglebondForces(&bridge, positions, utility_buffer, &potE_sum);
		force += computeDihedralForces(&bridge, positions, utility_buffer, &potE_sum);
	}
	__syncthreads();
	// --------------------------------------------------------------------------------------------------------------------------------------------------- //

	if (particle_id_bridge < bridge.n_particles) {
		ParticleRefCompact* p_ref = &bridge.particle_refs[particle_id_bridge];
		//box->compounds[p_ref->compound_id].forces[p_ref->local_id_compound] = force;
		box->compounds[p_ref->compound_id].forces[p_ref->local_id_compound] = 0;

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

#pragma warning (pop)