#include "Engine.cuh"

#include "ForceComputations.cuh"

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



__device__ Float3 computerIntercompoundLJForces(Float3* self_pos, uint8_t atomtype_self, float* potE_sum, uint32_t global_id_self, float* data_ptr,
	Compound* neighbor_compound, Float3* neighbor_positions, int neighborcompound_id, BondedParticlesLUT& bonded_particles_lut) {
	Float3 force(0.f);
	
	for (int neighborparticle_id = 0; neighborparticle_id < neighbor_compound->n_particles; neighborparticle_id++) {

		// If thread's assoc. particle is bonded to the particle in neighborcompound, continue
		if (*bonded_particles_lut.get(threadIdx.x, neighborparticle_id)) { continue; }

		const int neighborparticle_atomtype = neighbor_compound->atom_types[neighborparticle_id];	//// TEMPORARY, this is waaaay to many global mem accesses

		force += LimaForcecalc::calcLJForce(self_pos, &neighbor_positions[neighborparticle_id], data_ptr, potE_sum,
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
			//force += calcLJForce(&compound_state->positions[threadIdx.x], &compound_state->positions[i], data_ptr, potE_sum,

			force += LimaForcecalc::calcLJForce(&compound_state->positions[threadIdx.x], &compound_state->positions[i], data_ptr, potE_sum,
				calcSigma(compound->atom_types[threadIdx.x], compound->atom_types[i]), calcEpsilon(compound->atom_types[threadIdx.x], compound->atom_types[i]),
				compound->particle_global_ids[threadIdx.x], compound->particle_global_ids[i]
			);// *24.f * 1e-9;
		}
	}
	return force;
}

__device__ Float3 computeSolventToSolventLJForces(const SolventCoord& coord_self, const NeighborList& nlist, const SolventCoord* solvent_coords, float* data_ptr, float* potE_sum) {	// Specific to solvent kernel
	Float3 force{};
	const Float3 pos_self = coord_self.rel_position.toFloat3();

	for (int i = 0; i < nlist.n_solvent_neighbors; i++) {
		const SolventCoord& coord_neighbor = solvent_coords[nlist.neighborsolvent_ids[i]];

		// Check if the two solvents are too far away for Coord's to represent them
		if (!LIMAPOSITIONSYSTEM::canRepresentRelativeDist(coord_self.origo, coord_neighbor.origo)) { continue; }

		// Take hyperposition here. Both positions are now relative.
		const Float3 pos_neighbor = LIMAPOSITIONSYSTEM::getRelativeHyperposition(coord_self, coord_neighbor).toFloat3();

		if (threadIdx.x == 0 && (pos_self - pos_neighbor).len() == 0) {
			(pos_self - pos_neighbor).print('d');
			printf("n id %d id self %d\n\n\n", nlist.neighborsolvent_ids[i], threadIdx.x + blockIdx.x * THREADS_PER_SOLVENTBLOCK);
		}
		//continue;
		force += LimaForcecalc::calcLJForce(&pos_self, &pos_neighbor, data_ptr, potE_sum,
			forcefield_device.particle_parameters[ATOMTYPE_SOL].sigma,
			forcefield_device.particle_parameters[ATOMTYPE_SOL].epsilon);
	}
	return force;// *24.f * 1e-9;
}
__device__ Float3 computeSolventToCompoundLJForces(Float3* self_pos, int n_particles, Float3* positions, float* data_ptr, float* potE_sum, uint8_t atomtype_self) {	// Specific to solvent kernel
	Float3 force{};
	for (int i = 0; i < n_particles; i++) {
		force += LimaForcecalc::calcLJForce(self_pos, &positions[i], data_ptr, potE_sum,
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
		force += LimaForcecalc::calcLJForce(self_pos, &positions[i], data_ptr, potE_sum,
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

			LimaForcecalc::calcPairbondForces(
				&positions[pb->atom_indexes[0]],
				&positions[pb->atom_indexes[1]],
				pb,
				forces, potE
			);
		}
		if (blockIdx.x == 0 && threadIdx.x == 0) {
	/*		printf("\nbond_index %d n_bonds %d\n", bond_index, entity->n_singlebonds);
			forces[0].print('0');
			forces[1].print('1');*/
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


			LimaForcecalc::calcAnglebondForces(
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
			LimaForcecalc::calcDihedralbondForces(
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
//
//__device__ void integratePosition(Coord& coord, Coord& coord_tsub1, Float3* force, const float mass, const double dt, const float thermostat_scalar) {
//	//EngineUtils::applyHyperpos(pos, pos_tsub1);
//
//	//Float3 temp = *pos;
//	//float prev_vel = (*pos - *pos_tsub1).len();
//	Coord prev_vel = coord - coord_tsub1;
//
//	// [nm] - [nm] + [kg/mol*m*/s ^ 2] / [kg / mol] * [s] ^ 2 * (1e-9) ^ 2 = > [nm] - [nm] + []
//	//int32_t x = coord.x;
//	//int32_t dx = coord.x - coord_tsub1.x;
//	//int32_t ddx = static_cast<int32_t>(force->x * dt * dt / static_cast<double>(mass));
//	//*pos += (*pos - *pos_tsub1) + *force * dt * (dt / static_cast<double>(mass));
//
//	coord += (coord - coord_tsub1) * thermostat_scalar + *force * dt * dt / mass;
//	//pos->x = x + dx + ddx;
//	//coord.x = x + dx + ddx;
//
//	//Double3 pos_d{ *pos };
//
//	if (threadIdx.x == 0 && blockIdx.x == 0) {
//		//force->print('f');
//		//printf("dt %f mass %f\n", dt, mass);
//		//(*force*dt*dt / mass).print('F');
//		//uint32_t diff = coord.x - x - dx;
//		//printf("x  %d  dx %d  force %.10f ddx %d    x_ %d   dif %d\n", x, dx, force->x, ddx, dx + ddx, diff);
//	}
//}

// This function assumes that coord_tsub1 is already hyperpositioned to coord.
__device__ Coord integratePosition(const Coord& coord, const Coord& coord_tsub1, const Float3* force, const float mass, const double dt, const float thermostat_scalar) {

	return coord + (coord - coord_tsub1) * thermostat_scalar + *force * dt * dt / mass;
	//return coord + Coord{ 1000000, 0, 0 };
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

__device__ inline void LogSolventData(Box* box, const int solvent_id, const float& potE, const SolventCoord& coord) {
	const uint32_t index = EngineUtils::getLoggingIndexOfParticle(box->step, box->total_particles_upperbound, box->n_compounds, solvent_id);
	box->potE_buffer[index] = potE;
	box->traj_buffer[index] = coord.getAbsolutePositionLM();
}

__device__ void getCompoundHyperpositionsAsFloat3(const Coord& origo_self, const CompoundCoords* querycompound, Float3* output_buffer, Coord* utility_coord, int step) { 
	if (threadIdx.x == 0) {
		Coord querycompound_hyperorigo = LIMAPOSITIONSYSTEM::getHyperOrigo(origo_self, querycompound->origo);

		// calc Relative Position Shift from the origo-shift
		*utility_coord = LIMAPOSITIONSYSTEM::getRelShiftFromOrigoShift(querycompound_hyperorigo, origo_self);
	}
	__syncthreads();

	//Coord prev_rel_pos = box->compound_coord_array_prev[blockIdx.x].rel_positions[threadIdx.x] - rel_pos_shift;
	Coord queryparticle_coord = querycompound->rel_positions[threadIdx.x] + *utility_coord;
	output_buffer[threadIdx.x] = queryparticle_coord.toFloat3();
}



/// <summary>
/// 
/// </summary>
/// <param name="onehot_remainers">Must be size of MAX_SOLVENTS_IN_BLOCK</param>
/// <param name="utility">Must be large enough to store temporary sums for blelloch sum</param>
/// <returns></returns>
__device__ void doBlellochPrefixSum(uint8_t* onehot_remainers, uint8_t* utility) {
	// Forward scan
	for (int leap = 1; leap < MAX_SOLVENTS_IN_BLOCK - 1; leap *= 2) {
		if (threadIdx.x % leap == 0) {
			int index = threadIdx.x + leap;
		}
	}

}

// SLOW
__device__ void doSequentialPrefixSum(uint8_t* onehot_remainers, int n_elements) {
	for (int i = 1; i < n_elements; i++) {
		if (threadIdx.x == i) {
			onehot_remainers[i] = onehot_remainers[i - 1];
		}
		__syncthreads();
	}
}

__device__ uint8_t computePrefixSum(const bool remain, uint8_t* utility_buffer, int n_elements) {
	if (remain) {
		utility_buffer[threadIdx.x] = 1;
	}
	__syncthreads();
	doSequentialPrefixSum(utility_buffer, n_elements);
	//doBlellochPrefixSum
	return utility_buffer[threadIdx.x];
}













/// <summary>
/// </summary>
/// <param name="solventblock">In shared memory</param>
/// <param name="transferqueues">In shared memory</param>
/// <param name="relpos_next">Register</param>
/// <param name="transfermodules">Global memory</param>
__device__ void doSequentialForwarding(const SolventBlock& solventblock, STransferQueue* transferqueues, const Coord& relpos_next, SolventBlockTransfermodule* transfermodules) {
	const Coord transfer_direction = EngineUtils::getTransferDirection(relpos_next);
	const Coord blockId3d = SolventBlockHelpers::get3dIndex(blockIdx.x);
	const int new_blockid = EngineUtils::getNewBlockId(transfer_direction, blockId3d);

	// Sequential insertion in shared memory
	for (int i = 0; i < solventblock.n_solvents; i++) {
		if (threadIdx.x == i && new_blockid != blockIdx.x) {	// Only handle non-remain solvents
			const int queue_index = SolventBlockTransfermodule::getQueueIndex(transfer_direction);
			transferqueues[queue_index].addElement(relpos_next, solventblock.rel_pos[threadIdx.x]);
		}
		__syncthreads();
	}

	// Coaslescing copying to global memory
	for (int queue_index = 0; queue_index < 6; queue_index++) {
		STransferQueue& queue_local = transferqueues[queue_index];
		if (threadIdx.x < queue_local.n_elements) {
			const auto transferdir_queue = EngineUtils::getTransferDirection(queue_local.rel_positions[0]);		// Maybe use a utility-coord a precompute by thread0, or simply hardcode...
			const int blockid_global = EngineUtils::getNewBlockId(transfer_direction, blockId3d);

			STransferQueue* queue_global = &transfermodules[blockid_global].transfer_queues[queue_index];
			queue_global->rel_positions[threadIdx.x] = queue_local.rel_positions[threadIdx.x];
			queue_global->rel_positions_prev[threadIdx.x] = queue_local.rel_positions_prev[threadIdx.x];
		}
	}
	__syncthreads();
}

/// <summary>
/// Must be run AFTER doSequentialForwarding, as it erases all information about the transferring solvents
/// </summary>
/// <param name="solventblock_current">Solventblock in shared memory</param>
/// <param name="solventblock_next">Solventblock belong to cudablock at next step in global memory</param>
/// <param name="relpos_next"></param>
/// <param name="utility_buffer">Buffer of min size MAX_SOLVENTS_IN_BLOCK, maybe more for computing prefix sum</param>
/// <param name="remain_transfermodule">Transfermodule belonging to cudablock</param>
__device__ void purgeTransfersAndCompressRemaining(const SolventBlock& solventblock_current_local, SolventBlock* solventblock_next_global, 
	const Coord& relpos_next, uint8_t* utility_buffer, SolventBlockTransfermodule* remain_transfermodule) {
	// Clear buffer
	utility_buffer[threadIdx.x] = 0;

	// Find out which threads want to transfer their respective solvent
	const Coord transfer_direction = EngineUtils::getTransferDirection(relpos_next);
	const Coord blockId3d = SolventBlockHelpers::get3dIndex(blockIdx.x);
	const int new_blockid = EngineUtils::getNewBlockId(transfer_direction, blockId3d);
	const bool remain = new_blockid == blockIdx.x;

	// Compute prefix sum to find new index of solvent belonging to thread
	const uint8_t solventindex_new = computePrefixSum(remain, utility_buffer, solventblock_current_local.n_solvents);

	if (remain) {
		// Send current pos at threadindex to prevpos at the new index
		remain_transfermodule->remain_relpos_prev[solventindex_new] = solventblock_current_local.rel_pos[threadIdx.x];
		// Send the next pos 
		solventblock_next_global->rel_pos[solventindex_new] = relpos_next;
	}

	
	if (threadIdx.x == 0) {
		const int nsolventsinblock_next = __syncthreads_count(remain);
		remain_transfermodule->n_remain = nsolventsinblock_next;
		solventblock_next_global->n_solvents = nsolventsinblock_next;
	}
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

	// TODO: these compound_coords needs to be made const somehow. Changing these coords during the kernel is way too dangerous.
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
		//compound_state.positions[0].print('S');
		//force.print('F');
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
		// Dont think we need this anymore, if we bootstrap index N-1 with the positions aswell
		const int actual_stepindex_of_prev = box->step == 0 ? STEPS_PER_LOGTRANSFER - 1 : box->step - 1;
		const auto* coordarray_prev_ptr = CoordArrayQueueHelpers::getCoordarrayPtr(box->coordarray_circular_queue, actual_stepindex_of_prev, blockIdx.x);
		if (threadIdx.x == 0) {
			Coord prev_hyper_origo = LIMAPOSITIONSYSTEM::getHyperOrigo(compound_coords.origo, coordarray_prev_ptr->origo);
			rel_pos_shift = LIMAPOSITIONSYSTEM::getRelShiftFromOrigoShift(prev_hyper_origo, compound_coords.origo);
		}
		// Im not sure if this is correct..
		__syncthreads();

		Coord prev_rel_pos = coordarray_prev_ptr->rel_positions[threadIdx.x] - rel_pos_shift;


		if (threadIdx.x < compound.n_particles) {
			// Change this so it handles the return value??!?
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
#define thread_active (solvent_index < box->n_solvents)	// Remove this once shift is complete.
#define solvent_active (threadIdx.x < solventblock.n_solvents)
#define solvent_mass (forcefield_device.particle_parameters[ATOMTYPE_SOL].mass)
#define solventblock_ptr (CoordArrayQueueHelpers::getSolventBlockPtr(box->solventblockgrid_circurlar_queue, box->step, blockIdx.x))
//#define istransferstep (std::is_same<T, )

__global__ void solventForceKernel(Box* box) {
	__shared__ Float3 utility_buffer[MAX_COMPOUND_PARTICLES];	
	//__shared__ uint8_t utility_buffer_small[MAX_COMPOUND_PARTICLES];
	__shared__ uint8_t utility_buffer_small[MAX_SOLVENTS_IN_BLOCK];
	__shared__ CompoundCoords compound_coords;	// For neighboring compounds
	__shared__ SolventBlock solventblock;
	__shared__ SolventTransferqueue<SolventBlockTransfermodule::max_queue_size> transferqueues[6];		// TODO: Use template to make identical kernel, so the kernel with transfer is slower and larger, and the rest remain fast!!!!
	//__shared__ Coord coord_utility_buffer[MAX_SOLVENTS_IN_BLOCK + 6 * SolventBlockTransfermodule::max_queue_size];


	float potE_sum = 0;
	float data_ptr[4];	// Pot, force, closest particle, ?
	for (int i = 0; i < 4; i++)
		data_ptr[i] = 0;
	data_ptr[2] = 9999.f;




	if (threadIdx.x == 0) {
		solventblock.loadMeta(*solventblock_ptr);
		//solventblock.origo.print();
	}
	__syncthreads();
	solventblock.loadData(*solventblock_ptr);
	__syncthreads();

	//if (threadIdx.x == 0 && blockIdx.x == SolventBlockHelpers::get1dIndex({ 3, 3, 3 })) {
	//	printf("\n1d index: %d", SolventBlockHelpers::get1dIndex({ 3, 3, 3 }));
	//	solventblock.origo.print('o');
	//	solventblock.rel_pos[0].print('r');
	//}

	Float3 force(0.f);
	//const SolventCoord coord_self = thread_active ? *CoordArrayQueueHelpers::getSolventcoordPtr(box->solventcoordarray_circular_queue, box->step, solvent_index) : SolventCoord{};

	

	if (blockIdx.x == 0 && threadIdx.x == 0) {
		//coord.origo.print();
		//coord.rel_position.print();
	}

	// --------------------------------------------------------------- Molecule Interactions --------------------------------------------------------------- //
	for (int i = 0; i < box->n_compounds; i++) {
		//int n_compound_particles = box->compounds[i].n_particles;

		//// First all threads help loading the molecule
		//if (threadIdx.x < n_compound_particles) {
		//	// First load particles of neighboring compound
		//	auto* coordarray_ptr = CoordArrayQueueHelpers::getCoordarrayPtr(box->coordarray_circular_queue, box->step, blockIdx.x);
		//	compound_coords.loadData(*coordarray_ptr);
		//	LIMAPOSITIONSYSTEM::getRelativePositions(compound_coords.rel_positions, utility_buffer);

		//	// Then load atomtypes of neighboring compound
		//	utility_buffer_small[threadIdx.x] = box->compounds[i].atom_types[threadIdx.x];
		//}
		//__syncthreads();

	//  We can optimize here by loading and calculate the paired sigma and eps, jsut remember to loop threads, if there are many aomttypes.


		// Fuck me this is tricky. If we are too far away, we can complete skip� this calculation i guess?
		// Otherwise, we need to first compute a hyperpos, being careful not to change our original position.
		// 
		//Float3 solvent_pos{coord.}
		if (thread_active) {
			//EngineUtils::applyHyperpos(&utility_buffer[0], &solvent.pos);									// Move own particle in relation to compound-key-position
			//force += computeCompoundToSolventLJForces(&solvent.pos, n_compound_particles, utility_buffer, data_ptr, &potE_sum, ATOMTYPE_SOL, utility_buffer_small);
		}
		__syncthreads();
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //



	//// --------------------------------------------------------------- Solvent Interactions --------------------------------------------------------------- //
	if (thread_active) {
		//const auto* solventcoords_ptr = CoordArrayQueueHelpers::getSolventcoordPtr(box->solventcoordarray_circular_queue, box->step, 0);
		//force += computeSolventToSolventLJForces(coord_self, box->solvent_neighborlists[solvent_index], solventcoords_ptr, data_ptr, &potE_sum);
	}
	//// ----------------------------------------------------------------------------------------------------------------------------------------------------- //

	force = Float3{};

	//LogSolventData(box, solvent_index, potE_sum, coord_self);

	//SolventCoord coord_self_next{ coord_self };
	Coord relpos_next{ solventblock.rel_pos[threadIdx.x] };
	if (solvent_active) {

		//Coord relpos_prev = SolventBlockHelpers::isFirstStepAfterTransfer(box->step)
		//	? box->transfermodule_array[blockIdx.x].remain_queue.rel_positions_prev[threadIdx.x]
		//	: LIMAPOSITIONSYSTEM::getRelposPrev(box->solventblockgrid_circurlar_queue, blockIdx.x, box->step);

		Coord relpos_prev = LIMAPOSITIONSYSTEM::getRelposPrev(box->solventblockgrid_circurlar_queue, blockIdx.x, box->step);
		

		//// Make the above const, after we get this to work!
		if (box->step == 0 && blockIdx.x == 0 && threadIdx.x == 0) {
			relpos_prev.x -= 1000000;
			relpos_prev.z -= 1000000;
		}

		relpos_next = integratePosition(solventblock.rel_pos[threadIdx.x], relpos_prev, &force, solvent_mass, box->dt, box->thermostat_scalar);




		//-------

		//int step_prev = box->step == 0 ? 0 : box->step - 1;
		//const SolventCoord* coord_tsub1 = CoordArrayQueueHelpers::getSolventcoordPtr(box->solventcoordarray_circular_queue, step_prev, solvent_index);
		//Coord rel_hypercoord_tsub1 = LIMAPOSITIONSYSTEM::getRelativeHyperposition(coord_self, *coord_tsub1);
				//// Make the above const, after we get this to work!

		////if (threadIdx.x == 0) {
		////	force.print('f');
		////}

		//if (box->step == 0 && blockIdx.x == 0 && threadIdx.x == 0) { 
		//	rel_hypercoord_tsub1.x -= 1000000; 
		//	rel_hypercoord_tsub1.z -= 1000000;
		//}


		//const Coord newRelPos = integratePosition(coord_self.rel_position, rel_hypercoord_tsub1, &force, solvent_mass, box->dt, box->thermostat_scalar);
		//coord_self_next.rel_position = newRelPos;


		//LIMAPOSITIONSYSTEM::updateSolventcoord(coord_self_next);
		//LIMAPOSITIONSYSTEM::applyPBC(coord_self_next);
	}



	// Push new SolventCoord to global mem
	if (thread_active) {
		//*CoordArrayQueueHelpers::getSolventcoordPtr(box->solventcoordarray_circular_queue, box->step + 1, solvent_index) = coord_self_next;
	}

	if (solvent_active) {
		auto solventblock_next_ptr = CoordArrayQueueHelpers::getSolventBlockPtr(box->solventblockgrid_circurlar_queue, box->step + 1, blockIdx.x);

		solventblock_next_ptr->rel_pos[threadIdx.x] = relpos_next;
		if (threadIdx.x == 0) {
			solventblock_next_ptr->n_solvents = solventblock.n_solvents;
		}

		//if (box->step % STEPS_PER_SOLVENTBLOCKTRANSFER == SOLVENTBLOCK_TRANSFERSTEP) {
		//	//EngineUtils::doSolventTransfer(relpos_next, solventblock.rel_pos[threadIdx.x], box->transfermodule_array);

		//	solventblock_next_ptr->rel_pos[threadIdx.x] = relpos_next;
		//	if (threadIdx.x == 0) {
		//		solventblock_next_ptr->n_solvents = solventblock.n_solvents;
		//	}

		//	//doSequentialForwarding(solventblock, transferqueues, relpos_next, box->transfermodule_array);
		//	//purgeTransfersAndCompressRemaining(solventblock, solventblock_next_ptr,	relpos_next, utility_buffer_small, &box->transfermodule_array[blockIdx.x]);
		//}
		//else {
		//	solventblock_next_ptr->rel_pos[threadIdx.x] = relpos_next;
		//	if (threadIdx.x == 0) {
		//		solventblock_next_ptr->n_solvents = solventblock.n_solvents;
		//	}
		//}
	}
}
#undef solvent_index
#undef thread_active		// ALSO remove this
#undef solvent_mass
#undef solvent_active
#undef solventblock_ptr




// This is run before step.inc(), but will always publish results to the first array in grid!
__global__ void solventTransferKernel(Box* box) {
	SolventBlockTransfermodule* transfermodule_ptr = &box->transfermodule_array[blockIdx.x];
	
	SolventBlock* solventblock_ptr = CoordArrayQueueHelpers::getSolventBlockPtr(box->solventblockgrid_circurlar_queue, box->step + 1, blockIdx.x);


	// Temp implementation:
	//solventblock_ptr->rel_pos[threadIdx.x] = transfermodule_ptr->remain_queue.rel_positions[threadIdx.x];
	solventblock_ptr->rel_pos[threadIdx.x] = transfermodule_ptr->remain_queue.rel_positions[threadIdx.x];
	//transfermodule_ptr->remain_queue.rel_positions[threadIdx.x] = Coord{1,0,1};
}





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
		particle_coords[particle_id_bridge] = CoordArrayQueueHelpers::getCoordarrayPtr(box->coordarray_circular_queue, box->step, p_ref.compound_id)->rel_positions[p_ref.local_id_compound];

		// If we are on the right side, we need to shift 
		if (p_ref.compound_id == bridge.compound_id_right) { particle_coords[particle_id_bridge] += utility_coord; }
	}
	__syncthreads();

	// Load the now shifted relative coords into float3 positions for force calcs.
	LIMAPOSITIONSYSTEM::getRelativePositions(particle_coords, positions);

	float potE_sum = 0;
	Float3 force{};

	// ------------------------------------------------------------ Intercompund Operations ------------------------------------------------------------ //
	{											// So for the very first step, these �should all be 0, but they are not??										TODO: Look into this at some point!!!! 
		force += computePairbondForces(&bridge, positions, utility_buffer, &potE_sum);
		force += computeAnglebondForces(&bridge, positions, utility_buffer, &potE_sum);
		force += computeDihedralForces(&bridge, positions, utility_buffer, &potE_sum);
	}
	__syncthreads();
	// --------------------------------------------------------------------------------------------------------------------------------------------------- //

	if (particle_id_bridge < bridge.n_particles) {
		if (blockIdx.x == 0 && threadIdx.x == 0) {
			//force.print('f');			
		}
		ParticleRefCompact* p_ref = &bridge.particle_refs[particle_id_bridge];
		box->compounds[p_ref->compound_id].forces[p_ref->local_id_compound] = force;
		//box->compounds[p_ref->compound_id].forces[p_ref->local_id_compound] = 0;

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