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


// Must be called with ALL threads in a block
__device__ Coord getRandomCoord(int lcg_seed) {
	Coord randcoord;
	for (int i = 0; i < blockDim.x; i++) {
		if (threadIdx.x == i) {
			randcoord = Coord{
				EngineUtils::genPseudoRandomNum(lcg_seed),
				EngineUtils::genPseudoRandomNum(lcg_seed),
				EngineUtils::genPseudoRandomNum(lcg_seed)
			};
		}
		__syncthreads();
	}
	return randcoord;
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

__device__ Float3 computeSolventToSolventLJForces(const Float3& relpos_self, const Float3* relpos_others, int n_elements, bool exclude_own_index, float* data_ptr, float& potE_sum) {	// Specific to solvent kernel
	Float3 force{};

	for (int i = 0; i < n_elements; i++) {

		// If computing within block, dont compute force against thread's solvent
		if (exclude_own_index && threadIdx.x == i) { continue; }

		force += LimaForcecalc::calcLJForce(&relpos_self, &relpos_others[i], data_ptr, &potE_sum,
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

__device__ inline void LogSolventData(Box* box, const float& potE, const SolventBlock& solventblock, bool solvent_active, const Float3& force) {
	if (solvent_active) {
		const uint32_t index = EngineUtils::getLoggingIndexOfParticle(box->step, box->total_particles_upperbound, box->n_compounds, blockIdx.x * blockDim.x + threadIdx.x);
		box->potE_buffer[index] = potE;
		box->traj_buffer[index] = SolventBlockHelpers::extractAbsolutePositionLM(solventblock);

		if (solventblock.ids[threadIdx.x] > 13000) printf("\nhiewr: %u\n", solventblock.ids[threadIdx.x]);
		const auto debug_index = (box->step * box->total_particles_upperbound + box->n_compounds * MAX_COMPOUND_PARTICLES + solventblock.ids[threadIdx.x]) * DEBUGDATAF3_NVARS;
		box->debugdataf3[debug_index] = Float3(500);// Float3(solventblock.ids[threadIdx.x] * 10 + 1.f, solventblock.ids[threadIdx.x] * 10 + 2.f, solventblock.ids[threadIdx.x] * 10 + 3.f);
	}
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

// SLOW - Returns sum of actives before, thus must be -1 for 0-based index :)
__device__ void doSequentialPrefixSum(uint8_t* onehot_remainers, int n_elements) {
	for (int i = 1; i < n_elements; i++) {
		if (threadIdx.x == i) {
			onehot_remainers[i] += onehot_remainers[i - 1];
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

	const uint8_t solventindex_new = utility_buffer[threadIdx.x] - 1; // Underflow here doesn't matter, as the underflowing threads wont remain anyways :)
	return solventindex_new;	
}







/// <summary></summary>
/// <param name="solventblock">In shared memory</param>
/// <param name="transferqueues">In shared memory</param>
/// <param name="relpos_next">Register</param>
/// <param name="transfermodules">Global memory</param>
__device__ void transferOut(const Coord& transfer_dir, const SolventBlock& solventblock_current_local, const int new_blockid, 
	const Coord& relpos_next, STransferQueue* transferqueues, SolventBlockTransfermodule* transfermodules, const Coord& blockId3d) {

	// Sequential insertion in shared memory
	for (int i = 0; i < solventblock_current_local.n_solvents; i++) {
		if (threadIdx.x == i && new_blockid != blockIdx.x) {	// Only handle non-remain solvents
			const int queue_index = SolventBlockTransfermodule::getQueueIndex(transfer_dir);

			if (queue_index < 0 || queue_index > 5) { printf("\nGot unexpected queue index %d\n", queue_index); transfer_dir.print(); }
			if (!transferqueues[queue_index].addElement(relpos_next, solventblock_current_local.rel_pos[threadIdx.x], solventblock_current_local.ids[threadIdx.x]))
			{
				printf("What the fuck\n");
			}
		}
		__syncthreads();
	}

	// Coaslescing copying to global memory
	for (int queue_index = 0; queue_index < 6; queue_index++) {
		const STransferQueue& queue_local = transferqueues[queue_index];

		if (threadIdx.x < queue_local.n_elements) {

			const Coord transferdir_queue = EngineUtils::getTransferDirection(queue_local.rel_positions[0]);		// Maybe use a utility-coord a precompute by thread0, or simply hardcode...
			const int blockid_global = EngineUtils::getNewBlockId(transferdir_queue, blockId3d);
			if (blockid_global < 0 || blockid_global >= SolventBlockGrid::blocks_total) { printf("\nGot unexpected Block id index %d\n"); }
			STransferQueue* queue_global = &transfermodules[blockid_global].transfer_queues[queue_index];

			queue_global->fastInsert(
				queue_local.rel_positions[threadIdx.x], 
				queue_local.rel_positions_prev[threadIdx.x], 
				queue_local.ids[threadIdx.x], 
				transferdir_queue);




			// Debugging
			if (queue_global->rel_positions[threadIdx.x].x < -2 * static_cast<int32_t>(NANO_TO_LIMA) || queue_global->rel_positions[threadIdx.x].x > 2 * static_cast<int32_t>(NANO_TO_LIMA)
				|| queue_global->rel_positions[threadIdx.x].y < -2 * static_cast<int32_t>(NANO_TO_LIMA) || queue_global->rel_positions[threadIdx.x].y > 2 * static_cast<int32_t>(NANO_TO_LIMA)
				|| queue_global->rel_positions[threadIdx.x].z < -2 * static_cast<int32_t>(NANO_TO_LIMA) || queue_global->rel_positions[threadIdx.x].z > 2 * static_cast<int32_t>(NANO_TO_LIMA)
				) {
				printf("\n");
				transferdir_queue.print('t');
				queue_local.rel_positions[threadIdx.x].print('q');
				queue_global->rel_positions[threadIdx.x].print('Q');
			}

			if (threadIdx.x == 0) {
				if (queue_global->n_elements != 0) {
					printf("\nN elements was: %d in queue %d\n", queue_global->n_elements, queue_index);
					transferdir_queue.print('d');
				}

				queue_global->n_elements = queue_local.n_elements;
				if (queue_local.n_elements > 15) {
					printf("\nTransferring %d elements\n", queue_local.n_elements);
				}
			}
			



			// Only set n_elements if we get here, meaning atleast 1 new element. Otherwise it will just remain 0
			if (threadIdx.x == 0) { queue_global->n_elements = queue_local.n_elements; }
		}
	}
	__syncthreads();
}


/// <summary>
/// Must be run AFTER transferOut, as it erases all information about the transferring solvents
/// </summary>
/// <param name="solventblock_current">Solventblock in shared memory</param>
/// <param name="solventblock_next">Solventblock belong to cudablock at next step in global memory</param>
/// <param name="relpos_next"></param>
/// <param name="utility_buffer">Buffer of min size MAX_SOLVENTS_IN_BLOCK, maybe more for computing prefix sum</param>
/// <param name="remain_transfermodule">Transfermodule belonging to cudablock</param>
__device__ void compressRemainers(const SolventBlock& solventblock_current_local, SolventBlock* solventblock_next_global,
	const Coord& relpos_next, uint8_t* utility_buffer, SolventBlockTransfermodule* remain_transfermodule, const bool remain) {

	// Compute prefix sum to find new index of solvent belonging to thread
	const uint8_t solventindex_new = computePrefixSum(remain, utility_buffer, solventblock_current_local.n_solvents);


	if (remain) {
		// Send current pos at threadindex to prevpos at the new index
		remain_transfermodule->remain_relpos_prev[solventindex_new] = solventblock_current_local.rel_pos[threadIdx.x];
		// Send the next pos 
		solventblock_next_global->rel_pos[solventindex_new] = relpos_next;
		solventblock_next_global->ids[solventindex_new] = solventblock_current_local.ids[threadIdx.x];
	}

	const int nsolventsinblock_next = __syncthreads_count(remain);
	if (threadIdx.x == 0) {
		remain_transfermodule->n_remain = nsolventsinblock_next;
		solventblock_next_global->n_solvents = nsolventsinblock_next;	// Doesn't matter, since the transfer kernel handles this. Enabled for debugging now..
	}

	if (threadIdx.x >= nsolventsinblock_next) {
		solventblock_next_global->rel_pos[threadIdx.x] = Coord{42, 42, 42};
	}

}

__device__ void transferOutAndCompressRemainders(const SolventBlock& solventblock_current_local, SolventBlock* solventblock_next_global,
	const Coord& relpos_next, uint8_t* utility_buffer, SolventBlockTransfermodule* transfermodules_global, STransferQueue* transferqueues_local) {

	const Coord blockId3d = SolventBlockHelpers::get3dIndex(blockIdx.x);
	const Coord transfer_dir = threadIdx.x < solventblock_current_local.n_solvents ? EngineUtils::getTransferDirection(relpos_next) : Coord{ 0 };
	const int new_blockid = EngineUtils::getNewBlockId(transfer_dir, blockId3d);
	const bool remain = (blockIdx.x == new_blockid) && threadIdx.x < solventblock_current_local.n_solvents;

	
	transferOut(transfer_dir, solventblock_current_local, new_blockid, relpos_next, transferqueues_local, transfermodules_global, blockId3d);

	SolventBlockTransfermodule* remain_transfermodule = &transfermodules_global[blockIdx.x];
	compressRemainers(solventblock_current_local, solventblock_next_global, relpos_next, utility_buffer, remain_transfermodule, remain);
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





#include "cuda/std/cmath"
#include "math.h"


#include <curand.h>
#include <cuda.h>
#include <cuda_runtime.h>

#define solvent_active (threadIdx.x < solventblock.n_solvents)
#define solvent_mass (forcefield_device.particle_parameters[ATOMTYPE_SOL].mass)
#define solventblock_ptr (CoordArrayQueueHelpers::getSolventBlockPtr(box->solventblockgrid_circular_queue, box->step, blockIdx.x))
static_assert(MAX_SOLVENTS_IN_BLOCK > MAX_COMPOUND_PARTICLES, "solventForceKernel was about to reserve an insufficient amount of memory");
__global__ void solventForceKernel(Box* box) {
	__shared__ Float3 utility_buffer[MAX_SOLVENTS_IN_BLOCK];
	//__shared__ uint8_t utility_buffer_small[MAX_COMPOUND_PARTICLES];
	__shared__ uint8_t utility_buffer_small[MAX_SOLVENTS_IN_BLOCK];
	__shared__ CompoundCoords compound_coords;	// For neighboring compounds
	__shared__ SolventBlock solventblock;
	__shared__ SolventTransferqueue<SolventBlockTransfermodule::max_queue_size> transferqueues[6];		// TODO: Use template to make identical kernel, so the kernel with transfer is slower and larger, and the rest remain fast!!!!
	//__shared__ Coord coord_utility_buffer[MAX_SOLVENTS_IN_BLOCK + 6 * SolventBlockTransfermodule::max_queue_size];

	__shared__ int lcg_seed;
	lcg_seed = 12345;

	// Doubles as block_index_3d!
	const Coord block_origo = SolventBlockHelpers::get3dIndex(blockIdx.x);


	// Init queue, otherwise it will contain wierd values
	if (threadIdx.x < 6) { 
		transferqueues[threadIdx.x] = SolventTransferqueue<SolventBlockTransfermodule::max_queue_size>{};
	}



	float potE_sum = 0;
	float data_ptr[4];	// Pot, force, closest particle, ?
	for (int i = 0; i < 4; i++)
		data_ptr[i] = 0;
	data_ptr[2] = 9999.f;

	// temp
	utility_buffer[threadIdx.x] = Float3{0};
	utility_buffer_small[threadIdx.x] = 0;


	if (threadIdx.x == 0) {
		solventblock.loadMeta(*solventblock_ptr);
		//solventblock.origo.print();
	}
	__syncthreads();
	solventblock.loadData(*solventblock_ptr);
	__syncthreads();


	Float3 force(0.f);
	const Float3 relpos_self = solventblock.rel_pos[threadIdx.x].toFloat3();

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


		// Fuck me this is tricky. If we are too far away, we can complete skipå this calculation i guess?
		// Otherwise, we need to first compute a hyperpos, being careful not to change our original position.
		// 
		//Float3 solvent_pos{coord.}
		//if (thread_active) {
		//	//EngineUtils::applyHyperpos(&utility_buffer[0], &solvent.pos);									// Move own particle in relation to compound-key-position
		//	//force += computeCompoundToSolventLJForces(&solvent.pos, n_compound_particles, utility_buffer, data_ptr, &potE_sum, ATOMTYPE_SOL, utility_buffer_small);
		//}
		__syncthreads();
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //



	// --------------------------------------------------------------- Intrablock Solvent Interactions ----------------------------------------------------- //
	if (solvent_active) {
		utility_buffer[threadIdx.x] = solventblock.rel_pos[threadIdx.x].toFloat3();
	}
	__syncthreads();
	if (solvent_active) {		
		force += computeSolventToSolventLJForces(relpos_self, utility_buffer, solventblock.n_solvents, true, data_ptr, potE_sum);
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //

	// --------------------------------------------------------------- Interblock Solvent Interactions ----------------------------------------------------- //
	const int query_range = 3;
	for (int x = -query_range; x <= query_range; x++) {
		for (int y = -query_range; y <= query_range; y++) {
			for (int z = -query_range; z <= query_range; z++) {
				const Coord dir{ x,y,z };
				if (dir.isZero()) { continue; }

				const int blockindex_neighbor = EngineUtils::getNewBlockId(dir, block_origo);
				if (blockindex_neighbor < 0 || blockindex_neighbor >= SolventBlockGrid::blocks_total) { printf("\nWhat the fuck\n"); }

				const SolventBlock* solventblock_neighbor = CoordArrayQueueHelpers::getSolventBlockPtr(box->solventblockgrid_circular_queue, box->step, blockindex_neighbor);
				const int nsolvents_neighbor = solventblock_neighbor->n_solvents;

				if (threadIdx.x < nsolvents_neighbor) {
					utility_buffer[threadIdx.x] = (solventblock_neighbor->rel_pos[threadIdx.x] + (dir * static_cast<int32_t>(NANO_TO_LIMA))).toFloat3();
				}
				__syncthreads();

				if (solvent_active) {
					force += computeSolventToSolventLJForces(relpos_self, utility_buffer, nsolvents_neighbor, false, data_ptr, potE_sum);
				}
			}
		}
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //
	//if (box->step > 200)
	//force = Float3{};

	LogSolventData(box, potE_sum, solventblock, solvent_active, force);

	const Coord randcoord = getRandomCoord(lcg_seed);

	Coord relpos_next{};
	if (solvent_active) {
		Coord relpos_prev = LIMAPOSITIONSYSTEM::getRelposPrev(box->solventblockgrid_circular_queue, blockIdx.x, box->step);

		//// Make the above const, after we get this to work!
		const int scalar = 100000;
		//if (box->step == 0) { relpos_prev -= randcoord * scalar; }

		relpos_next = integratePosition(solventblock.rel_pos[threadIdx.x], relpos_prev, &force, solvent_mass, box->dt, box->thermostat_scalar);
		auto dif = (relpos_next - relpos_prev);
		if (std::abs(dif.x) > 15000000 || std::abs(dif.y) > 15000000 || std::abs(dif.z) > 15000000) { 
			printf("\nFuck mee %d\n", solventblock.ids[threadIdx.x]);
			dif.printS('D'); 
			force.print('F');
			relpos_prev.printS('p');
			relpos_next.printS('n');
			box->critical_error_encountered = true;
		}
	}



	// Push new SolventCoord to global mem
	auto solventblock_next_ptr = CoordArrayQueueHelpers::getSolventBlockPtr(box->solventblockgrid_circular_queue, box->step + 1, blockIdx.x);

	if (SolventBlockHelpers::isTransferStep(box->step)) {
		transferOutAndCompressRemainders(solventblock, solventblock_next_ptr, relpos_next, utility_buffer_small, box->transfermodule_array, transferqueues);
	}
	else {
		solventblock_next_ptr->rel_pos[threadIdx.x] = relpos_next;
		solventblock_next_ptr->ids[threadIdx.x] = solventblock.ids[threadIdx.x];
		if (threadIdx.x == 0) {
			solventblock_next_ptr->n_solvents = solventblock.n_solvents;
		}
	}
	
}
#undef solvent_index
#undef solvent_mass
#undef solvent_active
#undef solventblock_ptr




// This is run before step.inc(), but will always publish results to the first array in grid!
__global__ void solventTransferKernel(Box* box) {
	SolventBlockTransfermodule* transfermodule = &box->transfermodule_array[blockIdx.x];
	
	SolventBlock* solventblock_current = CoordArrayQueueHelpers::getSolventBlockPtr(box->solventblockgrid_circular_queue, box->step, blockIdx.x);
	SolventBlock* solventblock_next = CoordArrayQueueHelpers::getSolventBlockPtr(box->solventblockgrid_circular_queue, box->step + 1, blockIdx.x);

	// First overload what will become posrel_prev for the next step. This is needed since compaction has already been applied
	for (int index = threadIdx.x; index < transfermodule->n_remain; index += blockDim.x) {
		if (index > 255) { printf("\nWhat the fuck\n"); }
		solventblock_current->rel_pos[index] = transfermodule->remain_relpos_prev[index];
	}

	if (solventblock_next->n_solvents != transfermodule->n_remain) {
		printf("YOooooooooooooooooooooooooooooooooooooooooooo %d %d\n", solventblock_next->n_solvents, transfermodule->n_remain);
	}

	// Handling incoming transferring solvents
	int n_solvents_next = transfermodule->n_remain;
	for (int queue_index = 0; queue_index < SolventBlockTransfermodule::n_queues; queue_index++) {
		auto* queue = &transfermodule->transfer_queues[queue_index];
		if (threadIdx.x < queue->n_elements) {
			const int incoming_index = n_solvents_next + threadIdx.x;

			if (queue->rel_positions[threadIdx.x].x == 0) { queue->rel_positions[threadIdx.x].print('I'); }

			solventblock_next->rel_pos[incoming_index] = queue->rel_positions[threadIdx.x];
			solventblock_next->ids[incoming_index] = queue->ids[threadIdx.x];
			solventblock_current->rel_pos[incoming_index] = queue->rel_positions_prev[threadIdx.x];
		}
		n_solvents_next += queue->n_elements;

		// Signal that all elements of the queues have been moved
		__syncthreads();
		if (threadIdx.x == 0) {
			queue->n_elements = 0;
		}
	}

	// Finally update the solventblock_next with how many solvents it now contains
	if (threadIdx.x == 0) {
		solventblock_next->n_solvents = n_solvents_next;
	}
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
	{											// So for the very first step, these ´should all be 0, but they are not??										TODO: Look into this at some point!!!! 
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

#pragma warning (pop)