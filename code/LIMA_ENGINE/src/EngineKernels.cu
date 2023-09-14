
#include "Engine.cuh"
#include "ForceComputations.cuh"
#include "Utilities.h"
#include "KernelWarnings.cuh"

#pragma warning(push)
#pragma warning(disable: E0020)
#pragma warning(push)
#pragma warning(disable: 20054)

#pragma diag_suppress 20054

// ----------------------------------------------------------------------------------- FILE-SPECIFIC FORCEFIELD -------------------------------------------------------------------------------------------//

// Pre-calculate a solvent-X paired forcefield, to save ALOT of calc in kernels
__constant__ ForceField_NB forcefield_device;




void Engine::setDeviceConstantMemory() {
	//const int forcefield_bytes = sizeof(ForceField_NB);
	cudaMemcpyToSymbol(forcefield_device, &forcefield_host, sizeof(ForceField_NB), 0, cudaMemcpyHostToDevice);	// So there should not be a & before the device __constant__
	cudaDeviceSynchronize();
	LIMA_UTILS::genericErrorCheck("Error while moving forcefield to device\n");
}

__device__ inline float calcSigma(uint8_t atomtype1, uint8_t atomtype2) {
	return (forcefield_device.particle_parameters[atomtype1].sigma + forcefield_device.particle_parameters[atomtype2].sigma) * 0.5;
}
__device__ inline float calcEpsilon(uint8_t atomtype1, uint8_t atomtype2) {
	return sqrtf(forcefield_device.particle_parameters[atomtype1].epsilon * forcefield_device.particle_parameters[atomtype2].epsilon);
}


// ------------------------------------------------------------------------------------------- DEVICE FUNCTIONS -------------------------------------------------------------------------------------------//

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

// Two variants of this exists, with and without lut
__device__ Float3 computeIntercompoundLJForces(const Float3& self_pos, uint8_t atomtype_self, float& potE_sum, uint32_t global_id_self, float* data_ptr,
	Compound* neighbor_compound, Float3* neighbor_positions, int neighborcompound_id, BondedParticlesLUT& bonded_particles_lut) {
	Float3 force(0.f);
	
	for (int neighborparticle_id = 0; neighborparticle_id < neighbor_compound->n_particles; neighborparticle_id++) {

		// If thread's assoc. particle is bonded to the particle in neighborcompound, continue
		if (bonded_particles_lut.get(threadIdx.x, neighborparticle_id)) { continue; }

		const int neighborparticle_atomtype = neighbor_compound->atom_types[neighborparticle_id];	//// TEMPORARY, this is waaaay to many global mem accesses

		force += LimaForcecalc::calcLJForceOptim(self_pos, neighbor_positions[neighborparticle_id], data_ptr, potE_sum,
			calcSigma(atomtype_self, neighborparticle_atomtype), calcEpsilon(atomtype_self, neighborparticle_atomtype),
			LimaForcecalc::CalcLJOrigin::ComComInter,
			global_id_self, neighbor_compound->particle_global_ids[neighborparticle_id]
		);
	}
	return force * 24.f;
}

__device__ Float3 computeIntercompoundLJForces(const Float3& self_pos, uint8_t atomtype_self, float& potE_sum, uint32_t global_id_self, float* data_ptr,
	Compound* neighbor_compound, Float3* neighbor_positions, int neighborcompound_id) {
	Float3 force(0.f);

	for (int neighborparticle_id = 0; neighborparticle_id < neighbor_compound->n_particles; neighborparticle_id++) {
		const int neighborparticle_atomtype = neighbor_compound->atom_types[neighborparticle_id];	//// TEMPORARY, this is waaaay to many global mem accesses

		force += LimaForcecalc::calcLJForceOptim(self_pos, neighbor_positions[neighborparticle_id], data_ptr, potE_sum,
			calcSigma(atomtype_self, neighborparticle_atomtype), calcEpsilon(atomtype_self, neighborparticle_atomtype),
			LimaForcecalc::CalcLJOrigin::ComComInter,
			global_id_self, neighbor_compound->particle_global_ids[neighborparticle_id]
		);
	}
	return force * 24.f;
}

__device__ Float3 computeIntracompoundLJForces(Compound* compound, CompoundState* compound_state, float& potE_sum, float* data_ptr, BondedParticlesLUT* bonded_particles_lut) {
	Float3 force(0.f);
	if (threadIdx.x >= compound->n_particles) { return force; }
	for (int i = 0; i < compound->n_particles; i++) {

		// Skip if particle is self or bonded
		if (i == threadIdx.x || (bonded_particles_lut->get(threadIdx.x, i))) { continue; }

		force += LimaForcecalc::calcLJForceOptim(compound_state->positions[threadIdx.x], compound_state->positions[i], data_ptr, potE_sum,
			calcSigma(compound->atom_types[threadIdx.x], compound->atom_types[i]),
			calcEpsilon(compound->atom_types[threadIdx.x], compound->atom_types[i]),
			LimaForcecalc::CalcLJOrigin::ComComIntra,
			compound->particle_global_ids[threadIdx.x], compound->particle_global_ids[i]			
		);
	}
	return force * 24.f;
}

__device__ Float3 computeSolventToSolventLJForces(const Float3& relpos_self, const Float3* relpos_others, int n_elements, bool exclude_own_index, float* data_ptr, float& potE_sum) {	// Specific to solvent kernel
	Float3 force{};

	for (int i = 0; i < n_elements; i++) {

		// If computing within block, dont compute force against thread's solvent
		if (exclude_own_index && threadIdx.x == i) { continue; }

		force += LimaForcecalc::calcLJForceOptim(relpos_self, relpos_others[i], data_ptr, potE_sum,
			forcefield_device.particle_parameters[ATOMTYPE_SOLVENT].sigma,
			forcefield_device.particle_parameters[ATOMTYPE_SOLVENT].epsilon,
			exclude_own_index ? LimaForcecalc::CalcLJOrigin::SolSolIntra : LimaForcecalc::CalcLJOrigin::SolSolInter,
			threadIdx.x, i
		);
	}
	return force * 24.f;
}

__device__ Float3 computeSolventToCompoundLJForces(const Float3& self_pos, const int n_particles, Float3* positions, float* data_ptr, float& potE_sum, const uint8_t atomtype_self) {	// Specific to solvent kernel
	Float3 force{};
	for (int i = 0; i < n_particles; i++) {
		force += LimaForcecalc::calcLJForceOptim(self_pos, positions[i], data_ptr, potE_sum,
			calcSigma(atomtype_self, ATOMTYPE_SOLVENT), 
			calcEpsilon(atomtype_self, ATOMTYPE_SOLVENT),
			LimaForcecalc::CalcLJOrigin::SolCom,
			atomtype_self, ATOMTYPE_SOLVENT
		);
	}
	return force * 24.f;
}

__device__ Float3 computeCompoundToSolventLJForces(const Float3& self_pos, const int n_particles, const Float3* positions, float* data_ptr, float& potE_sum, const uint8_t* atomtypes_others, const int sol_id) {	// Assumes all positions are 
	Float3 force(0.f);
	for (int i = 0; i < n_particles; i++) {
		force += LimaForcecalc::calcLJForceOptim(self_pos, positions[i], data_ptr, potE_sum,
			calcSigma(ATOMTYPE_SOLVENT, atomtypes_others[i]),
			calcEpsilon(ATOMTYPE_SOLVENT, atomtypes_others[i]),
			LimaForcecalc::CalcLJOrigin::ComSol,
			sol_id, -1
		);
	}
	return force * 24.f;
}

__device__ Float3 computeSinglebondForces(SingleBond* singlebonds, int n_singlebonds, Float3* positions, Float3* utility_buffer, float* utilitybuffer_f, float* potE) {	// only works if n threads >= n bonds

	// First clear the buffer which will store the forces.
	utility_buffer[threadIdx.x] = Float3(0.f);
	utilitybuffer_f[threadIdx.x] = 0.f;
	__syncthreads();

	for (int bond_offset = 0; (bond_offset * blockDim.x) < n_singlebonds; bond_offset++) {
		SingleBond* pb = nullptr;
		Float3 forces[2] = { Float3{}, Float3{}};
		float potential = 0.f;
		const int bond_index = threadIdx.x + bond_offset * blockDim.x;
		
		if (bond_index < n_singlebonds) {
			pb = &singlebonds[bond_index];

			LimaForcecalc::calcSinglebondForces(
				positions[pb->atom_indexes[0]],
				positions[pb->atom_indexes[1]],
				*pb,
				forces, 
				potential
			);
		}

		for (int i = 0; i < blockDim.x; i++) {
			if (threadIdx.x == i && pb != nullptr) {
				for (int i = 0; i < 2; i++) {
					utility_buffer[pb->atom_indexes[i]] += forces[i];
					utilitybuffer_f[pb->atom_indexes[i]] += potential * 0.5f;
				}				
			}
			__syncthreads();
		}
	}

	*potE += utilitybuffer_f[threadIdx.x];
	return utility_buffer[threadIdx.x];
}

template <typename T>	// Can either be Compound or CompoundBridge
__device__ Float3 computeAnglebondForces(T* entity, Float3* positions, Float3* forces_interim, float* potentials_interim, float* potE) {
	
	// First clear the buffer which will store the forces.
	forces_interim[threadIdx.x] = Float3(0.f);
	potentials_interim[threadIdx.x] = 0.f;
	__syncthreads();

	for (int bond_offset = 0; (bond_offset * blockDim.x) < entity->n_anglebonds; bond_offset++) {
		AngleBond* ab = nullptr;
		Float3 forces[3] = { Float3{}, Float3{}, Float3{} };
		float potential = 0.f;
		const int bond_index = threadIdx.x + bond_offset * blockDim.x;

		if (bond_index < entity->n_anglebonds) {
			ab = &entity->anglebonds[bond_index];

			LimaForcecalc::calcAnglebondForces(
				positions[ab->atom_indexes[0]],
				positions[ab->atom_indexes[1]],
				positions[ab->atom_indexes[2]],
				*ab,
				forces, 
				potential
			);
		}


		for (int i = 0; i < blockDim.x; i++) {
			if (threadIdx.x == i && ab != nullptr) {
				for (int i = 0; i < ab->n_atoms; i++) {
					forces_interim[ab->atom_indexes[i]] += forces[i];
					potentials_interim[ab->atom_indexes[i]] += potential / 3.f;
				}
			}
			__syncthreads();
		}
	}

	*potE += potentials_interim[threadIdx.x];
	return forces_interim[threadIdx.x];
}

template <typename T>	// Can either be Compound or CompoundBridge
__device__ Float3 computeDihedralForces(T* entity, Float3* positions, Float3* forces_interim, float* potentials_interim, float* potE) {
	
	// First clear the buffer which will store the forces.
	forces_interim[threadIdx.x] = Float3(0.f);
	potentials_interim[threadIdx.x] = 0.f;
	__syncthreads();

	for (int bond_offset = 0; (bond_offset * blockDim.x) < entity->n_dihedrals; bond_offset++) {
		DihedralBond* db = nullptr;
		Float3 forces[4] = { Float3{}, Float3{}, Float3{}, Float3{}};
		float potential = 0.f;
		const int bond_index = threadIdx.x + bond_offset * blockDim.x;

		if (bond_index < entity->n_dihedrals) {
			db = &entity->dihedrals[bond_index];
			LimaForcecalc::calcDihedralbondForces(
				positions[db->atom_indexes[0]] / NANO_TO_LIMA,
				positions[db->atom_indexes[1]] / NANO_TO_LIMA,
				positions[db->atom_indexes[2]] / NANO_TO_LIMA,
				positions[db->atom_indexes[3]] / NANO_TO_LIMA,
				*db,
				forces,
				potential
			);
		}

		for (int i = 0; i < blockDim.x; i++) {
			if (threadIdx.x == i && db != nullptr) {
				for (int i = 0; i < 4; i++) {
					forces_interim[db->atom_indexes[i]] += forces[i];
					potentials_interim[db->atom_indexes[i]] += potential * 0.25f;
				}						
			}	
			__syncthreads();
		}
	}

	*potE += potentials_interim[threadIdx.x];
	return forces_interim[threadIdx.x];
}

__device__ Float3 computeImproperdihedralForces(ImproperDihedralBond* impropers, int n_impropers, Float3* positions, Float3* utility_buffer, float* potentials_interim, float* potE) {

	// First clear the buffer which will store the forces.
	utility_buffer[threadIdx.x] = Float3(0.f);
	potentials_interim[threadIdx.x] = 0.f;
	__syncthreads();

	for (int bond_offset = 0; (bond_offset * blockDim.x) < n_impropers; bond_offset++) {
		ImproperDihedralBond* db = nullptr;
		Float3 forces[4] = { Float3{}, Float3{}, Float3{}, Float3{} };
		float potential = 0.f;
		const int bond_index = threadIdx.x + bond_offset * blockDim.x;

		if (bond_index < n_impropers) {
			db = &impropers[bond_index];
			//printf("Firing %d of %d\n", bond_index, entity);
			LimaForcecalc::calcImproperdihedralbondForces(
				positions[db->atom_indexes[0]] / NANO_TO_LIMA,
				positions[db->atom_indexes[1]] / NANO_TO_LIMA,
				positions[db->atom_indexes[2]] / NANO_TO_LIMA,
				positions[db->atom_indexes[3]] / NANO_TO_LIMA,
				*db,
				forces,
				potential
			);
		}

		for (int i = 0; i < blockDim.x; i++) {
			if (threadIdx.x == i && db != nullptr) {
				for (int i = 0; i < db->n_atoms; i++) {
					utility_buffer[db->atom_indexes[i]] += forces[i];
					potentials_interim[db->atom_indexes[i]] += potential * 0.25f;
				}
			}
			__syncthreads();
		}
	}
	*potE += potentials_interim[threadIdx.x];
	return utility_buffer[threadIdx.x];
}


__device__ void getCompoundHyperpositionsAsFloat3(const NodeIndex& origo_self, const CompoundCoords* querycompound, Float3* output_buffer, Float3* utility_coord) { 
	if (threadIdx.x == 0) {
		const NodeIndex querycompound_hyperorigo = LIMAPOSITIONSYSTEM::getHyperNodeIndex(origo_self, querycompound->origo);

		KernelHelpersWarnings::assertHyperorigoIsValid(querycompound_hyperorigo, origo_self);

		// calc Relative LimaPosition Shift from the origo-shift
		*utility_coord = LIMAPOSITIONSYSTEM_HACK::getRelShiftFromOrigoShift(querycompound_hyperorigo, origo_self).toFloat3();
	}
	__syncthreads();

	// Eventually i could make it so i only copy the active particles in the compound
	if (threadIdx.x < MAX_COMPOUND_PARTICLES) {
		const Coord queryparticle_coord = querycompound->rel_positions[threadIdx.x];// +*utility_coord;
		output_buffer[threadIdx.x] = queryparticle_coord.toFloat3() + *utility_coord;
	}
	__syncthreads();
}



/// <summary>
/// </summary>
/// <param name="onehot_remainers">Must be size of MAX_SOLVENTS_IN_BLOCK</param>
/// <param name="utility">Must be large enough to store temporary sums for blelloch sum</param>
/// <returns></returns>
__device__ void doBlellochPrefixSum(uint8_t* onehot_remainers, uint8_t* utility) {
	// Forward scan
	//for (int leap = 1; leap < MAX_SOLVENTS_IN_BLOCK - 1; leap *= 2) {
	//	if (threadIdx.x % leap == 0) {
	//		int index = threadIdx.x + leap;
	//	}
	//}
}

// SLOW - Returns sum of actives before, thus must be -1 for 0-based index :)
__device__ void doSequentialPrefixSum(uint8_t* onehot_remainers, int n_elements) {
	for (int i = 1; i < n_elements; i++) {
		if (threadIdx.x == i) {
			onehot_remainers[i] += onehot_remainers[i - 1];
			KernelHelpersWarnings::verifyOnehotRemaindersIsValid(onehot_remainers, i);
		}
		__syncthreads();
	}
}

__device__ uint8_t computePrefixSum(const bool remain, uint8_t* utility_buffer, int n_elements) {
	utility_buffer[threadIdx.x] = static_cast<uint8_t>(remain);
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
__device__ void transferOut(const NodeIndex& transfer_dir, const SolventBlock& solventblock_current_local, const int new_blockid, 
	const Coord& relpos_next, STransferQueue* transferqueues, SolventBlockTransfermodule* transfermodules, const NodeIndex& blockId3d) {

	// Sequential insertion in shared memory
	for (int i = 0; i < solventblock_current_local.n_solvents; i++) {
		if (threadIdx.x == i && new_blockid != blockIdx.x) {	// Only handle non-remain solvents
			const int queue_index = SolventBlockTransfermodule::getQueueIndex(transfer_dir);

			KernelHelpersWarnings::transferoutVerifyQueueIndex(queue_index, transfer_dir);

			const bool insertionSuccess = transferqueues[queue_index].addElement(relpos_next, solventblock_current_local.rel_pos[threadIdx.x], solventblock_current_local.ids[threadIdx.x]);
			KernelHelpersWarnings::transferoutVerifyInsertion(insertionSuccess);
		}
		__syncthreads();
	}

	// Coaslescing copying to global memory
	for (int queue_index = 0; queue_index < 6; queue_index++) {
		const STransferQueue& queue_local = transferqueues[queue_index];

		if (threadIdx.x < queue_local.n_elements) {

			const NodeIndex transferdir_queue = LIMAPOSITIONSYSTEM::getTransferDirection(queue_local.rel_positions[0]);		// Maybe use a utility-coord a precompute by thread0, or simply hardcode...
			const int blockid_global = EngineUtils::getNewBlockId(transferdir_queue, blockId3d);
			KernelHelpersWarnings::assertValidBlockId(blockid_global);

			STransferQueue* queue_global = &transfermodules[blockid_global].transfer_queues[queue_index];

			queue_global->fastInsert(
				queue_local.rel_positions[threadIdx.x] - LIMAPOSITIONSYSTEM::nodeIndexToCoord(transferdir_queue),
				queue_local.ids[threadIdx.x]);

			KernelHelpersWarnings::transferOutDebug(queue_global, queue_local, transferdir_queue, queue_index);

			queue_global->n_elements = queue_local.n_elements;

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
		// solventindex_new is only valid for those who remain, the rest *might* have an index of -1
		solventblock_next_global->rel_pos[solventindex_new] = relpos_next;
		solventblock_next_global->ids[solventindex_new] = solventblock_current_local.ids[threadIdx.x];
	}


	const int nsolventsinblock_next = __syncthreads_count(remain);
	if (threadIdx.x == 0) {
		remain_transfermodule->n_remain = nsolventsinblock_next;
		solventblock_next_global->n_solvents = nsolventsinblock_next;	// Doesn't matter, since the transfer kernel handles this. Enabled for debugging now..
	}
}

__device__ void transferOutAndCompressRemainders(const SolventBlock& solventblock_current_local, SolventBlock* solventblock_next_global,
	const Coord& relpos_next, uint8_t* utility_buffer, SolventBlockTransfermodule* transfermodules_global, STransferQueue* transferqueues_local) {

	const NodeIndex blockId3d = SolventBlocksCircularQueue::get3dIndex(blockIdx.x);
	const NodeIndex transfer_dir = threadIdx.x < solventblock_current_local.n_solvents ? LIMAPOSITIONSYSTEM::getTransferDirection(relpos_next) : NodeIndex{};
	const int new_blockid = EngineUtils::getNewBlockId(transfer_dir, blockId3d);

	const bool remain = (blockIdx.x == new_blockid) && threadIdx.x < solventblock_current_local.n_solvents;

	
	transferOut(transfer_dir, solventblock_current_local, new_blockid, relpos_next, transferqueues_local, transfermodules_global, blockId3d);

	SolventBlockTransfermodule* remain_transfermodule = &transfermodules_global[blockIdx.x];
	compressRemainers(solventblock_current_local, solventblock_next_global, relpos_next, utility_buffer, remain_transfermodule, remain);
}






// ------------------------------------------------------------------------------------------- KERNELS -------------------------------------------------------------------------------------------//




#define compound_index blockIdx.x
__global__ void compoundKernel(SimulationDevice* sim) {
	__shared__ Compound compound;				// Mostly bond information
	__shared__ CompoundState compound_state;	// Relative position in [lm]
	__shared__ CompoundCoords compound_coords;	// Global positions in [lm]
	__shared__ NeighborList neighborlist;		
	__shared__ BondedParticlesLUT bonded_particles_lut;
	__shared__ Float3 utility_buffer[THREADS_PER_COMPOUNDBLOCK];
	utility_buffer[threadIdx.x] = Float3{};
	__shared__ float utility_buffer_f[THREADS_PER_COMPOUNDBLOCK];
	__shared__ Coord utility_coord;
	__shared__ Float3 utility_float3;
	__shared__ Coord rel_pos_shift;	// Maybe not needed, jsut use the utility one above?



	Box* box = sim->box;
	SimParams& simparams = *sim->params;

	if (threadIdx.x == 0) {
		compound.loadMeta(&box->compounds[blockIdx.x]);
		compound_state.setMeta(compound.n_particles);
		neighborlist.loadMeta(&box->compound_neighborlists[blockIdx.x]);
	}
	__syncthreads();


	compound.loadData(&box->compounds[blockIdx.x]);
	auto* coordarray_ptr = CoordArrayQueueHelpers::getCoordarrayRef(box->coordarray_circular_queue, simparams.step, blockIdx.x);

	// TODO: these compound_coords needs to be made const somehow. Changing these coords during the kernel is way too dangerous.
	compound_coords.loadData(*coordarray_ptr);
	neighborlist.loadData(&box->compound_neighborlists[blockIdx.x]);
	__syncthreads();

	compound_state.loadData(compound_coords);
	__syncthreads();

	float data_ptr[4]{};


	// TODO: Make this only if particle is part of bridge, otherwise skip the fetch and just use 0
	float potE_sum = compound.potE_interim[threadIdx.x];

	Float3 force = compound.forces[threadIdx.x];
	Float3 force_LJ_sol(0.f);

	// ------------------------------------------------------------ Intracompound Operations ------------------------------------------------------------ //
	{
		__syncthreads();
		force += computeSinglebondForces(compound.singlebonds, compound.n_singlebonds, compound_state.positions, utility_buffer, utility_buffer_f, &potE_sum);
		force += computeAnglebondForces(&compound, compound_state.positions, utility_buffer, utility_buffer_f, &potE_sum);
		force += computeDihedralForces(&compound, compound_state.positions, utility_buffer, utility_buffer_f, &potE_sum);
		force += computeImproperdihedralForces(compound.impropers, compound.n_improperdihedrals, compound_state.positions, utility_buffer, utility_buffer_f, &potE_sum);

		bonded_particles_lut.load(*box->bonded_particles_lut_manager->get(compound_index, compound_index));	// A lut always exists within a compound
		__syncthreads();
		force += computeIntracompoundLJForces(&compound, &compound_state, potE_sum, data_ptr, &bonded_particles_lut);
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------- //

	// --------------------------------------------------------------- Intercompound forces --------------------------------------------------------------- //	
	for (int i = 0; i < neighborlist.n_compound_neighbors; i++) {
		const uint16_t neighborcompound_id = neighborlist.neighborcompound_ids[i];
		const auto coords_ptr = CoordArrayQueueHelpers::getCoordarrayRef(box->coordarray_circular_queue, simparams.step, neighborcompound_id);
		getCompoundHyperpositionsAsFloat3(compound_coords.origo, coords_ptr, utility_buffer, &utility_float3);

		BondedParticlesLUT* compoundpair_lut = box->bonded_particles_lut_manager->get(compound_index, neighborcompound_id);
		const bool compounds_are_bonded = compoundpair_lut != nullptr;
		if (compounds_are_bonded) {
			bonded_particles_lut.load(*compoundpair_lut);
			__syncthreads();

			if (threadIdx.x < compound.n_particles) {
				force += computeIntercompoundLJForces(compound_state.positions[threadIdx.x], compound.atom_types[threadIdx.x], potE_sum, compound.particle_global_ids[threadIdx.x], data_ptr,
					&box->compounds[neighborcompound_id], utility_buffer, neighborcompound_id, bonded_particles_lut);
			}
		}
		else {
			if (threadIdx.x < compound.n_particles) {
				force += computeIntercompoundLJForces(compound_state.positions[threadIdx.x], compound.atom_types[threadIdx.x], potE_sum, compound.particle_global_ids[threadIdx.x], data_ptr,
					&box->compounds[neighborcompound_id], utility_buffer, neighborcompound_id);
			}
		}
		__syncthreads();
	}
	// ------------------------------------------------------------------------------------------------------------------------------------------------------ //

	// --------------------------------------------------------------- Solvation forces --------------------------------------------------------------- //
#ifdef ENABLE_SOLVENTS
	for (int i = 0; i < neighborlist.n_gridnodes; i++) {
		const int solventblock_id = neighborlist.gridnode_ids[i];
		const NodeIndex solventblock_hyperorigo = LIMAPOSITIONSYSTEM::getHyperNodeIndex(compound_coords.origo, SolventBlocksCircularQueue::get3dIndex(solventblock_id));

		const Float3 relpos_shift = LIMAPOSITIONSYSTEM_HACK::getRelShiftFromOrigoShift(solventblock_hyperorigo, compound_coords.origo).toFloat3();	// TODO: Only t0 needs to do this

		const SolventBlock* solventblock = box->solventblockgrid_circularqueue->getBlockPtr(solventblock_id, simparams.step);
		const int nsolvents_neighbor = solventblock->n_solvents;

		

		// There are many more solvents in a block, than threads in this kernel, so we need to loop in strides of blockdim.x
		__syncthreads();	// Dont load buffer before all are finished with the previous iteration
		for (uint32_t offset = 0; offset < nsolvents_neighbor; offset += blockDim.x) {
			const uint32_t solvent_index = offset + threadIdx.x;
			const int n_elements_this_stride = CPPD::min(nsolvents_neighbor - offset, blockDim.x);

			// Load the positions and add rel shift
			if (solvent_index < nsolvents_neighbor) {
				utility_buffer[threadIdx.x] = solventblock->rel_pos[solvent_index].toFloat3() + relpos_shift;
			}
			__syncthreads();

			if (threadIdx.x < compound.n_particles) {
				force += computeSolventToCompoundLJForces(compound_state.positions[threadIdx.x], n_elements_this_stride, utility_buffer, data_ptr, potE_sum, compound.atom_types[threadIdx.x]);
			}
			__syncthreads();
		}

	}
#endif
	// ------------------------------------------------------------------------------------------------------------------------------------------------ //

	// ------------------------------------------------------------ Integration ------------------------------------------------------------ //
	// From this point on, the origonal relpos is no longer acessible 
	{
		if (threadIdx.x < compound.n_particles) {	
			const float mass = forcefield_device.particle_parameters[compound.atom_types[threadIdx.x]].mass;

			if (blockIdx.x == 37 && threadIdx.x == 0 && simparams.step == 0 || true) {
				//printf("here %f %f %f %d\n", potE_sum, compound.vels_prev[threadIdx.x].len(), mass, compound.atom_types[threadIdx.x]);
			}

			const Float3 vel_now = EngineUtils::integrateVelocityVVS(compound.vels_prev[threadIdx.x], compound.forces_prev[threadIdx.x], force, simparams.constparams.dt, mass);
			const Coord pos_now = EngineUtils::integratePositionVVS(compound_coords.rel_positions[threadIdx.x], vel_now, force, mass, simparams.constparams.dt);

			compound.vels_prev[threadIdx.x] = vel_now * simparams.thermostat_scalar;
			compound.forces_prev[threadIdx.x] = force;

			// Save pos locally, but only push to box as this kernel ends
			compound_coords.rel_positions[threadIdx.x] = pos_now;
		}
	}
	// ------------------------------------------------------------------------------------------------------------------------------------- //


	// ------------------------------------------------------ PERIODIC BOUNDARY CONDITION ------------------------------------------------------------------------------- // 	
	{
		__shared__ Coord shift_lm;	// Use utility coord for this?
		if (threadIdx.x == 0) {
			shift_lm = LIMAPOSITIONSYSTEM::shiftOrigo(compound_coords);
		}
		__syncthreads();

		LIMAPOSITIONSYSTEM_HACK::shiftRelPos(compound_coords, shift_lm);
		__syncthreads();

		
	}
	LIMAPOSITIONSYSTEM_HACK::applyPBC(compound_coords);
	__syncthreads();
	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------ //

	EngineUtils::LogCompoundData(compound, box, compound_coords, &potE_sum, force, force_LJ_sol, simparams, sim->databuffers);

	// Push positions for next step
	auto* coordarray_next_ptr = CoordArrayQueueHelpers::getCoordarrayRef(box->coordarray_circular_queue, simparams.step + 1, blockIdx.x);
	coordarray_next_ptr->loadData(compound_coords);

	// Push vel and force for current step, for VelocityVS
	box->compounds[blockIdx.x].vels_prev[threadIdx.x] = compound.vels_prev[threadIdx.x];
	box->compounds[blockIdx.x].forces_prev[threadIdx.x] = compound.forces_prev[threadIdx.x];
}
#undef compound_index



#define solvent_active (threadIdx.x < solventblock.n_solvents)
#define solvent_mass (forcefield_device.particle_parameters[ATOMTYPE_SOLVENT].mass)
static_assert(MAX_SOLVENTS_IN_BLOCK > MAX_COMPOUND_PARTICLES, "solventForceKernel was about to reserve an insufficient amount of memory");
__global__ void solventForceKernel(SimulationDevice* sim) {
	__shared__ Float3 utility_buffer[MAX_SOLVENTS_IN_BLOCK];
	__shared__ uint8_t utility_buffer_small[MAX_SOLVENTS_IN_BLOCK];
	__shared__ SolventBlock solventblock;
	__shared__ SolventTransferqueue<SolventBlockTransfermodule::max_queue_size> transferqueues[6];		// TODO: Use template to make identical kernel, so the kernel with transfer is slower and larger, and the rest remain fast!!!!
	__shared__ int utility_int;
	__shared__ Coord utility_coord;
	__shared__ Float3 utility_float3;

	// Doubles as block_index_3d!
	const NodeIndex block_origo = SolventBlocksCircularQueue::get3dIndex(blockIdx.x);

	Box* box = sim->box;
	SimParams& simparams = *sim->params;
	SolventBlock* solventblock_ptr = box->solventblockgrid_circularqueue->getBlockPtr(blockIdx.x, simparams.step);

	// Init queue, otherwise it will contain wierd values // TODO: only do this on transferstep?
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
	}
	__syncthreads();
	solventblock.loadData(*solventblock_ptr);
	__syncthreads();	


	Float3 force(0.f);
	const Float3 relpos_self = solventblock.rel_pos[threadIdx.x].toFloat3();


	// --------------------------------------------------------------- Molecule Interactions --------------------------------------------------------------- //	
	{
		// Thread 0 finds n nearby compounds
		const CompoundGridNode* compoundgridnode = box->compound_grid->getBlockPtr(blockIdx.x);
		if (threadIdx.x == 0) { utility_int = compoundgridnode->n_nearby_compounds; }
		__syncthreads();



		for (int i = 0; i < utility_int; i++) {
			const int16_t neighborcompound_index = compoundgridnode->nearby_compound_ids[i];
			const Compound* neighborcompound = &box->compounds[neighborcompound_index];
			const int n_compound_particles = neighborcompound->n_particles;

			// All threads help loading the molecule
			// First load particles of neighboring compound
			const CompoundCoords* coordarray_ptr = CoordArrayQueueHelpers::getCoordarrayRef(box->coordarray_circular_queue, simparams.step, neighborcompound_index);
			getCompoundHyperpositionsAsFloat3(solventblock.origo, coordarray_ptr, utility_buffer, &utility_float3);	// This should take n_comp_part aswell!


			// Then load atomtypes of neighboring compound
			if (threadIdx.x < n_compound_particles) {
				utility_buffer_small[threadIdx.x] = neighborcompound->atom_types[threadIdx.x];
			}
			__syncthreads();

			//  We can optimize here by loading and calculate the paired sigma and eps, jsut remember to loop threads, if there are many aomttypes.


			if (solvent_active) {
				force += computeCompoundToSolventLJForces(relpos_self, n_compound_particles, utility_buffer, data_ptr, potE_sum, utility_buffer_small, solventblock.ids[threadIdx.x]);
			}
			__syncthreads();
		}
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //



	// --------------------------------------------------------------- Intrablock Solvent Interactions ----------------------------------------------------- //
	{
		__syncthreads(); // Sync since use of utility
		if (solvent_active) {
			utility_buffer[threadIdx.x] = relpos_self;
		}
		__syncthreads();
		if (solvent_active) {
			force += computeSolventToSolventLJForces(relpos_self, utility_buffer, solventblock.n_solvents, true, data_ptr, potE_sum);
		}
		__syncthreads(); // Sync since use of utility
	}	
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //

	// --------------------------------------------------------------- Interblock Solvent Interactions ----------------------------------------------------- //
	const int query_range = 1;
	for (int x = -query_range; x <= query_range; x++) {
		for (int y = -query_range; y <= query_range; y++) {
			for (int z = -query_range; z <= query_range; z++) {
				const NodeIndex dir{ x,y,z };
				if (dir.isZero()) { continue; }

				const int blockindex_neighbor = EngineUtils::getNewBlockId(dir, block_origo);
				KernelHelpersWarnings::assertValidBlockId(blockindex_neighbor);

				const SolventBlock* solventblock_neighbor = box->solventblockgrid_circularqueue->getBlockPtr(blockindex_neighbor, simparams.step);
				const int nsolvents_neighbor = solventblock_neighbor->n_solvents;
				const Float3 origoshift_offset = LIMAPOSITIONSYSTEM::nodeIndexToCoord(dir).toFloat3();

				// All threads help loading the solvent, and shifting it's relative position reletive to this solventblock
				__syncthreads();
				if (threadIdx.x < nsolvents_neighbor) {
					utility_buffer[threadIdx.x] = solventblock_neighbor->rel_pos[threadIdx.x].toFloat3() + origoshift_offset;
				}
				__syncthreads();

				if (solvent_active) {
					force += computeSolventToSolventLJForces(relpos_self, utility_buffer, nsolvents_neighbor, false, data_ptr, potE_sum);
				}
				__syncthreads();
			}
		}
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //

	

	Coord relpos_next{};

	if (solvent_active) {
		const auto mass = forcefield_device.particle_parameters[ATOMTYPE_SOLVENT].mass;
		Solvent& solventdata_ref = box->solvents[solventblock.ids[threadIdx.x]];	// Solvent private data, for VVS

		const Float3 vel_now = EngineUtils::integrateVelocityVVS(solventdata_ref.vel_prev, solventdata_ref.force_prev, force, simparams.constparams.dt, mass);
		const Coord pos_now = EngineUtils::integratePositionVVS(solventblock.rel_pos[threadIdx.x], vel_now, force, mass, simparams.constparams.dt);

		solventdata_ref.vel_prev = vel_now * simparams.thermostat_scalar;
		solventdata_ref.force_prev = force;

		// Save pos locally, but only push to box as this kernel ends
		relpos_next = pos_now;

		EngineUtils::LogSolventData(box, potE_sum, solventblock, solvent_active, force, vel_now, simparams.step, sim->databuffers);
	}



	// Push new SolventCoord to global mem
	SolventBlock* solventblock_next_ptr = box->solventblockgrid_circularqueue->getBlockPtr(blockIdx.x, simparams.step + 1);

	if (SolventBlocksCircularQueue::isTransferStep(simparams.step)) {
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
__global__ void solventTransferKernel(SimulationDevice* sim) {
	Box* box = sim->box;
	SimParams& simparams = *sim->params;

	SolventBlockTransfermodule* transfermodule = &box->transfermodule_array[blockIdx.x];
	
	SolventBlock* solventblock_current = box->solventblockgrid_circularqueue->getBlockPtr(blockIdx.x, simparams.step);
	SolventBlock* solventblock_next = box->solventblockgrid_circularqueue->getBlockPtr(blockIdx.x, simparams.step + 1);

	SolventTransferWarnings::assertSolventsEqualNRemain(*solventblock_next, *transfermodule);

	// Handling incoming transferring solvents
	int n_solvents_next = transfermodule->n_remain;
	for (int queue_index = 0; queue_index < SolventBlockTransfermodule::n_queues; queue_index++) {
		auto* queue = &transfermodule->transfer_queues[queue_index];
		if (threadIdx.x < queue->n_elements) {
			const int incoming_index = n_solvents_next + threadIdx.x;

			solventblock_next->rel_pos[incoming_index] = queue->rel_positions[threadIdx.x];
			solventblock_next->ids[incoming_index] = queue->ids[threadIdx.x];
		}
		n_solvents_next += queue->n_elements;

		// Signal that all elements of the queues have been moved
		__syncthreads();
		if (threadIdx.x == 0) {
			queue->n_elements = 0;
		}
	}

	SolventTransferWarnings::assertMaxPlacedSolventsIsWithinLimits(n_solvents_next, simparams.critical_error_encountered);

	// Finally update the solventblock_next with how many solvents it now contains
	if (threadIdx.x == 0) {
		solventblock_next->n_solvents = n_solvents_next;
	}
}





#define particle_id_bridge threadIdx.x
__global__ void compoundBridgeKernel(SimulationDevice* sim) {
	__shared__ CompoundBridge bridge;
	__shared__ Float3 positions[MAX_PARTICLES_IN_BRIDGE];
	__shared__ Float3 utility_buffer[MAX_PARTICLES_IN_BRIDGE];
	__shared__ float utility_buffer_f[MAX_PARTICLES_IN_BRIDGE];
	__shared__ Coord utility_coord[MAX_COMPOUNDS_IN_BRIDGE];

	SimParams& simparams = *sim->params;
	Box* box = sim->box;

	if (threadIdx.x == 0) {
		bridge.loadMeta(&box->bridge_bundle->compound_bridges[blockIdx.x]);

		
	}
	__syncthreads();

	// TODO: we dont need to do this for the first compound, as it will always be 0,0,0
	if (threadIdx.x < bridge.n_compounds) {
		// Calculate necessary shift in relative positions for right, so right share the origo with left.
		utility_coord[threadIdx.x] = LIMAPOSITIONSYSTEM_HACK::getRelativeShiftBetweenCoordarrays(box->coordarray_circular_queue, simparams.step, bridge.compound_ids[0], bridge.compound_ids[threadIdx.x]);
	}


	bridge.loadData(&box->bridge_bundle->compound_bridges[blockIdx.x]);
	__syncthreads();

	if (particle_id_bridge < bridge.n_particles) {
		ParticleReference& p_ref = bridge.particle_refs[particle_id_bridge];

		BridgeWarnings::verifyPRefValid(p_ref, bridge);

		const CompoundCoords* coordarray = CoordArrayQueueHelpers::getCoordarrayRef(box->coordarray_circular_queue, simparams.step, p_ref.compound_id);

		Coord relpos = coordarray->rel_positions[p_ref.local_id_compound];
		relpos += utility_coord[p_ref.compoundid_local_to_bridge];
		positions[threadIdx.x] = relpos.toFloat3();
	}
	__syncthreads();

	float potE_sum = 0;
	Float3 force{};

	// ------------------------------------------------------------ Intercompund Operations ------------------------------------------------------------ //
	{											// So for the very first step, these ´should all be 0, but they are not??										TODO: Look into this at some point!!!! 
		force += computeSinglebondForces(bridge.singlebonds, bridge.n_singlebonds, positions, utility_buffer, utility_buffer_f, &potE_sum);
		force += computeAnglebondForces(&bridge, positions, utility_buffer, utility_buffer_f, &potE_sum);
		force += computeDihedralForces(&bridge, positions, utility_buffer, utility_buffer_f, &potE_sum);
		force += computeImproperdihedralForces(bridge.impropers, bridge.n_improperdihedrals, positions, utility_buffer, utility_buffer_f, &potE_sum);
	}
	__syncthreads();
	// --------------------------------------------------------------------------------------------------------------------------------------------------- //

	if (particle_id_bridge < bridge.n_particles) {
		ParticleReference* p_ref = &bridge.particle_refs[particle_id_bridge];
		box->compounds[p_ref->compound_id].forces[p_ref->local_id_compound] = force;

		sim->box->compounds[p_ref->compound_id].potE_interim[p_ref->local_id_compound] = potE_sum;
	}
}

#pragma warning (pop)
#pragma warning (pop)
