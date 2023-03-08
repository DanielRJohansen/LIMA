#pragma once

#include<iostream>
#include "LimaTypes.cuh"
#include "Constants.cuh"

#include "Simulation.cuh"
#include "Forcefield.cuh"
#include <cmath>

#include <cuda.h>
#include <device_launch_parameters.h>
#include <cuda_runtime_api.h>

namespace ForceCalc {
//	-
};

#include "cuda/std/cmath"
//#include "cuda/std//utility"

namespace CPPD {
	__device__ __host__ constexpr int32_t ceil(float num) {
		return (static_cast<float>(static_cast<int32_t>(num)) == num)
			? static_cast<int32_t>(num)
			: static_cast<int32_t>(num) + ((num > 0) ? 1 : 0);
	}

	template <typename T>
	__device__ __host__ static T max(const T l, const T r) {
		return r > l ? r : l;
	}

	__device__ __host__ static int32_t abs(const int32_t val) {
		return val < 0 ? -val : val;
	}
}




namespace LIMAPOSITIONSYSTEM {
	//static Coord coordFromAbsPos(Float3 pos)
	static CompoundCoords positionCompound(CompoundState& state,  int key_particle_index=0) {
		CompoundCoords compoundcoords{};
		//Float3& key_pos = state.positions[key_particle_index];
		//compoundcoords.origo = Float3(key_pos);	
		compoundcoords.origo = Float3(0);	// Temp, use the one above in future

		double default_norm_dist = 1.0;	// By default, 2 nm has a relative distance of 1.0 (int float) and 2^29 (uint32_t

		for (int i = 0; i < state.n_particles; i++) {
			double x = (static_cast<double>(state.positions[i].x) - static_cast<double>(compoundcoords.origo.x)) / default_norm_dist;
			double y = (static_cast<double>(state.positions[i].y) - static_cast<double>(compoundcoords.origo.y)) / default_norm_dist;
			double z = (static_cast<double>(state.positions[i].z) - static_cast<double>(compoundcoords.origo.z)) / default_norm_dist;

			Float3 rel_pos_nm{ x, y, z };
			Coord rel_coord{ rel_pos_nm * NANO_TO_LIMA };
			compoundcoords.rel_positions[i] = rel_coord;


		}
		return compoundcoords;
	}

	__device__ static Float3 getGlobalPositionNM(CompoundCoords& coords) {
		return coords.origo.toFloat3() + coords.rel_positions[threadIdx.x].toFloat3() / NANO_TO_LIMA;
	}

	__device__ static Float3 getGlobalPositionFM(CompoundCoords& coords) {
		return coords.origo.toFloat3() * NANO_TO_FEMTO + coords.rel_positions[threadIdx.x].toFloat3() * LIMA_TO_FEMTO;
	}

	// Returns position in LimaMetres
	__device__ static Float3 getGlobalPosition(const CompoundCoords& coords) {
		return coords.origo.toFloat3() * NANO_TO_LIMA + coords.rel_positions[threadIdx.x].toFloat3();
	}

	// Returns positions in LimaMetres
	__device__ static void getGlobalPositions(CompoundCoords& coords, CompoundState& state) {
		state.positions[threadIdx.x] = getGlobalPosition(coords);
	}

	__device__ __host__ static bool canRepresentRelativeDist(const Coord& origo_a, const Coord& origo_b) {
		const auto diff = origo_a - origo_b;
		return std::abs(diff.x) < MAX_REPRESENTABLE_DIFF_NM && std::abs(diff.y) < MAX_REPRESENTABLE_DIFF_NM && std::abs(diff.z) < MAX_REPRESENTABLE_DIFF_NM;
	}

	__device__ static void getRelativePositions(Coord* coords, Float3* positions) {
		positions[threadIdx.x] = coords[threadIdx.x].toFloat3();
	}

	__device__ static void getRelativePositions(CompoundCoords& coords, CompoundState& state) {
		//state.positions[threadIdx.x] = coords.rel_positions[threadIdx.x].toFloat3();
		//or 
		getRelativePositions(coords.rel_positions, state.positions);
	}

	//static void applyHyperpos(CompoundCoords& lhs, CompoundCoords& rhs) {
	//	// TODO: IMplement
	//}
	
	// Calculate the necessary shift in LM of all elements of FROM, assuming the origo has been shifted to TO
	__device__ static Coord getRelShiftFromOrigoShift(const Coord& origo_from, const Coord& origo_to) {
		return (origo_from - origo_to) * static_cast<int32_t>(NANO_TO_LIMA);
	}

	__device__ static void applyHyperpos(const Coord& static_coord, Coord& movable_coord) {
		movable_coord.x += static_cast<int32_t>(BOX_LEN_NM) * ((static_coord.x - movable_coord.x) > static_cast<int32_t>( BOX_LEN_HALF_NM));
		movable_coord.x -= static_cast<int32_t>(BOX_LEN_NM) * ((static_coord.x - movable_coord.x) < static_cast<int32_t>(-BOX_LEN_HALF_NM));

		movable_coord.y += static_cast<int32_t>(BOX_LEN_NM) * ((static_coord.y - movable_coord.y) > static_cast<int32_t>(BOX_LEN_HALF_NM));
		movable_coord.y -= static_cast<int32_t>(BOX_LEN_NM) * ((static_coord.y - movable_coord.y) < static_cast<int32_t>(-BOX_LEN_HALF_NM));

		movable_coord.z += static_cast<int32_t>(BOX_LEN_NM) * ((static_coord.z - movable_coord.z) > static_cast<int32_t>(BOX_LEN_HALF_NM));
		movable_coord.z -= static_cast<int32_t>(BOX_LEN_NM) * ((static_coord.z - movable_coord.z) < static_cast<int32_t>(-BOX_LEN_HALF_NM));
	}


	//static void alignCoordinates(CompoundCoords& lhs, CompoundCoords& rhs) {
	//	applyHyperpos(lhs, rhs);

	//	Coord common_origo = {
	//		std::min(lhs.origo.x, rhs.origo.x) + ((std::max(lhs.origo.x, rhs.origo.x) - std::min(lhs.origo.x, rhs.origo.x)) / 2),
	//		std::min(lhs.origo.y, rhs.origo.y) + ((std::max(lhs.origo.y, rhs.origo.y) - std::min(lhs.origo.y, rhs.origo.y)) / 2),
	//		std::min(lhs.origo.z, rhs.origo.z) + ((std::max(lhs.origo.z, rhs.origo.z) - std::min(lhs.origo.z, rhs.origo.z)) / 2)
	//	};
	//}

	__device__ static Coord getHyperOrigo(const Coord& self, const Coord& other) {
		Coord temp = other;
		applyHyperpos(self, temp);
		return temp;
	}

	// The following two functions MUST ALWAYS be used together
	// Shift refers to the wanted difference in the relative positions, thus origo must move -shift.
	// Keep in mind that origo is in nm and rel_pos are in lm
	// ONLY CALL FROM THREAD 0
	__device__ static Coord shiftOrigo(CompoundCoords& coords, const int keyparticle_index=0) {
		//Coord shift_nm = -coords.rel_positions[keyparticle_index] / static_cast<uint32_t>(NANO_TO_LIMA);	// OPTIM. If LIMA wasn't 100 femto, but rather a power of 2, we could do this much better!

		Coord shift_nm = coords.rel_positions[keyparticle_index] / static_cast<int32_t>(NANO_TO_LIMA);	// OPTIM. If LIMA wasn't 100 femto, but rather a power of 2, we could do this much better!
		
		if (threadIdx.x == 0 && blockIdx.x == 1) {
			//printf("x %d shift.x %d  shift_out%d\n", coords.rel_positions[keyparticle_index].x, shift_nm.x, shift_nm.x * static_cast<int32_t>(NANO_TO_LIMA));
			//printf("%d %d %d\n", shift_nm.x, shift_nm.y, shift_nm.z);
		}
		
		coords.origo += shift_nm;
		return -shift_nm * static_cast<int32_t>(NANO_TO_LIMA);
	}


	__device__ static void shiftRelPos(CompoundCoords& coords, const Coord& shift_lm) {
		coords.rel_positions[threadIdx.x] += shift_lm;
	}


	static void moveCoordinate(Coord& coord, Float3 delta /*[nm]*/) {	// TODO: some checks to make this safe
		Coord delta_c{ delta / LIMA_TO_FEMTO };
		coord += delta_c;
	}

	__device__ static void applyPBC(CompoundCoords& coords) {
		if (threadIdx.x != 0) { return; }
		coords.origo.x += static_cast<int32_t>(BOX_LEN_NM) * (coords.origo.x < 0);
		coords.origo.x -= static_cast<int32_t>(BOX_LEN_NM) * (coords.origo.x > static_cast<int32_t>(BOX_LEN_NM));
		coords.origo.y += static_cast<int32_t>(BOX_LEN_NM) * (coords.origo.y < 0);
		coords.origo.y -= static_cast<int32_t>(BOX_LEN_NM) * (coords.origo.y > static_cast<int32_t>(BOX_LEN_NM));
		coords.origo.z += static_cast<int32_t>(BOX_LEN_NM) * (coords.origo.z < 0);
		coords.origo.z -= static_cast<int32_t>(BOX_LEN_NM) * (coords.origo.z > static_cast<int32_t>(BOX_LEN_NM));
	}

	// DANGEROUS: Only wraps around once! Min is assumed to be 0
	__device__ static int getIntegarWraparound(int num, const int max) {
		num += max * (num < 0);
		num -= max * (num > max);
		return num;
	}

	__device__ static Coord getRelativeShiftBetweenCoordarrays(CompoundCoords* coordarray_circular_queue, int step, int compound_index_left, int compound_index_right) {
		Coord& coord_origo_left = CoordArrayQueueHelpers::getCoordarrayPtr(coordarray_circular_queue, step, compound_index_left)->origo;
		Coord& coord_origo_right = CoordArrayQueueHelpers::getCoordarrayPtr(coordarray_circular_queue, step, compound_index_right)->origo;

		Coord hyperorigo_right = LIMAPOSITIONSYSTEM::getHyperOrigo(coord_origo_left, coord_origo_right);

		// Calculate necessary shift in relative position for all particles in other, so they share origo with left
		return (coord_origo_left - hyperorigo_right) * static_cast<uint32_t>(NANO_TO_LIMA);	// This fucks up when the diff is > ~20
	}

	// Calculates the relative position of movable_solvent, relative to another staticcoords's origo.
	// Returns false if it is not possible to represent the position as coord. In that case, we should avoid
	// following computations..
	__device__ static Coord getRelativeHyperposition(const SolventCoord& static_solvent, const SolventCoord& movable_solvent) {
		const Coord hyperorigo_other = LIMAPOSITIONSYSTEM::getHyperOrigo(static_solvent.origo, movable_solvent.origo);

		// calc Relative Position Shift from the origo-shift
		//const Coord relPosShiftOfMovable = LIMAPOSITIONSYSTEM::getRelShift(static_solvent.origo, hyperorigo_other);
		// 
		//const Coord relPosShiftOfMovable = LIMAPOSITIONSYSTEM::getRelShiftFromOrigoShift(static_solvent.origo, hyperorigo_movable);
		const Coord relPosShiftOfMovable = LIMAPOSITIONSYSTEM::getRelShiftFromOrigoShift(hyperorigo_other, static_solvent.origo);
		return movable_solvent.rel_position + relPosShiftOfMovable;
	}

	__device__ __host__ static void applyPBC(Coord& origo) {
		origo.x += BOX_LEN_NM_INT * (origo.x < 0);
		origo.x -= BOX_LEN_NM_INT * (origo.x >= BOX_LEN_NM_INT);
		origo.y += BOX_LEN_NM_INT * (origo.y < 0);
		origo.y -= BOX_LEN_NM_INT * (origo.y >= BOX_LEN_NM_INT);
		origo.z += BOX_LEN_NM_INT * (origo.z < 0);
		origo.z -= BOX_LEN_NM_INT * (origo.z >= BOX_LEN_NM_INT);
	}

	__device__ __host__ static void applyPBC(SolventCoord& coord) {	// Only changes position if position is outside of box;
		Coord& pos = coord.origo;
		applyPBC(pos);
	}

	// Moves origo when necessary, so all relpos dimensions are positive
	__host__ static void forceRelposPositive(SolventCoord& coord) {
		for (int i = 0; i < 3; i++) {
			if (*coord.rel_position.get(i) < 0) {
				(*coord.origo.get(i))--;
				(*coord.rel_position.get(i)) += static_cast<int32_t>(NANO_TO_LIMA);
			}
		}
	}



	// ReCenter origo of solvent, and the relative pos around said origo
	__device__ static void updateSolventcoord(SolventCoord& coord) {
		const int shift_at = static_cast<int32_t>(NANO_TO_LIMA) / 2;
		Coord shift_nm = coord.rel_position / static_cast<int32_t>(shift_at);	// OPTIM. If LIMA wasn't 100 femto, but rather a power of 2, we could do this much better! 

		SolventCoord tmp = coord;
		coord.origo += shift_nm;
		coord.rel_position -= shift_nm * static_cast<int32_t>(NANO_TO_LIMA);
		
		if (blockIdx.x + threadIdx.x == 0 && shift_nm.x != 0) {
			tmp.origo.print('o');
			tmp.rel_position.print('r');
			shift_nm.print('s');
			coord.origo.print('O');
			coord.rel_position.print('R');
		}
	}

	// Get the relpos_prev, if the solvent was in the same solventblock last step
	__device__ static Coord getRelposPrev(SolventBlockGrid* solventblockgrid_circularqueue, const int solventblock_id, const int current_step) {
		const int step_prev = current_step == 0 ? 0 : current_step - 1;	// Unnecessary if we use transfermodule's prevpos for first step!!!!!!!!! TODOTODO
		auto blockPtr = CoordArrayQueueHelpers::getSolventBlockPtr(solventblockgrid_circularqueue, step_prev, solventblock_id);
		return blockPtr->rel_pos[threadIdx.x];
	}

	//__device__ static bool willSolventRe

	// Get the relpos_prev, if the solvent was NOT in the same solventblock last step
	//__device__ static Coord getRelposPrevAfterTransfer(SolventBlockGrid* solventblockgrid_circularqueue, const int solventblock_id, const int step) {


	//__device__ static void applyPBC(Compound* compound);
};








namespace EngineUtils {
	__device__ __host__ static inline void applyHyperpos(const Float3* static_particle, Float3* movable_particle) {
		for (int i = 0; i < 3; i++) {
			//*movable_particle->placeAt(i) += BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) > BOX_LEN_HALF);
			//*movable_particle->placeAt(i) -= BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) < -BOX_LEN_HALF);	// use at not X!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			*movable_particle->placeAt(i) += BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) > BOX_LEN_HALF);
			*movable_particle->placeAt(i) -= BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) < -BOX_LEN_HALF);
		}
	}

	__device__ __host__ static float calcHyperDist(const Float3* p1, const Float3* p2) {
		Float3 temp = *p2;
		applyHyperpos(p1, &temp);
		return (*p1 - temp).len();
	}

	__device__ __host__ static float calcHyperDistDenormalized(const Float3* p1, const Float3* p2) {
		Float3 temp = *p2;
		applyHyperpos(p1, &temp);
		return (*p1 - temp).len() * NORMALIZER;
	}

	__device__ __host__ static float calcKineticEnergy(const Float3* pos1, const Float3* pos2, const float mass, const float elapsed_time) {
		const float vel = calcHyperDist(pos1, pos2) / elapsed_time * NORMALIZER;
		const float kinE = 0.5f * mass * vel * vel;
		return kinE;
	}

	__device__ __host__ static void applyPBC(Float3* current_position) {	// Only changes position if position is outside of box;
		for (int dim = 0; dim < 3; dim++) {
			/**current_position->placeAt(dim) += BOX_LEN * (current_position->at(dim) < 0.f);
			*current_position->placeAt(dim) -= BOX_LEN * (current_position->at(dim) > BOX_LEN);*/
			*current_position->placeAt(dim) += BOX_LEN * (current_position->at(dim) < 0.f);
			*current_position->placeAt(dim) -= BOX_LEN * (current_position->at(dim) > BOX_LEN);
		}
	}

	static void __host__ genericErrorCheck(const char* text) {
		cudaDeviceSynchronize();
		cudaError_t cuda_status = cudaGetLastError();
		if (cuda_status != cudaSuccess) {
			std::cout << "\nCuda error code: " << cuda_status << " - " << cudaGetErrorString(cuda_status) << std::endl;
			fprintf(stderr, text);
			exit(1);
		}
	}

	static float calcSpeedOfParticle(const float mass /*[kg]*/, const float temperature /*[K]*/) { // 
		const float R = 8.3144f;								// Ideal gas constants - J/(Kelvin*mol)
		const float v_rms = static_cast<float>(sqrt(3.f * R * temperature / mass));
		return v_rms;	// [m/s]
	}

	// For solvents, compound_id = n_compounds and particle_id = solvent_index
	__device__ static uint32_t getLoggingIndexOfParticle(uint32_t step, uint32_t total_particles_upperbound, uint32_t compound_id, uint32_t particle_id_local) {
		const uint32_t steps_since_transfer = (step % STEPS_PER_LOGTRANSFER);
		const uint32_t step_offset = steps_since_transfer * total_particles_upperbound;
		const uint32_t compound_offset = compound_id * MAX_COMPOUND_PARTICLES;
		return step_offset + compound_offset + particle_id_local;
	}

	__host__ static size_t getAlltimeIndexOfParticle(uint64_t step, uint32_t total_particles_upperbound, uint32_t compound_id, uint32_t particle_id_local) {
		const uint32_t step_offset = static_cast<uint32_t>(step) * total_particles_upperbound;
		const uint32_t compound_offset = compound_id * MAX_COMPOUND_PARTICLES;
		return step_offset + compound_offset + particle_id_local;
	}

	__device__ int static getNewBlockId(const Coord& transfer_direction, const Coord& origo) {
		auto newCoord3d = transfer_direction + origo;
		LIMAPOSITIONSYSTEM::applyPBC(newCoord3d);
		return SolventBlockGrid::get1dIndex(newCoord3d);
	}

	__device__ static SolventBlockTransfermodule* getTransfermoduleTargetPtr(SolventBlockTransfermodule* transfermodule_array, int blockId, const Coord& transfer_direction) {
		Coord new3dIndex = SolventBlockGrid::get3dIndex(blockId) + transfer_direction;
		LIMAPOSITIONSYSTEM::applyPBC(new3dIndex);
		auto index = SolventBlockGrid::get1dIndex(new3dIndex);
		return &transfermodule_array[index];
	}

	__device__ static Coord getOnehotDirection(const Coord coord, int32_t threshold) {
		const int32_t magnitude_x = std::abs(coord.x);
		const int32_t magnitude_y = std::abs(coord.y);
		const int32_t magnitude_z = std::abs(coord.z);

		if (magnitude_x < threshold && magnitude_y < threshold && magnitude_z < threshold) { return Coord{ 0,0,0 }; }

		// Determine which magnitude is the largest
		if (magnitude_x >= magnitude_y && magnitude_x >= magnitude_z) {
			// The x component has the largest magnitude
			return Coord{ coord.x < 0 ? -1 : 1, 0, 0 };
		}
		else if (magnitude_y >= magnitude_x && magnitude_y >= magnitude_z) {
			// The y component has the largest magnitude
			return Coord{ 0, coord.y < 0 ? -1 : 1, 0 };
		}
		else {
			// The z component has the largest magnitude
			return Coord{ 0, 0, coord.z < 0 ? -1 : 1 };
		}
	}
	__device__ static Coord getOnehotDirection1(const Coord coord, int32_t threshold) {
		
		//if (coord.x > 0 && coord.x != std::abs(coord.x)) {
		//	printf("\n %d %d\n", coord.x, std::abs(coord.x));
		//}
		//if (coord.y > 0 && coord.y != std::abs(coord.y)) {
		//	printf("\n %d %d\n", coord.y, std::abs(coord.y));
		//}
		//if (coord.z > 0 && coord.z != std::abs(coord.z)) {
		//	printf("\n %d %d\n", coord.z, std::abs(coord.z));
		//}


		const int32_t max_val{ 
			CPPD::max(
			int32_t{CPPD::max(std::abs(coord.x), std::abs(coord.y))},
			std::abs(coord.z))
		};

		int max2 = std::abs(coord.x);
		if (std::abs(coord.y) > max2) { max2 = std::abs(coord.y); }
		if (std::abs(coord.z) > max2) { max2 = std::abs(coord.z); }

		int max3 = CPPD::max(coord.x, coord.y);
		max3 = CPPD::max(coord.z, max3);

		if (max_val > threshold) {
			// Could be optimized by using threshold as a template argument
			Coord onehot =  coord / max_val;
			if (onehot.x + onehot.y + onehot.z > 1) {
				//printf("\nMax val %d max2 %d max3 %d\n", max_val, max2, max3);
				onehot.print('O');
				coord.print('C');
				//printf("\nABS: Max  %d  %d  %d\n", std::abs(coord.x), std::abs(coord.y), std::abs(coord.z));
				//printf("issame %d %d %d\n", coord.x == std::abs(coord.x), coord.y == std::abs(coord.y), coord.z == std::abs(coord.z));
				printf("x %d %d %d\n", coord.x, CPPD::abs(coord.x), coord.x < 0);
				printf("y %d %d %d\n", coord.y, CPPD::abs(coord.y), coord.x < 0);
				printf("z %d %d %d\n", coord.z, CPPD::abs(coord.z), coord.x < 0);
				printf("\n");
				//__nvvm_atom_max()
			}
			// Handle the case where some coordinates are identical
			onehot.y = onehot.x != 0 ? 0 : onehot.y;
			onehot.z = onehot.x != onehot.y ? 0 : onehot.z;	// x and y can't both be 1, since x already has priority over x
			return onehot;
		}

		return Coord{};
	}

	// returns an int between -21 and 21
	__device__ static int genPseudoRandomNum(int& seed) {
		const unsigned int a = 1664525;
		const unsigned int c = 1013904223;
		seed = a * seed + c;
		return seed / 100000000;
	}
	//__device__ static Coord getOnehotDirection(const Coord& coord, int32_t threshold) {
	//	int32_t max_val{ CPPD::max(
	//		CPPD::max(std::abs(coord.x), std::abs(coord.y)),
	//		std::abs(coord.z))
	//	};
	//	if (max_val > threshold) {
	//		// Could be optimized by using threshold as a template argument
	//		return coord / max_val;
	//	}
	//		
	//	return Coord{};
	//}

	// Since coord is rel to 0,0,0 of a block, we need to offset the positions so they are scattered around the origo instead of above it
	// We also need a threshold of half a blocklen, otherwise we should not transfer, and return{0,0,0}
	__device__ static Coord getTransferDirection(const Coord relpos) {
		const int32_t blocklen_half = static_cast<int32_t>(NANO_TO_LIMA) / 2;
		const Coord rel_blockcenter{ blocklen_half };
		if (relpos.x < INT32_MIN + blocklen_half || relpos.y < INT32_MIN + blocklen_half || relpos.z < INT32_MIN + blocklen_half ) {
			printf("\nWe have underflow!\n");
			relpos.print('R');
		}
		if (relpos.x > INT32_MAX - blocklen_half || relpos.y > INT32_MAX - blocklen_half || relpos.z > INT32_MAX - blocklen_half) {
			printf("\nWe have overflow!\n");
			relpos.print('R');
		}
		return EngineUtils::getOnehotDirection(relpos - rel_blockcenter, blocklen_half);
	}

	//__device__ static void doSolventTransfer(const Coord& relpos, const Coord& relpos_prev, SolventBlockTransfermodule* transfermodule_array) {
	//	//const Coord transfer_direction = getTransferDirection(relpos);
	//	//const Coord blockId3d = SolventBlockHelpers::get3dIndex(blockIdx.x);

	//	//int new_blockid = getNewBlockId(transfer_direction, blockId3d);
	//	//if (new_blockid == blockIdx.x) {
	//	//	transfermodule_array[blockIdx.x].remain_queue.addElement(threadIdx.x, relpos, relpos_prev);
	//	//}
	//	//else {
	//	//	// Coord 
	//	//	// TEMP:
	//	//	transfermodule_array[blockIdx.x].remain_queue.addElement(threadIdx.x, relpos, relpos_prev);
	//	//	transfermodule_array[new_blockid].getQueuePtr()/
	//	//}
	//}

	/*__device__ static void compressTransferqueue(SolventTransferqueue<SolventBlockTransfermodule::max_queue_size>* transferqueues) {

		for (int i = 0; i < transferqueue.used_size; i += blockDim.x) {

		}
	}*/
};


struct CompoundGridNode {
public:
	__host__ void addCompound(int16_t compound_id);
	__device__ int16_t getNElements() { return n_nearby_compounds; }
	__device__ int16_t getElement(int index) { return nearby_compound_ids[index]; }

private:
	static const int max_elements = 64;
	// A particle belonging to this node coord, can iterate through this list
	// to find all appropriate nearby compounds;
	int16_t nearby_compound_ids[64];	// MAX_COMPOUNDS HARD LIMIT
	int16_t n_nearby_compounds = 0;
};

// Class for signaling compound origo's and quickly searching nearby compounds using a coordinate on the grid
class CompoundGrid : private BoxGrid<CompoundGridNode> {

	// Each compound in kernel will transmit their origo. Will be transferring from device to host by nlist
	Coord compound_origos[MAX_COMPOUNDS];




};

using CompoundGrid = BoxGrid<CompoundGridNode>;





// LIMA algorithm Library
namespace LAL {
	__device__ constexpr int getBlellochTablesize(int n) {
		float nf = static_cast<float>(n);
		return CPPD::ceil(nf * log2f(nf) * 2.f);
	}
}