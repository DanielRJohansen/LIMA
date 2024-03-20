#include "ChargeOcttree.cuh"
#include "EngineUtils.cuh"
#include "BoundaryCondition.cuh"

#include "Engine.cuh"
#include "SimulationDevice.cuh"
#include "ForceComputations.cuh"



using namespace SCA;

template<typename T>
__device__ void ParallelSum(T* arrayptr, int array_len) {				// Places the result at pos 0 of input_array
	T temp;			// This is a lazy solution, but maybe it is also fast? Definitely simple..
	for (int i = 1; i < array_len; i *= 2) {	// Distributed averaging							// Make a generic and SAFER function for this, PLEASE OK??
		if ((threadIdx.x + i) < array_len) {
			temp = arrayptr[threadIdx.x] + arrayptr[threadIdx.x + i];
		}
		__syncthreads();
		arrayptr[threadIdx.x] = temp;
		__syncthreads();
	}
}

__device__ void ParallelSum(float* arrayPtr, int arrayLen, int stride) {
	const int index = threadIdx.x * stride;
	float temp{};
	for (int i = 1; i < arrayLen; i *= 2) {			
		if ((index + i) < arrayLen) {
			temp = arrayPtr[index] + arrayPtr[index + i];
		}
		__syncthreads();
		if ((index + i) < arrayLen) {
			arrayPtr[index] = temp;
		}
		__syncthreads();
	}
}

// diff = self-other
__device__ Float3 CalcCoulumbHalfforce(const Float3& diff, float charge) {
	const float dist = diff.len();	// [nm]
	const float magnitude = coulumb_constant * charge / (dist * dist);
	return diff.norm() * magnitude * LIMA_TO_NANO;
}

namespace OcttreeHelpers {
	// Where n is the difference in depth for the parent.
	// Please dont call with negative n..
	__device__ static Int3 GetNParentNodeId3D(Int3 id3D, int n) {
		return Int3{ id3D.x >> n, id3D.y >> n, id3D.z >> n };
	}

	/*__device__ static size_t getAbsParentIndex(const Int3& index_3d, const int n_nodes_per_dim, size_t depth_offset) {
		return GetAbsIndex(index_3d / 2, n_nodes_per_dim / 2, depth_offset / 8);
	}*/

	__device__ static void applyPBC(Int3& index3d, int dim) {
		index3d.x += dim * (index3d.x < 0);
		index3d.x -= dim * (index3d.x >= dim);
		index3d.y += dim * (index3d.y < 0);
		index3d.y -= dim * (index3d.y >= dim);
		index3d.z += dim * (index3d.z < 0);
		index3d.z -= dim * (index3d.z >= dim);
	}

	// Get abs pos in nm
	__device__ static Float3 getAbsPos(const Int3& index3d, float nodelen_nm) {
		return Float3{
			static_cast<float>(index3d.x) * nodelen_nm + nodelen_nm * 0.5f,
			static_cast<float>(index3d.y) * nodelen_nm + nodelen_nm * 0.5f,
			static_cast<float>(index3d.z) * nodelen_nm + nodelen_nm * 0.5f
		};
	}
}


// Specialized to handle ALOT of data
__global__ void propagateChargesUp(SimulationDevice* simdev, int depth, int nNodesPerDim, size_t depthOffset) {
	float2 buffer;

	float chargeSum = 0.f;

	for (int zOff = 0; zOff < 2; zOff++) {
		for (int yOff = 0; yOff < 2; yOff++) {
			const Int3 queryIndex3D{
				threadIdx.x * 2,
				blockIdx.y * 2 + yOff,
				blockIdx.z * 2 + zOff};

			const size_t query_index = OcttreeHelpers::GetAbsIndex(queryIndex3D, nNodesPerDim * 2, depthOffset); // Adjust depth_offset as needed
			buffer = *reinterpret_cast<float2*>(&simdev->charge_octtree->chargenodes[query_index]);
			chargeSum += buffer.x;
			chargeSum += buffer.y;
		}
	}


	const Int3 index3D{ threadIdx.x, blockIdx.y, blockIdx.z };
	const size_t node_index = OcttreeHelpers::GetAbsIndex(index3D, nNodesPerDim, simdev->charge_octtree->depthOffsets[depth]);
	simdev->charge_octtree->chargenodes[node_index] = chargeSum;
}


// 32 threads per block. 1 block per node in the force-tree
__global__ void DownwardSweep(SimulationDevice* simdev) {	
	__shared__ Float3 forces[32];

	forces[threadIdx.x] = Float3{};

	// At leaf level in the forcetree
	const Int3 forcenodeId3D{ blockIdx.x, blockIdx.y, blockIdx.z};
	const Float3 forcenodeAbsPos = OcttreeHelpers::getAbsPos(forcenodeId3D, ChargeOctTree::forcetree_leafnodelen);


	for (int depth = 3; depth < ChargeOctTree::depthForcetree; depth++)
	{
		const int range = depth == 3 ? 3 : 2 + 5;

		Int3 currentParentId3D = OcttreeHelpers::GetNParentNodeId3D(forcenodeId3D, ChargeOctTree::depthForcetree - depth);
		const int nNodesPerDim = 1 << depth;
		const float nodelenNM = BOX_LEN_NM / static_cast<float>(nNodesPerDim);

		// An odd side will always have to break into an extra cube "left", and an even side must break down an extra cube "right"
		const int x_odd = (currentParentId3D.x % 2) == 1;
		const int y_odd = (currentParentId3D.y % 2) == 1;
		const int z_odd = (currentParentId3D.z % 2) == 1;


		const size_t depthOffset = simdev->charge_octtree->depthOffsets[depth];


		const int surroundingNodesTotal = range * 2 * range * 2 * range * 2;
		for (int index = threadIdx.x; index < surroundingNodesTotal; index += blockDim.x) {

			const int x_off = (index / (range * 2 * range * 2)) - range - x_odd;
			const int y_off = ((index / (range * 2)) % (range * 2)) - range - y_odd;
			const int z_off = (index % (range * 2)) - range - z_odd;


			// If we are inside this border, we need to go to next level of precision
			if (std::abs(x_off) < 2 && std::abs(y_off) < 2 && std::abs(z_off) < 2) {	// These 3^3 out of 6^3 threads never have to do any work
				continue;
			}

			Int3 query_index3d = currentParentId3D + Int3{ x_off, y_off, z_off };
			OcttreeHelpers::applyPBC(query_index3d, nNodesPerDim);

			const size_t query_index = OcttreeHelpers::GetAbsIndex(query_index3d, nNodesPerDim, depthOffset);
			const float query_chargenode = simdev->charge_octtree->chargenodes[query_index];

			Float3 pos_abs_query = OcttreeHelpers::getAbsPos(query_index3d, nodelenNM);
			LIMAPOSITIONSYSTEM::applyHyperposNM<PeriodicBoundaryCondition>(&forcenodeAbsPos, &pos_abs_query);

			forces[threadIdx.x] += CalcCoulumbHalfforce(forcenodeAbsPos - pos_abs_query, query_chargenode);
		}
	}
	__syncthreads();

	ParallelSum(forces, 32);

	// Push the force to the tree
	if (threadIdx.x == 0) {
		//const size_t forcenodeIndex = OcttreeHelpers::GetAbsIndex(forcenodeId3D, ChargeOctTree::nLeafnodesPerDim_Forcetree, simdev->charge_octtree->depthOffsets[ChargeOctTree::depthForcetree]);

		const size_t forcenodeIndex = OcttreeHelpers::GetIndexAtDepth(forcenodeId3D, ChargeOctTree::nLeafnodesPerDim_Forcetree);
		simdev->charge_octtree->halfforceNodes[forcenodeIndex] = forces[0];
	}
}


__global__ void CompoundsFetchChargeforce(SimulationDevice* simdev) {
	__shared__ int n_particles;
	__shared__ CompoundCoords* compoundcoords;
	__shared__ Compound* compound;

	if (threadIdx.x == 0) {
		compound = &simdev->box->compounds[blockIdx.x];
		n_particles = compound->n_particles;
		compoundcoords = simdev->box->compoundcoordsCircularQueue->getCoordarrayRef(simdev->signals->step, blockIdx.x);
	}
	__syncthreads();

	if (threadIdx.x >= n_particles) {
		return;
	}
	Float3 abspos = LIMAPOSITIONSYSTEM::getAbsolutePositionNM(compoundcoords->origo, compoundcoords->rel_positions[threadIdx.x]);
	LIMAPOSITIONSYSTEM::applyBCNM<PeriodicBoundaryCondition>(&abspos);
	
	const Float3 charge_force = simdev->charge_octtree->getForceAtLeaf(abspos, compound->atom_charges[threadIdx.x]);
	simdev->box->compounds[blockIdx.x].forces_interim[threadIdx.x] += charge_force;
}





namespace SCA {
	float* chargenodesPtr;	// For quickly clearing all charges

	void wipeCharges() {
		cudaMemset(chargenodesPtr, 0, ChargeOctTree::nNodesChargetree * sizeof(float));
	}

	ChargeOctTree::ChargeOctTree() {
		cudaMalloc(&chargenodes, nNodesChargetree * sizeof(float));
		//cudaMalloc(&halfcharge_forces, n_nonleaf_nodes * sizeof(Float3));
		cudaMalloc(&halfforceNodes, force_n_leafnodes * sizeof(Float3));
		//cudaMalloc(&nodelists, n_nonleaf_nodes * sizeof(NodeList));
		chargenodesPtr = chargenodes;

		static const std::array<size_t, ChargeOctTree::depthChargetree+1> depthOffsetsTemp = ChargeOctTree::getDepthOffsets();
		cudaMalloc(&depthOffsets, sizeof(depthOffsetsTemp));
		cudaMemcpy(depthOffsets, depthOffsetsTemp.data(), sizeof(depthOffsetsTemp), cudaMemcpyHostToDevice);
	}
	// TODO: Make destructor


	int handleElectrostatics(SimulationDevice* sim_dev, const BoxParams& boxparams) {
		const auto t0 = std::chrono::high_resolution_clock::now();

		for (int depth = ChargeOctTree::depthChargetree - 1; depth >= 3; depth--) {
			const uint32_t parentnodesPerDim = OcttreeHelpers::getNNodesPerDim(depth);
			const uint3 blockDims{ 1u, parentnodesPerDim, parentnodesPerDim };
			propagateChargesUp << < blockDims, parentnodesPerDim >> > (sim_dev, depth, parentnodesPerDim, OcttreeHelpers::getDepthOffset(depth + 1));
			//LIMA_UTILS::genericErrorCheck("SCA error during propagate up");
		}

		{
			const uint32_t dim = static_cast<uint32_t>(ChargeOctTree::nLeafnodesPerDim_Forcetree);
			const uint3 gridDims{ dim, dim, dim };
			DownwardSweep << <gridDims, 32 >> > (sim_dev);
			//LIMA_UTILS::genericErrorCheck("SCA error during propagate down");
		}	
		
		const auto t1 = std::chrono::high_resolution_clock::now();


		CompoundsFetchChargeforce<<<boxparams.n_compounds, MAX_COMPOUND_PARTICLES>>>(sim_dev);

		LIMA_UTILS::genericErrorCheck("SCA error during distributeForces");

		wipeCharges();


		const auto t2 = std::chrono::high_resolution_clock::now();
		return static_cast<int>(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t0).count());
	}
}