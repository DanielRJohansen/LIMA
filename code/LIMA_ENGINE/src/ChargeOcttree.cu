#include "ChargeOcttree.cuh"
//#include "Simulation.cuh"
#include "EngineUtils.cuh"
//#include "Engine.cuh"
#include "BoundaryCondition.cuh"


#include "Engine.cuh"
#include "SimulationDevice.cuh"

using namespace SCA;

template<typename T>
__device__ void ParallelSum(T* arrayptr, int array_len) {				// Places the result at pos 0 of input_array
	T temp;			// This is a lazy soluation, but maybe it is also fast? Definitely simple..
	for (int i = 1; i < array_len; i *= 2) {	// Distributed averaging							// Make a generic and SAFER function for this, PLEASE OK??
		if ((threadIdx.x + i) < array_len) {
			temp = arrayptr[threadIdx.x] + arrayptr[threadIdx.x + i];
		}
		__syncthreads();
		arrayptr[threadIdx.x] = temp;
		__syncthreads();
	}
}

// diff = self-other
__device__ Float3 CalcCoulumbHalfforce(const Float3& diff, float charge) {
	const float dist = diff.len();	// [nm]
	const float magnitude = coulumb_constant * charge / (dist * dist);
	if (dist == 0.f) { printf("Yo why is dist 0"); }
	return diff.norm() * magnitude * LIMA_TO_NANO;
}

namespace OcttreeHelpers {
	// Where n is the difference in depth for the parent.
	// Please dont call with negative n..
	__device__ static Int3 GetNParentNodeId3D(Int3 id3D, int n) {
		return Int3{ id3D.x >> n, id3D.y >> n, id3D.z >> n };
	}

	__device__ static size_t getAbsParentIndex(const Int3& index_3d, const int n_nodes_per_dim, size_t depth_offset) {
		return GetAbsIndex(index_3d / 2, n_nodes_per_dim / 2, depth_offset / 8);
	}

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






// Loop over all 8 children
//for (int x_off = 0; x_off < 2; x_off++) {
//	for (int y_off = 0; y_off < 2; y_off++) {
//		for (int z_off = 0; z_off < 2; z_off++) {



//__global__ void propagateChargesUp(ChargeOctTree* tree, int depth, size_t depth_offset, size_t depth_offset_child, int n_nodes_per_dim, int n_nodes_per_dim_child, float nodelen_nm) {
//
//	__shared__ float charges[64][8];
//
//	charges[threadIdx.y][threadIdx.x] = 0;
//
//	// Threads distributed along x axis to maximize coalescence
//	const Int3 index_3d{ blockIdx.x, blockIdx.y, threadIdx.y };
//	if (index_3d.x >= n_nodes_per_dim) { return; }
//
//
//	const size_t node_index = OcttreeHelpers::GetAbsIndex(index_3d, n_nodes_per_dim, depth_offset);
//	const Float3 pos_abs_self = OcttreeHelpers::getAbsPos(index_3d, nodelen_nm);
//
//	float charge = 0.f;
//
//				
//	const int x_off = threadIdx.x % 2;
//	const int y_off = (threadIdx.x / 2) % 2;
//	const int z_off = (threadIdx.x / 4) % 2;
//
//	// Dont need PBC, since we are only looking at direct children
//	const Int3 query_index3d = index_3d * 2 + Int3{ x_off, y_off, z_off };
//
//	const size_t query_index = OcttreeHelpers::GetAbsIndex(query_index3d, n_nodes_per_dim_child, depth_offset_child);
//	const float query_chargenode = tree->chargenodes[query_index];
//	charges[threadIdx.y][threadIdx.x] += query_chargenode;
//
//	ParallelSum<float>(charges[threadIdx.y], 8);
//
//	if (threadIdx.x == 0) {
//		tree->chargenodes[node_index] = charges[threadIdx.y][0];
//	}
//}

__global__ void propagateChargesUp(ChargeOctTree* tree) {
	__shared__ float charges[64][8];


	charges[threadIdx.y][threadIdx.x] = 0;

	//int depth = ChargeOctTree::depthForcetree-1;
	//int n_nodes_per_dim = 1 << depth; // Assuming the tree is complete and balanced

	//// Calculate initial 3D index and node index
	const Int3 index_3d{ blockIdx.x, blockIdx.y, threadIdx.y };


	for (int depth = ChargeOctTree::depthForcetree - 1; depth >= 3; depth--) {
		const int n_nodes_per_dim = 1 << depth;

		if (index_3d.x >= n_nodes_per_dim || index_3d.y >= n_nodes_per_dim) return;	// But we still need all threads



		const size_t depth_offset = tree->depthOffsets[depth];

		float charge = 0.f;
		const int x_off = threadIdx.x % 2;
		const int y_off = (threadIdx.x / 2) % 2;
		const int z_off = (threadIdx.x / 4) % 2;

		Int3 query_index3d = index_3d * 2 + Int3{ x_off, y_off, z_off };
		size_t query_index = OcttreeHelpers::GetAbsIndex(query_index3d, n_nodes_per_dim * 2, depth_offset); // Adjust depth_offset as needed
		float query_charge = tree->chargenodes[query_index];
		charges[threadIdx.y][threadIdx.x] += query_charge;

		ParallelSum<float>(charges[threadIdx.y], 8);

		const size_t node_index = OcttreeHelpers::GetAbsIndex(index_3d, n_nodes_per_dim, depth_offset);

		if (threadIdx.x == 0) {
			tree->chargenodes[node_index] = charges[threadIdx.y][0];
		}

		__syncthreads(); // Ensure all threads have written their results

		// Prepare for the next iteration (move up the tree)
		if (index_3d.x % 2 == 0 && index_3d.y % 2 == 0 && index_3d.z % 2 == 0) {
			// This block will participate in the next level
			index_3d = index_3d / 2;
			//node_index = OcttreeHelpers::GetAbsIndex(index_3d, n_nodes_per_dim / 2, depth_offset); // Adjust depth_offset for the parent level
		}
		else {
			break; // This block does not participate in computing charges at the next level
		}

		depth--;
		n_nodes_per_dim /= 2;
		//depth_offset = ...; // Update depth_offset to point to the start of the block at the new depth
		__syncthreads(); // Ensure all threads are synchronized before potentially looping again
	}
}


// 32 threads per block. 1 block per node in the force-tree
__global__ void DownwardSweep(ChargeOctTree* tree) {	
	__shared__ Float3 forces[32];

	forces[threadIdx.x] = Float3{};

	// At leaf level in the forcetree
	const Int3 forcenodeId3D{ blockIdx.x, blockIdx.y, blockIdx.z};
	const Float3 forcenodeAbsPos = OcttreeHelpers::getAbsPos(forcenodeId3D, ChargeOctTree::forcetree_leafnodelen);
	

	//const float forcenodeCharge = tree->chargenodes[OcttreeHelpers::GetAbsIndex(forcenodeId3D, ChargeOctTree::force_n_leafnodes_per_dim, forcenodes_depthoffset)];
	const float forcenodeCharge = tree->chargenodes[OcttreeHelpers::GetAbsIndex(forcenodeId3D, ChargeOctTree::nLeafnodesPerDim_Forcetree, tree->depthOffsets[ChargeOctTree::depthForcetree])];

	if (std::abs(forcenodeCharge) < 1e-6) {
		return;
	}

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

		//const size_t depthOffset = ChargeOctTree::getDepthOffset(depth);
		//const size_t depthOffset = depthOffsets[depth];
		const size_t depthOffset = tree->depthOffsets[depth];


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
			const float query_chargenode = tree->chargenodes[query_index];



			Float3 pos_abs_query = OcttreeHelpers::getAbsPos(query_index3d, nodelenNM);
			LIMAPOSITIONSYSTEM::applyHyperposNM<PeriodicBoundaryCondition>(&forcenodeAbsPos, &pos_abs_query);

			forces[threadIdx.x] += CalcCoulumbHalfforce(forcenodeAbsPos - pos_abs_query, query_chargenode);
		}
	}
	__syncthreads();

	//EngineUtils::ParallelSum(forces, 6 * 6 * 6);
	ParallelSum(forces, 32);


	// Push the force to the tree
	if (threadIdx.x == 0) {
		//if (forces[0].len() != 0.f)
		//	forces[0].print('F', true);

		const size_t forcenodeIndex = OcttreeHelpers::GetIndexAtDepth(forcenodeId3D, ChargeOctTree::nLeafnodesPerDim_Forcetree);
		/*forcenodeId3D.print('W');
		printf("Sending to index %d\n", (int)forcenodeIndex);*/
		tree->halfforceNodes[forcenodeIndex] = forces[0];
	}
}


__global__ void CompoundsFetchChargeforce(SimulationDevice* simdev) {
	__shared__ int n_particles;
	__shared__ CompoundCoords* compoundcoords;
	__shared__ Compound* compound;

	if (threadIdx.x == 0) {
		compound = &simdev->box->compounds[blockIdx.x];
		n_particles = compound->n_particles;
		compoundcoords = CoordArrayQueueHelpers::getCoordarrayRef(simdev->box->coordarray_circular_queue, simdev->signals->step, blockIdx.x);
	}
	__syncthreads();

	if (threadIdx.x >= n_particles) {
		return;
	}
	Float3 abspos = LIMAPOSITIONSYSTEM::getAbsolutePositionNM(compoundcoords->origo, compoundcoords->rel_positions[threadIdx.x]);
	LIMAPOSITIONSYSTEM::applyBCNM<PeriodicBoundaryCondition>(&abspos);
	
	const Float3 charge_force = simdev->charge_octtree->getForceAtLeaf(abspos, compound->atom_charges[threadIdx.x]);
	//charge_force.print('f', true);

	//printf("my charge %f\n", compound->atom_charges[threadIdx.x]);

	simdev->box->compounds[blockIdx.x].forces_interim[threadIdx.x] += charge_force;

	//printf("Block %d force.x %f\n", blockIdx.x, charge_force.x);
}



void wipeCharges(ChargeOctTree* tree_dev) {
	cudaMemset(tree_dev->chargenodes, 0, ChargeOctTree::nNodesChargetree * sizeof(float));
}


namespace SCA {
	ChargeOctTree::ChargeOctTree() {
		cudaMalloc(&chargenodes, nNodesChargetree * sizeof(float));
		//cudaMalloc(&halfcharge_forces, n_nonleaf_nodes * sizeof(Float3));
		cudaMalloc(&halfforceNodes, force_n_leafnodes * sizeof(Float3));
		//cudaMalloc(&nodelists, n_nonleaf_nodes * sizeof(NodeList));

		static const std::array<size_t, ChargeOctTree::depthChargetree> depthOffsetsTemp = ChargeOctTree::getDepthOffsets();
		cudaMalloc(&depthOffsets, sizeof(depthOffsetsTemp));
		const auto a = sizeof(depthOffsetsTemp);
		cudaMemcpy(depthOffsets, depthOffsetsTemp.data(), sizeof(depthOffsetsTemp), cudaMemcpyHostToDevice);
	}
	// TODO: Make destructor



	int handleElectrostatics(SimulationDevice* sim_dev, const BoxParams& boxparams) {
		const auto t0 = std::chrono::high_resolution_clock::now();


		for (int depth = ChargeOctTree::depthChargetree - 1; depth >= 3; depth--) {
			// Threads distributed along x axis to maximize coalescence
			const int n_nodes_per_dim = OcttreeHelpers::getNNodesPerDim(depth);
			const uint3 gridDims{ n_nodes_per_dim, n_nodes_per_dim, 1 };
			const uint3 blockDims{ 8, n_nodes_per_dim, 1 };


			if (blockDims.y > 64) throw std::runtime_error("Shared mem can't handle more right now");
			propagateChargesUp<<<gridDims, blockDims >>>(sim_dev->charge_octtree, depth, OcttreeHelpers::getDepthOffset(depth), OcttreeHelpers::getDepthOffset(depth + 1),
				n_nodes_per_dim, OcttreeHelpers::getNNodesPerDim(depth + 1), OcttreeHelpers::getNodelenAtDepth(depth));
			LIMA_UTILS::genericErrorCheck("SCA error during propagate up");
		}

		propagateChargesUp << <gridDims, blockDims >> > (sim_dev->charge_octtree, depth, OcttreeHelpers::getDepthOffset(depth), OcttreeHelpers::getDepthOffset(depth + 1),
			n_nodes_per_dim, OcttreeHelpers::getNNodesPerDim(depth + 1), OcttreeHelpers::getNodelenAtDepth(depth));
		LIMA_UTILS::genericErrorCheck("SCA error during propagate up");

		const auto t1 = std::chrono::high_resolution_clock::now();
		{
			const uint3 gridDims{ ChargeOctTree::nLeafnodesPerDim_Forcetree, ChargeOctTree::nLeafnodesPerDim_Forcetree, ChargeOctTree::nLeafnodesPerDim_Forcetree };
			DownwardSweep << <gridDims, 32 >> > (sim_dev->charge_octtree);
			LIMA_UTILS::genericErrorCheck("SCA error during propagate down");
		}	
		
		CompoundsFetchChargeforce<<<boxparams.n_compounds, MAX_COMPOUND_PARTICLES>>>(sim_dev);

		//distributeChargesCompounds << < boxparams.n_compounds, MAX_COMPOUND_PARTICLES >> > (sim_dev);
		LIMA_UTILS::genericErrorCheck("SCA error during distributeForces");

		wipeCharges(sim_dev->charge_octtree);




		const auto t2 = std::chrono::high_resolution_clock::now();
		//const int compounds_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
		return static_cast<int>(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count());
	}
}