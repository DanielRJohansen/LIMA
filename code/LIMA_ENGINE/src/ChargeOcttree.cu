#include "ChargeOcttree.cuh"
//#include "Simulation.cuh"
#include "EngineUtils.cuh"
//#include "Engine.cuh"
#include "BoundaryCondition.cuh"


#include "Engine.cuh"
#include "SimulationDevice.cuh"

using namespace SCA;


// diff = self-other
__device__ HalfChargeForce calcForce(const Float3& diff, float charge_other) {
	const float dist = diff.len();	// [nm]
	const float magnitude = coulumb_constant * charge_other / (dist * dist);
	if (dist == 0.f) { printf("Yo why is dist 0"); }
	return HalfChargeForce{ diff * magnitude};
}

__global__ void propagateChargesUp(ChargeOctTree* tree, int depth, size_t depth_offset, size_t depth_offset_child, int n_nodes_per_dim, int n_nodes_per_dim_child, float nodelen_nm) {
	// Threads distributed along x axis to maximize coalescence
	const Int3 index_3d{ blockIdx.x * blockDim.x + threadIdx.x, blockIdx.y, blockIdx.z };
	if (index_3d.x > n_nodes_per_dim) { return; }

	const size_t node_index = ChargeOctTree::getAbsIndex(index_3d, n_nodes_per_dim, depth_offset);
	const Float3 pos_abs_self = ChargeOctTree::getAbsPos(index_3d, nodelen_nm);

	float charge = 0.f;

	// Loop over all 8 children
	for (int x_off = 0; x_off < 2; x_off++) {
		for (int y_off = 0; y_off < 2; y_off++) {
			for (int z_off = 0; z_off < 2; z_off++) {
				Int3 query_index3d = index_3d * 2 + Int3{ x_off, y_off, z_off };
				ChargeOctTree::applyPBC(query_index3d, n_nodes_per_dim_child);

				const size_t query_index = ChargeOctTree::getAbsIndex(query_index3d, n_nodes_per_dim_child, depth_offset_child);
				const ChargeNode query_chargenode = tree->nodes[query_index];

				charge += query_chargenode.charge;
			}
		}
	}

	//if (charge > 0.f)
	//	printf("depth %d charge %f\n", depth, charge);

	tree->nodes[node_index].charge = charge;
}


__global__ void propagateChargesDown(ChargeOctTree* tree, int depth, size_t depth_offset, int n_nodes_per_dim, float nodelen_nm)
{
	// Threads distributed along x axis to maximize coalescence
	const Int3 index_3d{ blockIdx.x * blockDim.x + threadIdx.x, blockIdx.y, blockIdx.z};
	if (index_3d.x > n_nodes_per_dim) { return; }

	const size_t node_index = ChargeOctTree::getAbsIndex(index_3d, n_nodes_per_dim, depth_offset);
	const Float3 pos_abs_self = ChargeOctTree::getAbsPos(index_3d, nodelen_nm);

	HalfChargeForce halfchargeforce_local{};

	// An odd side will always have to break into an extra cube "left", and an even side must break down an extra cube "right"
	const int x_even = (blockIdx.x % 2) == 0;
	const int y_even = (blockIdx.y % 2) == 0;
	const int z_even = (blockIdx.z % 2) == 0;

	const int range = depth == 3 ? 3 : 6;

	for (int x_off = -range * 2 + 1 - x_even; x_off <= range * 2 - x_even; x_off++) {
		for (int y_off = -range * 2 + 1 - y_even; y_off <= range * 2 - y_even; y_off++) {
			for (int z_off = -range * 2 + 1 - z_even; z_off <= range * 2 - z_even; z_off++) {

				// If we are inside this border, we need to go to next level of precision
				if (abs(x_off) >= range || abs(y_off) >= range || abs(z_off) >= range) 
				{
					Int3 query_index3d = index_3d + Int3{ x_off, y_off, z_off };
					ChargeOctTree::applyPBC(query_index3d, n_nodes_per_dim);

					const size_t query_index = ChargeOctTree::getAbsIndex(query_index3d, n_nodes_per_dim, depth_offset);
					const ChargeNode query_chargenode = tree->nodes[query_index];
					Float3 pos_abs_query = ChargeOctTree::getAbsPos(query_index3d, nodelen_nm);
					//LIMAPOSITIONSYSTEM::applyHyperpos(pos_abs_self, pos_abs_query);
					LIMAPOSITIONSYSTEM::applyHyperposNM<PeriodicBoundaryCondition>(&pos_abs_self, &pos_abs_query);
					halfchargeforce_local += calcForce(pos_abs_self - pos_abs_query, query_chargenode.charge);
				}


				
			}
		}
	}

	if (halfchargeforce_local.len() > 0.f && depth > 5) {
		//printf("depth %d index %llu\n", depth, node_index);
		//halfchargeforce_local.print('h');
	}

	const HalfChargeForce halfcharge_force_parent = tree->halfcharge_forces[ChargeOctTree::getAbsParentIndex(index_3d, n_nodes_per_dim, depth_offset)];
	// Do a correction of the parent force here, seing as we now know which quadrant of the parent node box we are in

	if (node_index == 165924) {
		//(halfcharge_force_parent).print('O');
	}

	// Push the halfchargeforce to device
	tree->halfcharge_forces[node_index] = halfchargeforce_local + halfcharge_force_parent;
	//tree->halfcharge_forces[node_index] = Float3{ 1.f };
}

void wipeCharges(ChargeOctTree* tree_dev) {
	cudaMemset(tree_dev->nodes, 0, ChargeOctTree::n_nodes * sizeof(float));
}

__global__ void distributeChargesCompounds(SimulationDevice* simdev, ChargeOctTree* chargetree) {

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
	//LIMAPOSITIONSYSTEM::applyPBCNM(&abspos);	// Since any particle might be slightly outside box
	LIMAPOSITIONSYSTEM::applyBCNM<PeriodicBoundaryCondition>(&abspos);

	//abspos.print('p');
	const Float3 charge_force = chargetree->getChargeforceAtLeaf(abspos, compound->atom_charges[threadIdx.x]);
	charge_force.print('f');

	compound->forces_interim[threadIdx.x] += charge_force;
}

namespace SCA {
	ChargeOctTree::ChargeOctTree() {
		cudaMalloc(&nodes, n_nodes * sizeof(ChargeNode));
		cudaMalloc(&halfcharge_forces, n_nonleaf_nodes * sizeof(HalfChargeForce));
		//cudaMalloc(&nodelists, n_nonleaf_nodes * sizeof(NodeList));
	}







	void handleElectrostatics(ChargeOctTree* tree_dev, const BoxParams& boxparams, SimulationDevice* sim_dev) {
		const auto t0 = std::chrono::high_resolution_clock::now();

		const int n_threads = 32;
		for (int depth = ChargeOctTree::tree_depth - 1; depth > 1; depth--) {
			// Threads distributed along x axis to maximize coalescence
			const int n_nodes_per_dim = ChargeOctTree::getNNodesPerDim(depth);
			const uint3 blockdims{ (n_nodes_per_dim + n_threads - 1) / n_threads, n_nodes_per_dim, n_nodes_per_dim };

			propagateChargesUp<<<blockdims, n_threads>>>(tree_dev, depth, ChargeOctTree::getDepthOffset(depth), ChargeOctTree::getDepthOffset(depth + 1),
				n_nodes_per_dim, ChargeOctTree::getNNodesPerDim(depth + 1), ChargeOctTree::getNodelenAtDepth(depth));
			LIMA_UTILS::genericErrorCheck("SCA error during propagate up");
		}

		const auto t1 = std::chrono::high_resolution_clock::now();

		for (int depth = 3; depth < ChargeOctTree::tree_depth; depth++) {
			// Threads distributed along x axis to maximize coalescence
			const int n_nodes_per_dim = ChargeOctTree::getNNodesPerDim(depth);
			const uint3 blockdims{ (n_nodes_per_dim + n_threads - 1) / n_threads, n_nodes_per_dim, n_nodes_per_dim };

			propagateChargesDown << <blockdims, n_threads >> > (tree_dev, depth, ChargeOctTree::getDepthOffset(depth),
				n_nodes_per_dim, ChargeOctTree::getNodelenAtDepth(depth));
			LIMA_UTILS::genericErrorCheck("SCA error during propagate down");
		}

		
		
		distributeChargesCompounds << < boxparams.n_compounds, MAX_COMPOUND_PARTICLES >> > (sim_dev, tree_dev);
		LIMA_UTILS::genericErrorCheck("SCA error during distributeForces");

		wipeCharges(tree_dev);

		const auto t2 = std::chrono::high_resolution_clock::now();





		const int compounds_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
	}
}