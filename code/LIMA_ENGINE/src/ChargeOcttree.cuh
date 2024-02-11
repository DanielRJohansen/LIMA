#pragma once 

#include "LimaTypes.cuh"
#include "Constants.h"
#include "Simulation.cuh"

//#include <cufft.h>
#include <math.h>
#include "Utilities.h"

#define ENABLE_SCA


// Superior coulumb algorithm
namespace SCA {
	// The distance between nodes at the lowest level will be no greater than this value
	constexpr float min_grid_spacing = 0.08;	// [nm]

	constexpr float coulumb_constant = 8.99e+9 * 1e-18; // N m^2/C^2 -> [N nm^2 /C^2]

	// Contains direction and magnitude
	using HalfChargeForce = Float3;



	struct ChargeNode {
		// Pos is implicit from the position of the node
		float charge{};
	};

	// All 1d indices are absolute, while the 3d indices are only relative to the current depth!
	class ChargeOctTree {
	public:


		struct NodeList {
			static const int max_nodes = 32;
			size_t n_nodes;
			size_t node_indices[max_nodes];
		};


		ChargeOctTree();

		__device__ void pushChargeToLeaf(const Float3& abspos_of_charge, const float charge) {
			const Int3 index3d = getClosestChargeLeafnodeIndex3d(abspos_of_charge);
			const size_t abs_index = getAbsIndex(index3d, n_leafnodes_per_dim, depth_offset_leafs);
			nodes[abs_index].charge = charge;
			//printf("%f\n", nodes[abs_index].charge);
		}

		__device__ Float3 getChargeforceAtLeaf(const Float3& abspos_of_charge, const float charge) {
			const Int3 index3d = getClosestChargeLeafnodeIndex3d(abspos_of_charge);
			//index3d.print('i');
			const size_t abs_index = getAbsParentIndex(index3d, n_leafnodes_per_dim, depth_offset_leafs);
			//printf("index %llu\n", abs_index);
			const HalfChargeForce& halfchargeforce = halfcharge_forces[abs_index];
			//halfchargeforce.print('h');
			return halfchargeforce * charge;
		}

		//static Int3 get3dIndex(size_t index, int n_nodes_per_dim) {
		//	int z = index / (n_nodes_per_dim * n_nodes_per_dim);
		//	index -= (z * n_nodes_per_dim * n_nodes_per_dim);
		//	int y = index / n_nodes_per_dim;
		//	int x = index % n_nodes_per_dim;
		//	return Int3(x, y, z);
		//}

		__device__ static size_t getAbsIndex(const Int3& index_3d, int n_nodes_per_dim, size_t depth_offset) {
			// Calc index relative to current depth
			const size_t rel_index = static_cast<size_t>(index_3d.z) * n_nodes_per_dim * n_nodes_per_dim
				+ index_3d.y * n_nodes_per_dim
				+ index_3d.x;
			return rel_index + depth_offset;
		}

		__device__ static size_t getAbsParentIndex(const Int3& index_3d, const int n_nodes_per_dim, size_t depth_offset) {
			return getAbsIndex(index_3d / 2, n_nodes_per_dim / 2, depth_offset / 8);
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

		static constexpr size_t getDepthOffset(int depth) {
			return LAL::powi(8, depth - 1);
		}

		static constexpr int getNNodesPerDim(const int depth) {
			return LAL::powi(2, depth);
		}
		static constexpr size_t getNNodesAtDepth(const int depth) {
			return (1 - LAL::powi(8, depth + 1)) / 7;
		}
		static constexpr float getNodelenAtDepth(const int depth) {
			return static_cast<float>(static_cast<double>(BOX_LEN_NM) / LAL::powi(2, depth));
		}


		//static constexpr int tree_depth = static_cast<int>(LAL::ceilf_constexpr(LAL::log2f(BOX_LEN_NM / min_grid_spacing)));
		static constexpr int tree_depth = LAL::ceilFloatToInt(
			LAL::log2f_constexpr(BOX_LEN_NM / min_grid_spacing));
		static constexpr size_t n_nodes = (LAL::powi(8, tree_depth + 1) - 1) / 7;

		ChargeNode* nodes = nullptr;	// Always on device

		// A vector with magnitude carrying only the charge of other.
		// Multiply with charge_self to get actual force on a paticle
		HalfChargeForce* halfcharge_forces = nullptr;
	private:

		__device__ constexpr Int3 getClosestChargeLeafnodeIndex3d(const Float3& abs_pos) {
			return Int3{
				static_cast<int>(abs_pos.x * abspos_to_leafnumber),
				static_cast<int>(abs_pos.y * abspos_to_leafnumber),
				static_cast<int>(abs_pos.z * abspos_to_leafnumber)
			};
		}

		//NodeList* nodelists;

	

		//constexpr int depthFn() { return static_cast<int>(ceilf(-log2f(min_grid_spacing / BOX_LEN_NM))); }

		// Not counting the root node
		//const int depth = static_cast<int>(ceilf(-LAL::log4f(min_grid_spacing / BOX_LEN_NM)));
		static constexpr int n_leafnodes_per_dim = LAL::powi(2, tree_depth);
		static constexpr float abspos_to_leafnumber = static_cast<float>(n_leafnodes_per_dim) / BOX_LEN_NM;

		// Root also gets charge, although it is not relevant
		static constexpr size_t n_nonleaf_nodes = (LAL::powi(8, tree_depth)-1) / 7;
		static constexpr size_t n_leafnodes = (LAL::powi(8, tree_depth));

		const size_t depth_offset_leafs = getDepthOffset(tree_depth);


	};



	void handleElectrostatics(ChargeOctTree* tree, Simulation& sim);
} // End of namespace SPA