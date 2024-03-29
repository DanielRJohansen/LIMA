#pragma once 

#include "LimaTypes.cuh"
#include "Constants.h"
#include "Simulation.cuh"
#include <math.h>
#include "Utilities.h"


class SimulationDevice;
struct BoxParams;


// Important!
// There is no root node in the charge tree, so the first node of depth=1is at index 0
// The reason for this, is that the root note containing 1 float (charge) fucks up the alignment of 
// all subsequent nodes, since we want each node to align with a 64bit cacheline



namespace OcttreeHelpers {

	// Calc index relative to current depth
	__device__ static size_t GetIndexAtDepth(const Int3& index3D, int nNodesPerDim) {
		return static_cast<size_t>(index3D.z) * nNodesPerDim * nNodesPerDim
			+ index3D.y * nNodesPerDim
			+ index3D.x;
	}
	__device__ static size_t GetAbsIndex(const Int3& index_3d, int n_nodes_per_dim, size_t depth_offset) {
		return depth_offset + GetIndexAtDepth(index_3d, n_nodes_per_dim);
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
	static constexpr size_t getDepthOffset(int depth) {
		// No offset at depth 1!
		return depth == 0 
			? 0 
			: (LAL::powi(8, depth) - 1) / 7 - 1; // -1 to account for the lack of a root node
	}
}

// Superior coulumb algorithm
namespace SCA {
	

	__device__ __host__ inline bool DoRecalc(int step) {
		const int STEPS_PER_RECALC = 10;
		return step % STEPS_PER_RECALC == 0;
	}

	// The distance between nodes at the lowest level will be no greater than this value
	constexpr float min_grid_spacing = 0.08f;	// [nm]

	constexpr float coulumb_constant = 8.99e+9 * 1e-18; // N m^2/C^2 -> [N nm^2 /C^2]



	// All 1d indices are absolute, while the 3d indices are only relative to the current depth!
	class ChargeOctTree {
	public:
		ChargeOctTree();

		__device__ void pushChargeToLeaf(const Float3& abspos_of_charge, const float charge) {
			const Int3 index3d = getClosestChargeLeafnodeIndex3d(abspos_of_charge);
			const size_t abs_index = OcttreeHelpers::GetAbsIndex(
				index3d, nLeafnodesPerDim_Chargetree, depth_offset_leafs);
			chargenodes[abs_index] = charge;
		}
		__device__ Float3 getForceAtLeaf(const Float3& abspos_of_charge, const float charge) {
			const Int3 index3d = getClosestForceLeafnodeIndex3d(abspos_of_charge);
			const size_t absIndex = OcttreeHelpers::GetIndexAtDepth(index3d, nLeafnodesPerDim_Forcetree);
			const Float3 halfChargeForce = halfforceNodes[absIndex];
			return halfChargeForce * charge;
		}




		static constexpr int depthChargetree = LAL::ceilFloatToInt(LAL::log2f(BOX_LEN_NM / min_grid_spacing));
		static constexpr int depthForcetree = depthChargetree - 2;	// So ~< 64 particles per leafnode, maybe too many


		//static constexpr size_t nNodesChargetree = (LAL::powi(8, depthChargetree + 1) - 1) / 7-1;	// -1 because there is no root node
		static constexpr size_t nNodesChargetree = OcttreeHelpers::getDepthOffset(depthChargetree + 1);

		// Always on device
		float* chargenodes = nullptr;		// Full tree
		Float3* halfforceNodes = nullptr;	// Only leafnodes

		size_t* depthOffsets;// TODO: move this array to constant memory?

		static constexpr int nLeafnodesPerDim_Forcetree = LAL::powi(2, depthForcetree);
		static constexpr int nLeafnodesPerDim_Chargetree = LAL::powi(2, depthChargetree);


		static constexpr std::array <size_t, depthChargetree+1> getDepthOffsets() {
			std::array<size_t, depthChargetree+1> depth_offsets{};

			for (int i = 0; i < depthChargetree+1; i++) {
				depth_offsets[i] = OcttreeHelpers::getDepthOffset(i);
			}
			return depth_offsets;
		}

		static constexpr float forcetree_leafnodelen = BOX_LEN_NM / nLeafnodesPerDim_Forcetree;



	private:

		__device__ constexpr Int3 getClosestChargeLeafnodeIndex3d(const Float3& abs_pos) {
			return Int3{
				static_cast<int>(abs_pos.x * abspos_to_chargeleafnumber),
				static_cast<int>(abs_pos.y * abspos_to_chargeleafnumber),
				static_cast<int>(abs_pos.z * abspos_to_chargeleafnumber)
			};
		}

		__device__ constexpr Int3 getClosestForceLeafnodeIndex3d(const Float3& abs_pos) {
			return Int3{
				static_cast<int>(abs_pos.x * abspos_to_forceleafnumber),
				static_cast<int>(abs_pos.y * abspos_to_forceleafnumber),
				static_cast<int>(abs_pos.z * abspos_to_forceleafnumber)
			};
		}

		static constexpr float abspos_to_chargeleafnumber = 
			static_cast<float>(nLeafnodesPerDim_Chargetree) / BOX_LEN_NM;
		static constexpr float abspos_to_forceleafnumber = 
			static_cast<float>(nLeafnodesPerDim_Forcetree) / BOX_LEN_NM;

		const size_t depth_offset_leafs = OcttreeHelpers::getDepthOffset(depthChargetree);

		static constexpr size_t force_n_leafnodes = (LAL::powi(8, depthForcetree));

	};
	
	// Instantiates a chargeocttree on the device, and returns the device ptr to the tree
	//ChargeOctTree* MakeChargeOcttree();

	// returns timing [ys]
	int handleElectrostatics(SimulationDevice* sim_dev, const BoxParams& boxparams);
} // End of namespace SPA