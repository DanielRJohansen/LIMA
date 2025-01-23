#pragma once

// TODO?: Pre-calculate a solvent-X paired forcefield, to save ALOT of calc in kernels
static_assert(sizeof(ForceField_NB) < 64000, "ForceFieldNB too large for constant memory");



struct BoxSize {
	void Set(int boxSizeNM) {
		//assert(NANO_TO_LIMA_i * boxSizeNM < INT32_MAX);
		boxSizeNM_i = boxSizeNM;
		boxSizeNM_f = static_cast<float>(boxSizeNM);
		blocksPerDim = BoxGrid::NodesPerDim(boxSizeNM);
                blocksPerDimHalf = blocksPerDim/2;
	}

	int boxSizeNM_i = 0;
	float boxSizeNM_f = 0;
	int blocksPerDim = 0;	// for boxGrid
        int blocksPerDimHalf = 0;
};

namespace DeviceConstants {

	__constant__ ForceField_NB forcefield;
	__constant__ ForcefieldTinymol tinymolForcefield;
	__constant__ BoxSize boxSize;
	__constant__ float cutoffNM;
	__constant__ float cutoffNmReciprocal;
	__constant__ float cutoffNmSquaredReciprocal;
	__constant__ float ewaldKappa;

	__constant__ float thermostatScalar;

	__constant__ NonbondedInteractionParams nonbondedinteractionParams[ForceField_NB::MAX_TYPES * ForceField_NB::MAX_TYPES];


	// Precomputed values

	static constexpr int BSPLINE_LUT_SIZE = 128;
	// We store only w0, w1 => total size = 2*N.
	__constant__ float bsplineTable[2 * (BSPLINE_LUT_SIZE)]; // precomputed 4th order bsplines [0,1]


        static constexpr int ERFC_LUT_SIZE = 32;
	__constant__ float erfcForcescalarTable[ERFC_LUT_SIZE]; // precomputed scalers [0, 1], where 1=cutoffNM
	__constant__ float erfcPotentialscalarTable[ERFC_LUT_SIZE];
}
