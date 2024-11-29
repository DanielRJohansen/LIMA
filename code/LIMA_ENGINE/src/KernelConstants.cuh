#pragma once

// TODO?: Pre-calculate a solvent-X paired forcefield, to save ALOT of calc in kernels
static_assert(sizeof(ForceField_NB) < 64000, "ForceFieldNB too large for constant memory");



struct BoxSize {
	void Set(int boxSizeNM) {
		//assert(NANO_TO_LIMA_i * boxSizeNM < INT32_MAX);
		boxSizeNM_i = boxSizeNM;
		boxSizeNM_f = static_cast<float>(boxSizeNM);
		blocksPerDim = BoxGrid::NodesPerDim(boxSizeNM);
	}

	int boxSizeNM_i = 0;
	float boxSizeNM_f = 0;
	int blocksPerDim = 0;	// for boxGrid
};


__constant__ ForceField_NB forcefield_device;
__constant__ ForcefieldTinymol tinymolForcefield_device;
__constant__ BoxSize boxSize_device;
__constant__ float cutoffNm_device;
__constant__ float cutoffNmSquaredReciprocal_device;
__constant__ float thermostatScalar_device;

__constant__ NonbondedInteractionParams nonbondedInteractionParams_device[ForceField_NB::MAX_TYPES* ForceField_NB::MAX_TYPES];



//
//inline void SetConstantMem(int boxSizeNM) {
//	BoxSize boxSize;
//	boxSize.Set(boxSizeNM);
//	cudaMemcpyToSymbol(boxSize_device, &boxSize, sizeof(BoxSize));	
//}