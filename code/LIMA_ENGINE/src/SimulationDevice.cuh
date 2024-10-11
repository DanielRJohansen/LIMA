#pragma once

// TOdo move impl to .cpp file

#include "Bodies.cuh"
#include "Simulation.cuh"
#include "EngineBodies.cuh"
#include "ChargeOcttree.cuh"

//using namespace SCA;




struct BoxDevice {
	BoxDevice() {}
	void CopyDataToHost(Box& boxDev);
	void DeleteBox();

	BoxParams boxparams;

	Compound* compounds = nullptr;
	CompoundCoords* compoundcoordsCircularQueue = nullptr;

	Solvent* solvents = nullptr;
	SolventBlocksCircularQueue* solventblockgrid_circularqueue = nullptr;

	CompoundBridgeBundleCompact* bridge_bundle = nullptr;

	// BondedParticlesLUT data - NEVER access directly, use the bpLUTHelpers namespace
	BondedParticlesLUT* bpLUTs = nullptr;

	UniformElectricField uniformElectricField;	
};



/// <summary>
/// All members of this will only ever exist on device. Immediately after creating this class
/// it must also be moved to device.
/// </summary>
struct SimulationDevice {
	SimulationDevice(const SimulationDevice&) = delete;

	SimulationDevice(const SimParams& params_host, Box* box_host);

	// Recursively free members. Use cudaFree on *this immediately after
	void deleteMembers();

	// Compounds signal where they are on a grid, handled by NLists. Used by solvents to load nearby compounds.
	CompoundGridNode* compound_grid = nullptr;

	// Compounds can see which compounds are near them
	NeighborList* compound_neighborlists = nullptr;

	// Module used to move solvents to a new block, in parallel
	SolventBlockTransfermodule* transfermodule_array = nullptr;

	SimParams* params;
	SimSignals* signals;

	BoxDevice* box;
	DatabuffersDevice* databuffers;

	//ChargeOctTree* charge_octtree;
	Electrostatics::ChargeNode* chargeGrid = nullptr;
	float* chargeGridChargeSums = nullptr;

	// potE should be divided equally between all the particles in the node
	ForceAndPotential* chargeGridOutputForceAndPot = nullptr; // {Float3 force [1/l N/mol], float potE [J/mol]}

};