#pragma once
#include "Bodies.cuh"
#include "Simulation.cuh"
#include "CompoundBuilder.h"

#include <vector>

class BoxBuilder
{
public:
	BoxBuilder(std::unique_ptr<LimaLogger> logger) :
		m_logger(std::move(logger))
	{
		srand(290128309);
	};
	void buildBox(Simulation* simulation);
	void addCompoundCollection(Simulation* simulation, CompoundCollection& coll);		// Can only use a single "add" function per Simulation for now!!!!!!!!!!!!!
	//void addScatteredMolecules(Simulation* simulation, Compound* molecule, int n_copies);
	//void addDoubleMembrane(Simulation* simulation, Compound* molecule);
	void finishBox(Simulation* simulation);
	int solvateBox(Simulation* simulation);					// Returns # of solvate compounds placed
	int solvateBox(Simulation* simulation, const std::vector<Float3>& solvate_positions);	// Returns # of solvate compounds placed

	// This function expects all ptr's of simulation->box to be pre-allocated on host
	void copyBoxState(Simulation* simulation, Box* boxsrc, const SimParams& simparams_src, uint32_t boxsrc_current_step);

private:
	void insertCompoundInBox(const CompoundFactory& compound, Simulation* simulation);
	
	void setupDataBuffers(Simulation& simulation, const uint64_t n_steps);
	void setupTrainingdataBuffers(Simulation& simulation, const uint64_t n_steps);

	// -------------- Functions for compound manipulation BEFORE integration -------------- //
	//void placeMultipleCompoundsRandomly(Simulation* simulation, Compound* template_compound, int n_copies);
	Compound* randomizeCompound(Compound* template_compound);
	void moveCompound(Compound* compound, Float3 vector);

	//Float3 calcCompoundCom(Compound* compound);
	void rotateCompound(Compound* compound, Float3 xyz_rot);
	//BoundingBox calcCompoundBoundingBox(Compound* compound);
	bool spaceAvailable(const Box& box, Compound* compound);
	bool spaceAvailable(const Box& box, Float3 particle_center, bool verbose=true);	// Ignore radius for this, as it just check against bounding boxes. 
	//bool verifyPairwiseParticleMindist(Compound* a, Compound* b);
	//What about other solvents then? Not problem now while solvents are placed on a grid, but what about later?



	// ------------------------------------------------------------------------------------ //




	// ---------------------------------------------------- Variables ---------------------------------------------------- //
	const float box_len = BOX_LEN;
	const float box_base = 0;

	const std::unique_ptr<LimaLogger> m_logger;

	const float MIN_NONBONDED_DIST = 0.2f;


	// If molecule is offset, each solvent from .gro file must be aswell
	Float3 most_recent_offset_applied = Float3(0.f);	

	
	

	// ---------------------------------------------------- Helper functions ---------------------------------------------------- //
	Float3 get3Random() {	// Returns 3 numbers between 0-1
		return Float3(
			(float) (rand() % RAND_MAX / (double) RAND_MAX),
			(float) (rand() % RAND_MAX / (double) RAND_MAX),
			(float) (rand() % RAND_MAX / (double) RAND_MAX)
		);
	}
	Float3 get3RandomSigned() {	// returns 3 numbers between -0.5->0.5
		return get3Random() - Float3(0.5f);
	}

	float random() {
		return static_cast<float>(rand() % 10000) / 10000.f * 2.f - 1.f;
	}
};

