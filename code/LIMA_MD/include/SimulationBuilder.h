/// This file is used to build simulations, that can be saved to .gro and .top files.
/// This is NOT used for loading a simulation into LIMA in any ways
#pragma once

#include "MDFiles.h"
#include "Lipids.h"
#include "MembraneGeometry.h"

#include <algorithm>
#include <functional>


struct AtomtypeSelect {
	const TopologyFile::AtomsEntry atomtype;
	const float percentage;
};
using AtomsSelection = std::vector<AtomtypeSelect>;

namespace SimulationBuilder {
	using namespace MDFiles;
	void DistributeParticlesInBox(GroFile& grofile, TopologyFile& topfile, const AtomsSelection& particles,
		float minDistBetweenAnyParticle=0.1f, float particlesPerNm3=32.f);


	// 33.4 is the density of water at 300K, but in some nodes we may have less solvents due to collisions, so we aim a bit higher
	const int defaultSolventsPerNm3 = 34;
	const std::string defaultWaterModel = "spce";
	// waterModel is a path to an .itp with a single moleculetype, or a model name (tip3p, tips3p, spc, spce, ...) resolved
	// to that model's .itp in the system's forcefield dir, or CHARMM27's. The water must be 3 atoms, made rigid by
	// [ settles ] or given by bonds and an angle with explicit parameters, since its geometry is taken from those
	void SolvateGrofile(GroFile& grofile, TopologyFile& topfile, int desiredSolventsPerNm3 = defaultSolventsPerNm3,
		const std::string& waterModel = defaultWaterModel);

	void InsertSubmoleculeInSimulation(GroFile& targetGrofile, TopologyFile& targetTopol,
		GroFile& submolGro, const std::shared_ptr<TopologyFile>& submolTop, Float3 targetCenter);

	void InsertSubmoleculesInSimulation(GroFile& targetGrofile, TopologyFile& targetTopol,
		GroFile& submolGro, const std::shared_ptr<TopologyFile>& submolTop, int nMoleculesToInsert,
		bool rotateRandomly);

	// Copies the molecules of a whole system, centered at targetCenter, keeping the molecule types that keep() accepts.
	// The source forcefield's parameters are appended to the target's. Both keep the target's [ defaults ], so a source
	// with different 1-4 scaling gets the target's (see TODO.md)
	void InsertSystem(GroFile& targetGrofile, TopologyFile& targetTopol, const GroFile& sourceGrofile,
		const TopologyFile& sourceTopol, Float3 targetCenter,
		const std::function<bool(const TopologyFile::Moleculetype&)>& keep);
	
	void InsertSubmoleculesOnSphere(
		GroFile& targetGrofile,
		TopologyFile& targetTopol,
		Lipids::Selection,
		int nMoleculesToInsert,
		float sphereRadius,
		const Float3& sphereCenter
	);


	void CreateMembrane(GroFile& grofile, TopologyFile& topfile, const Lipids::Selection& lipidselection,
		const MembraneGeometry::Figure& geometry, int randomSeed = 0);
	void CreateMembrane(GroFile& grofile, TopologyFile& topfile, const Lipids::Selection& lipidselection, 
		float membraneCenter, int randomSeed = 0);

	// The minimum is derived from the lipid length and enough inner-leaflet area
	// to pack a small, but meaningful, closed surface.
	float MinimumSphereRadius(const Lipids::Selection& lipidselection);
};
