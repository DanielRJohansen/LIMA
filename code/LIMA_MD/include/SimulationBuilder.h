/// This file is used to build simulations, that can be saved to .gro and .top files.
/// This is NOT used for loading a simulation into LIMA in any ways
#pragma once

#include "MDFiles.h"
#include "Geometry.cuh"
#include "Lipids.h"
#include "MembraneGeometry.h"

#include <algorithm>


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
	void SolvateGrofile(GroFile& grofile, TopologyFile& topfile, int desiredSolventsPerNm3 = defaultSolventsPerNm3);

	void InsertSubmoleculeInSimulation(GroFile& targetGrofile, TopologyFile& targetTopol,
		GroFile& submolGro, const std::shared_ptr<TopologyFile>& submolTop, Float3 targetCenter);

	void InsertSubmoleculesInSimulation(GroFile& targetGrofile, TopologyFile& targetTopol,
		GroFile& submolGro, const std::shared_ptr<TopologyFile>& submolTop, int nMoleculesToInsert, 
		bool rotateRandomly);
	
	void InsertSubmoleculesOnSphere(
		GroFile& targetGrofile,
		TopologyFile& targetTopol,
		Lipids::Selection,
		int nMoleculesToInsert,
		float sphereRadius,
		const Float3& sphereCenter
	);


	void CreateMembrane(GroFile& grofile, TopologyFile& topfile, const Lipids::Selection& lipidselection,
		const MembraneGeometry::Figure& geometry);
	void CreateMembrane(GroFile& grofile, TopologyFile& topfile, const Lipids::Selection& lipidselection, 
		float membraneCenter);

	// The minimum is derived from the lipid length and enough inner-leaflet area
	// to pack a small, but meaningful, closed surface.
	float MinimumSphereRadius(const Lipids::Selection& lipidselection);
};
