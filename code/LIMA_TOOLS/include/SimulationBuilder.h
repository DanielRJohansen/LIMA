/// This file is used to build simulations, that can be saved to .gro and .top files.
/// This is NOT used for loading a simulation into LIMA in any ways
#pragma once

#include "MDFiles.h"
#include "Geometry.cuh"
#include "Lipids.h"

#include <algorithm>


struct AtomtypeSelect {
	const TopologyFile::AtomsEntry atomtype;
	const float percentage;
};
using AtomsSelection = std::vector<AtomtypeSelect>;

namespace SimulationBuilder {
	using namespace MDFiles;
	using Geometry::Plane;

	void DistributeParticlesInBox(GroFile& grofile, TopologyFile& topfile, const AtomsSelection& particles,
		float minDistBetweenAnyParticle=0.1f, float particlesPerNm3=32.f);


	void SolvateGrofile(GroFile& grofile);

	/// <summary>
	/// Insert submolecules into a possibly empty simulation. The molecules will be spread randomly around in the box
	/// The input grofile does NOT need to be centered around origo, but the grofile must have the atoms
	/// in the correct position relative to eachother, as it uses raw positions, instead of hyperpos relative to index 0
	/// </summary>
	void InsertSubmoleculesInSimulation(GroFile& targetGrofile, TopologyFile& targetTopol,
		const GroFile& submolGro, const std::shared_ptr<TopologyFile>& submolTop, int nMoleculesToInsert);
	
	void InsertSubmoleculesOnSphere(
		GroFile& targetGrofile,
		TopologyFile& targetTopol,
		Lipids::Selection,
		int nMoleculesToInsert,
		float sphereRadius,
		const Float3& sphereCenter
	);


	FilePair CreateMembrane(const Lipids::Selection& lipidselection, Float3 boxSize, float membraneCenter);
	void CreateMembrane(GroFile& grofile, TopologyFile& topfile, const Lipids::Selection& lipidselection, Float3 boxSize, float membraneCenter);
};
