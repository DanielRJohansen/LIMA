#pragma once

#include "Display.h"
#include "Programs.h"
#include "SimulationBuilder.h"
#include "Forcefield.h"
#include "ConvexHullEngine.cuh"
#include "BoxImageBuilder.h"
#include "TimeIt.h"

#include <glm.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <gtx/rotate_vector.hpp>
#undef GLM_ENABLE_EXPERIMENTAL

void Programs::GetForcefieldParams(const GroFile& grofile, const TopologyFile& topfile, const fs::path& workdir) {
	LIMAForcefield forcefield{topfile.forcefieldInclude->contents};
	
	std::vector<int> ljtypeIndices;
	for (const auto& atom : topfile.GetAllElements<TopologyFile::AtomsEntry>()) {
		ljtypeIndices.push_back(forcefield.GetActiveLjParameterIndex(atom.type));
	}
	ForceField_NB forcefieldNB = forcefield.GetActiveLjParameters();

	// Open csv file for write
	std::ofstream file;
	file.open(workdir / "appliedForcefield.itp");
	if (!file.is_open()) {
		std::cerr << "Could not open file for writing forcefield parameters\n";
		return;
	}

	{
		file << "[ atoms ]\n";
		file << "; type mass sigma[nm] epsilon[J/mol] \n";
		int atomIndex = 0;
		for (auto atom : topfile.GetAllElements<TopologyFile::AtomsEntry>()) {
			const int ljtypeIndex = ljtypeIndices[atomIndex++];
			file << atom.type << " "
				<< forcefieldNB.particle_parameters[ljtypeIndex].sigmaHalf * 2.
				<< " " << forcefieldNB.particle_parameters[ljtypeIndex].epsilonSqrt * forcefieldNB.particle_parameters[ljtypeIndex].epsilonSqrt << "\n";
		}
		file << "\n";
	}

	SimParams params;
	params.em_variant = true;
	auto boximage = LIMA_MOLECULEBUILD::buildMolecules(grofile,	topfile, V1, {}, false, params);

	std::vector<std::string> atomNames;
	for (auto atom : topfile.GetAllElements<TopologyFile::AtomsEntry>()) {
		atomNames.emplace_back(atom.type);
	}

	{
		file << "[ bondtypes ]\n";
		file << "; name_i name_j  b0[nm] kb[J/(mol*nm^2)]\n";
		for (const auto& bond : boximage->topology.singlebonds) {
			for (int i = 0; i < bond.nAtoms; i++)
				file << atomNames[bond.global_atom_indexes[i]] << " ";
			file << bond.params.b0 << " " << bond.params.kb << "\n";
		}
		file << "\n";
	}

	{
		file << "[ angletypes ]\n";
		file << "; name_i name_j name_k theta0[rad] ktheta[J/(mol*rad^2)]\n";
		for (const auto& angle : boximage->topology.anglebonds) {
			for (int i = 0; i < angle.nAtoms; i++)
				file << atomNames[angle.global_atom_indexes[i]] << " ";
			file << angle.params.theta0 << " " << angle.params.kTheta << " " << angle.params.ub0 << " " << angle.params.kUB << "\n";
		}
		file << "\n";
	}

	{
		file << "[ dihedraltypes ]\n";
		file << "; name_i name_j name_k name_l phi0[rad] kphi[J/(mol*rad^2)] multiplicity\n";
		for (const auto& dihedral : boximage->topology.dihedralbonds) {
			for (int i = 0; i < dihedral.nAtoms; i++)
				file << atomNames[dihedral.global_atom_indexes[i]] << " ";
			file << static_cast<float>(dihedral.params.phi_0) << " " << static_cast<float>(dihedral.params.k_phi) << " " << static_cast<float>(dihedral.params.n) << "\n";
		}
		file << "\n";
	}

	{
		file << "[ dihedraltypes ]\n";
		file << "; name_i name_j name_k name_l psi0[rad] kpsi[J/(mol*rad^2)]\n";
		for (const auto& improper : boximage->topology.improperdihedralbonds) {
			for (int i = 0; i < improper.nAtoms; i++)
				file << atomNames[improper.global_atom_indexes[i]] << " ";
			file << improper.params.psi_0 << " " << improper.params.k_psi << "\n";
		}
		file << "\n";
	}

	file.close();
}


void Programs::MoveMoleculesUntillNoOverlap(MoleculeHullCollection& mhCol, Float3 boxSize, bool renderProgress) {

	ConvexHullEngine chEngine{};

	auto d = renderProgress ? std::make_shared<Display>() : nullptr;
	auto renderCallback = [&d, &mhCol, &boxSize]() mutable {
		if (d != nullptr)
			d->Render(std::make_unique<Rendering::MoleculehullTask>(mhCol, boxSize));
	};
	chEngine.MoveMoleculesUntillNoOverlap(mhCol, boxSize, std::ref(renderCallback));

	
	if (renderProgress) {
		//TimeIt::PrintTaskStats("FindIntersect");
		TimeIt::PrintTaskStats("FindIntersectIteration");
	}
}







MoleculeHullCollection Programs::MakeLipidVesicle(GroFile& grofile, TopologyFile& topfile, Lipids::Selection lipidsSelection, float vesicleRadius, Float3 vesicleCenter, std::optional<int> numLipids) {

	const float area = 4.f * PI * vesicleRadius * vesicleRadius;
	const int nLipids = numLipids.value_or(static_cast<int>(area * 0.9f));		

	SimulationBuilder::InsertSubmoleculesOnSphere(grofile, topfile,
		lipidsSelection,
		nLipids, vesicleRadius, vesicleCenter
	);

	std::vector<MoleculeHullFactory> moleculeContainers;

	//for (const auto& molecule : topfile.GetAllSubMolecules()) {
	//	moleculeContainers.push_back({});

	//	for (int globalparticleIndex = molecule.globalIndexOfFirstParticle; globalparticleIndex <= molecule.GlobalIndexOfFinalParticle(); globalparticleIndex++) {
	//		moleculeContainers.back().AddParticle(grofile.atoms[globalparticleIndex].position, grofile.atoms[globalparticleIndex].atomName[0]);
	//	}
	//	
	//	moleculeContainers.back().CreateConvexHull();
	//}


	MoleculeHullCollection mhCol{ moleculeContainers, grofile.box_size };

	return mhCol;
}

void Programs::StaticbodyEnergyMinimize(GroFile& grofile, const TopologyFile& topfile, bool render) {
	std::vector<MoleculeHullFactory> moleculeContainers;
	int globalParticleIndex = 0;

	for (const auto& molecule : topfile.GetSystem().molecules) {
		moleculeContainers.push_back({});

		for (const auto& atom : molecule.moleculetype->atoms) {			
			moleculeContainers.back().AddParticle(grofile.atoms[globalParticleIndex].position, atom.atomname[0]);
			globalParticleIndex++;
		}

		moleculeContainers.back().CreateConvexHull();
	}

	MoleculeHullCollection mhCol{ moleculeContainers, grofile.box_size };

	MoveMoleculesUntillNoOverlap(mhCol, grofile.box_size, render);
}
