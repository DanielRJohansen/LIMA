#pragma once

#include "Programs.h"
#include "SimulationBuilder.h"
#include "Environment.h"
#include "BoxBuilder.cuh"
#include "Forcefield.h"
#include "Statistics.h"

namespace lfs = Filehandler;

MDFiles::FilePair Programs::CreateMembrane(Environment& env, LipidsSelection& lipidselection, bool carryout_em, float centerCoordinate, bool writeFiles) {
	// Load the files for each lipid
	for (auto& lipid : lipidselection) {
		const std::string lipid_path = env.main_dir + "/resources/Lipids/" + lipid.lipidname + "/";
		lipid.grofile = std::make_unique<GroFile>(lipid_path + lipid.lipidname + ".gro");
		lipid.topfile = std::make_unique<TopologyFile>(lipid_path + lipid.lipidname + ".itp");
	}

	BoxBuilder boxbuilder( std::make_unique<LimaLogger>());


	// Insert the x lipids with plenty of distance in a non-pbc box
	auto [monolayerGro, monolayerTop] = SimulationBuilder::buildMembrane(lipidselection, Float3{ env.getSimPtr()->box_host->boxparams.boxSize });

	// Create simulation and run on the newly created files in the workfolder
	monolayerGro->printToFile(env.work_dir / "molecule/monolayer.gro");
	monolayerTop->printToFile(env.work_dir / "molecule/monolayer.top");

	// Monolayer energy Minimization NoBC
	if (carryout_em) {
		SimParams ip{};
		ip.bc_select = NoBC;
		ip.n_steps = carryout_em ? 20000 : 0;
		ip.snf_select = HorizontalSqueeze;
		ip.em_variant = true;
		env.CreateSimulation(*monolayerGro, *monolayerTop, ip);

		// Draw each lipid towards the center - no pbc
		env.run();
		
		if (!boxbuilder.verifyAllParticlesIsInsideBox(*env.getSimPtr(), 0.06f)) { return { {},{} }; }	// FAIL
		*monolayerGro = env.writeBoxCoordinatesToFile();	
	}

	// Copy each particle, and flip them around the xy plane, so the monolayer becomes a bilayer
	auto [bilayerGro, bilayerTop] = SimulationBuilder::makeBilayerFromMonolayer({ std::move(monolayerGro), std::move(monolayerTop) }, Float3{ env.getSimPtr()->box_host->boxparams.boxSize });

	bilayerTop->printToFile(env.work_dir / "molecule/bilayer.top");

	// Run EM for a while - with pbc
	if (carryout_em) {
		SimParams ip{};
		ip.n_steps = carryout_em ? 3000 : 0;
		ip.dt = 50.f;
		ip.bc_select = carryout_em ? PBC : NoBC;	// Cannot insert compounds with PBC, if they are not in box
		env.CreateSimulation(*bilayerGro, *bilayerTop, ip);

		// Draw each lipid towards the center - no pbc
		env.run();

		*bilayerGro = env.writeBoxCoordinatesToFile();
	}


	// Ensure the membrane is centered around the centerCoordinate
	{
		double sum = 0;
		for (auto& particle : bilayerGro->atoms) {
			sum += particle.position.z;
		}
		const float diff = centerCoordinate - (sum / static_cast<double>(bilayerGro->atoms.size()));
		for (auto& particle : bilayerGro->atoms) {
			particle.position.z += diff;
		}
	}

	// Save box to .gro and .top file
	if (writeFiles) {
		bilayerGro->printToFile(env.work_dir / "molecule/membrane.gro");
		bilayerTop->printToFile(env.work_dir / "molecule/membrane.top");
	}

	return { std::move(bilayerGro), std::move(bilayerTop) };
}

void Programs::SetMoleculeCenter(GroFile& grofile, Float3 targetCenter) {
	
	// First make sure the molecule is not split due to PBC
	for (auto& particle : grofile.atoms) {
		const Float3 center{ 0, 0, 0 };
		BoundaryConditionPublic::applyHyperposNM(center, particle.position, grofile.box_size.x, BoundaryConditionSelect::PBC);
	}

	Double3 sum = { 0,0,0 };
	for (auto& particle : grofile.atoms) 
		sum += particle.position;

	const Double3 currentCenter = sum / static_cast<double>(grofile.atoms.size());
	const Float3 diff = targetCenter - Float3{ currentCenter.x, currentCenter.y, currentCenter.z };

	for (auto& particle : grofile.atoms) {
		particle.position += diff;
	}
}

void Programs::EnergyMinimize(Environment& env, GroFile& grofile, const TopologyFile& topfile, bool solvate, float boxlen_nm) 
{
	SimParams simparams;
	simparams.em_variant = true;
	simparams.n_steps = 4000;
	simparams.dt = 10.f;	// 0.5 ls

	grofile.box_size = Float3{ boxlen_nm };

	env.CreateSimulation(grofile, topfile, simparams);
	env.run();

	GroFile EnergyMinimizedGro = env.writeBoxCoordinatesToFile();
	
	// Now save the new positions to the original gro file, sans any eventually newly added solvents
	assert(grofile.atoms.size() <= EnergyMinimizedGro.atoms.size());
	for (int i = 0; i < grofile.atoms.size(); i++) {
		assert(grofile.atoms[i].atomName == EnergyMinimizedGro.atoms[i].atomName);	// Sanity check the name is unchanged
		grofile.atoms[i].position = EnergyMinimizedGro.atoms[i].position;
	}
}



void Programs::GetForcefieldParams(const GroFile& grofile, const TopologyFile& topfile, const fs::path& workdir) {
	LIMAForcefield forcefield{};

	std::vector<int> ljtypeIndices;
	for (auto& atom : topfile.GetAllAtoms()) {
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
		for (auto atom : topfile.GetAllAtoms()) {
			const int ljtypeIndex = ljtypeIndices[atomIndex++];
			file << atom.type << " "
				<< forcefieldNB.particle_parameters[ljtypeIndex].sigma * LIMA_TO_NANO
				<< " " << forcefieldNB.particle_parameters[ljtypeIndex].epsilon << "\n";
		}
		file << "\n";
	}

	SimParams params;
	params.em_variant = true;
	auto boximage = LIMA_MOLECULEBUILD::buildMolecules(grofile,	topfile, V1, {}, false, params);

	std::vector<std::string> atomNames;
	for (auto atom : topfile.GetAllAtoms()) {
		atomNames.emplace_back(atom.type);
	}

	{
		file << "[ bondtypes ]\n";
		file << "; name_i name_j  b0[nm] kb[J/(mol*nm^2)]\n";
		for (const auto& bond : boximage->topology.singlebonds) {
			for (int i = 0; i < bond.n_atoms; i++)
				file << atomNames[bond.global_atom_indexes[i]] << " ";
			file << bond.params.b0 * LIMA_TO_NANO << " " << bond.params.kb / LIMA_TO_NANO / LIMA_TO_NANO << "\n";
		}
		file << "\n";
	}

	{
		file << "[ angletypes ]\n";
		file << "; name_i name_j name_k theta0[rad] ktheta[J/(mol*rad^2)]\n";
		for (const auto& angle : boximage->topology.anglebonds) {
			for (int i = 0; i < angle.n_atoms; i++)
				file << atomNames[angle.global_atom_indexes[i]] << " ";
			file << angle.params.theta_0 << " " << angle.params.k_theta << "\n";
		}
		file << "\n";
	}

	{
		file << "[ dihedraltypes ]\n";
		file << "; name_i name_j name_k name_l phi0[rad] kphi[J/(mol*rad^2)] multiplicity\n";
		for (const auto& dihedral : boximage->topology.dihedralbonds) {
			for (int i = 0; i < dihedral.n_atoms; i++)
				file << atomNames[dihedral.global_atom_indexes[i]] << " ";
			file << static_cast<float>(dihedral.params.phi_0) << " " << static_cast<float>(dihedral.params.k_phi) << " " << static_cast<float>(dihedral.params.n) << "\n";
		}
		file << "\n";
	}

	{
		file << "[ dihedraltypes ]\n";
		file << "; name_i name_j name_k name_l psi0[rad] kpsi[J/(mol*rad^2)]\n";
		for (const auto& improper : boximage->topology.improperdihedralbonds) {
			for (int i = 0; i < improper.n_atoms; i++)
				file << atomNames[improper.global_atom_indexes[i]] << " ";
			file << improper.params.psi_0 << " " << improper.params.k_psi << "\n";
		}
		file << "\n";
	}

	file.close();
}














void MoveMoleculesUntillNoOverlap(GroFile& grofile, const TopologyFile& topfile) {

}


void Programs::MakeLipidVesicle(GroFile& grofile, TopologyFile& topfile) {

	SimulationBuilder::InsertSubmoleculesInSimulation(grofile, topfile,
		GroFile{ Filehandler::GetLimaDir() / "resources/lipids/POPC/POPC.gro" },
		std::make_shared<TopologyFile>(TopologyFile{ Filehandler::GetLimaDir() / "resources/lipids/POPC/POPC.itp" }),
		grofile.box_size.x);

	std::vector<MoleculeHullFactory> moleculeContainers;

	for (const auto& molecule : topfile.GetAllSubMolecules()) {
		moleculeContainers.push_back({});
	

		for (int globalparticleIndex = molecule.globalIndexOfFirstParticle; globalparticleIndex <= molecule.GlobalIndexOfFinalParticle(); globalparticleIndex++) {
			moleculeContainers.back().AddParticle(grofile.atoms[globalparticleIndex].position, grofile.atoms[globalparticleIndex].atomName[0]);
		}
		
				

		std::vector<Float3> positions{ moleculeContainers.back().GetParticles().begin(), moleculeContainers.back().GetParticles().end()};


		moleculeContainers.back().convexHull = {positions};
	}

	Display d(Full);

	MoleculeHullCollection mhCol{ moleculeContainers, grofile.box_size };

	while (true) {
		d.checkWindowStatus();
		d.Render(mhCol, grofile.box_size);
	}


}