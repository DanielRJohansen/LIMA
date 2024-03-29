#pragma once

#include "Programs.h"
#include "SimulationBuilder.h"
#include "Environment.h"
#include "BoxBuilder.cuh"

namespace lfs = Filehandler;

void Programs::CreateMembrane(Environment& env, LipidsSelection& lipidselection, bool carryout_em, float centerCoordinate) {
	// Load the files for each lipid
	for (auto& lipid : lipidselection) {
		const std::string lipid_path = env.main_dir + "/resources/Lipids/" + lipid.lipidname + "/";
		lipid.grofile = std::make_unique<ParsedGroFile>(lipid_path + lipid.lipidname + ".gro");
		lipid.topfile = MDFiles::loadTopologyFile(lipid_path + lipid.lipidname + ".itp");
	}

	BoxBuilder boxbuilder( std::make_unique<LimaLogger>());


	// Insert the x lipids with plenty of distance in a non-pbc box
	auto [monolayerGro, monolayerTop] = SimulationBuilder::buildMembrane(lipidselection, env.getSimPtr()->box_host->boxparams.dims);

	// Create simulation and run on the newly created files in the workfolder
	monolayerGro.printToFile(lfs::pathJoin(env.work_dir, "/molecule/monolayer.gro"));
	monolayerTop.printToFile(lfs::pathJoin(env.work_dir, "/molecule/monolayer.top"));

	// Monolayer energy Minimization NoBC
	{
		SimParams ip{};
		ip.bc_select = NoBC;
		ip.n_steps = carryout_em ? 15000 : 0;
		ip.snf_select = HorizontalSqueeze;
		ip.em_variant = true;
		//CreateSimulation(lfs::pathJoin(work_dir, "/molecule/membrane.gro"), lfs::pathJoin(work_dir, "/molecule/membrane.top"), ip);
		env.CreateSimulation(monolayerGro, monolayerTop, ip);

		// Draw each lipid towards the center - no pbc
		env.run();

		if (carryout_em) {
			if (!boxbuilder.verifyAllParticlesIsInsideBox(*env.getSimPtr(), 0.06f)) { return; }	// FAIL
			monolayerGro = env.writeBoxCoordinatesToFile();
		}
	}

	//const ParsedGroFile monolayer_grofile_em = carryout_em ? env.writeBoxCoordinatesToFile() : monolayerGro;
	//auto em_monolayer_grofile = MDFiles::loadGroFile(work_dir + "/molecule/membrane.gro");


	// Copy each particle, and flip them around the xy plane, so the monolayer becomes a bilayer
	auto [bilayerGro, bilayerTop] = SimulationBuilder::makeBilayerFromMonolayer({ monolayerGro, monolayerTop}, env.getSimPtr()->box_host->boxparams.dims);

	bilayerTop.printToFile(lfs::pathJoin(env.work_dir, "/molecule/bilayer.top"));

	// Run EM for a while - with pbc
	{
		SimParams ip{};
		ip.n_steps = carryout_em ? 3000 : 0;
		ip.dt = 50.f;
		ip.bc_select = carryout_em ? PBC : NoBC;	// Cannot insert compounds with PBC, if they are not in box
		env.CreateSimulation(bilayerGro, bilayerTop, ip);

		// Draw each lipid towards the center - no pbc
		env.run();

		if (carryout_em) {
			bilayerGro = env.writeBoxCoordinatesToFile();
		}
	}
	//ParsedGroFile bilayer_grofile_em = carryout_em ? env.writeBoxCoordinatesToFile() : bilayerGro;


	// Ensure the membrane is centered around the centerCoordinate
	{
		double sum = 0;
		for (auto& particle : bilayerGro.atoms) {
			sum += particle.position.z;
		}
		const float diff = centerCoordinate - (sum / static_cast<double>(bilayerGro.atoms.size()));
		for (auto& particle : bilayerGro.atoms) {
			particle.position.z += diff;
		}
	}

	// Save box to .gro and .top file
	bilayerGro.printToFile(lfs::pathJoin(env.work_dir, "/molecule/membrane.gro"));
	//bilayerfiles.first.printToFile(lfs::pathJoin(work_dir, "/molecule/membrane.gro"));	// TEMP
	bilayerTop.printToFile(lfs::pathJoin(env.work_dir, "/molecule/membrane.top"));
}

