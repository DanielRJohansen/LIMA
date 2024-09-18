#pragma once

#include "Programs.h"
#include "SimulationBuilder.h"
#include "Environment.h"
#include "BoxBuilder.cuh"
#include "Forcefield.h"
#include "Statistics.h"
#include "MoleculeGraph.h"
#include "ConvexHullEngine.cuh"

#include <glm.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <gtx/rotate_vector.hpp>
#undef GLM_ENABLE_EXPERIMENTAL

void Programs::SetMoleculeCenter(GroFile& grofile, Float3 targetCenter) {
	
	// First make sure the molecule is not split due to PBC
	for (auto& particle : grofile.atoms) {
		//const Float3 center{ 0, 0, 0 };
		const Float3 center = grofile.atoms[0].position;
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
	ForcefieldManager forcefield{};
	
	std::vector<int> ljtypeIndices;
	for (auto& atom : topfile.GetAllAtoms()) {
		ljtypeIndices.push_back(forcefield.GetActiveLjParameterIndex(topfile.GetForcefieldPaths(), atom.type));
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
			for (int i = 0; i < bond.nAtoms; i++)
				file << atomNames[bond.global_atom_indexes[i]] << " ";
			file << bond.params.b0 * LIMA_TO_NANO << " " << bond.params.kb / LIMA_TO_NANO / LIMA_TO_NANO << "\n";
		}
		file << "\n";
	}

	{
		file << "[ angletypes ]\n";
		file << "; name_i name_j name_k theta0[rad] ktheta[J/(mol*rad^2)]\n";
		for (const auto& angle : boximage->topology.anglebonds) {
			for (int i = 0; i < angle.nAtoms; i++)
				file << atomNames[angle.global_atom_indexes[i]] << " ";
			file << angle.params.theta0 << " " << angle.params.kTheta << " " << angle.params.ub0 * LIMA_TO_NANO << " " << angle.params.kUB / LIMA_TO_NANO / LIMA_TO_NANO << "\n";
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


void Programs::MoveMoleculesUntillNoOverlap(MoleculeHullCollection& mhCol, Float3 boxSize) {
	Display d(Full, boxSize);



	d.Render(std::make_unique<Rendering::MoleculehullTask>(mhCol, boxSize));

	ConvexHullEngine chEngine{};

	auto renderCallback = [&d, &mhCol, &boxSize]() {
		d.Render(std::make_unique<Rendering::MoleculehullTask>(mhCol, boxSize));
	};

	chEngine.MoveMoleculesUntillNoOverlap(mhCol, boxSize, renderCallback);
	
	TimeIt::PrintTaskStats("FindIntersect");
	TimeIt::PrintTaskStats("FindIntersectIteration");

	
	d.Render(std::make_unique<Rendering::MoleculehullTask>(mhCol, boxSize));
}







MoleculeHullCollection Programs::MakeLipidVesicle(GroFile& grofile, TopologyFile& topfile, LipidsSelection lipidsSelection, float vesicleRadius, Float3 vesicleCenter, std::optional<int> numLipids) {

	const float area = 4.f * PI * vesicleRadius * vesicleRadius;
	const int nLipids = numLipids.value_or(static_cast<int>(area * 0.9f));		

	SimulationBuilder::InsertSubmoleculesOnSphere(grofile, topfile,
		lipidsSelection,
		nLipids, vesicleRadius, vesicleCenter
	);

	std::vector<MoleculeHullFactory> moleculeContainers;

	for (const auto& molecule : topfile.GetAllSubMolecules()) {
		moleculeContainers.push_back({});

		for (int globalparticleIndex = molecule.globalIndexOfFirstParticle; globalparticleIndex <= molecule.GlobalIndexOfFinalParticle(); globalparticleIndex++) {
			moleculeContainers.back().AddParticle(grofile.atoms[globalparticleIndex].position, grofile.atoms[globalparticleIndex].atomName[0]);
		}
		
		moleculeContainers.back().CreateConvexHull();
	}


	MoleculeHullCollection mhCol{ moleculeContainers, grofile.box_size };

	return mhCol;
}



MDFiles::FilePair Programs::CreateMembrane(const fs::path& workDir, LipidsSelection& lipidsSelection, Float3 boxSize, float membraneCenterZ, EnvMode envmode) {

	auto [grofile, topfile] = SimulationBuilder::CreateMembrane(lipidsSelection, boxSize, membraneCenterZ);

	Environment env{ workDir, envmode, false};
	SimParams params;
	params.em_variant = true;
	params.bc_select = BoundaryConditionSelect::NoBC;
	params.dt = 1.f;
	params.n_steps = 1000;
	params.snf_select = BoxEdgePotential;
	env.CreateSimulation(*grofile, *topfile, params);
	env.run(false);
	grofile = std::make_shared<GroFile>(env.writeBoxCoordinatesToFile(std::nullopt));

	params.dt = 20;
	env.CreateSimulation(*env.getSim(), params);
	env.run(false);
	grofile = std::make_shared<GroFile>(env.writeBoxCoordinatesToFile(std::nullopt));

	params.snf_select = None;
	params.bc_select = BoundaryConditionSelect::PBC;
	env.CreateSimulation(*env.getSim(), params);
	env.run(false);

	return {grofile, topfile};
}


float CalculateRadius(const std::span<const Float3>& positions) {
	const Float3 mean = Statistics::Mean(positions);

	float maxDistSquared = 0;
	for (const Float3& pos : positions) {
		const float dist = (pos - mean).lenSquared();
		if (dist > maxDistSquared)
			maxDistSquared = dist;
	}
	return std::sqrt(maxDistSquared);
}

float CalculateTotalRadius(const std::vector<Float3>& positions, const std::vector<int>& partition) {
	std::vector<std::span<const Float3>> sections(partition.size());

	int head = 0;

	for (int i = 0; i < partition.size(); i++) {
		sections[i] = std::span<const Float3>(positions.data() + head, partition[i]);
		head += partition[i];
	}
	

	float radiusSum = 0;
	for (const auto& section : sections) {
		radiusSum += CalculateRadius(section);
	}
	return radiusSum;
}

void BruteForcePartition(const std::vector<int>& acceptableSectionSizes, const std::vector<Float3>& positions, int index, std::vector<int>& currentPartition, 
	float& minRadius, std::vector<int>& bestCombination) 
{
	if (index == positions.size()) {
		// Calculate radii for this partition
		double currentRadius = CalculateTotalRadius(positions, currentPartition);

		// Keep track of the minimum radius configuration
		if (currentRadius < minRadius) {
			minRadius = currentRadius;
			bestCombination = currentPartition;
		}
		return;
	}

	// Try each section size from acceptableSectionSizes
	for (int size : acceptableSectionSizes) {
		if (index + size <= positions.size()) {
			currentPartition.push_back(size);
			BruteForcePartition(acceptableSectionSizes, positions, index + size, currentPartition, minRadius, bestCombination);
			currentPartition.pop_back();
		}
	}
}

std::vector<int> DetermineAcceptableSectionsizesBasedOnAtomcount(int nAtoms) {
	int maxSize = 200;

	// SectionSizes based in 32 threads/warp
	std::vector<int> idealSectionsizes {
		64, 63, 62, 61, 60, 59, 58,
		32, 31, 30, 29, 28, 27, 26, 
	};

	// Figure out if the atomcount  can be represented by these ideal sectionsizes
	{
		std::vector<bool> dp(nAtoms + 1, false);
		dp[0] = true; // base case: we can always make the sum 0

		for (int num : idealSectionsizes) {
			for (int i = nAtoms; i >= num; --i) {
				if (dp[i - num]) {
					dp[i] = true;
				}
			}
		}

		if (dp[nAtoms]) {
			return idealSectionsizes;
		}
	}

	// Otherwise, return this list that can represent all atomcounts
	{
		
		return std::vector<int> {
			64, 63, 62, 61, 60, 59, 58, 57, 56,
			32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21
		};
	}	
}


void Programs::ReorderLipidAndDivideIntoCompoundsizedSections(GroFile& grofile, TopologyFile& topfile) {
	std::cout << grofile.m_path << "\n";
	grofile.box_size = Float3{ 6.f };
	SetMoleculeCenter(grofile, grofile.box_size / 2.f);

	Display d(Full, grofile.box_size);	
	d.Render( std::make_unique<Rendering::GrofileTask>( grofile, ColoringMethod::GradientFromAtomid ), false);
	//std::this_thread::sleep_for(std::chrono::milliseconds(500));
	LimaMoleculeGraph::reorderoleculeParticlesAccoringingToSubchains(grofile, topfile);

	topfile.printToFile(std::string("out.itp"));

	d.Render(std::make_unique<Rendering::GrofileTask>(grofile, ColoringMethod::GradientFromAtomid), false);

	// Figure out how to prioritise only using the better sizes, IF that combination can be made
	std::vector<int> sectionSIzes = DetermineAcceptableSectionsizesBasedOnAtomcount(grofile.atoms.size());

	std::vector<Float3> positions(grofile.atoms.size());
	for (int i = 0; i < grofile.atoms.size(); i++) {
		positions[i] = grofile.atoms[i].position;
	}

	float minRadius = FLT_MAX;
	std::vector<int> bestPartition;
	std::vector<int> currentPartition;
	BruteForcePartition(sectionSIzes, positions, 0, currentPartition, minRadius, bestPartition);

	int cummulativeIndex = 0;
	for (int index : bestPartition) {
		topfile.GetLocalAtoms()[cummulativeIndex].section_name = ";lipid_section";
		cummulativeIndex += index;
	}

	Environment env{ grofile.m_path.parent_path(), Headless, false };
	SimParams params;
	params.n_steps = 2;
	params.dt = 0;
	params.data_logging_interval = 1;
	params.em_variant = true;
	env.CreateSimulation(grofile, topfile, params);
	env.run(false);
	auto sim = env.getSim();

	std::vector<Compound> compoundsVec(sim->box_host->boxparams.n_compounds);
	memcpy(compoundsVec.data(), sim->box_host->compounds, sim->box_host->boxparams.n_compounds * sizeof(Compound));
	d.Render(std::make_unique<Rendering::SimulationTask>(sim->traj_buffer->GetBufferAtStep(0), compoundsVec, sim->box_host->boxparams, 0, 0.f, ColoringMethod::GradientFromCompoundId), false);
	//std::this_thread::sleep_for(std::chrono::milliseconds(500));

	//grofile.printToFile();
	//topfile.printToFile();
}

