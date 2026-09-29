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
			d->Submit(0, std::make_unique<Rendering::MoleculehullTask>(mhCol, boxSize));
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

SimulationJob Programs::MakeMembraneJob(fs::path workDir, Lipids::Selection composition,
	Float3 boxSize, MembraneGeometry::Figure geometry, int seed, SimParams params, EnvMode mode) {
	SimulationJob job;
	job.workDir = std::move(workDir);
	job.grofile.emplace();
	job.grofile->box_size = boxSize;
	job.topfile.emplace();
	job.simParams = std::move(params);
	job.mode = mode;
	job.preprocess = [composition = std::move(composition), geometry = std::move(geometry), seed](
		GroFile& coordinates, TopologyFile& topology, SimParams&) {
		SimulationBuilder::CreateMembrane(coordinates, topology, composition, geometry, seed);
	};
	return job;
}

SimulationJob Programs::MakeSimulationJob(fs::path workDir, MolecularSystem system,
	SimParams params, EnvMode mode) {
	SimulationJob job;
	job.workDir = std::move(workDir);
	job.grofile = std::move(system.coordinates);
	job.topfile = std::move(system.topology);
	job.simParams = std::move(params);
	job.mode = mode;
	return job;
}

std::vector<Programs::WorkflowInput> Programs::MakeMembraneInputs(
	const std::vector<Lipids::Selection>& compositions, const std::vector<int>& seeds,
	Float3 boxSize, MembraneGeometry::Figure geometry) {
	std::vector<WorkflowInput> inputs;
	inputs.reserve(compositions.size() * seeds.size());
	for (const Lipids::Selection& composition : compositions) {
		const std::string compositionName = Lipids::NameSelection(composition);
		for (const int seed : seeds) {
			const std::string name = compositionName + std::format("_seed{}", seed);
			inputs.push_back({
				.name = name,
				.tags = { { "composition", compositionName }, { "seed", std::to_string(seed) } },
				.MakeJob = [composition, boxSize, geometry, seed](fs::path runDir, SimParams params, EnvMode mode) {
					return MakeMembraneJob(std::move(runDir), composition, boxSize, geometry,
						seed, std::move(params), mode);
				}
			});
		}
	}
	return inputs;
}

Programs::SimulationWorkflow::SimulationWorkflow(fs::path workDir, EnvMode mode)
	: workDir(std::move(workDir)), mode(mode) {}

void Programs::SimulationWorkflow::AddInputs(std::vector<WorkflowInput> newInputs) {
	inputs.insert(inputs.end(), std::make_move_iterator(newInputs.begin()),
		std::make_move_iterator(newInputs.end()));
}

void Programs::SimulationWorkflow::AddStage(WorkflowStage stage) {
	stages.push_back(std::move(stage));
}

void Programs::SimulationWorkflow::CompareDensityProfiles(std::string compositionTag,
	std::string temperatureTag, fs::path output) {
	densityProfileComparison.emplace(
		std::move(compositionTag), std::move(temperatureTag), std::move(output));
}

void Programs::SimulationWorkflow::Run(Environment& environment) {
	if (inputs.empty()) throw std::invalid_argument("A simulation workflow requires at least one input");
	if (stages.empty()) throw std::invalid_argument("A simulation workflow requires at least one stage");

	struct PendingRun {
		std::string name;
		std::map<std::string, std::string> tags;
		fs::path runDir;
		SimulationHandle handle;
	};
	std::vector<PendingRun> pending;
	const WorkflowStage& firstStage = stages.front();
	if (!firstStage.variants.empty())
		throw std::invalid_argument("The first workflow stage cannot have variants");
	pending.reserve(inputs.size());
	for (const WorkflowInput& input : inputs) {
		const fs::path runDir = workDir / input.name;
		SimulationJob job = input.MakeJob(runDir, firstStage.params, mode);
		job.name = input.name;
		job.outputs = firstStage.outputs;
		pending.push_back({ input.name, input.tags, runDir, environment.Submit(std::move(job)) });
	}

	for (size_t stageId = 1; stageId < stages.size(); ++stageId) {
		const WorkflowStage& stage = stages[stageId];
		std::vector<WorkflowVariant> variants = stage.variants;
		if (variants.empty()) variants.push_back({ .name = stage.name });
		std::vector<PendingRun> next;
		next.reserve(pending.size() * variants.size());
		for (PendingRun& parent : pending) {
			SimulationResult completed = parent.handle.Get();
			MolecularSystem parentSystem = completed.FinalSystem();
			for (size_t variantId = 0; variantId < variants.size(); ++variantId) {
				const WorkflowVariant& variant = variants[variantId];
				SimParams params = stage.params;
				if (variant.Configure) variant.Configure(params);
				MolecularSystem system = variantId + 1 == variants.size()
					? std::move(parentSystem) : parentSystem;
				const std::string runName = variant.name.empty() ? stage.name : variant.name;
				const fs::path runDir = parent.runDir / runName;
				SimulationJob job = MakeSimulationJob(runDir, std::move(system), std::move(params), mode);
				job.name = parent.name + " " + runName;
				job.outputs = stage.outputs;
				std::map<std::string, std::string> tags = parent.tags;
				for (const auto& [key, value] : variant.tags)
					tags.insert_or_assign(key, value);
				next.push_back({ job.name, std::move(tags), runDir,
					environment.Submit(std::move(job)) });
			}
		}
		pending = std::move(next);
	}

	for (PendingRun& run : pending) run.handle.Get();
	if (!densityProfileComparison) return;
	if (!stages.back().outputs.contains(OutputSelect::DensityProfile))
		throw std::invalid_argument("Density-profile comparison requires density-profile output");

	std::map<std::pair<std::string, std::string>, std::vector<fs::path>> groupedProfiles;
	for (const PendingRun& run : pending) {
		const auto composition = run.tags.find(densityProfileComparison->compositionTag);
		const auto temperature = run.tags.find(densityProfileComparison->temperatureTag);
		if (composition == run.tags.end() || temperature == run.tags.end())
			throw std::invalid_argument("Density-profile comparison references a missing workflow tag");
		groupedProfiles[{ composition->second, temperature->second }].push_back(
			run.runDir / "density_profile.csv");
	}
	std::vector<SimAnalysis::DensityProfileGroup> groups;
	groups.reserve(groupedProfiles.size());
	for (auto& [key, profiles] : groupedProfiles)
		groups.push_back({ key.first, std::stof(key.second), std::move(profiles) });
	SimAnalysis::CompareDensityProfiles(groups, workDir / densityProfileComparison->output);
}
