#include <chrono>
#include <filesystem>
#include <string>
#include <optional>
#include <numeric>

#include "Environment.h"
#include "MDFiles.h"
#include "CompoundBuilder.h"
#include "Display.h"
#include "BoxBuilder.cuh"
#include "Engine.cuh"
#include "UpgradeableFileFormat.h"
#include "SimulationBuilder.h"
#include "MoleculeUtils.h"

void Environment::InsertMolecule(GroFile& grofile, TopologyFile& topfile, LiveEdit::InsertMolecule& insertionCmd, SimParams simparams) {
	insertionCmd.groPath = FixPath(insertionCmd.groPath);
	insertionCmd.topPath = FixPath(insertionCmd.topPath);

	// First load the new data
	GroFile newmolGro(insertionCmd.groPath);
	auto newmolTop = std::make_shared<TopologyFile>(insertionCmd.topPath);
	MoleculeUtils::CenterMolecule(newmolGro, newmolTop->GetMoleculeType());	// Make molecule whole

	BoundingBox newmolBb(newmolGro.atoms | std::views::transform([](auto& a) { return a.position; }));


	// Save the current state to current files
	WriteBoxCoordinatesToFile(grofile);
	Float3 defaultInsertSite = Float3{ grofile.box_size.x / 2, grofile.box_size.y / 2, grofile.box_size.z - (newmolBb.Dimensions().z / 2.f) };
	Float3 insertionPosition = insertionCmd.position.value_or(defaultInsertSite);

	engine.reset(); // Kill the engine, since it holds a pointer to the sim which we will now change under it. We will make a new engine after creating the new sim
	SimulationBuilder::InsertSubmoleculeInSimulation(grofile, topfile, newmolGro, newmolTop, insertionPosition);
	if (!topfile.forcefieldInclude.has_value())
		topfile.forcefieldInclude = TopologyFile::ForcefieldInclude("charmm27.ff/forcefield.itp");
	CreateSimulation(grofile, topfile, simparams);
	//SimulationBuilder
	//FileUtils::mer
	//MDFiles::
}

void GatherPositionsIntoVector(std::vector<Float3>& dst, CudaBuffer<PersistentCluster>& pcBuffer, int nPclusters) {
	dst.resize(nPclusters * PersistentCluster::maxParticles);
	std::vector<PersistentCluster> pcHost = GenericCopyToHost(pcBuffer.Get(), nPclusters); // TODO: Reuse mem here somehow, this'll be slow..
	for (int pcId = 0; pcId < nPclusters; pcId++) {
		const PersistentCluster& pc = pcHost[pcId];
		for (int pid = 0; pid < PersistentCluster::maxParticles; pid++) {
			dst[pcId * PersistentCluster::maxParticles + pid] = pc.pqd[pid].position;
		}
	}
}

void Environment::HandleMoveMoleculeCommand(const LiveEdit::MoveMolecule& newDragCommand, const LiveEdit::MoveMolecule& prevDragCommand, const std::set<int>& activeSelection, std::vector<Float3>& fixedMovements, std::vector<Rotation>& fixedRotations) {
	if (newDragCommand.draggingForce == prevDragCommand.draggingForce && newDragCommand.rotation == prevDragCommand.rotation) {
		return;
	}

	fixedMovements.clear();
	fixedRotations.clear();
	if (newDragCommand.draggingForce.len() > 0) {
		fixedMovements.resize(simulation->box_host->boxparams.totalParticles);
		for (const auto& id : activeSelection) {
			fixedMovements[id] = newDragCommand.draggingForce * .05f;
		}
	}
	else if (newDragCommand.rotation.len() > 0) {
		Float3 rotationCenter = Float3{ 5.f };
		fixedRotations.resize(simulation->box_host->boxparams.totalParticles);
		for (const auto& id : activeSelection) {
			fixedRotations[id] = Rotation{ rotationCenter, newDragCommand.rotation * 0.01f};
		}
	}

	engine->SetFixedParticleMovementBuffer(fixedMovements);
	engine->SetFixedParticleRotationBuffer(fixedRotations);
}

void UpdateSelection(std::set<int>& selection, LimaMoleculeGraph::MoleculeGraph& molGraph, int pid) {
	if (selection.contains(pid)) {
		return; // This operation will just yield the same set
	}

	selection.clear();
	for (const auto& node : molGraph.BFS(pid)) {
		selection.insert(node.atomid);
	}
}

void UpdateSelection(std::set<int>& selection, const LiveEdit::SelectAtomsBasedOnQualifier& cmd, const Simulation& sim) {
	selection.clear();
	for (int pcid = 0; pcid < sim.box_host->boxparams.totalParticles / PersistentCluster::maxParticles; pcid++) {
		for (int pid = 0; pid < PersistentCluster::maxParticles; pid++) {
			const int gpid = sim.box_host->persistentClustersMetadata[pcid].particleIdsGlobal[pid];
			if (gpid == -1)
				continue;

			switch (cmd.qualifier)
			{
				case LiveEdit::SelectAtomsBasedOnQualifier::Qualifier::All:
					selection.insert(gpid);
					break;
				case LiveEdit::SelectAtomsBasedOnQualifier::Qualifier::Solvent:
					if (sim.box_host->persistentClustersMetadata[pcid].isSolvent) {
						selection.insert(gpid);
					}
					break;
				case LiveEdit::SelectAtomsBasedOnQualifier::Qualifier::Nonsolvent:
					if (!sim.box_host->persistentClustersMetadata[pcid].isSolvent) {
						selection.insert(gpid);
					}
					break;
			default:
				break;
			}
		}
	}
}

void Environment::BuildMembrane(const LiveEdit::BuildMembrane& cmd, GroFile& grofile, TopologyFile& topfile) {
	Lipids::Selection lipidselection;
	for (const auto [name, percentage] : cmd.lipids) {
		lipidselection.emplace_back(Lipids::Select(name, work_dir, percentage));
	}
	float membraneCenterZ = cmd.membraneCenterZ.value_or(grofile.box_size.z / 2.f);
	SimulationBuilder::CreateMembrane(grofile, topfile, lipidselection, membraneCenterZ);
	SimParams simparams = simulation->simparams_host;
	CreateSimulation(grofile, topfile, simparams);
}

void Environment::UpdateForcemask(std::vector<Float3>& forceMaskVec, const std::set<int>& activeSelection, const Float3& newForcemask) {
	forceMaskVec.clear();
	if (newForcemask != Float3{ 0.f } && !activeSelection.empty()) {
		forceMaskVec.resize(simulation->box_host->boxparams.totalParticles, Float3{ 1.f });
		for (const int& pid : activeSelection) {
			forceMaskVec[pid] = newForcemask;
		}
	}
	
	engine->SetForceMask(forceMaskVec);
}

void Environment::LiveEdit(GroFile& grofile, TopologyFile& topfile) {
	simulation->simparams_host.n_steps = 0;
	simulation->simparams_host.data_logging_interval = 0;
	simulation->simparams_host.em_variant = false;

	display = std::make_unique<Display>();
	display->WaitForDisplayReady();
	display->Render(std::make_unique<Rendering::SimulationTask>(
		simulation->box_host->persistentClusters, simulation->box_host->persistentClustersMetadata, simulation->box_host->boxparams, coloringMethod, simStatus
	), false);
	display->allowUserInputs = true;


	bool shouldExit = false;
	std::vector<Float3> positionData;
	bool shouldUpdateRender = true;
	bool canAcceptNewCommand = true;

	std::set<int> activeSelection{};
	//std::optional<Float3> centerOfActiveSelection;

	// MoleculeDragging
	LiveEdit::MoveMolecule prevDragmoleculeCmd{};
	std::vector<int> affectedParticleIds; // TODO: Remove this
	std::vector<Float3> fixedMovements;
	std::vector<Rotation> fixedRotations;
	std::vector<Float3> forceMask;
	

	// Control stepping
	int remainingStepsCount = 0;
	bool runContinous = false;

	auto GetNextCommand = [&]() -> std::optional<LiveEdit::Command> {
		if (!liveEditCommandsQueue.empty()) {
			LiveEdit::Command cmd = liveEditCommandsQueue.front();
			liveEditCommandsQueue.pop_front();
			return cmd;
		}
		return display->GetLiveEditCommand();
		};

	while (true) {
		if (shouldExit) {
			break;
		}

		if (display->DisplaySelfTerminated()) {
			break;
		}

		// Poll interface for new commands, and execute if any		
		if (canAcceptNewCommand) {
			if (auto newCmd = GetNextCommand()) {
				std::visit(
					[&](auto&& cmd) {
						using T = std::decay_t<decltype(cmd)>;

						if constexpr (std::is_same_v<T, LiveEdit::Invalid>) {
							return;
						}
						else if constexpr (std::is_same_v<T, LiveEdit::InsertMolecule>) {
							InsertMolecule(grofile, topfile, cmd, simulation->simparams_host);
							fixedMovements.resize(simulation->box_host->boxparams.totalParticles, Float3{ 0 });
							prevDragmoleculeCmd = LiveEdit::MoveMolecule{};
							forceMask.resize(simulation->box_host->boxparams.totalParticles, Float3{ 1.f }); // expand the forcemask, leaving the existing mask untouched
							display->Render(std::make_unique<Rendering::SimulationTask>(
								simulation->box_host->persistentClusters, simulation->box_host->persistentClustersMetadata, simulation->box_host->boxparams, coloringMethod, simStatus
							));
						}
						else if constexpr (std::is_same_v<T, LiveEdit::MoveMolecule>) {
							HandleMoveMoleculeCommand(cmd, prevDragmoleculeCmd, activeSelection, fixedMovements, fixedRotations);
							prevDragmoleculeCmd = cmd;
							
							if (cmd.draggingForce.len() > 0 || cmd.rotation.len() > 0)
								remainingStepsCount = 50;
							simulation->simparams_host.em_variant = false;
						}
						else if constexpr (std::is_same_v<T, LiveEdit::BuildMembrane>) {
							assert(simulation->box_host->boxparams.totalParticles == 0); // TODO: Change this to a user warning msg or something, and bail
							BuildMembrane(cmd, grofile, topfile);
							simulation->simparams_host.em_variant = true;
							remainingStepsCount = 4000;
							canAcceptNewCommand = false;
							display->Render(std::make_unique<Rendering::SimulationTask>(
								simulation->box_host->persistentClusters, simulation->box_host->persistentClustersMetadata, simulation->box_host->boxparams, coloringMethod, simStatus
							));
						}
						else if constexpr (std::is_same_v<T, LiveEdit::TogglePause>) {
							runContinous = !runContinous;
						}
						else if constexpr (std::is_same_v<T, LiveEdit::AtomSelected>) {
							UpdateSelection(activeSelection, *boximage->systemGraph, cmd.particleId);
							display->UpdateSelection(activeSelection);
						}
						else if constexpr (std::is_same_v<T, LiveEdit::SelectAtomsBasedOnQualifier>) {
							UpdateSelection(activeSelection, cmd, *simulation);
							display->UpdateSelection(activeSelection);
						}
						else if constexpr (std::is_same_v<T, LiveEdit::AddForcemaskToSelection>) {
							UpdateForcemask(forceMask, activeSelection, cmd.forcemask);
						}
						else {
							//static_assert(always_false<T>, "Non-exhaustive visitor!");
						}
					},
					*newCmd
				);
			}
		}



		// Run engine
		if (!engine && simulation->box_host->boxparams.totalParticles > 0) {
			engine = std::make_unique<Engine>(
				simulation.get(),
				simulation->simparams_host.bc_select,
				std::make_unique<LimaLogger>(LimaLogger::compact, m_mode, "engine", work_dir));

			auto& pcBuffer = engine->OffloadPclusterState();
			GatherPositionsIntoVector(positionData, pcBuffer, simulation->box_host->persistentClusters.size());
			shouldUpdateRender = true;
		}
		//printf("Step count %d\n", remainingStepsCount);
		if (engine && (remainingStepsCount > 0 || runContinous)) {
			// Add step logic here
			//shouldUpdateRender = true;
			engine->step();
			auto& pcBuffer = engine->OffloadPclusterState();
			GatherPositionsIntoVector(positionData, pcBuffer, simulation->box_host->persistentClusters.size());
			shouldUpdateRender = true;
			remainingStepsCount--;
			if (remainingStepsCount == 0) {
				// check engine if we should continue..
			}
			UpdateSimstatus(false);
		}

		if (shouldUpdateRender) {
			display->Render(std::make_unique<Rendering::SimulationTaskUpdate>(
				positionData.data(), simStatus
			), false);
			shouldUpdateRender = false;
		}

		if (remainingStepsCount == 0) {
			canAcceptNewCommand = true;
		}
	}
}