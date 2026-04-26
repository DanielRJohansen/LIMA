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


struct LiveEditData {
	// Selection
	std::set<int> activeSelection{};
	std::optional<int> selectedParticleId = std::nullopt;

	// MoleculeDragging
	LiveEdit::MoveMolecule prevDragmoleculeCmd{};
	std::vector<Float3> fixedMovements;
	std::vector<Rotation> fixedRotations;
	std::vector<Float3> forceMask;
	std::vector<Float3> elasticPositions; // Atom will experience a SNF force towards the non-nan components of its elastic position

	// etc
	std::vector<Float3> positionData;

	// Control stepping
	int remainingStepsCount = 0;
	bool runContinous = false;
};


void Environment::InsertMolecule(LiveEditData* liveeditData, GroFile& grofile, TopologyFile& topfile, LiveEdit::InsertMolecule& insertionCmd, SimParams simparams) {
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

	
	if (!liveeditData->fixedMovements.empty())
		liveeditData->fixedMovements.resize(simulation->box->boxparams.totalParticles, Float3{ 0 });
	if (!liveeditData->fixedRotations.empty())
		liveeditData->fixedRotations.resize(simulation->box->boxparams.totalParticles, Rotation{});
	if (!liveeditData->forceMask.empty())
		liveeditData->forceMask.resize(simulation->box->boxparams.totalParticles, Float3{ 1.f }); // expand the forcemask, leaving the existing mask untouched
	if (!liveeditData->elasticPositions.empty())
		liveeditData->elasticPositions.resize(simulation->box->boxparams.totalParticles, Float3(NAN));

	liveeditData->prevDragmoleculeCmd = LiveEdit::MoveMolecule{};
	
	display->Render(std::make_unique<Rendering::SimulationTask>(
		simulation->box->persistentClusters, simulation->box->persistentClustersMetadata, simulation->box->boxparams, simStatus
	));
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

void Environment::HandleMoveMoleculeCommand(LiveEditData* liveeditData, const LiveEdit::MoveMolecule& cmd) {
	if (cmd.draggingForce.len() > 0 || cmd.rotation.len() > 0)
		liveeditData->remainingStepsCount = 1;

	if (cmd.draggingForce == liveeditData->prevDragmoleculeCmd.draggingForce && cmd.rotation == liveeditData->prevDragmoleculeCmd.rotation) {
		return;
	}

	liveeditData->fixedMovements.clear();
	liveeditData->fixedRotations.clear();
	if (cmd.draggingForce.len() > 0) {
		liveeditData->fixedMovements.resize(simulation->box->boxparams.totalParticles);
		for (const auto& id : liveeditData->activeSelection) {
			liveeditData->fixedMovements[id] = cmd.draggingForce * .05f;
		}
	}
	else if (cmd.rotation.len() > 0) {
		Float3 rotationCenter{};
		const int gpid = liveeditData->selectedParticleId.value_or(-1);
		if (gpid != -1 && gpid < liveeditData->positionData.size()) {
			auto [pcid, pid] = boximage->gpidToPcidAndPid[gpid];
			rotationCenter = liveeditData->positionData[pcid * PersistentCluster::maxParticles + pid];
		}

		liveeditData->fixedRotations.resize(simulation->box->boxparams.totalParticles);
		for (const auto& id : liveeditData->activeSelection) {
			liveeditData->fixedRotations[id] = Rotation{ rotationCenter, cmd.rotation * 0.03f};
		}
	}

	simulation->simParams.em_variant = false;
	engine->SetFixedParticleMovementBuffer(liveeditData->fixedMovements);
	engine->SetFixedParticleRotationBuffer(liveeditData->fixedRotations);
	liveeditData->prevDragmoleculeCmd = cmd;
}

void Environment::UpdateSelection(LiveEditData* liveeditData, const LiveEdit::AtomSelected& cmd) {
	liveeditData->selectedParticleId = cmd.particleId;
	if (liveeditData->activeSelection.contains(cmd.particleId)) {
		return; // This operation will just yield the same set
	}

	liveeditData->activeSelection.clear();	
	for (const auto& node : boximage->systemGraph->BFS(cmd.particleId)) {
		liveeditData->activeSelection.insert(node.atomid);
	}
	display->UpdateSelection(liveeditData->activeSelection);
}

void Environment::UpdateSelection(LiveEditData* liveeditData, const LiveEdit::SelectAtomsBasedOnQualifier& cmd) {
	liveeditData->selectedParticleId = std::nullopt;
	liveeditData->activeSelection.clear();
	for (int pcid = 0; pcid < simulation->box->persistentClusters.size(); pcid++) {
		for (int pid = 0; pid < PersistentCluster::maxParticles; pid++) {
			const int gpid = simulation->box->persistentClustersMetadata[pcid].particleIdsGlobal[pid];
			if (gpid == -1)
				continue;

			switch (cmd.qualifier)
			{
				case LiveEdit::SelectAtomsBasedOnQualifier::Qualifier::All:
					liveeditData->activeSelection.insert(gpid);
					break;
				case LiveEdit::SelectAtomsBasedOnQualifier::Qualifier::Solvent:
					if (simulation->box->persistentClustersMetadata[pcid].isSolvent) {
						liveeditData->activeSelection.insert(gpid);
					}
					break;
				case LiveEdit::SelectAtomsBasedOnQualifier::Qualifier::Nonsolvent:
					if (!simulation->box->persistentClustersMetadata[pcid].isSolvent) {
						liveeditData->activeSelection.insert(gpid);
					}
					break;
			default:
				break;
			}
		}
	}
	display->UpdateSelection(liveeditData->activeSelection);
}

void Environment::BuildMembrane(LiveEditData* liveeditData, const LiveEdit::BuildMembrane& cmd, GroFile& grofile, TopologyFile& topfile) {
	Lipids::Selection lipidselection;
	for (const auto [name, percentage] : cmd.lipids) {
		lipidselection.emplace_back(Lipids::Select(name, workDir, percentage));
	}
	float membraneCenterZ = cmd.membraneCenterZ.value_or(grofile.box_size.z / 2.f);
	SimulationBuilder::CreateMembrane(grofile, topfile, lipidselection, membraneCenterZ);
	SimParams simparams = simulation->simParams;	
	CreateSimulation(grofile, topfile, simparams);

	simulation->simParams.em_variant = true;
	liveeditData->remainingStepsCount = 4000;
	display->Render(std::make_unique<Rendering::SimulationTask>(
		simulation->box->persistentClusters, simulation->box->persistentClustersMetadata, simulation->box->boxparams, simStatus
	));

	EM(liveeditData);
}

void Environment::UpdateForcemask(LiveEditData* liveeditData, const LiveEdit::AddForcemaskToSelection& cmd) {
	liveeditData->forceMask.clear();
	if (cmd.forcemask != Float3{ 0.f } && !liveeditData->activeSelection.empty()) {
		liveeditData->forceMask.resize(simulation->box->boxparams.totalParticles, Float3{ 1.f });
		for (const int& pid : liveeditData->activeSelection) {
			liveeditData->forceMask[pid] = cmd.forcemask;
		}
	}
	
	engine->SetForceMask(liveeditData->forceMask);
}

void Environment::UpdateElasticPosition(LiveEditData* liveeditData, const LiveEdit::ElasticPosition& cmd) {
	liveeditData->elasticPositions.clear();
	simulation->simParams.snf_select.erase(SupernaturalForcesSelect::ElasticPosition);
	const bool anyComponentActive = cmd.x || cmd.y || cmd.z;
	if (anyComponentActive && !liveeditData->activeSelection.empty()) {
		liveeditData->elasticPositions.resize(simulation->box->boxparams.totalParticles, Float3{ NAN, NAN, NAN});
		for (const int& pid : liveeditData->activeSelection) {
			Float3 currentPosition = liveeditData->positionData[pid];
			liveeditData->elasticPositions[pid] = Float3 {
				cmd.x ? currentPosition.x : NAN,
				cmd.y ? currentPosition.y : NAN,
				cmd.z ? currentPosition.z : NAN
			};
		}
		simulation->simParams.snf_select.insert(SupernaturalForcesSelect::ElasticPosition);
	}

	engine->SetElasticPositions(liveeditData->elasticPositions);
}

void Environment::EM(LiveEditData* liveeditData) {
	simulation->simParams.em_variant = true;
	liveeditData->remainingStepsCount = 4000;
}

void Environment::LiveEdit(GroFile& grofile, TopologyFile& topfile) {
	simulation->simParams.n_steps = 0;
	simulation->simParams.data_logging_interval = 0;
	simulation->simParams.em_variant = false;

	display = std::make_unique<Display>();
	display->WaitForDisplayReady();
	display->Render(std::make_unique<Rendering::SimulationTask>(
		simulation->box->persistentClusters, simulation->box->persistentClustersMetadata, simulation->box->boxparams, simStatus
	), false);
	display->allowUserInputs = true;

	LiveEditData liveeditData{};

	bool shouldExit = false;
	bool shouldUpdateRender = true;
	

	//bool canAcceptNewCommand = true;
	auto CanAcceptNewCommand = [&]() -> bool {
		return !(simulation->simParams.em_variant && liveeditData.remainingStepsCount > 0);
		};

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
		if (CanAcceptNewCommand()) {
			if (auto newCmd = GetNextCommand()) {
				std::visit(
					[&](auto&& cmd) {
						using T = std::decay_t<decltype(cmd)>;

						if constexpr (std::is_same_v<T, LiveEdit::Invalid>) {
							return;
						}
						else if constexpr (std::is_same_v<T, LiveEdit::InsertMolecule>) {
							InsertMolecule(&liveeditData, grofile, topfile, cmd, simulation->simParams);
						}
						else if constexpr (std::is_same_v<T, LiveEdit::MoveMolecule>) {
							HandleMoveMoleculeCommand(&liveeditData, cmd);
						}
						else if constexpr (std::is_same_v<T, LiveEdit::BuildMembrane>) {
							assert(simulation->box->boxparams.totalParticles == 0); // TODO: Change this to a user warning msg or something, and bail
							BuildMembrane(&liveeditData, cmd, grofile, topfile);
						}
						else if constexpr (std::is_same_v<T, LiveEdit::TogglePause>) {
							liveeditData.runContinous = !liveeditData.runContinous;
							simulation->simParams.em_variant = false;
						}
						else if constexpr (std::is_same_v<T, LiveEdit::AtomSelected>) {
							UpdateSelection(&liveeditData, cmd);
						}
						else if constexpr (std::is_same_v<T, LiveEdit::SelectAtomsBasedOnQualifier>) {
							UpdateSelection(&liveeditData, cmd);
						}
						else if constexpr (std::is_same_v<T, LiveEdit::AddForcemaskToSelection>) {
							UpdateForcemask(&liveeditData, cmd);
						}
						else if constexpr (std::is_same_v<T, LiveEdit::ElasticPosition>) {
							UpdateElasticPosition(&liveeditData, cmd);
						}
						else if constexpr (std::is_same_v<T, LiveEdit::EnergyMinimize>) {
							EM(&liveeditData);
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
		if (!engine && simulation->box->boxparams.totalParticles > 0) {
			engine = std::make_unique<Engine>(
				simulation.get(),
				simulation->simParams.bc_select);

			engine->SetFixedParticleMovementBuffer(liveeditData.fixedMovements);
			engine->SetFixedParticleRotationBuffer(liveeditData.fixedRotations);
			engine->SetForceMask(liveeditData.forceMask);
			engine->SetElasticPositions(liveeditData.elasticPositions);

			auto& pcBuffer = engine->OffloadPclusterState();
			GatherPositionsIntoVector(liveeditData.positionData, pcBuffer, simulation->box->persistentClusters.size());
			shouldUpdateRender = true;
		}
		//printf("Step count %d\n", remainingStepsCount);
		if (engine && (liveeditData.remainingStepsCount > 0 || liveeditData.runContinous)) {
			// Add step logic here
			//shouldUpdateRender = true;
			engine->step();
			UpdateSimstatus(false, true);

			auto& pcBuffer = engine->OffloadPclusterState();
			GatherPositionsIntoVector(liveeditData.positionData, pcBuffer, simulation->box->persistentClusters.size());
			shouldUpdateRender = true;
			liveeditData.remainingStepsCount--;
			if (liveeditData.remainingStepsCount == 0) {
				// check engine if we should continue..
			}
			if (simulation->simParams.em_variant && engine->runstatus.greatestForce < simulation->simParams.em_force_tolerance) {
				simulation->simParams.em_variant = false;
				liveeditData.remainingStepsCount = 0;
			}
		}

		if (shouldUpdateRender) {
			display->Render(std::make_unique<Rendering::SimulationTaskUpdate>(
				liveeditData.positionData.data(), simStatus
			), false);
			shouldUpdateRender = false;
		}
	}
}