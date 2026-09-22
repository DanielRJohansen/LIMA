#pragma once

#include "SimulationData.h"
#include "BatchCompatibility.h"
#include <algorithm>
#include <type_traits>

namespace EngineBatch {
	template<typename T>
	void Append(std::vector<T>& dst, const std::vector<T>& src) {
		dst.insert(dst.end(), src.begin(), src.end());
	}

	inline int CheckedCount(size_t count) {
		if (count > INT_MAX) throw std::invalid_argument("Engine batch exceeds 32-bit work indices");
		return static_cast<int>(count);
	}

	inline void Validate(const std::vector<Simulation*>& simulations) {
		if (simulations.empty()) throw std::invalid_argument("Engine requires a nonempty batch");
		for (size_t i = 0; i < simulations.size(); ++i) {
			const auto* sim = simulations[i];
			if (!sim || !sim->box) throw std::invalid_argument("Null simulation in batch");
			if (std::find(simulations.begin(), simulations.begin() + i, sim) != simulations.begin() + i)
				throw std::invalid_argument("Duplicate simulation in batch");
			if (sim->finished || sim->getStep() != 0)
				throw std::invalid_argument("Engine requires fresh simulations at step zero");
			const auto& params = sim->simParams;
			if (params.stepsPerNlistupdate <= 0 || params.steps_per_temperature_measurement <= 0 || params.data_logging_interval < 0)
				throw std::invalid_argument("Invalid engine measurement or update interval");
			if (params.n_steps.value > INT64_MAX)
				throw std::invalid_argument("Simulation step count exceeds engine limit");
			const auto& box = *sim->box;
			const auto dim = box.boxparams.boxSize;
			if (dim.x <= 0 || dim.y <= 0 || dim.z <= 0 || dim.x >= 1024 || dim.y >= 1024 || dim.z >= 1024)
				throw std::invalid_argument("Unsupported engine box dimensions");
			if (box.persistentClusters.size() != box.persistentClustersMetadata.size() || box.persistentClusters.size() != box.pclusterInterimStates.size()
				|| box.particlesBondedToParticle.size() != box.boxparams.totalParticles || box.pclustersBondedToPcluster.size() != box.persistentClusters.size())
				throw std::invalid_argument("Incomplete simulation particle metadata");
			if (params.data_logging_interval > 0 && (!sim->traj_buffer || !sim->potE_buffer || !sim->vel_buffer || !sim->forceBuffer))
				throw std::invalid_argument("Simulation logging buffers have not been prepared");
			if (sim->box->persistentClusters.empty()) throw std::invalid_argument("Cannot simulate an empty box");
			if (auto mismatch = FindIncompatibility(*simulations.front(), *sim))
				throw std::invalid_argument("Incompatible batch parameter: " + std::string(*mismatch));
		}
	}

	inline void Pack(EngineBatchData& batch, const std::vector<Simulation*>& simulations) {
		Validate(simulations);
		batch.params = simulations.front()->simParams;
		batch.boxSize = simulations.front()->box->boxparams.boxSize;
		batch.ewaldKappa = PhysicsUtils::CalcEwaldkappa(batch.params.cutoff_nm);
		std::vector<PersistentCluster> pclusters;
		std::vector<PersistentClusterMeta> metadata;
		std::vector<PersistentclusterInterimState> states;
		std::vector<int> owners;
		BondGroups bonds;
		for (auto* simulation : simulations) {
			const auto& box = *simulation->box;
			const int simulationId = CheckedCount(batch.simulations.size());
			EngineSimulationData sim;
			sim.simulation = simulation;
			auto& layout = sim.device;
			layout.particles = { batch.nParticles, box.boxparams.totalParticles };
			layout.pclusters = { CheckedCount(pclusters.size()), CheckedCount(box.persistentClusters.size()) };
			layout.bondgroups = { CheckedCount(bonds.groups.size()), CheckedCount(box.bondgroups.size()) };
			layout.gridnodes = { batch.nGridnodes, batch.boxSize.InnerProduct() };
			layout.logOffset = pclusters.size() * PersistentCluster::maxParticles * DatabuffersDeviceController::nStepsInBuffer;
			layout.dt = simulation->simParams.dt;
			layout.uniformElectricField = box.uniformElectricField;
			layout.active = simulation->simParams.n_steps != 0;
			sim.runstatus.simulation_finished = !layout.active;
			batch.nParticles = CheckedCount(size_t(batch.nParticles) + layout.particles.count);
			batch.nGridnodes = CheckedCount(size_t(batch.nGridnodes) + layout.gridnodes.count);
			CheckedCount(size_t(batch.nGridnodes) * PClusterTransfermodule::maxClustersPerBlock);
			Append(pclusters, box.persistentClusters);
			Append(states, box.pclusterInterimStates);
			owners.insert(owners.end(), layout.pclusters.count, simulationId);
			for (auto meta : box.persistentClustersMetadata) {
				for (int lane = 0; lane < PersistentCluster::maxParticles; ++lane) {
					if (meta.particleIdsGlobal[lane] >= 0) meta.particleIdsGlobal[lane] += layout.particles.offset;
					auto& refs = meta.bondgroupReferences[lane];
					for (int i = 0; i < refs.nBondgroupApperances; ++i) {
						refs.bondgroupApperances[i].indexInForceEnergiesBondgroups += CheckedCount(bonds.particles.size());
						refs.bondgroupApperances[i].bondgroupId += layout.bondgroups.offset;
					}
				}
				metadata.push_back(meta);
			}
			for (auto group : box.bondgroups.groups) {
				group.indexOfFirstParticle += CheckedCount(bonds.particles.size());
				group.indexOfFirstSinglebond += CheckedCount(bonds.singlebonds.size());
				group.indexOfFirstPairbond += CheckedCount(bonds.pairbonds.size());
				group.indexOfFirstAnglebond += CheckedCount(bonds.anglebonds.size());
				group.indexOfFirstDihedralbond += CheckedCount(bonds.dihedralbonds.size());
				group.indexOfFirstImproperdihedralbond += CheckedCount(bonds.improperdihedralbonds.size());
				bonds.groups.push_back(group);
			}
			for (auto ref : box.bondgroups.particles) {
				ref.pcid += layout.pclusters.offset;
				bonds.particles.push_back(ref);
			}
			Append(bonds.singlebonds, box.bondgroups.singlebonds);
			Append(bonds.pairbonds, box.bondgroups.pairbonds);
			Append(bonds.anglebonds, box.bondgroups.anglebonds);
			Append(bonds.dihedralbonds, box.bondgroups.dihedralbonds);
			Append(bonds.improperdihedralbonds, box.bondgroups.improperdihedralbonds);
			for (auto neighbors : box.particlesBondedToParticle) {
				neighbors.AddOffset(layout.particles.offset);
				batch.particlesBondedToParticle.push_back(neighbors);
			}
			for (auto neighbors : box.pclustersBondedToPcluster) {
				neighbors.AddOffset(layout.pclusters.offset);
				batch.pclustersBondedToPcluster.push_back(neighbors);
			}
			batch.simulations.push_back(sim);
		}
		batch.nPclusters = CheckedCount(pclusters.size());
		CheckedCount(pclusters.size() * PersistentCluster::maxParticles);
		batch.pClusterDevice.SetData(pclusters);
		batch.pClusterMetaDevice.SetData(metadata);
		batch.pclusterSimulationIds.SetData(owners);
		batch.boxState.pclusterInterimStates = GenericCopyToDevice(states);
		batch.bondgroupDescriptors.SetData(bonds.groups);
		batch.bondgroupParticles.SetData(bonds.particles);
		batch.bondgroupSinglebonds.SetData(bonds.singlebonds);
		batch.bondgroupPairbonds.SetData(bonds.pairbonds);
		batch.bondgroupAnglebonds.SetData(bonds.anglebonds);
		batch.bondgroupDihedralbonds.SetData(bonds.dihedralbonds);
		batch.bondgroupImproperdihedralbonds.SetData(bonds.improperdihedralbonds);
		batch.forceEnergyInterims = std::make_unique<ForceEnergyInterims>(CheckedCount(bonds.particles.size()), batch.nParticles, batch.nPclusters);
		batch.forcesMagnitudeSquareDevice.Expand(batch.nParticles);
		cudaMemset(batch.forcesMagnitudeSquareDevice.Get(), 0, sizeof(float) * batch.nParticles);
		cudaMalloc(&batch.adamState, sizeof(AdamState) * batch.nPclusters * PersistentCluster::maxParticles);
		cudaMemset(batch.adamState, 0, sizeof(AdamState) * batch.nPclusters * PersistentCluster::maxParticles);
	}
}
