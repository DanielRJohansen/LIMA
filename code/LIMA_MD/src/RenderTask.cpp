#include "RenderTask.h"

#include <stdexcept>
#include <utility>

Rendering::AtomRenderTask::AtomRenderTask(const GroFile& grofile, bool shouldShowSolvents)
	: positions(grofile.atoms.size())
	, atoms(grofile.atoms.size())
	, packedPositionIndices(grofile.atoms.size())
	, boxSize(grofile.box_size)
	, showSolvents(shouldShowSolvents)
{
	for (std::size_t atomId = 0; atomId < grofile.atoms.size(); ++atomId) {
		const GroRecord& atom = grofile.atoms[atomId];
		positions[atomId] = atom.position;
		packedPositionIndices[atomId] = static_cast<int>(atomId);
		atoms[atomId].atomLetter = atom.atomName[0];
		atoms[atomId].isSolvent = atom.residueName == "SOL" || atom.residueName == "TIP3";
	}
}

Rendering::AtomRenderTask::AtomRenderTask(
	const std::vector<PersistentCluster>& pclusters,
	const std::vector<PersistentClusterMeta>& pcMeta,
	const BoxParams& boxparams,
	SimStatus initialSimStatus,
	BackboneChains initialBackboneChains)
	: positions(boxparams.totalParticles)
	, atoms(boxparams.totalParticles)
	, packedPositionIndices(boxparams.totalParticles, -1)
	, boxSize(boxparams.BoxSizeFloat())
	, simStatus(initialSimStatus)
	, backboneChains(std::move(initialBackboneChains))
{
	for (std::size_t pcid = 0; pcid < pcMeta.size(); ++pcid) {
		for (int pid = 0; pid < PersistentCluster::maxParticles; ++pid) {
			const int globalParticleId = pcMeta[pcid].particleIdsGlobal[pid];
			if (globalParticleId < 0)
				continue;

			const std::size_t atomId = static_cast<std::size_t>(globalParticleId);
			if (atomId >= atoms.size() || pcid >= pclusters.size())
				throw std::runtime_error("Invalid simulation particle mapping in render task");

			positions[atomId] = pclusters[pcid].pqd[pid].position;
			packedPositionIndices[atomId] = static_cast<int>(pcid * PersistentCluster::maxParticles + pid);
			atoms[atomId] = {
				pcMeta[pcid].atomLetter[pid],
				pclusters[pcid].pqd[pid].params.charge,
				static_cast<int>(pcid),
				pcMeta[pcid].isSolvent
			};
		}
	}
}
