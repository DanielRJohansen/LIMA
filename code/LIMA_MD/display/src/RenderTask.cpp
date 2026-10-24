#include "RenderTask.h"

#include "BoxBuilder.cuh"
#include "BoxImageBuilder.h"
#include "Format.h"
#include "RenderUtilities.cuh"

#include <algorithm>
#include <stdexcept>
#include <array>
#include <utility>
#include <numeric>
#include <ranges>
#include <unordered_map>

std::vector<Rendering::MoleculeInfo> Rendering::GetMoleculeInfo(const BoxImage& boxImage) {
	std::unordered_map<std::string, size_t> typeCounts;
	std::vector<MoleculeInfo> molecules;
	for (const auto& instance : boxImage.topology.moleculeInstances) {
		MoleculeInfo molecule;
		molecule.name = instance.type->name;
		molecule.number = ++typeCounts[molecule.name];
		molecule.atomIds.resize(instance.type->atoms.size());
		std::iota(molecule.atomIds.begin(), molecule.atomIds.end(), instance.particleOffset);
		molecules.push_back(std::move(molecule));
	}
	for (auto& molecule : molecules)
		molecule.typeCount = typeCounts[molecule.name];
	return molecules;
}

std::vector<Rendering::MoleculeInfo> Rendering::GetMoleculeInfo(const Simulation& simulation) {
	return simulation.boxImage ? GetMoleculeInfo(*simulation.boxImage) : std::vector<MoleculeInfo>{};
}

std::vector<Rendering::MoleculeInfo> Rendering::GetMoleculeInfo(const GroFile& grofile) {
	std::unordered_map<std::string, size_t> typeCounts;
	std::vector<MoleculeInfo> molecules;
	std::size_t begin = 0;
	while (begin < grofile.atoms.size()) {
		std::size_t end = begin + 1;
		while (end < grofile.atoms.size()
			&& grofile.atoms[end].residue_number == grofile.atoms[begin].residue_number
			&& grofile.atoms[end].residueName == grofile.atoms[begin].residueName)
			++end;
		MoleculeInfo molecule;
		molecule.name = std::string{ grofile.atoms[begin].residueName.View() };
		molecule.number = ++typeCounts[molecule.name];
		molecule.atomIds.resize(end - begin);
		std::iota(molecule.atomIds.begin(), molecule.atomIds.end(), static_cast<int>(begin));
		molecules.push_back(std::move(molecule));
		begin = end;
	}
	for (auto& molecule : molecules)
		molecule.typeCount = typeCounts[molecule.name];
	return molecules;
}

std::optional<Rendering::MoleculeInfo> Rendering::GetMoleculeInfo(const BoxImage& boxImage, int atomId) {
	const auto molecules = GetMoleculeInfo(boxImage);
	const auto molecule = std::ranges::find_if(molecules, [atomId](const MoleculeInfo& candidate) {
		return std::ranges::find(candidate.atomIds, atomId) != candidate.atomIds.end();
	});
	if (molecule == molecules.end())
		return std::nullopt;
	return *molecule;
}

Rendering::AtomRenderTask::AtomRenderTask(const GroFile& grofile, bool shouldShowSolvents)
	: positions(grofile.atoms.size())
	, atoms(grofile.atoms.size())
	, packedPositionIndices(grofile.atoms.size())
	, boxSize(grofile.box_size)
	, lipidModels(BoxBuilder::InterpretLipids(grofile))
	, particleCategories(BoxBuilder::InterpretParticleCategories(grofile, lipidModels))
	, molecules(GetMoleculeInfo(grofile))
	, showSolvents(shouldShowSolvents)
{
	for (std::size_t atomId = 0; atomId < grofile.atoms.size(); ++atomId) {
		const GroRecord& atom = grofile.atoms[atomId];
		positions[atomId] = atom.position;
		packedPositionIndices[atomId] = static_cast<int>(atomId);
		atoms[atomId].atomLetter = RenderUtilities::AtomLetter(atom.atomName.View(), atom.residueName.View());
		atoms[atomId].isSolvent = atom.residueName == "SOL" || atom.residueName == "TIP3";
	}
}

Rendering::AtomRenderTask::AtomRenderTask(
	const std::vector<PersistentCluster>& pclusters,
	const std::vector<PersistentClusterMeta>& pcMeta,
	const BoxParams& boxparams,
	SimStatus initialSimStatus,
	BackboneChains initialBackboneChains,
	std::vector<MoleculeInfo> initialMolecules,
	LipidModels initialLipidModels,
	ParticleCategories initialParticleCategories)
	: positions(boxparams.totalParticles)
	, atoms(boxparams.totalParticles)
	, packedPositionIndices(boxparams.totalParticles, -1)
	, boxSize(boxparams.BoxSizeFloat())
	, simStatus(initialSimStatus)
	, backboneChains(std::move(initialBackboneChains))
	, lipidModels(std::move(initialLipidModels))
	, particleCategories(std::move(initialParticleCategories))
	, molecules(std::move(initialMolecules))
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

namespace {
	std::string CountLabel(std::size_t count, const char* unit) {
		return Lima::Format("{} {}{}", count, unit, count == 1 ? "" : "s");
	}

	const char* CategoryName(ParticleCategory category) {
		switch (category) {
		case ParticleCategory::Protein: return "protein";
		case ParticleCategory::Lipid: return "lipids";
		case ParticleCategory::Water: return "water";
		case ParticleCategory::Ion: return "ions";
		case ParticleCategory::Other: return "other";
		}
		return "other";
	}
}

std::vector<Rendering::SelectionOption> Rendering::SelectionOptions(const AtomRenderTask& task, int atomId, bool hasSelection)
{
	std::vector<SelectionOption> options;
	const auto AddAtoms = [&options](std::string label, std::vector<int> atomIds, std::size_t moleculeCount = 0) {
		std::string count = moleculeCount > 0 ? CountLabel(moleculeCount, "molecule") : CountLabel(atomIds.size(), "atom");
		MoleculeInfo selection{ label, 0, moleculeCount, std::move(atomIds) };
		options.push_back({ std::move(label), std::move(count), std::move(selection) });
	};

	if (atomId >= 0 && atomId < static_cast<int>(task.atoms.size())) {
		const char element = task.atoms[atomId].atomLetter;
		options.push_back({ Lima::Format("This atom (#{})", atomId), std::string(1, element),
			MoleculeInfo{ Lima::Format("Atom #{}", atomId), 0, 0, { atomId } } });

		// Ions are lettered I, which is no element
		if (element != 'I' && element != ' ') {
			std::vector<int> sameElement;
			for (int id = 0; id < static_cast<int>(task.atoms.size()); ++id)
				if (task.atoms[id].atomLetter == element)
					sameElement.push_back(id);
			AddAtoms(Lima::Format("All {} atoms", element), std::move(sameElement));
		}

		const auto molecule = std::ranges::find_if(task.molecules, [atomId](const MoleculeInfo& candidate) {
			return std::ranges::find(candidate.atomIds, atomId) != candidate.atomIds.end();
		});
		if (molecule != task.molecules.end()) {
			options.push_back({ "This molecule", Lima::Format("{} #{}", molecule->name, molecule->number), *molecule });
			std::vector<int> sameType;
			for (const MoleculeInfo& candidate : task.molecules)
				if (candidate.name == molecule->name)
					sameType.insert(sameType.end(), candidate.atomIds.begin(), candidate.atomIds.end());
			AddAtoms(Lima::Format("All {}", molecule->name), std::move(sameType), molecule->typeCount);
		}
	}

	// Only categories present in the system. Kept even when one is everything, as "All lipids" is what you look for after clicking a lipid
	constexpr std::array categories{ ParticleCategory::Protein, ParticleCategory::Lipid,
		ParticleCategory::Water, ParticleCategory::Ion, ParticleCategory::Other };
	std::array<std::vector<int>, categories.size()> categoryAtoms;
	for (int id = 0; id < static_cast<int>(task.particleCategories.size()) && id < static_cast<int>(task.atoms.size()); ++id)
		categoryAtoms[static_cast<std::size_t>(std::ranges::find(categories, task.particleCategories[id]) - categories.begin())].push_back(id);
	for (std::size_t i = 0; i < categories.size(); ++i) {
		if (categoryAtoms[i].empty())
			continue;
		const std::size_t lipidCount = categories[i] == ParticleCategory::Lipid ? task.lipidModels.size() : 0;
		AddAtoms(Lima::Format("All {}", CategoryName(categories[i])), std::move(categoryAtoms[i]), lipidCount);
	}

	std::vector<int> all(task.atoms.size());
	std::iota(all.begin(), all.end(), 0);
	AddAtoms("All atoms", std::move(all));

	if (hasSelection)
		options.push_back({ "Clear selection", "", MoleculeInfo{} });
	return options;
}
