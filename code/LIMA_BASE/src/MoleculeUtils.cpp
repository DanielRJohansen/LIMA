#include "MoleculeUtils.h"
#include "BoundaryConditionPublic.h"
#include "MoleculeGraph.h"
#include <unordered_set>
#include <numeric>

Float3 MoleculeUtils::GeometricCenter(const GroFile& grofile) {
	Float3 bbMin{ FLT_MAX }, bbMax{ -FLT_MAX };

	for (const auto& atom : grofile.atoms) {
		bbMin = Float3::ElementwiseMin(bbMin, atom.position);
		bbMax = Float3::ElementwiseMax(bbMax, atom.position);
	}

	return (bbMin + bbMax) / 2;
}

float MoleculeUtils::Radius(const GroFile& grofile, const Float3& center) {
	auto maxDistSq = std::transform_reduce(
		grofile.atoms.begin(), grofile.atoms.end(), 0.f,
		[](float a, float b) { return std::max(a, b); },
		[&center](const auto& atom) { return (atom.position - center).lenSquared(); }
	);

	return std::sqrtf(maxDistSq);
}


void MoleculeUtils::MakeMoleculeWholeAfterPBCFragmentation(GroFile& grofile, const TopologyFile& topfile) {
	LimaMoleculeGraph::MoleculeGraph graph = LimaMoleculeGraph::createGraph(topfile);

	std::unordered_set<int> visited{};

	for (auto& node : graph.BFS(1)) {
		const Float3 nodePostion = grofile.atoms[node.atomid].position;

		for (const auto& neighbor : node.getNeighbors()) {
			if (visited.contains(neighbor->atomid))
				continue;

			Float3& neighborPosition = grofile.atoms[neighbor->atomid].position;
			BoundaryConditionPublic::applyHyperposNM(nodePostion, neighborPosition, grofile.box_size.x, BoundaryConditionSelect::PBC);

			visited.insert(neighbor->atomid);
		}

		visited.insert(node.atomid);
	}
}


void MoleculeUtils::CenterMolecule(GroFile& grofile, const TopologyFile& topfile, std::optional<Float3> targetCenter) {
	MakeMoleculeWholeAfterPBCFragmentation(grofile, topfile);

	const Float3 currentCenter = GeometricCenter(grofile);
	const Float3 diff = targetCenter.value_or(grofile.box_size/2.f) - Float3{currentCenter.x, currentCenter.y, currentCenter.z};

	for (auto& particle : grofile.atoms) {
		particle.position += diff;
	}
}