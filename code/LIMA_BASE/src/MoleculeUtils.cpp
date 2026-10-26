#include "MoleculeUtils.h"
#include "BoundaryConditionPublic.h"
#include "MoleculeGraph.h"

#include <unordered_set>
#include <numeric>
#include <algorithm>
#include <functional>
#include <cfloat>
#include <queue>
#include <stdexcept>

namespace {
	void MakeMoleculeWholeAtOffset(
		GroFile& grofile,
		const TopologyFile::Moleculetype& moltype,
		const std::size_t atomOffset
	) {
		const std::size_t atomCount = moltype.atoms.size();
		if (atomOffset + atomCount > grofile.atoms.size()) {
			throw std::runtime_error("Topology contains more atoms than the GRO file");
		}

		const LimaMoleculeGraph::MoleculeGraph graph{ moltype };
		std::unordered_set<int> visited;

		// A molecule type can contain more than one connected component. Unwrap each
		// component independently instead of assuming atom 1 belongs to the only one.
		for (const auto& [rootId, root] : graph.nodes) {
			if (!visited.insert(rootId).second) continue;

			std::queue<int> pending;
			pending.push(rootId);
			while (!pending.empty()) {
				const int nodeId = pending.front();
				pending.pop();
				if (nodeId < 0 || static_cast<std::size_t>(nodeId) >= atomCount) {
					throw std::runtime_error("Topology atom id is outside its molecule's atom range");
				}

				const Float3 nodePosition = grofile.atoms[atomOffset + nodeId].position;
				for (const int neighborId : graph.nodes.at(nodeId).getNeighbors()) {
					if (!visited.insert(neighborId).second) continue;
					if (neighborId < 0 || static_cast<std::size_t>(neighborId) >= atomCount) {
						throw std::runtime_error("Topology atom id is outside its molecule's atom range");
					}

					Float3& neighborPosition = grofile.atoms[atomOffset + neighborId].position;
					BoundaryConditionPublic::applyHyperposNM(
						nodePosition, neighborPosition, grofile.box_size, BoundaryConditionSelect::PBC);
					pending.push(neighborId);
				}
			}
		}
	}
}

Float3 MoleculeUtils::GeometricCenter(const GroFile& grofile) {
	Float3 bbMin{ FLT_MAX };
	Float3 bbMax{ -FLT_MAX };

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


float MoleculeUtils::Diameter(const GroFile& grofile) {
	if (grofile.atoms.empty())
		return 0.f;
	// Exact, but prunes the O(n^2) search: a pair is at most ri + rj apart, with r the distance to the centroid. With the
	// atoms sorted by r, both loops stop once no remaining pair can beat the longest distance found
	Float3 centroid{};
	for (const auto& atom : grofile.atoms)
		centroid += atom.position;
	centroid = centroid / static_cast<float>(grofile.atoms.size());

	struct Atom { Float3 position; float r; };
	std::vector<Atom> atoms;
	atoms.reserve(grofile.atoms.size());
	for (const auto& atom : grofile.atoms)
		atoms.push_back({ atom.position, (atom.position - centroid).len() });
	std::ranges::sort(atoms, std::greater{}, &Atom::r);

	float diameterSq = 0.f;
	for (size_t i = 0; i < atoms.size(); i++) {
		if ((atoms[i].r + atoms[0].r) * (atoms[i].r + atoms[0].r) <= diameterSq)
			break;
		for (size_t j = i + 1; j < atoms.size(); j++) {
			if ((atoms[i].r + atoms[j].r) * (atoms[i].r + atoms[j].r) <= diameterSq)
				break;
			diameterSq = std::max(diameterSq, (atoms[i].position - atoms[j].position).lenSquared());
		}
	}
	return std::sqrt(diameterSq);
}

void MoleculeUtils::SetBoxAndCenter(GroFile& grofile, Float3 boxSize) {
	if (boxSize.x <= 0.f || boxSize.y <= 0.f || boxSize.z <= 0.f)
		throw std::invalid_argument("Box size must be positive");
	grofile.box_size = boxSize;
	const Float3 translation = boxSize / 2.f - GeometricCenter(grofile);
	for (auto& atom : grofile.atoms)
		atom.position += translation;
}

void MoleculeUtils::MakeMoleculeWholeAfterPBCFragmentation(GroFile& grofile, const TopologyFile::Moleculetype& moltype) {
	MakeMoleculeWholeAtOffset(grofile, moltype, 0);
}

void MoleculeUtils::MakeMoleculeWholeAfterPBCFragmentation(GroFile& grofile, const TopologyFile& topfile) {
	std::size_t atomOffset = 0;
	for (const auto& molecule : topfile.GetSystem().Instances()) {
		MakeMoleculeWholeAtOffset(grofile, *molecule.moleculetype, atomOffset);
		atomOffset += molecule.moleculetype->atoms.size();
	}

	if (atomOffset != grofile.atoms.size()) {
		throw std::runtime_error("GRO file contains atoms not described by the topology");
	}
}

void MoleculeUtils::FitMoleculeInBox(GroFile& grofile, const float padding) {
	if (padding < 0.f) throw std::invalid_argument("Box padding cannot be negative");
	if (grofile.atoms.empty()) throw std::runtime_error("Cannot fit an empty GRO file in a box");

	Float3 bbMin{ FLT_MAX };
	Float3 bbMax{ -FLT_MAX };
	for (const auto& atom : grofile.atoms) {
		bbMin = Float3::ElementwiseMin(bbMin, atom.position);
		bbMax = Float3::ElementwiseMax(bbMax, atom.position);
	}

	const Float3 translation = Float3{ padding } - bbMin;
	for (auto& atom : grofile.atoms) {
		atom.position += translation;
	}
	grofile.box_size = bbMax - bbMin + Float3{ 2.f * padding };
}


void MoleculeUtils::CenterMolecule(GroFile& grofile, const TopologyFile::Moleculetype& topfile, std::optional<Float3> targetCenter) {
	MakeMoleculeWholeAfterPBCFragmentation(grofile, topfile);

	const Float3 currentCenter = GeometricCenter(grofile);
	const Float3 diff = targetCenter.value_or(grofile.box_size/2.f) - Float3{currentCenter.x, currentCenter.y, currentCenter.z};

	for (auto& particle : grofile.atoms) {
		particle.position += diff;
	}
}

void MoleculeUtils::RotateMolecule(GroFile& grofile, Float3 rotation) {
	const Float3 center = GeometricCenter(grofile);

	std::function<void(Float3&)> position_transform = [&](Float3& pos) {
		pos -= center;
		pos = Float3::rodriguesRotatation(pos, Float3(0, 0, 1), rotation.z);
		pos = Float3::rodriguesRotatation(pos, Float3(0, 1, 0), rotation.y);
		pos = Float3::rodriguesRotatation(pos, Float3(1, 0, 0), rotation.x);
		pos += center;
	};

	std::for_each(grofile.atoms.begin(), grofile.atoms.end(), [&](auto& atom) { position_transform(atom.position); });
}
