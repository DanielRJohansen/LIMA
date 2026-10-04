#include "BoxBuilder.cuh"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <string_view>
#include <unordered_set>
#include <unordered_map>
#include <cstdint>

#include <glm/glm.hpp>

namespace {

constexpr float pi = 3.14159265358979323846f;

const std::unordered_set<std::string_view> proteinResidues = {
	"ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYM", "CYX", "GLN", "GLU", "GLH",
	"GLY", "HIS", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE",
	"PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE", "ACE", "NME"
};

struct InterpretedResidue {
	static constexpr std::size_t noAtom = std::numeric_limits<std::size_t>::max();

	std::size_t nAtom = noAtom;
	std::size_t caAtom = noAtom;
	std::size_t cAtom = noAtom;
	std::size_t oAtom = noAtom;
	std::size_t hAtom = noAtom;
	glm::vec3 n{};
	glm::vec3 ca{};
	glm::vec3 c{};
	glm::vec3 o{};
	glm::vec3 h{};
	SecondaryStructure secondaryStructure = SecondaryStructure::Coil;
};

using InterpretedChain = std::vector<InterpretedResidue>;
using InterpretedChains = std::vector<InterpretedChain>;

glm::vec3 ToGlmVec3(const Float3& value)
{
	return { value.x, value.y, value.z };
}

glm::vec3 MinimumImage(glm::vec3 delta, const Float3& box)
{
	const std::array<float, 3> lengths{ box.x, box.y, box.z };
	for (int axis = 0; axis < 3; ++axis) {
		if (lengths[axis] > 0.f)
			delta[axis] -= std::round(delta[axis] / lengths[axis]) * lengths[axis];
	}
	return delta;
}

float DihedralDegrees(const glm::vec3& a, const glm::vec3& b, const glm::vec3& c, const glm::vec3& d)
{
	const glm::vec3 b0 = a - b;
	const glm::vec3 b1 = c - b;
	const glm::vec3 b2 = d - c;
	if (glm::dot(b0, b0) < 1e-10f || glm::dot(b1, b1) < 1e-10f || glm::dot(b2, b2) < 1e-10f)
		return 0.f;

	const glm::vec3 axis = glm::normalize(b1);
	const glm::vec3 v = b0 - glm::dot(b0, axis) * axis;
	const glm::vec3 w = b2 - glm::dot(b2, axis) * axis;
	if (glm::dot(v, v) < 1e-10f || glm::dot(w, w) < 1e-10f)
		return 0.f;

	return std::atan2(glm::dot(glm::cross(axis, v), w), glm::dot(v, w)) * 180.f / pi;
}

bool IsHelicalAngle(float phi, float psi)
{
	return phi >= -105.f && phi <= -30.f && psi >= -90.f && psi <= 20.f;
}

bool IsExtendedAngle(float phi, float psi)
{
	return phi >= -180.f && phi <= -65.f
		&& ((psi >= 55.f && psi <= 180.f) || (psi >= -180.f && psi <= -125.f));
}

glm::vec3 AmideHydrogen(const InterpretedChain& chain, std::size_t residueIndex)
{
	const InterpretedResidue& residue = chain[residueIndex];
	if (residue.hAtom != InterpretedResidue::noAtom)
		return residue.h;

	if (residueIndex > 0) {
		const glm::vec3 towardCa = glm::normalize(residue.ca - residue.n);
		const glm::vec3 towardPreviousC = glm::normalize(chain[residueIndex - 1].c - residue.n);
		const glm::vec3 direction = -(towardCa + towardPreviousC);
		if (glm::dot(direction, direction) > 1e-8f)
			return residue.n + glm::normalize(direction) * 0.10f;
	}
	return residue.n;
}

bool IsBackboneHydrogenBond(
	const InterpretedResidue& acceptor,
	const InterpretedChain& donorChain,
	std::size_t donorIndex)
{
	const InterpretedResidue& donor = donorChain[donorIndex];
	if (glm::length(acceptor.o - donor.n) > 0.36f)
		return false;

	const glm::vec3 hydrogen = AmideHydrogen(donorChain, donorIndex);
	const glm::vec3 donorDirection = donor.n - hydrogen;
	const glm::vec3 acceptorDirection = acceptor.o - hydrogen;
	if (glm::dot(donorDirection, donorDirection) < 1e-8f || glm::dot(acceptorDirection, acceptorDirection) < 1e-8f)
		return true;

	const float cosine = glm::clamp(
		glm::dot(glm::normalize(donorDirection), glm::normalize(acceptorDirection)), -1.f, 1.f);
	return std::acos(cosine) * 180.f / pi >= 115.f;
}

void MarkRuns(
	InterpretedChain& chain,
	const std::vector<bool>& candidates,
	std::size_t minimumLength,
	SecondaryStructure assignment,
	bool preserveHelices)
{
	std::size_t begin = 0;
	while (begin < candidates.size()) {
		while (begin < candidates.size() && !candidates[begin])
			++begin;
		std::size_t end = begin;
		while (end < candidates.size() && candidates[end])
			++end;

		if (end - begin >= minimumLength) {
			for (std::size_t i = begin; i < end; ++i) {
				if (!preserveHelices || chain[i].secondaryStructure != SecondaryStructure::Helix)
					chain[i].secondaryStructure = assignment;
			}
		}
		begin = end;
	}
}

std::optional<InterpretedResidue> MakeResidue(
	const GroFile& grofile,
	std::size_t begin,
	std::size_t end)
{
	if (begin == end || !proteinResidues.contains(grofile.atoms[begin].residueName.View()))
		return std::nullopt;

	InterpretedResidue residue{};
	for (std::size_t atomIndex = begin; atomIndex < end; ++atomIndex) {
		const std::string_view atomName = grofile.atoms[atomIndex].atomName.View();
		if (atomName == "N") residue.nAtom = atomIndex;
		else if (atomName == "CA") residue.caAtom = atomIndex;
		else if (atomName == "C") residue.cAtom = atomIndex;
		else if (atomName == "O" || atomName == "O1" || atomName == "OT1") residue.oAtom = atomIndex;
		else if (atomName == "H" || atomName == "HN" || atomName == "H1") residue.hAtom = atomIndex;
	}

	if (residue.nAtom == InterpretedResidue::noAtom || residue.caAtom == InterpretedResidue::noAtom
		|| residue.cAtom == InterpretedResidue::noAtom || residue.oAtom == InterpretedResidue::noAtom)
		return std::nullopt;

	residue.n = ToGlmVec3(grofile.atoms[residue.nAtom].position);
	residue.ca = ToGlmVec3(grofile.atoms[residue.caAtom].position);
	residue.c = ToGlmVec3(grofile.atoms[residue.cAtom].position);
	residue.o = ToGlmVec3(grofile.atoms[residue.oAtom].position);
	if (residue.hAtom != InterpretedResidue::noAtom)
		residue.h = ToGlmVec3(grofile.atoms[residue.hAtom].position);
	return residue;
}

void Translate(InterpretedResidue& residue, const glm::vec3& shift)
{
	residue.n += shift;
	residue.ca += shift;
	residue.c += shift;
	residue.o += shift;
	if (residue.hAtom != InterpretedResidue::noAtom)
		residue.h += shift;
}

void AssignSecondaryStructure(InterpretedChains& structure)
{
	std::vector<std::vector<bool>> extendedCandidates(structure.size());
	for (std::size_t chainIndex = 0; chainIndex < structure.size(); ++chainIndex) {
		InterpretedChain& chain = structure[chainIndex];
		std::vector<bool> helixAngles(chain.size(), false);
		std::vector<bool> sheetAngles(chain.size(), false);
		for (std::size_t i = 1; i + 1 < chain.size(); ++i) {
			const float phi = DihedralDegrees(chain[i - 1].c, chain[i].n, chain[i].ca, chain[i].c);
			const float psi = DihedralDegrees(chain[i].n, chain[i].ca, chain[i].c, chain[i + 1].n);
			helixAngles[i] = IsHelicalAngle(phi, psi);
			sheetAngles[i] = IsExtendedAngle(phi, psi);
		}

		MarkRuns(chain, helixAngles, 4, SecondaryStructure::Helix, false);

		std::vector<bool> alphaBonds(chain.size(), false);
		for (std::size_t i = 0; i + 4 < chain.size(); ++i)
			alphaBonds[i] = IsBackboneHydrogenBond(chain[i], chain, i + 4);
		for (std::size_t i = 0; i + 1 < alphaBonds.size(); ++i) {
			if (alphaBonds[i] && alphaBonds[i + 1]) {
				const std::size_t last = std::min(i + 5, chain.size() - 1);
				for (std::size_t residueIndex = i; residueIndex <= last; ++residueIndex)
					chain[residueIndex].secondaryStructure = SecondaryStructure::Helix;
			}
		}
		extendedCandidates[chainIndex] = std::move(sheetAngles);
	}

	std::vector<std::vector<bool>> sheetContacts;
	sheetContacts.reserve(structure.size());
	for (const InterpretedChain& chain : structure)
		sheetContacts.emplace_back(chain.size(), false);

	// Residues can only be in contact if their CAs are within maxCaDistance, so instead of comparing all pairs,
	// we put the residues in a grid, and only compare with residues in the neighboring cells. The criteria are symmetric
	// and only ever set contacts to true, so checking each pair once gives the same result as checking all ordered pairs
	constexpr float maxCaDistance = 0.75f;
	constexpr float cellSize = 0.8f; // Slightly larger than maxCaDistance, so rounding can never hide a neighbor
	struct ResidueRef { std::size_t chain; std::size_t index; };
	auto CellOf = [](const glm::vec3& position) { return glm::ivec3(glm::floor(position / cellSize)); };
	auto CellKey = [](const glm::ivec3& cell) {
		constexpr uint64_t mask = (1ull << 21) - 1;
		return ((static_cast<uint64_t>(cell.x) & mask) << 42) | ((static_cast<uint64_t>(cell.y) & mask) << 21) | (static_cast<uint64_t>(cell.z) & mask);
		};

	// Helices are never part of a sheet contact
	std::unordered_map<uint64_t, std::vector<ResidueRef>> grid;
	for (std::size_t chain = 0; chain < structure.size(); ++chain)
		for (std::size_t i = 0; i < structure[chain].size(); ++i)
			if (structure[chain][i].secondaryStructure != SecondaryStructure::Helix)
				grid[CellKey(CellOf(structure[chain][i].ca))].push_back({ chain, i });

	for (std::size_t chainA = 0; chainA < structure.size(); ++chainA) {
		for (std::size_t i = 0; i < structure[chainA].size(); ++i) {
			InterpretedResidue& a = structure[chainA][i];
			if (a.secondaryStructure == SecondaryStructure::Helix)
				continue;

			const glm::ivec3 cellA = CellOf(a.ca);
			for (int dx = -1; dx <= 1; dx++) for (int dy = -1; dy <= 1; dy++) for (int dz = -1; dz <= 1; dz++) {
				const auto cell = grid.find(CellKey(cellA + glm::ivec3{ dx, dy, dz }));
				if (cell == grid.end())
					continue;

				for (const auto [chainB, j] : cell->second) {
					// Each unordered pair once
					if (chainB < chainA || (chainB == chainA && j <= i))
						continue;
					if (chainA == chainB && j - i < 3)
						continue;
					InterpretedResidue& b = structure[chainB][j];
					if (glm::length(a.ca - b.ca) > maxCaDistance
						|| (!extendedCandidates[chainA][i] && !extendedCandidates[chainB][j]))
						continue;

					if (IsBackboneHydrogenBond(a, structure[chainB], j)
						|| IsBackboneHydrogenBond(b, structure[chainA], i)) {
						sheetContacts[chainA][i] = true;
						sheetContacts[chainB][j] = true;
					}
				}
			}
		}
	}

	for (std::size_t chainIndex = 0; chainIndex < structure.size(); ++chainIndex) {
		std::vector<bool> sheetResidues(structure[chainIndex].size(), false);
		for (std::size_t i = 0; i < sheetContacts[chainIndex].size(); ++i) {
			if (!sheetContacts[chainIndex][i])
				continue;
			const std::size_t begin = i > 0 ? i - 1 : i;
			const std::size_t end = std::min(i + 1, sheetResidues.size() - 1);
			for (std::size_t neighbor = begin; neighbor <= end; ++neighbor)
				sheetResidues[neighbor] = extendedCandidates[chainIndex][neighbor];
		}
		MarkRuns(structure[chainIndex], sheetResidues, 2, SecondaryStructure::Sheet, true);
	}
}

BackboneChains ToBackboneChains(const InterpretedChains& interpretedChains)
{
	BackboneChains result;
	result.reserve(interpretedChains.size());
	for (const InterpretedChain& interpretedChain : interpretedChains) {
		BackboneChain chain;
		chain.points.reserve(interpretedChain.size());
		for (const InterpretedResidue& residue : interpretedChain) {
			chain.points.push_back({ static_cast<int>(residue.caAtom), residue.secondaryStructure });
		}
		if (chain.points.size() >= 2)
			result.push_back(std::move(chain));
	}
	return result;
}

} // namespace

BackboneChains BoxBuilder::InterpretBackboneChains(const GroFile& grofile)
{
	InterpretedChains interpretedChains;
	InterpretedChain currentChain;

	std::size_t begin = 0;
	while (begin < grofile.atoms.size()) {
		std::size_t end = begin + 1;
		while (end < grofile.atoms.size()
			&& grofile.atoms[end].residue_number == grofile.atoms[begin].residue_number
			&& grofile.atoms[end].residueName == grofile.atoms[begin].residueName)
			++end;

		std::optional<InterpretedResidue> residue = MakeResidue(grofile, begin, end);
		if (!residue.has_value()) {
			if (currentChain.size() >= 2)
				interpretedChains.push_back(std::move(currentChain));
			currentChain.clear();
			begin = end;
			continue;
		}

		if (!currentChain.empty()) {
			const glm::vec3 rawPeptideBond = residue->n - currentChain.back().c;
			const glm::vec3 peptideBond = MinimumImage(rawPeptideBond, grofile.box_size);
			if (glm::length(peptideBond) > 0.22f) {
				if (currentChain.size() >= 2)
					interpretedChains.push_back(std::move(currentChain));
				currentChain.clear();
			}
			else {
				Translate(*residue, peptideBond - rawPeptideBond);
			}
		}

		currentChain.push_back(*residue);
		begin = end;
	}

	if (currentChain.size() >= 2)
		interpretedChains.push_back(std::move(currentChain));

	AssignSecondaryStructure(interpretedChains);
	return ToBackboneChains(interpretedChains);
}
