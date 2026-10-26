#include "BoxBuilder.cuh"
#include "RenderUtilities.cuh"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <queue>
#include <set>

namespace {

// Longest heavy-atom bond in a lipid is P-O, ~0.16 nm. The shortest non-bonded heavy-atom distance in a residue is ~0.24 nm
constexpr float maxBondLength = 0.2f;
// Choline methyls and protein side chains have short carbon runs, acyl and sterol chains have long ones
constexpr int minTailCarbonRun = 6;
// A tail that shares all but a few atoms with a longer one, is a methyl branch (eg. on a sterol ring) and not a chain
constexpr std::size_t minTailOwnAtoms = 3;
constexpr std::size_t minHeavyAtoms = 10;
// Bonds are found pairwise, so skip residues far larger than any lipid (eg. a whole molecule in one residue)
constexpr std::size_t maxHeavyAtoms = 200;

float DistanceMinimumImage(const Float3& a, const Float3& b, const Float3& box)
{
	const std::array<float, 3> delta{ b.x - a.x, b.y - a.y, b.z - a.z };
	const std::array<float, 3> lengths{ box.x, box.y, box.z };
	float squared = 0.f;
	for (int axis = 0; axis < 3; ++axis) {
		float d = delta[axis];
		if (lengths[axis] > 0.f)
			d -= std::round(d / lengths[axis]) * lengths[axis];
		squared += d * d;
	}
	return std::sqrt(squared);
}

struct HeavyAtomGraph {
	std::vector<int> atomIds;	// Index into the grofile
	std::vector<char> letters;
	std::vector<std::vector<int>> neighbors;	// Local indices

	bool IsTerminal(int local) const { return neighbors[local].size() == 1; }
};

HeavyAtomGraph BuildGraph(const GroFile& grofile, std::size_t begin, std::size_t end)
{
	HeavyAtomGraph graph;
	for (std::size_t atomId = begin; atomId < end; ++atomId) {
		const GroRecord& atom = grofile.atoms[atomId];
		const char letter = RenderUtilities::AtomLetter(atom.atomName.View(), atom.residueName.View());
		if (letter == 'H')
			continue;
		graph.atomIds.push_back(static_cast<int>(atomId));
		graph.letters.push_back(letter);
	}
	if (graph.atomIds.size() < minHeavyAtoms || graph.atomIds.size() > maxHeavyAtoms)
		return graph;

	graph.neighbors.resize(graph.atomIds.size());
	for (std::size_t i = 0; i < graph.atomIds.size(); ++i) {
		for (std::size_t j = i + 1; j < graph.atomIds.size(); ++j) {
			const float distance = DistanceMinimumImage(grofile.atoms[graph.atomIds[i]].position,
				grofile.atoms[graph.atomIds[j]].position, grofile.box_size);
			if (distance < maxBondLength) {
				graph.neighbors[i].push_back(static_cast<int>(j));
				graph.neighbors[j].push_back(static_cast<int>(i));
			}
		}
	}
	return graph;
}

LipidHeadgroup ClassifyHeadgroup(const HeavyAtomGraph& graph, bool hasPhosphate)
{
	if (!hasPhosphate)
		return LipidHeadgroup::Sterol;

	bool hasTerminalAmine = false;
	for (std::size_t i = 0; i < graph.letters.size(); ++i) {
		if (graph.letters[i] != 'N')
			continue;
		const auto terminalCarbons = std::ranges::count_if(graph.neighbors[i], [&graph](int n) {
			return graph.letters[n] == 'C' && graph.IsTerminal(n);
		});
		if (terminalCarbons >= 3)
			return LipidHeadgroup::Choline;
		hasTerminalAmine |= graph.IsTerminal(static_cast<int>(i));
	}
	if (!hasTerminalAmine)
		return LipidHeadgroup::Anionic;

	// Serine carries a carboxyl group, a carbon with two terminal oxygens
	for (std::size_t i = 0; i < graph.letters.size(); ++i) {
		if (graph.letters[i] != 'C')
			continue;
		const auto terminalOxygens = std::ranges::count_if(graph.neighbors[i], [&graph](int n) {
			return graph.letters[n] == 'O' && graph.IsTerminal(n);
		});
		if (terminalOxygens >= 2)
			return LipidHeadgroup::Serine;
	}
	return LipidHeadgroup::Ethanolamine;
}

std::optional<LipidModel> InterpretResidue(const GroFile& grofile, std::size_t begin, std::size_t end)
{
	const HeavyAtomGraph graph = BuildGraph(grofile, begin, end);
	if (graph.neighbors.empty())
		return std::nullopt;

	auto head = std::ranges::find(graph.letters, 'P');
	const bool hasPhosphate = head != graph.letters.end();
	if (!hasPhosphate)
		head = std::ranges::find(graph.letters, 'O');
	if (head == graph.letters.end())
		return std::nullopt;
	const int headLocal = static_cast<int>(std::distance(graph.letters.begin(), head));

	// Breadth first from the head, so the path to every atom is a shortest path, and paths share their prefixes
	std::vector<int> parent(graph.atomIds.size(), -1);
	std::vector<bool> visited(graph.atomIds.size(), false);
	std::queue<int> queue;
	queue.push(headLocal);
	visited[headLocal] = true;
	while (!queue.empty()) {
		const int current = queue.front();
		queue.pop();
		for (int neighbor : graph.neighbors[current]) {
			if (visited[neighbor])
				continue;
			visited[neighbor] = true;
			parent[neighbor] = current;
			queue.push(neighbor);
		}
	}

	std::vector<std::vector<int>> candidates;
	for (std::size_t i = 0; i < graph.atomIds.size(); ++i) {
		if (graph.letters[i] != 'C' || !graph.IsTerminal(static_cast<int>(i)) || !visited[i])
			continue;
		std::vector<int> path;
		int carbonRun = 0;
		bool inCarbonRun = true;
		for (int current = static_cast<int>(i); current != -1; current = parent[current]) {
			path.push_back(current);
			inCarbonRun &= graph.letters[current] == 'C';
			carbonRun += inCarbonRun ? 1 : 0;
		}
		if (carbonRun >= minTailCarbonRun) {
			std::ranges::reverse(path);
			candidates.push_back(std::move(path));
		}
	}
	if (candidates.empty())
		return std::nullopt;

	// Longest first, so a methyl branch near the end of a chain is dropped rather than the chain itself
	std::ranges::stable_sort(candidates, [](const auto& a, const auto& b) { return a.size() > b.size(); });

	LipidModel lipid;
	lipid.headParticleId = graph.atomIds[headLocal];
	lipid.firstParticleId = static_cast<int>(begin);
	lipid.particleCount = static_cast<int>(end - begin);
	lipid.headgroup = ClassifyHeadgroup(graph, hasPhosphate);

	std::set<int> covered;
	for (const auto& path : candidates) {
		std::size_t branch = 0;
		for (std::size_t i = 0; i < path.size(); ++i)
			if (covered.contains(path[i]))
				branch = i;
		if (path.size() - branch - 1 < minTailOwnAtoms)
			continue;

		std::vector<int> tail;
		for (std::size_t i = branch; i < path.size(); ++i) {
			tail.push_back(graph.atomIds[path[i]]);
			covered.insert(path[i]);
		}
		lipid.tails.push_back(std::move(tail));
	}
	return lipid;
}

} // namespace

LipidModels BoxBuilder::InterpretLipids(const GroFile& grofile)
{
	LipidModels lipids;
	std::size_t begin = 0;
	while (begin < grofile.atoms.size()) {
		std::size_t end = begin + 1;
		while (end < grofile.atoms.size()
			&& grofile.atoms[end].residue_number == grofile.atoms[begin].residue_number
			&& grofile.atoms[end].residueName == grofile.atoms[begin].residueName)
			++end;

		if (std::optional<LipidModel> lipid = InterpretResidue(grofile, begin, end))
			lipids.push_back(std::move(*lipid));
		begin = end;
	}
	return lipids;
}

ParticleCategories BoxBuilder::InterpretParticleCategories(const GroFile& grofile, const LipidModels& lipids)
{
	static const std::set<std::string_view> waterResidues{ "SOL", "SPC", "SPCE", "TIP3", "TIP3P", "TIP4", "TIP4P", "TIP5P", "HOH", "WAT" };

	ParticleCategories categories(grofile.atoms.size(), ParticleCategory::Other);
	std::size_t begin = 0;
	while (begin < grofile.atoms.size()) {
		std::size_t end = begin + 1;
		while (end < grofile.atoms.size()
			&& grofile.atoms[end].residue_number == grofile.atoms[begin].residue_number
			&& grofile.atoms[end].residueName == grofile.atoms[begin].residueName)
			++end;

		const GroRecord& first = grofile.atoms[begin];
		ParticleCategory category = ParticleCategory::Other;
		if (waterResidues.contains(first.residueName.View()))
			category = ParticleCategory::Water;
		else if (end - begin == 1 && RenderUtilities::AtomLetter(first.atomName.View(), first.residueName.View()) == 'I')
			category = ParticleCategory::Ion;
		else if (IsProteinResidue(first.residueName.View()))
			category = ParticleCategory::Protein;
		std::fill(categories.begin() + begin, categories.begin() + end, category);
		begin = end;
	}

	for (const LipidModel& lipid : lipids)
		std::fill_n(categories.begin() + lipid.firstParticleId, lipid.particleCount, ParticleCategory::Lipid);
	return categories;
}
