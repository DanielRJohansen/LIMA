#include "MoleculeGraph.h"

#include "queue"
#include "unordered_set"
#include <stack>
#include <algorithm>

using namespace LimaMoleculeGraph;
using std::string;
using std::vector;

void MoleculeGraph::Node::addNeighbor(Node* neighbor) {
	if (!neighbor->isHydrogen())
		n_nonhydrogen_neighbors++;
	neighbors.push_back(neighbor);
}

void MoleculeGraph::connectNodes(int left_id, int right_id) {
	nodes[left_id].addNeighbor(&nodes[right_id]);
	nodes[right_id].addNeighbor(&nodes[left_id]);
}

MoleculeGraph LimaMoleculeGraph::createGraph(const TopologyFile& topolfile) {
	MoleculeGraph graph;

	for (const auto& atom : topolfile.GetAllAtoms()) {
		graph.addNode(atom.id, atom.atomname);
	}

	for (const auto& bond : topolfile.GetAllSinglebonds()) {
		graph.connectNodes(bond.ids[0], bond.ids[1]);
	}

	return graph;
}

int GetNumNonvisitedNonHydrogenNeighbors(const MoleculeGraph::Node* node, const std::unordered_set<int>& visitedNodes) {
	int nNonvisitedNonHydrogenNeighbors = 0;
	for (const MoleculeGraph::Node* neighbor : node->getNeighbors()) {
		if (!visitedNodes.contains(neighbor->atomid) && !neighbor->isHydrogen())
			nNonvisitedNonHydrogenNeighbors++;
	}
	return nNonvisitedNonHydrogenNeighbors;
}
std::vector<const MoleculeGraph::Node*> GetNonvisitedNonHydrogenNeighbors(const MoleculeGraph::Node* node, const std::unordered_set<int>& visitedNodes) {
	std::vector<const MoleculeGraph::Node*> neighbors;
	for (const MoleculeGraph::Node* neighbor : node->getNeighbors()) {
		if (!visitedNodes.contains(neighbor->atomid) && !neighbor->isHydrogen())
			neighbors.push_back(neighbor);
	}
	return neighbors;
}


std::pair<const MoleculeGraph::Node*, int> FindFurthestNode(
	const MoleculeGraph& molgraph,
	const MoleculeGraph::Node* start, 
	std::unordered_set<int> visitedNodes = {} // by default an empty set, but we can pass in a set of visited nodes if wanted
) {
	// Queue for BFS, storing pairs of node pointer and current depth
	std::queue<std::pair<const MoleculeGraph::Node*, int>> q;
	
	q.push({ start, 0 });
	visitedNodes.insert(start->atomid);

	// Track the furthest node and its distance
	const MoleculeGraph::Node* furthestNode = start;
	int maxDistance = 0;

	while (!q.empty()) {
		auto [currentNode, currentDistance] = q.front();
		q.pop();

		// Update the furthest node if the current distance is greater
		if (currentDistance > maxDistance) {
			furthestNode = currentNode;
			maxDistance = currentDistance;
		}

		// Get all non-visited, non-hydrogen neighbors
		auto neighbors = GetNonvisitedNonHydrogenNeighbors(currentNode, visitedNodes);

		// Add them to the queue and mark as visited
		for (const auto& neighbor : neighbors) {
			if (visitedNodes.find(neighbor->atomid) == visitedNodes.end()) {
				visitedNodes.insert(neighbor->atomid);
				q.push({ neighbor, currentDistance + 1 });
			}
		}
	}

	return { furthestNode, maxDistance };
}

const MoleculeGraph::Node* FindRootNodeInGraph(const MoleculeGraph& molgraph) {
	const MoleculeGraph::Node* const initialGuessForRootnode = &molgraph.nodes.find(molgraph.highestNodeId)->second;

	auto [secondGuessForRootnode, depth1] = FindFurthestNode(molgraph, initialGuessForRootnode);

	auto [thirdGuessForRootnode, depth2] = FindFurthestNode(molgraph, secondGuessForRootnode);
	
	assert(depth2 >= depth1-1); // Minus 1 because the initial guess might be a hydrogen, in which case we we move 1 step extra.

	return thirdGuessForRootnode;
}


void AddNodeToMapping(const MoleculeGraph::Node* node, std::vector<int>& mapping, int& next_new_id) {
	if (node->atomid >= mapping.size())
		throw std::runtime_error("Node id is too high, shoudnt happen as we keep track of highest id");

	if (mapping[node->atomid] != -1)
		int a = 90;
	mapping[node->atomid] = next_new_id++;
	for (const MoleculeGraph::Node* neighbor : node->getNeighbors()) {
		if (neighbor->isHydrogen()) {
			assert(mapping[neighbor->atomid] == -1);
			mapping[neighbor->atomid] = next_new_id++;
		}
	}
}

void SortSidechainsByDepth(const MoleculeGraph& molgraph, const std::unordered_set<int>& visitedNodes, std::vector<const MoleculeGraph::Node*>& neighbors) {
	// Measure the depths of each sidechain
	std::vector<int> sidechainDepths(neighbors.size());
	for (int i = 0; i < neighbors.size(); i++) {
		sidechainDepths[i] = FindFurthestNode(molgraph, neighbors[i], visitedNodes).second;
	}

	// Now sort the neighbors wrt to the depths, so we take the shallowest one first
	// Pair neighbors with their respective depths
	std::vector<std::pair<const MoleculeGraph::Node*, int>> neighborsWithDepths;
	for (int i = 0; i < neighbors.size(); i++) {
		neighborsWithDepths.emplace_back(neighbors[i], sidechainDepths[i]);
	}

	// Sort based on the second element of the pair (depths)
	std::sort(neighborsWithDepths.begin(), neighborsWithDepths.end(),
		[](const std::pair<const MoleculeGraph::Node*, int>& a, const std::pair<const MoleculeGraph::Node*, int>& b) {
			return a.second < b.second; // Ascending order of depth
		});

	// Extract sorted neighbors
	neighbors.clear();
	for (const auto& pair : neighborsWithDepths) {
		neighbors.push_back(pair.first);
	}
}

void AddSidechainToMapping(const MoleculeGraph& molgraph, std::unordered_set<int>& visitedNodes, std::vector<int>& mapping, int& next_new_id, const MoleculeGraph::Node* currentNode) {
	if (visitedNodes.contains(currentNode->atomid))
		return; // If we have loops, we might encounter a start of a chain already selected (i think)

	while (true) {
		// First add current node
		AddNodeToMapping(currentNode, mapping, next_new_id);
		visitedNodes.insert(currentNode->atomid);
	
		// Now find out if we need to go to a new node/nodes
		const int nNonvisitedNonhydrogenNeighbors = GetNumNonvisitedNonHydrogenNeighbors(currentNode, visitedNodes);

		if (nNonvisitedNonhydrogenNeighbors == 0)
			return;
		else if (nNonvisitedNonhydrogenNeighbors == 1) {
			// Only 1 nonvisited non-hydrogen neighbor, go to that one
			currentNode = GetNonvisitedNonHydrogenNeighbors(currentNode, visitedNodes)[0];
			continue;
		}
		else {
			// We need to search multiple sidechains
			std::vector<const MoleculeGraph::Node*> neighbors = GetNonvisitedNonHydrogenNeighbors(currentNode, visitedNodes);
			SortSidechainsByDepth(molgraph, visitedNodes, neighbors);

			for (const auto neighbor : neighbors) {
				AddSidechainToMapping(molgraph, visitedNodes, mapping, next_new_id, neighbor);
			}
			return;
		}
	}
}

std::vector<int> MakeParticleReorderMapping(const MoleculeGraph& molgraph) {
	const MoleculeGraph::Node* currentNode = FindRootNodeInGraph(molgraph);

	//std::vector<bool> visitedNodes(molgraph.highestNodeId+1, false);
	std::unordered_set<int> visitedNodes;
	std::vector<int> mapping(molgraph.highestNodeId+1, -1);
	int next_new_id = 0;

	AddSidechainToMapping(molgraph, visitedNodes, mapping, next_new_id, currentNode);

	return mapping;
}






template<typename T>
void overwriteParticleIds(std::vector<T>& bonds, const std::vector<int>& map) {
	for (auto& bond : bonds) {
		for (int i = 0; i < bond.n; i++) {
			bond.ids[i] = map[bond.ids[i]];
		}
	}	
}


void LimaMoleculeGraph::reorderoleculeParticlesAccoringingToSubchains(GroFile& grofile, TopologyFile& topfile) {

	std::vector<std::array<std::string, 2>> sbAtomtypesCheck;
	for (const auto& sb : topfile.GetLocalSinglebonds()) {
		sbAtomtypesCheck.push_back({ topfile.GetLocalAtoms()[sb.ids[0]].atomname, topfile.GetLocalAtoms()[sb.ids[1]].atomname });
	}

	const MoleculeGraph molgraph = createGraph(topfile);

	const std::vector<int> map = MakeParticleReorderMapping(molgraph);

	// Overwrite all references to gro_ids in the files	
	for (auto& atom : grofile.atoms) {
		if (map[atom.gro_id-1] < 0)
			throw std::runtime_error("Invalid gro_id in map");
		atom.gro_id = map[atom.gro_id-1] + 1;  // +1 to convert from lima id to groID
	}
	// We can only do this part if the topol does NOT have includes
	if (topfile.GetLocalMolecules().size() > 0)
		throw std::runtime_error("Cannot reorder topol with includes");
	for (auto& atom : topfile.GetLocalAtoms()) {
		atom.id = map[atom.id];
	}

	overwriteParticleIds<>(topfile.GetLocalSinglebonds(), map);
	overwriteParticleIds<>(topfile.GetLocalPairs(), map);
	overwriteParticleIds<>(topfile.GetLocalAnglebonds(), map);
	overwriteParticleIds<>(topfile.GetLocalDihedralbonds(), map);
	overwriteParticleIds<>(topfile.GetLocalImproperDihedralbonds(), map);





	// Re-sort all entries with the new groids
	std::sort(grofile.atoms.begin(), grofile.atoms.end(), [](const GroRecord& a, const GroRecord& b) {return a.gro_id < b.gro_id; });

	std::sort(topfile.GetLocalAtoms().begin(), topfile.GetLocalAtoms().end(), [](const auto& a, const auto& b) {return a.id < b.id; });

	// Check that we didn't mess up the singlebonds. If all these are identical, the other bonds must be aswell
	for (int i = 0; i < topfile.GetLocalSinglebonds().size(); i++) {
		if (
			sbAtomtypesCheck[i][0] != topfile.GetLocalAtoms()[topfile.GetLocalSinglebonds()[i].ids[0]].atomname ||
			sbAtomtypesCheck[i][1] != topfile.GetLocalAtoms()[topfile.GetLocalSinglebonds()[i].ids[1]].atomname
			)
		{
			auto a = topfile.GetLocalAtoms()[topfile.GetLocalSinglebonds()[i].ids[0]].atomname;
			auto b = topfile.GetLocalAtoms()[topfile.GetLocalSinglebonds()[i].ids[1]].atomname;
			throw std::runtime_error("Reordering of particles messed up singlebonds" + std::to_string(i));
		}
	}

	std::sort(topfile.GetLocalSinglebonds().begin(), topfile.GetLocalSinglebonds().end(), [](const auto& a, const auto& b) { return a.ids[0] < b.ids[0]; });
	std::sort(topfile.GetLocalPairs().begin(), topfile.GetLocalPairs().end(), [](const auto& a, const auto& b) { return a.ids[0] < b.ids[0]; });
	std::sort(topfile.GetLocalAnglebonds().begin(), topfile.GetLocalAnglebonds().end(), [](const auto& a, const auto& b) { return a.ids[0] < b.ids[0]; });
	std::sort(topfile.GetLocalDihedralbonds().begin(), topfile.GetLocalDihedralbonds().end(), [](const auto& a, const auto& b) { return a.ids[0] < b.ids[0]; });
	std::sort(topfile.GetLocalImproperDihedralbonds().begin(), topfile.GetLocalImproperDihedralbonds().end(), [](const auto& a, const auto& b) { return a.ids[0] < b.ids[0]; });


	for (int i = 0; i < topfile.GetLocalAtoms().size(); i++) {
		if (topfile.GetLocalAtoms()[i].id != i)
			throw std::runtime_error("Reordering of particles messed up atoms" + std::to_string(i));
	}



}