#include "MoleculeGraph.h"


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

int GetNumNonvisitedNonHydrogenNeighbors(const MoleculeGraph::Node* node, const std::vector<bool>& visitedNodes) {
	int nNonvisitedNonHydrogenNeighbors = 0;
	for (const MoleculeGraph::Node* neighbor : node->getNeighbors()) {
		if (!visitedNodes[neighbor->atomid] && !neighbor->isHydrogen())
			nNonvisitedNonHydrogenNeighbors++;
	}
	return nNonvisitedNonHydrogenNeighbors;
}

std::vector<const MoleculeGraph::Node*> GetNonvisitedNonHydrogenNeighbors(const MoleculeGraph::Node* node, const std::vector<bool>& visitedNodes) {
	std::vector<const MoleculeGraph::Node*> neighbors;
	for (const MoleculeGraph::Node* neighbor : node->getNeighbors()) {
		if (!visitedNodes[neighbor->atomid] && !neighbor->isHydrogen())
			neighbors.push_back(neighbor);
	}
	return neighbors;
}

const std::pair<const MoleculeGraph::Node*, int> FindFurthestNode(const MoleculeGraph& molgraph, std::vector<bool> visitedNodes, const MoleculeGraph::Node* currentNode, int currentDistance) {
	while (true) {
		int nNonvisitedNonHydrogenNeighbors = GetNumNonvisitedNonHydrogenNeighbors(currentNode, visitedNodes);

		if (nNonvisitedNonHydrogenNeighbors == 0)
			return std::make_pair(currentNode, currentDistance);
		if (nNonvisitedNonHydrogenNeighbors == 1) {
			const auto neighbor = GetNonvisitedNonHydrogenNeighbors(currentNode, visitedNodes)[0];
			visitedNodes[neighbor->atomid] = true;
			currentNode = neighbor;
			currentDistance++;
		}
		else {
			// We need to search multiple sidechains
			std::vector<std::pair<const MoleculeGraph::Node*, int>> furthestNodes;

			for (const auto neighbor : GetNonvisitedNonHydrogenNeighbors(currentNode, visitedNodes)) {
				visitedNodes[neighbor->atomid] = true;
				currentDistance++;
				furthestNodes.push_back(FindFurthestNode(molgraph, visitedNodes, neighbor, currentDistance));
			}


			const int furthestNodeIndex = std::max_element(furthestNodes.begin(), furthestNodes.end(), [](const auto& a, const auto& b) {return a.second < b.second; }) - furthestNodes.begin();
			return furthestNodes[furthestNodeIndex];
		}
	}
}

const MoleculeGraph::Node* FindRootNodeInGraph(const MoleculeGraph& molgraph) {
	const MoleculeGraph::Node* const initialGuessForRootnode = &molgraph.nodes.find(molgraph.highestNodeId)->second;

	std::vector<bool> visitedNodes1(molgraph.highestNodeId+2, false);
	visitedNodes1[initialGuessForRootnode->atomid] = true;
	auto [secondGuessForRootnode, depth1] = FindFurthestNode(molgraph, visitedNodes1, initialGuessForRootnode, 0);

	std::vector<bool> visitedNodes2(molgraph.highestNodeId+2, false);
	visitedNodes2[secondGuessForRootnode->atomid] = true;
	auto [thirdGuessForRootnode, depth2] = FindFurthestNode(molgraph, visitedNodes2, secondGuessForRootnode, 0);

	assert(depth2 >= depth1);

	return thirdGuessForRootnode;
}


void AddNodeToMapping(const MoleculeGraph::Node* node, std::vector<int>& mapping, int& next_new_id) {
	if (node->atomid >= mapping.size())
		throw std::runtime_error("Node id is too high, shoudnt happen as we keep track of highest id");

	mapping[node->atomid] = next_new_id++;
	for (const MoleculeGraph::Node* neighbor : node->getNeighbors()) {
		if (neighbor->isHydrogen()) {
			mapping[neighbor->atomid] = next_new_id++;
		}
	}
}

void SortSidechainsByDepth(const MoleculeGraph& molgraph, const std::vector<bool>& visitedNodes, std::vector<const MoleculeGraph::Node*>& neighbors) {
	// Measure the depths of each sidechain
	std::vector<int> sidechainDepths(neighbors.size());
	for (int i = 0; i < neighbors.size(); i++) {
		sidechainDepths[i] = FindFurthestNode(molgraph, visitedNodes, neighbors[i], 0).second;
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

void AddSidechainToMapping(const MoleculeGraph& molgraph, std::vector<bool>& visitedNodes, std::vector<int>& mapping, int& next_new_id, const MoleculeGraph::Node* currentNode) {
	if (visitedNodes[currentNode->atomid])
		return; // If we have loops, we might encounter a start of a chain already selected (i think)

	while (true) {
		// First add current node
		visitedNodes[currentNode->atomid] = true;
		AddNodeToMapping(currentNode, mapping, next_new_id);

		// Now find out if we need to go to a new node/nodes
		const int nNonvisitedNonhydrogenNeighbors = GetNumNonvisitedNonHydrogenNeighbors(currentNode, visitedNodes);

		if (nNonvisitedNonhydrogenNeighbors == 0)
			break;
		if (nNonvisitedNonhydrogenNeighbors == 1) {
			// Only 1 nonvisited non-hydrogen neighbor, go to that one
			currentNode = GetNonvisitedNonHydrogenNeighbors(currentNode, visitedNodes)[0];
		}
		else {
			// We need to search multiple sidechains
			std::vector<const MoleculeGraph::Node*> neighbors = GetNonvisitedNonHydrogenNeighbors(currentNode, visitedNodes);
			SortSidechainsByDepth(molgraph, visitedNodes, neighbors);

			for (const auto neighbor : neighbors) {
				AddSidechainToMapping(molgraph, visitedNodes, mapping, next_new_id, neighbor);
			}
		}
	}
}

std::vector<int> MakeParticleReorderMapping(const MoleculeGraph& molgraph) {
	const MoleculeGraph::Node* currentNode = FindRootNodeInGraph(molgraph);

	std::vector<bool> visitedNodes(molgraph.highestNodeId+2, false);

	std::vector<int> mapping(molgraph.highestNodeId+2, -1);
	int next_new_id = 1;

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
	const MoleculeGraph molgraph = createGraph(topfile);

	const std::vector<int> map = MakeParticleReorderMapping(molgraph);

	// Overwrite all references to gro_ids in the files	
	for (auto& atom : grofile.atoms) {
		if (map[atom.gro_id-1] == -1)
			throw std::runtime_error("Invalid gro_id in map");
		atom.gro_id = map[atom.gro_id-1];
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
	std::sort(topfile.GetLocalSinglebonds().begin(), topfile.GetLocalSinglebonds().end(), [](const auto& a, const auto& b) { return a.ids[0] < b.ids[0]; });
	std::sort(topfile.GetLocalPairs().begin(), topfile.GetLocalPairs().end(), [](const auto& a, const auto& b) { return a.ids[0] < b.ids[0]; });
	std::sort(topfile.GetLocalAnglebonds().begin(), topfile.GetLocalAnglebonds().end(), [](const auto& a, const auto& b) { return a.ids[0] < b.ids[0]; });
	std::sort(topfile.GetLocalDihedralbonds().begin(), topfile.GetLocalDihedralbonds().end(), [](const auto& a, const auto& b) { return a.ids[0] < b.ids[0]; });
	std::sort(topfile.GetLocalImproperDihedralbonds().begin(), topfile.GetLocalImproperDihedralbonds().end(), [](const auto& a, const auto& b) { return a.ids[0] < b.ids[0]; });
}