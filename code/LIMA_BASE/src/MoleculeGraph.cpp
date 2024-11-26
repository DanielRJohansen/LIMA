#include "MoleculeGraph.h"

#include "queue"
#include "unordered_set"
#include <stack>
#include <algorithm>

#include <unordered_map>

using namespace LimaMoleculeGraph;
using std::string;
using std::vector;



MoleculeTree::MoleculeTree(int rootId) {
	tree.insert({ rootId, {} });
}
std::vector<int> MoleculeTree::GetAllChildIdsAndSelf(int parentId, std::unordered_map<int, int> nodeIdsNumDownstreamNodes) const {
	std::vector<int> ids;
	ids.reserve(nodeIdsNumDownstreamNodes.at(parentId));

	std::stack<int> nodeStack;
	nodeStack.push(parentId);

	while (!nodeStack.empty()) {
		int currentId = nodeStack.top();
		nodeStack.pop();

		ids.push_back(currentId);

		for (const int childId : tree.at(currentId)) {
			nodeStack.push(childId);
		}
	}

	if (ids.size() != nodeIdsNumDownstreamNodes.at(parentId))
		throw std::runtime_error("Number of nodes in tree does not match number of nodes in numDownstreamNodes");

	return ids;
}


void MoleculeGraph::Node::addNeighbor(Node* neighbor) {
	if (!neighbor->isHydrogen())
		n_nonhydrogen_neighbors++;
	neighbors.push_back(neighbor);
}

void MoleculeGraph::connectNodes(int left_id, int right_id) {
	nodes.at(left_id).addNeighbor(&nodes.at(right_id));
	nodes.at(right_id).addNeighbor(&nodes.at(left_id));
}

MoleculeGraph::MoleculeGraph(const TopologyFile::Moleculetype& molecule, std::optional<const std::unordered_set<int>> allowedIds) {
	/*if (molecule.atoms.front().id != 0 || molecule.atoms.back() != molecule.atoms.size() - 1)
		throw std::runtime_error("Atoms must be ordered from 0 to n-1");*/

	//nodes.resize(molecule.atoms.size());

	std::vector<std::pair<int, Node>> temp_nodes;
	temp_nodes.reserve(molecule.atoms.size());
	for (const auto& atom : molecule.atoms) {
		if (allowedIds.has_value() && !allowedIds.value().contains(atom.id))
			continue;

		temp_nodes.emplace_back(atom.id, Node(atom.id, atom.atomname));
	}
	nodes.insert(temp_nodes.begin(), temp_nodes.end());
	
	root = &nodes.at(temp_nodes[0].second.atomid);

	for (const auto& bond : molecule.singlebonds) {
		connectNodes(bond.ids[0], bond.ids[1]);
	}
}

MoleculeGraph::MoleculeGraph(const std::vector<std::pair<int, std::string>>& atoms, const std::vector<std::array<int, 2>>& edges) {
	std::vector<std::pair<int, Node>> temp_nodes;
	temp_nodes.reserve(atoms.size());
	for (const auto& atom : atoms) {
		temp_nodes.emplace_back(atom.first, Node(atom.first, atom.second));
	}

	nodes.insert(temp_nodes.begin(), temp_nodes.end());
	root = &nodes.at(0);
	
	for (const auto& edge : edges) {
		connectNodes(edge[0], edge[1]);
	}
}


MoleculeTree MoleculeGraph::ConstructMoleculeTree() const {
	std::unordered_set<int> visited;
	MoleculeTree moleculeTree(root->atomid);

	std::queue<const Node*> nodeQueue;
	nodeQueue.push(root);
	visited.insert(root->atomid);
	

	while (!nodeQueue.empty()) {
		const Node* current = nodeQueue.front();
		nodeQueue.pop();

		for (const Node* neighbor : current->getNeighbors()) {
			const Node* neighbor_ptr = neighbor; // Assign Node* to NodePtr
			if (visited.insert(neighbor->atomid).second) {
				nodeQueue.push(neighbor_ptr);
				moleculeTree.AddChild(current->atomid, neighbor->atomid);
			}
		}
	}

	return moleculeTree;
}


std::unordered_map<int, int> MoleculeGraph::ComputeNumDownstreamNodes() const 
{
	MoleculeTree moleculeTree = ConstructMoleculeTree();

	std::stack<const Node*> processStack; // The order in which to actually process NumDownstreamNodes
	for (const Node& node : BFS(root->atomid)) {
		processStack.push(&node);
	}

	std::unordered_map<int, int> nodeIdToNumDownstreamnodes;

	while (!processStack.empty()) {
		const Node& current = *processStack.top();
		processStack.pop();
		int numDownstreamnodes = 1;
		for (const int& neighborId : moleculeTree.GetChildIds(current.atomid)) {
			numDownstreamnodes += nodeIdToNumDownstreamnodes.at(neighborId);
		}
		nodeIdToNumDownstreamnodes.insert({ current.atomid, numDownstreamnodes });
	}

	return nodeIdToNumDownstreamnodes;
}


//MoleculeGraph::MoleculeGraph(const Node* root) {
//	for (const auto node : BFSRange(root)) {
//		nodes.insert({ node.atomid, node });
//	}
//}

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
	//const MoleculeGraph::Node* const initialGuessForRootnode = &molgraph.nodes.find(molgraph.highestNodeId)->second;
	const MoleculeGraph::Node* const initialGuessForRootnode = &(molgraph.nodes.begin()->second);

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
	std::vector<int> mapping(molgraph.nodes.size(), -1);
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


void LimaMoleculeGraph::reorderoleculeParticlesAccoringingToSubchains(GroFile& grofile, TopologyFile::Moleculetype& molecule) {

	std::vector<std::array<std::string, 2>> sbAtomtypesCheck;
	for (const auto& sb : molecule.singlebonds) {
		//sbAtomtypesCheck.push_back({ topfile.GetLocalAtoms()[sb.ids[0]].atomname, topfile.GetLocalAtoms()[sb.ids[1]].atomname });
		sbAtomtypesCheck.push_back({ molecule.atoms[sb.ids[0]].atomname, molecule.atoms[sb.ids[1]].atomname });
	}

	const MoleculeGraph molgraph(molecule);

	const std::vector<int> map = MakeParticleReorderMapping(molgraph);

	// Overwrite all references to gro_ids in the files	
	for (auto& atom : grofile.atoms) {
		if (map[atom.gro_id-1] < 0)
			throw std::runtime_error("Invalid gro_id in map");
		atom.gro_id = map[atom.gro_id-1] + 1;  // +1 to convert from lima id to groID
	}
	// We can only do this part if the topol does NOT have includes
	/*if (molecule.GetLocalMolecules().size() > 0)
		throw std::runtime_error("Cannot reorder topol with includes");*/
	for (auto& atom : molecule.atoms) {
		atom.id = map[atom.id];
	}

	overwriteParticleIds<>(molecule.singlebonds, map);
	overwriteParticleIds<>(molecule.pairs, map);
	overwriteParticleIds<>(molecule.anglebonds, map);
	overwriteParticleIds<>(molecule.dihedralbonds, map);
	overwriteParticleIds<>(molecule.improperdihedralbonds, map);



	// Re-sort all entries with the new groids
	std::sort(grofile.atoms.begin(), grofile.atoms.end(), [](const GroRecord& a, const GroRecord& b) {return a.gro_id < b.gro_id; });

	std::sort(molecule.atoms.begin(), molecule.atoms.end(), [](const auto& a, const auto& b) {return a.id < b.id; });

	// Check that we didn't mess up the singlebonds. If all these are identical, the other bonds must be aswell
	for (int i = 0; i < molecule.singlebonds.size(); i++) {
		if (
			sbAtomtypesCheck[i][0] != molecule.atoms[molecule.singlebonds[i].ids[0]].atomname ||
			sbAtomtypesCheck[i][1] != molecule.atoms[molecule.singlebonds[i].ids[1]].atomname
			)
		{
			auto a = molecule.atoms[molecule.singlebonds[i].ids[0]].atomname;
			auto b = molecule.atoms[molecule.singlebonds[i].ids[1]].atomname;
			throw std::runtime_error("Reordering of particles messed up singlebonds" + std::to_string(i));
		}
	}

	std::sort(molecule.singlebonds.begin(), molecule.singlebonds.end(), [](const auto& a, const auto& b) { return a.ids[0] < b.ids[0]; });
	std::sort(molecule.pairs.begin(), molecule.pairs.end(), [](const auto& a, const auto& b) { return a.ids[0] < b.ids[0]; });
	std::sort(molecule.anglebonds.begin(), molecule.anglebonds.end(), [](const auto& a, const auto& b) { return a.ids[0] < b.ids[0]; });
	std::sort(molecule.dihedralbonds.begin(), molecule.dihedralbonds.end(), [](const auto& a, const auto& b) { return a.ids[0] < b.ids[0]; });
	std::sort(molecule.improperdihedralbonds.begin(), molecule.improperdihedralbonds.end(), [](const auto& a, const auto& b) { return a.ids[0] < b.ids[0]; });


	for (int i = 0; i < molecule.atoms.size(); i++) {
		if (molecule.atoms[i].id != i)
			throw std::runtime_error("Reordering of particles messed up atoms" + std::to_string(i));
	}



}




int MarkAllNodes(std::unordered_set<int>& visited, MoleculeGraph::BFSRange<const MoleculeGraph::Node*> nodes) {
	int count = 0;
	for (const auto& node : nodes) {
		visited.insert(node.atomid);
		count++;
	}
	return count;

}

bool MoleculeGraph::GraphIsDisconnected() const {
	std::unordered_set<int> visited;
	
	MarkAllNodes(visited, BFSRange(&nodes.at(0)));

	return visited.size() != nodes.size();
}

std::vector<std::vector<int>> MoleculeGraph::GetListOfListsofConnectedNodeids() const {
	std::vector<vector<int>> subGraphs;
	std::unordered_set<int> visited;

	for (int i = 0; i < nodes.size(); i++) {
		if (visited.contains(i))
			continue;
		// Start a new empty molecule
		subGraphs.push_back({});
		// Fill said molecule
		for (const auto& node : BFS(i)) {
			visited.insert(node.atomid);
			subGraphs.back().emplace_back(node.atomid);
		}
	}

	return subGraphs;
}

//std::vector<MoleculeGraph> MoleculeGraph::GetListOfConnectedGraphs() const {
//	std::vector<MoleculeGraph> subGraphs;
//	std::unordered_set<int> visited;
//
//	for (int i = 0; i < nodes.size(); i++) {
//		if (visited.contains(i))
//			continue;
//		
//		MarkAllNodes(visited, BFS(i));
//		subGraphs.emplace_back(MoleculeGraph(&nodes.at(i)));		
//	}
//
//	return subGraphs;
//}
