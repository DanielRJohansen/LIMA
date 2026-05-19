#include "MoleculeGraph.h"

#include <algorithm>
#include <cassert>
#include <queue>
#include <ranges>
#include <stack>
#include <unordered_map>
#include <unordered_set>

#include "TimeIt.h"

using namespace LimaMoleculeGraph;
using std::string;
using std::vector;


MoleculeTree::MoleculeTree(int rootId) {
	tree.insert({ rootId, {} });
}

void MoleculeTree::ForSelfAndAllChildrenIds(
	int parentId,
	const std::function<void(int)>& visitor
) const
{
	std::stack<int> nodeStack;
	nodeStack.push(parentId);

	while (!nodeStack.empty()) {
		const int currentId = nodeStack.top();
		nodeStack.pop();

		visitor(currentId);

		for (const int childId : tree.at(currentId)) {
			nodeStack.push(childId);
		}
	}
}


void MoleculeGraph::Node::addNeighbor(int neighborId, bool neighborIsHydrogen) {
	if (nNeighbors >= maxNeighbors)
		throw std::runtime_error("Exceeded maximum number of neighbors for node " + std::to_string(atomid));
	if (!neighborIsHydrogen)
		n_nonhydrogen_neighbors++;
	neighborIds[nNeighbors++] = neighborId;
}

void MoleculeGraph::connectNodes(int left_id, int right_id) {
	Node& left = nodes.at(left_id);
	Node& right = nodes.at(right_id);

	left.addNeighbor(right_id, right.isHydrogen());
	right.addNeighbor(left_id, left.isHydrogen());
}

MoleculeGraph::MoleculeGraph(const TopologyFile::Moleculetype& molecule, std::optional<const std::unordered_set<int>> allowedIds) {
	std::vector<std::pair<int, Node>> temp_nodes;
	temp_nodes.reserve(molecule.atoms.size());

	for (const auto& atom : molecule.atoms) {
		if (allowedIds.has_value() && !allowedIds.value().contains(atom.id))
			continue;

		temp_nodes.emplace_back(atom.id, Node(atom.id, atom.atomname));
	}

	nodes.insert(temp_nodes.begin(), temp_nodes.end());

	for (const auto& bond : molecule.singlebonds) {
		if (!nodes.contains(bond.ids[0]) || !nodes.contains(bond.ids[1]))
			continue;

		connectNodes(bond.ids[0], bond.ids[1]);
	}
}

MoleculeGraph::MoleculeGraph(
	const std::vector<std::pair<int, std::string>>& atoms,
	const std::vector<std::array<int, 2>>& edges
) {
	std::vector<std::pair<int, Node>> temp_nodes;
	temp_nodes.reserve(atoms.size());

	for (const auto& atom : atoms) {
		temp_nodes.emplace_back(atom.first, Node(atom.first, atom.second));
	}

	nodes.insert(temp_nodes.begin(), temp_nodes.end());

	for (const auto& edge : edges) {
		connectNodes(edge[0], edge[1]);
	}
}


MoleculeTree MoleculeGraph::ConstructMoleculeTree() const {
	if (nodes.empty())
		throw std::runtime_error("Cannot construct molecule tree from empty graph");

	std::unordered_set<int> visited;
	const int rootId = nodes.begin()->first;
	MoleculeTree moleculeTree(rootId);

	std::queue<int> nodeQueue;
	nodeQueue.push(rootId);
	visited.insert(rootId);

	while (!nodeQueue.empty()) {
		const int currentId = nodeQueue.front();
		nodeQueue.pop();

		for (const int neighborId : nodes.at(currentId).getNeighbors()) {
			if (visited.insert(neighborId).second) {
				nodeQueue.push(neighborId);
				moleculeTree.AddChild(currentId, neighborId);
			}
		}
	}

	return moleculeTree;
}


std::unordered_map<int, int> MoleculeGraph::ComputeNumDownstreamNodes(const MoleculeTree& moleculeTree) const
{
	if (nodes.empty())
		return {};

	std::stack<int> processStack;
	for (const Node& node : BFS(nodes.begin()->first)) {
		processStack.push(node.atomid);
	}

	std::unordered_map<int, int> nodeIdToNumDownstreamnodes;

	while (!processStack.empty()) {
		const int currentId = processStack.top();
		processStack.pop();

		int numDownstreamnodes = 1;
		for (const int neighborId : moleculeTree.GetChildIds(currentId)) {
			numDownstreamnodes += nodeIdToNumDownstreamnodes.at(neighborId);
		}

		nodeIdToNumDownstreamnodes.insert({ currentId, numDownstreamnodes });
	}

	return nodeIdToNumDownstreamnodes;
}


int GetNumNonvisitedNonHydrogenNeighbors(
	const MoleculeGraph& molgraph,
	int nodeId,
	const std::unordered_set<int>& visitedNodes
) {
	int nNonvisitedNonHydrogenNeighbors = 0;

	for (const int neighborId : molgraph.nodes.at(nodeId).getNeighbors()) {
		const MoleculeGraph::Node& neighbor = molgraph.nodes.at(neighborId);
		if (!visitedNodes.contains(neighborId) && !neighbor.isHydrogen())
			nNonvisitedNonHydrogenNeighbors++;
	}

	return nNonvisitedNonHydrogenNeighbors;
}

std::vector<int> GetNonvisitedNonHydrogenNeighbors(
	const MoleculeGraph& molgraph,
	int nodeId,
	const std::unordered_set<int>& visitedNodes
) {
	std::vector<int> neighbors;

	for (const int neighborId : molgraph.nodes.at(nodeId).getNeighbors()) {
		const MoleculeGraph::Node& neighbor = molgraph.nodes.at(neighborId);
		if (!visitedNodes.contains(neighborId) && !neighbor.isHydrogen())
			neighbors.push_back(neighborId);
	}

	return neighbors;
}


std::pair<int, int> FindFurthestNode(
	const MoleculeGraph& molgraph,
	int startNodeId,
	std::unordered_set<int> visitedNodes = {}
) {
	std::queue<std::pair<int, int>> q;

	q.push({ startNodeId, 0 });
	visitedNodes.insert(startNodeId);

	int furthestNodeId = startNodeId;
	int maxDistance = 0;

	while (!q.empty()) {
		const auto [currentNodeId, currentDistance] = q.front();
		q.pop();

		if (currentDistance > maxDistance) {
			furthestNodeId = currentNodeId;
			maxDistance = currentDistance;
		}

		const auto neighbors = GetNonvisitedNonHydrogenNeighbors(molgraph, currentNodeId, visitedNodes);

		for (const int neighborId : neighbors) {
			if (!visitedNodes.contains(neighborId)) {
				visitedNodes.insert(neighborId);
				q.push({ neighborId, currentDistance + 1 });
			}
		}
	}

	return { furthestNodeId, maxDistance };
}

int FindRootNodeInGraph(const MoleculeGraph& molgraph) {
	const int initialGuessForRootnode = molgraph.nodes.begin()->first;

	auto [secondGuessForRootnode, depth1] = FindFurthestNode(molgraph, initialGuessForRootnode);

	auto [thirdGuessForRootnode, depth2] = FindFurthestNode(molgraph, secondGuessForRootnode);

	assert(depth2 >= depth1 - 1);

	return thirdGuessForRootnode;
}


void AddNodeToMapping(const MoleculeGraph& molgraph, int nodeId, std::vector<int>& mapping, int& next_new_id) {
	if (nodeId >= mapping.size())
		throw std::runtime_error("Node id is too high, shoudnt happen as we keep track of highest id");

	if (mapping[nodeId] != -1)
		int a = 90;

	mapping[nodeId] = next_new_id++;

	for (const int neighborId : molgraph.nodes.at(nodeId).getNeighbors()) {
		const MoleculeGraph::Node& neighbor = molgraph.nodes.at(neighborId);

		if (neighbor.isHydrogen()) {
			assert(mapping[neighborId] == -1);
			mapping[neighborId] = next_new_id++;
		}
	}
}

void SortSidechainsByDepth(
	const MoleculeGraph& molgraph,
	const std::unordered_set<int>& visitedNodes,
	std::vector<int>& neighbors
) {
	std::vector<int> sidechainDepths(neighbors.size());

	for (int i = 0; i < neighbors.size(); i++) {
		sidechainDepths[i] = FindFurthestNode(molgraph, neighbors[i], visitedNodes).second;
	}

	std::vector<std::pair<int, int>> neighborsWithDepths;
	neighborsWithDepths.reserve(neighbors.size());

	for (int i = 0; i < neighbors.size(); i++) {
		neighborsWithDepths.emplace_back(neighbors[i], sidechainDepths[i]);
	}

	std::sort(
		neighborsWithDepths.begin(),
		neighborsWithDepths.end(),
		[](const std::pair<int, int>& a, const std::pair<int, int>& b) {
			return a.second < b.second;
		}
	);

	neighbors.clear();

	for (const auto& pair : neighborsWithDepths) {
		neighbors.push_back(pair.first);
	}
}

void AddSidechainToMapping(
	const MoleculeGraph& molgraph,
	std::unordered_set<int>& visitedNodes,
	std::vector<int>& mapping,
	int& next_new_id,
	int currentNodeId
) {
	if (visitedNodes.contains(currentNodeId))
		return;

	while (true) {
		AddNodeToMapping(molgraph, currentNodeId, mapping, next_new_id);
		visitedNodes.insert(currentNodeId);

		const int nNonvisitedNonhydrogenNeighbors = GetNumNonvisitedNonHydrogenNeighbors(molgraph, currentNodeId, visitedNodes);

		if (nNonvisitedNonhydrogenNeighbors == 0) {
			return;
		}
		else if (nNonvisitedNonhydrogenNeighbors == 1) {
			currentNodeId = GetNonvisitedNonHydrogenNeighbors(molgraph, currentNodeId, visitedNodes)[0];
			continue;
		}
		else {
			std::vector<int> neighbors = GetNonvisitedNonHydrogenNeighbors(molgraph, currentNodeId, visitedNodes);
			SortSidechainsByDepth(molgraph, visitedNodes, neighbors);

			for (const int neighborId : neighbors) {
				AddSidechainToMapping(molgraph, visitedNodes, mapping, next_new_id, neighborId);
			}

			return;
		}
	}
}

std::vector<int> MakeParticleReorderMapping(const MoleculeGraph& molgraph) {
	if (molgraph.nodes.empty())
		return {};

	const int currentNodeId = FindRootNodeInGraph(molgraph);

	std::unordered_set<int> visitedNodes;
	std::vector<int> mapping(molgraph.nodes.rbegin()->first + 1, -1);
	int next_new_id = 0;

	AddSidechainToMapping(molgraph, visitedNodes, mapping, next_new_id, currentNodeId);

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
		sbAtomtypesCheck.push_back({ molecule.atoms[sb.ids[0]].atomname, molecule.atoms[sb.ids[1]].atomname });
	}

	const MoleculeGraph molgraph(molecule);

	const std::vector<int> map = MakeParticleReorderMapping(molgraph);

	for (auto& atom : grofile.atoms) {
		if (map[atom.gro_id - 1] < 0)
			throw std::runtime_error("Invalid gro_id in map");

		atom.gro_id = map[atom.gro_id - 1] + 1;
	}

	for (auto& atom : molecule.atoms) {
		atom.id = map[atom.id];
	}

	overwriteParticleIds<>(molecule.singlebonds, map);
	overwriteParticleIds<>(molecule.pairbonds, map);
	overwriteParticleIds<>(molecule.anglebonds, map);
	overwriteParticleIds<>(molecule.dihedralbonds, map);
	overwriteParticleIds<>(molecule.improperdihedralbonds, map);

	std::sort(grofile.atoms.begin(), grofile.atoms.end(), [](const GroRecord& a, const GroRecord& b) { return a.gro_id < b.gro_id; });

	std::sort(molecule.atoms.begin(), molecule.atoms.end(), [](const auto& a, const auto& b) { return a.id < b.id; });

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
	std::sort(molecule.pairbonds.begin(), molecule.pairbonds.end(), [](const auto& a, const auto& b) { return a.ids[0] < b.ids[0]; });
	std::sort(molecule.anglebonds.begin(), molecule.anglebonds.end(), [](const auto& a, const auto& b) { return a.ids[0] < b.ids[0]; });
	std::sort(molecule.dihedralbonds.begin(), molecule.dihedralbonds.end(), [](const auto& a, const auto& b) { return a.ids[0] < b.ids[0]; });
	std::sort(molecule.improperdihedralbonds.begin(), molecule.improperdihedralbonds.end(), [](const auto& a, const auto& b) { return a.ids[0] < b.ids[0]; });

	for (int i = 0; i < molecule.atoms.size(); i++) {
		if (molecule.atoms[i].id != i)
			throw std::runtime_error("Reordering of particles messed up atoms" + std::to_string(i));
	}
}


int MarkAllNodes(std::unordered_set<int>& visited, MoleculeGraph::BFSRange<const MoleculeGraph> nodes) {
	int count = 0;

	for (const auto& node : nodes) {
		visited.insert(node.atomid);
		count++;
	}

	return count;
}

std::optional<int> MoleculeGraph::DistanceBetweenNodes(int id0, int id1, int maxSearchDepth) const {
	auto bfs = BFS(id0);

	for (auto it = bfs.begin(); it != bfs.end(); ++it) {
		const Node& node = *it;

		if (node.atomid == id1)
			return it.Depth();

		if (it.Depth() > maxSearchDepth)
			return std::nullopt;
	}

	return std::nullopt;
}


bool MoleculeGraph::GraphIsDisconnected() const {
	if (nodes.empty())
		return false;

	std::unordered_set<int> visited;

	MarkAllNodes(visited, BFS(nodes.begin()->first));

	return visited.size() != nodes.size();
}

std::vector<std::vector<int>> MoleculeGraph::GetListOfListsofConnectedNodeids() const {
	std::vector<vector<int>> subGraphs;
	std::unordered_set<int> visited;
	visited.reserve(nodes.size());

	for (const auto& [nodeId, node] : nodes) {
		if (visited.contains(nodeId))
			continue;

		subGraphs.push_back({});

		for (const auto& bfsNode : BFS(nodeId)) {
			visited.insert(bfsNode.atomid);
			subGraphs.back().emplace_back(bfsNode.atomid);
		}
	}

	return subGraphs;
}