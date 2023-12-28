#include "MoleculeGraph.h"


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

MoleculeGraph createGraph(const SimpleParsedFile& topolfile) {
	MoleculeGraph graph;

	for (auto& row : topolfile.rows) 
	{
		// Insert nodes
		if (row.section == "atoms") {
			const string atomname = row.words[1];
			const int id = std::stoi(row.words[2]);
			graph.addNode(id, atomname);
		}
		else if (row.section == "bonds") {
			graph.connectNodes(std::stoi(row.words[0]), std::stoi(row.words[1]));
		}
	}

	return graph;
}

int getIdOfValidChainStartingPoint(const MoleculeGraph& molgraph, const std::set<int>& ids_already_in_a_chain) {
	for (const auto& elem : molgraph.nodes) {
		const MoleculeGraph::Node& node = elem.second;

		if (
			!node.isHydrogen() && 
			!ids_already_in_a_chain.contains(node.atomid) &&
			node.getNNonhydrogenNeighbors() < 2
			)

			return node.atomid;
	}

	return -1;	// No valid chain starting point found
}

// Returns the next node that is NOT hydrogen, and not already in a node
int getIdOfNextValidNodeInChain(const std::set<int>& ids_already_in_a_chain, const std::vector<MoleculeGraph::Node*>& currentnode_neighbors)
{
	for (const MoleculeGraph::Node* neighbornode : currentnode_neighbors) {
		if (!neighbornode->isHydrogen() && !ids_already_in_a_chain.contains(neighbornode->atomid))
			return neighbornode->atomid;
	}

	return -1;
}

void addConnectedHydrogensToChain(std::set<int>& ids_already_in_a_chain, const std::vector<MoleculeGraph::Node*>& currentnode_neighbors, Chain& chain) {
	for (const MoleculeGraph::Node* neighbornode : currentnode_neighbors) {
		if (neighbornode->isHydrogen())
			chain.insert(neighbornode->atomid, ids_already_in_a_chain);
	}
}

vector<Chain> LimaMoleculeGraph::getAllSubchains(const MoleculeGraph& molgraph) {
	vector<Chain> chains{};
	std::set<int> ids_already_in_a_chain{};

	while (true)
	{
		const int id_of_starting_node = getIdOfValidChainStartingPoint(molgraph, ids_already_in_a_chain);
		if (id_of_starting_node == -1) { break; }


		// Create new chain
		chains.emplace_back(Chain{});
		const MoleculeGraph::Node* current_node = &molgraph.nodes.at(4);

		chains.back().insert(id_of_starting_node, ids_already_in_a_chain);
		ids_already_in_a_chain.insert(id_of_starting_node);

		while (true) {

			// is the current node an intersection of chains, then find a new chain
			if (current_node->getNNonhydrogenNeighbors() > 2)
				break;

			addConnectedHydrogensToChain(ids_already_in_a_chain, current_node->getNeighbors(), chains.back());

			const int next_id = getIdOfNextValidNodeInChain(ids_already_in_a_chain, current_node->getNeighbors());
			if (next_id == -1)
				break;

			chains.back().insert(next_id, ids_already_in_a_chain);
			current_node = &molgraph.nodes.at(next_id);
		}

	}

	return chains;
}

void LimaMoleculeGraph::reorderoleculeParticlesAccoringingToSubchains(fs::path gro_path, fs::path top_path) {

}