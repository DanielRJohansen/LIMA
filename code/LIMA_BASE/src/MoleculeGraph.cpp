#include "MoleculeGraph.h"
#include <stack>

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

MoleculeGraph LimaMoleculeGraph::createGraph(const ParsedTopologyFile& topolfile) {
	MoleculeGraph graph;

	for (const auto& atom : topolfile.atoms.entries) {
		graph.addNode(atom.nr, atom.atom);
	}

	for (const auto& bond : topolfile.singlebonds.entries) {
		graph.connectNodes(bond.atom_indexes[0], bond.atom_indexes[1]);
	}

	return graph;
}

int getIdOfValidChainStartingPoint(const MoleculeGraph& molgraph) {
	for (const auto& elem : molgraph.nodes) {
		const MoleculeGraph::Node& node = elem.second;

		if (!node.isHydrogen() && node.getNNonhydrogenNeighbors() < 2)
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
			chain.append(*neighbornode, ids_already_in_a_chain);
	}
}



vector<Chain> LimaMoleculeGraph::getAllSubchains(const MoleculeGraph& molgraph) {
	vector<Chain> chains{};
	std::set<int> ids_already_in_a_chain{};

	// Stack used to perform depth first search for chains, which ensure coherent ordering of chains
	std::stack<int> chainstartNodeIdsStack;
	chainstartNodeIdsStack.push(getIdOfValidChainStartingPoint(molgraph));

	// Assign all particles to chains
	while (true)
	{
		//const int id_of_starting_node = getIdOfValidChainStartingPoint(molgraph, ids_already_in_a_chain);
		if (chainstartNodeIdsStack.empty()) { break; }
		const int id_of_starting_node = chainstartNodeIdsStack.top();
		chainstartNodeIdsStack.pop();
		if (id_of_starting_node == -1) { break; }


		// Create new chain
		chains.emplace_back(Chain{});
		const MoleculeGraph::Node* current_node = &molgraph.nodes.at(4);

		int next_id = id_of_starting_node;
		while (true) {
			current_node = &molgraph.nodes.at(next_id);

			chains.back().append(molgraph.nodes.at(next_id), ids_already_in_a_chain);

			addConnectedHydrogensToChain(ids_already_in_a_chain, current_node->getNeighbors(), chains.back());


			// Is the current node an intersection of chains, then find the chainstart ids of adjecent chains, and terminate current chain
			if (current_node->getNNonhydrogenNeighbors() > 2) {
				for (int i = 0; i < current_node->getNNonhydrogenNeighbors(); i++) {

					// Push all nearby non-hydrogens to the stack of chain-starts
					const MoleculeGraph::Node& neighbor = *current_node->getNeighbors()[i];
					if (!neighbor.isHydrogen() && !ids_already_in_a_chain.contains(neighbor.atomid))
						chainstartNodeIdsStack.push(neighbor.atomid);
				}

				break;
			}


			next_id = getIdOfNextValidNodeInChain(ids_already_in_a_chain, current_node->getNeighbors());
			if (next_id == -1)
				break;

		}
	}

	// Finally find parent chain of each chain
	const int n_chains = chains.size();

	for (int rhs_chain_id = 1; rhs_chain_id < n_chains; rhs_chain_id++) {
		const int rhs_chain_start_nodeid = chains[rhs_chain_id].getNodeIds().at(0);
		const MoleculeGraph::Node& rhs_chain_startnode = molgraph.nodes.at(rhs_chain_start_nodeid);

		// Lambda to find parent chain, so i can throw if we fail.
		auto findParentchainId = [&]() {
			// Loop over all chains left of current chain
			for (int lhs_chain_id = 0; lhs_chain_id < rhs_chain_id; lhs_chain_id++) {
				const int lhs_chain_end_nodeid = chains[lhs_chain_id].getBack();

				if (rhs_chain_startnode.isConnected(lhs_chain_end_nodeid)) {
					return lhs_chain_id;
				}
			}
			throw std::runtime_error("Failed to find parent chain");
		};

		chains[rhs_chain_id].parentchain_id = findParentchainId();	
	}

	return chains;
}




// Returns a mapping where from=index in vector, and to= value on that index
// Remember that index 0 must not be used!
std::vector<int> makeParticleReorderMapping(const std::vector<Chain>& subchains) {
	std::vector<int> map(1);
	map[0] = -1;
	int next_new_id = 1;	// Sadly MD files are 1d indexes

	for (const Chain& chain : subchains) {
		for (const int id : chain.getNodeIds()) {
			if (id >= map.size())
				map.resize(id + 1);

			map[id] = next_new_id++;
		}
	}
	return map;
}

void mergeTinySidechainsIntoParentchains(std::vector<Chain>& chains) {
	const int cutoff_chain_size = 5;	// Including hydrogens

	for (int chain_id = 1; chain_id < chains.size(); chain_id++) {
		const Chain& chain = chains[chain_id];
		if (chain.getNNodesInChain() < cutoff_chain_size) {
			// Merge the tinychain into its parent
			chains[chain.parentchain_id].insert(chain);

			// All chains that has the current chain as parent, should update to the chains parent. 
			// Chains that has a parent id higher than the current chain, should dec that index
			for (int i = chain_id + 1; i < chains.size(); i++) {
				if (chains[i].parentchain_id == chain_id)
					chains[i].parentchain_id = chain.parentchain_id;
				else if (chains[i].parentchain_id > chain_id)
					chains[i].parentchain_id--;
			}

			// Move all subsequent chains 1 step forward
			for (int i = chain_id; i < chains.size() - 1; i++) {
				chains[i] = chains[i + 1];
			}
			chains.pop_back();

			// Now list is 1 shorter, so dec iterator
			chain_id--;
		}
	}
}

// Reorder so subsequent subchains always comes with the shortest one first
void reorderSubchains(std::vector<Chain>& subchains, const MoleculeGraph& molgraph)
{
	// Two chains that have the same parent but are out of order should be swapped
	auto sortCondition = [](const Chain& a, const Chain& b) {
		if (a.parentchain_id == b.parentchain_id && a.getNNodesInChain() > b.getNNodesInChain()) { return true; }
		return false;
		};
	std::sort(subchains.begin(), subchains.end(), sortCondition);
}

template<typename T>
void overwriteParticleIds(std::vector<T> bonds, const std::vector<int>& map) {
	for (auto bond : bonds) {
		for (int i = 0; i < bond.n; i++) {
			bond.atom_indexes[i] = map[bond.atom_indexes[i]];
		}
	}	
}

void LimaMoleculeGraph::reorderoleculeParticlesAccoringingToSubchains(fs::path gro_path, fs::path top_path) {
	ParsedGroFile grofile = MDFiles::loadGroFile(gro_path);
	ParsedTopologyFile topfile = MDFiles::loadTopologyFile(top_path);

	const MoleculeGraph molgraph = createGraph(topfile);
	std::vector<Chain> chains = getAllSubchains(molgraph);

	mergeTinySidechainsIntoParentchains(chains);
	reorderSubchains(chains, molgraph);

	// dont use index 0 of map
	const std::vector<int> map = makeParticleReorderMapping(chains);

	// Overwrite all references to gro_ids in the files
	for (auto& atom : grofile.atoms) {
		atom.gro_id = map[atom.gro_id];
	}

	for (auto atom : topfile.atoms.entries) {
		atom.nr = map[atom.nr];
	}
	overwriteParticleIds<>(topfile.singlebonds.entries, map);
	overwriteParticleIds<>(topfile.pairs.entries, map);
	overwriteParticleIds<>(topfile.anglebonds.entries, map);
	overwriteParticleIds<>(topfile.dihedralbonds.entries, map);
	overwriteParticleIds<>(topfile.improperdihedralbonds.entries, map);

	// Re-sort all entries with the new groids
	std::sort(grofile.atoms.begin(), grofile.atoms.end(), [](const GroRecord& a, const GroRecord& b) {return a.gro_id < b.gro_id; });

	std::sort(topfile.atoms.entries.begin(), topfile.atoms.entries.end(), [](const auto& a, const auto& b) {return a.nr > b.nr; });
	std::sort(topfile.singlebonds.entries.begin(), topfile.singlebonds.entries.end(), [](const auto& a, const auto& b) {return a.atom_indexes[0] > b.atom_indexes[0]; });
	std::sort(topfile.pairs.entries.begin(), topfile.pairs.entries.end(), [](const auto& a, const auto& b) { return a.atom_indexes[0] > b.atom_indexes[0]; });
	std::sort(topfile.anglebonds.entries.begin(), topfile.anglebonds.entries.end(), [](const auto& a, const auto& b) { return a.atom_indexes[0] > b.atom_indexes[0]; });
	std::sort(topfile.dihedralbonds.entries.begin(), topfile.dihedralbonds.entries.end(), [](const auto& a, const auto& b) { return a.atom_indexes[0] > b.atom_indexes[0]; });
	std::sort(topfile.improperdihedralbonds.entries.begin(), topfile.improperdihedralbonds.entries.end(), [](const auto& a, const auto& b) { return a.atom_indexes[0] > b.atom_indexes[0]; });

	gro_path.replace_filename(gro_path.stem().string() + "_new" + gro_path.extension().string());
	top_path.replace_filename(top_path.stem().string() + "_new" + top_path.extension().string());

	grofile.printToFile(gro_path);
	topfile.printToFile(top_path);
}