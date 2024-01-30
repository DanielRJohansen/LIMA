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

void Chain::append(const MoleculeGraph::Node& node, std::set<int>& ids_already_in_a_chain) {
	if (ids_already_in_a_chain.contains(node.atomid))
		throw std::runtime_error("Didn't expect this atom to already be in a chain");
	ids_already_in_a_chain.insert(node.atomid);

	nodeids.push_back(node.atomid);

	if (!node.isHydrogen()) {
		backbone_len++;
		height++;
	}
}

void Chain::append(std::unique_ptr<Chain> chain) {
	height = std::max(height, chain->height + backbone_len);
	subchains.push_back(std::move(chain));
}

void Chain::sort() {
	auto sort_condition = [](auto& a, auto& b) { return a->getHeight() < b->getHeight(); };
	std::sort(subchains.begin(), subchains.end(), sort_condition);

	for (auto& subchain : subchains)
		subchain->sort();
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

std::unique_ptr<Chain> makeSubChain(int start_id, const MoleculeGraph& molgraph, std::set<int>& ids_already_in_a_chain) {
	auto chain = std::make_unique<Chain>();

	int current_id = start_id;
	while (true) {
		const MoleculeGraph::Node* current_node = &molgraph.nodes.at(current_id);
		
		// Add the current node in backbone
		chain->append(molgraph.nodes.at(current_id), ids_already_in_a_chain);

		// Add hydrogens connected to node
		addConnectedHydrogensToChain(ids_already_in_a_chain, current_node->getNeighbors(), *chain);

		// Is the current node an intersection of chains, then find the chainstart ids of adjecent chains, and terminate current chain
		if (current_node->getNNonhydrogenNeighbors() > 2) {

			for (int i = 0; i < current_node->getNNeighbors(); i++) {

				// Push all nearby non-hydrogens to the stack of chain-starts
				const MoleculeGraph::Node& neighbor = *current_node->getNeighbors()[i];
				if (!neighbor.isHydrogen() && !ids_already_in_a_chain.contains(neighbor.atomid)) {
					chain->append(std::move(makeSubChain(neighbor.atomid, molgraph, ids_already_in_a_chain)));
				}
			}
			break;
		}

		// Go to next node in backbone of chain
		current_id = getIdOfNextValidNodeInChain(ids_already_in_a_chain, current_node->getNeighbors());
		if (current_id == -1)
			break;
	}
	
	return std::move(chain);
}

// Returns the root chain
std::unique_ptr<Chain> getChainTree(const MoleculeGraph& molgraph) {
	
	std::set<int> ids_already_in_a_chain{};

	const int start_id = getIdOfValidChainStartingPoint(molgraph);

	std::unique_ptr<Chain> root = makeSubChain(start_id, molgraph, ids_already_in_a_chain);

	return root;
}

void makeParticleReorderMapping(const Chain& chain, std::vector<int>& map, int& next_new_id) {
	for (const int id : chain.getNodeIds()) {
		if (id >= map.size())
			map.resize(id + 2);

		map[id] = next_new_id++;
	}

	for (auto& subchain : chain.getSubchains()) {
		makeParticleReorderMapping(*subchain, map, next_new_id);
	}
}

// Returns a mapping where from=index in vector, and to= value on that index
// Remember that index 0 must not be used!
std::vector<int> makeParticleReorderMapping(const Chain& root_chain) {
	std::vector<int> map(1);
	map[0] = -1;
	int next_new_id = 1;	// Sadly MD files are 1d indexes

	makeParticleReorderMapping(root_chain, map, next_new_id);

	return map;
}

template<typename T>
void overwriteParticleIds(std::vector<T>& bonds, const std::vector<int>& map) {
	for (auto& bond : bonds) {
		for (int i = 0; i < bond.n; i++) {
			bond.atom_indexes[i] = map[bond.atom_indexes[i]];
		}
	}	
}

void LimaMoleculeGraph::reorderoleculeParticlesAccoringingToSubchains(const fs::path& gro_path_in, const fs::path& top_path_in, const fs::path& gro_path_out, const fs::path& top_path_out) {
	auto grofile = MDFiles::loadGroFile(gro_path_in);
	auto topfile = MDFiles::loadTopologyFile(top_path_in);

	const MoleculeGraph molgraph = createGraph(*topfile);
	std::unique_ptr<Chain> root_chain = getChainTree(molgraph);

	root_chain->sort();

	// dont use index 0 of map
	const std::vector<int> map = makeParticleReorderMapping(*root_chain);

	// Overwrite all references to gro_ids in the files
	for (auto& atom : grofile->atoms) {
		atom.gro_id = map[atom.gro_id];
	}

	for (auto& atom : topfile->atoms.entries) {
		atom.nr = map[atom.nr];
	}
	overwriteParticleIds<>(topfile->singlebonds.entries, map);
	overwriteParticleIds<>(topfile->pairs.entries, map);
	overwriteParticleIds<>(topfile->anglebonds.entries, map);
	overwriteParticleIds<>(topfile->dihedralbonds.entries, map);
	overwriteParticleIds<>(topfile->improperdihedralbonds.entries, map);

	// Re-sort all entries with the new groids
	std::sort(grofile->atoms.begin(), grofile->atoms.end(), [](const GroRecord& a, const GroRecord& b) {return a.gro_id < b.gro_id; });

	std::sort(topfile->atoms.entries.begin(), topfile->atoms.entries.end(), [](const auto& a, const auto& b) {return a.nr < b.nr; });
	std::sort(topfile->singlebonds.entries.begin(), topfile->singlebonds.entries.end(), [](const auto& a, const auto& b) {return a.atom_indexes[0] < b.atom_indexes[0]; });
	std::sort(topfile->pairs.entries.begin(), topfile->pairs.entries.end(), [](const auto& a, const auto& b) { return a.atom_indexes[0] < b.atom_indexes[0]; });
	std::sort(topfile->anglebonds.entries.begin(), topfile->anglebonds.entries.end(), [](const auto& a, const auto& b) { return a.atom_indexes[0] < b.atom_indexes[0]; });
	std::sort(topfile->dihedralbonds.entries.begin(), topfile->dihedralbonds.entries.end(), [](const auto& a, const auto& b) { return a.atom_indexes[0] < b.atom_indexes[0]; });
	std::sort(topfile->improperdihedralbonds.entries.begin(), topfile->improperdihedralbonds.entries.end(), [](const auto& a, const auto& b) { return a.atom_indexes[0] < b.atom_indexes[0]; });

	grofile->printToFile(gro_path_out);
	topfile->printToFile(top_path_out);
}