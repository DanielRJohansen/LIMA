#pragma once

#include <vector>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <set>

#include "Filehandling.h"
#include <filesystem>
#include "MDFiles.h"

namespace LimaMoleculeGraph {
	namespace fs = std::filesystem;

	struct MoleculeGraph {
		struct Node {
			Node() {}	// Not sure why this is needed to compile, but dont use!
			Node(int id, const std::string& atomname) : atomid(id), atomname(atomname) { 
				if (atomname == "") throw std::runtime_error("Cannot add noname atom to graph"); 
			};
			void addNeighbor(Node* neighbor);

			bool isHydrogen() const { return atomname[0] == 'H'; }
			int getNNeighbors() const { return neighbors.size(); }
			int getNNonhydrogenNeighbors() const { return n_nonhydrogen_neighbors; }
			const std::vector<Node*> getNeighbors() const { return neighbors; }

			bool isConnected(int id) const {
				for (const auto neighbor : neighbors) {
					if (neighbor->atomid == id)
						return true;
				}
				return false;
			}

			const int atomid{-1};

		private:
			const std::string atomname{};
			int n_nonhydrogen_neighbors{};
			std::vector<Node*> neighbors;
		};


		void addNode(int node_id, const std::string& atomname) 
		{ nodes.insert({ node_id, Node{node_id, atomname} }); }
		void connectNodes(int left_id, int right_id);

		std::unordered_map<int, Node> nodes;
	};

	class Chain {
		// First node is guaranteed to be non-hydrogen. Last node is not.
		std::vector<int> nodeids{};	// All atoms, TRUE ids
		int terminusIndex = -1;  // Index of the last non-hydrogen atom

	public:
		void append(const MoleculeGraph::Node& node, std::set<int>& ids_already_in_a_chain) {
			if (ids_already_in_a_chain.contains(node.atomid))
				throw std::runtime_error("Didn't expect this atom to already be in a chain");
			ids_already_in_a_chain.insert(node.atomid);

			nodeids.push_back(node.atomid);

			if (!node.isHydrogen()) terminusIndex = nodeids.size() - 1;
		}

		// Inserts a chain just before the current terminus
		void insert(const Chain& chain) {
			//nodeids.insert(nodeids.end() - 1, chain.nodeids.begin(), chain.nodeids.end());
			
			// Calculate the insertion position based on the terminus index
			auto insertPos = nodeids.begin() + terminusIndex;
			nodeids.insert(insertPos, chain.nodeids.begin(), chain.nodeids.end());

			// Update the terminus index to reflect the new position
			terminusIndex += chain.nodeids.size();
		}

		// Returns terminus non-hydrogen member of the chain
		int getBack() {
			if (terminusIndex != -1) {
				return nodeids[terminusIndex];
			}
			throw std::runtime_error("Terminus is not set");
		}

		const std::vector<int>& getNodeIds() const { return nodeids; }
		int getNNodesInChain() const { return nodeids.size(); }
		int parentchain_id = -1;
	};

	MoleculeGraph createGraph(const ParsedTopologyFile& topolfile);

	// Selects a random terminus carbon, and starts making chains from that point.
	// 1-atom functional groups are treated as separate chains. Particles in chains are in order
	// of how they appear in the chain, and thus the ids are not guaranteed to be sequential.
	// The order of the chains follows a DepthFirst like structure
	std::vector<Chain> getAllSubchains(const MoleculeGraph& molgraph);

	/// <summary>
	/// Reorders the particles so they are sequential inside of the subchains, and that particles of
	/// touching subchains are in order aswell
	/// </summary>
	/// <param name="gro_path">Reads and overwrites indices& ids of particles</param>
	/// <param name="top_path">Reads and overwrites indices& ids of particles + bonds</param>
	void reorderoleculeParticlesAccoringingToSubchains(fs::path gro_path, fs::path top_path);

};

