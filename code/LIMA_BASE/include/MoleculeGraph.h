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

		std::vector<std::unique_ptr<Chain>> subchains; 

		int backbone_len = 0;
		int height = 0;	// Furthest len to a leaf in any subchain

	public:
		// Dont do this after using other append
		void append(const MoleculeGraph::Node& node, std::set<int>& ids_already_in_a_chain);
		void append(std::unique_ptr<Chain> chain);
		int getHeight() const { return height; }

		const std::vector<std::unique_ptr<Chain>>& getSubchains() const { return subchains; }

		void sort();



		const std::vector<int>& getNodeIds() const { return nodeids; }
	};

	MoleculeGraph createGraph(const ParsedTopologyFile& topolfile);

	/// <summary>
	/// Reorders the particles so they are sequential inside of the subchains, and that particles of
	/// touching subchains are in order aswell
	/// TODO: This function should throw if trying to work with any files containing ; Residue in the itp
	/// file, as that structure may not be kept.
	/// </summary>
	void reorderoleculeParticlesAccoringingToSubchains(const fs::path& gro_path_in, const fs::path& top_path_in, 
		const fs::path& gro_path_out, const fs::path& top_path_out);

};

