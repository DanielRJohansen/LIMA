#pragma once

#include <vector>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <set>

#include "Filehandling.h"
#include <filesystem>
namespace LimaMoleculeGraph {
	namespace fs = std::filesystem;

	struct MoleculeGraph {
		struct Node {
			Node() {}	// Not sure why this is needed to compile, but dont use!
			Node(int id, const std::string& atomname) : atomid(id), atomname(atomname) { if (atomname == "") throw std::runtime_error("Cannot add noname atom to graph"); };
			void addNeighbor(Node* neighbor);

			bool isHydrogen() const { return atomname[0] == 'H'; }
			int getNNeighbors() const { return neighbors.size(); }
			int getNNonhydrogenNeighbors() const { return n_nonhydrogen_neighbors; }
			const std::vector<Node*> getNeighbors() const { return neighbors; }


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
		std::vector<int> nodeids;

	public:
		void insert(int id, std::set<int>& ids_already_in_a_chain) {
			if (ids_already_in_a_chain.contains(id))
				throw std::runtime_error("Didn't expect this atom to already be in a chain");
			ids_already_in_a_chain.insert(id);

			nodeids.push_back(id);
		}
	};

	MoleculeGraph createGraph(const SimpleParsedFile& topolfile);
	std::vector<Chain> getAllSubchains(const MoleculeGraph& molgraph);

	/// <summary>
	/// Reorders the particles so they are sequential inside of the subchains, and that particles of
	/// touching subchains are in order aswell
	/// </summary>
	/// <param name="gro_path">Reads and overwrites indices& ids of particles</param>
	/// <param name="top_path">Reads and overwrites indices& ids of particles + bonds</param>
	void reorderoleculeParticlesAccoringingToSubchains(fs::path gro_path, fs::path top_path);

};

