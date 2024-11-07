#pragma once

#include <vector>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <set>

#include "Filehandling.h"
#include <filesystem>
#include "MDFiles.h"

#include <ranges>
#include <iterator>
#include <unordered_set>
#include <queue>
#include <optional>
#include <type_traits>


namespace LimaMoleculeGraph {
	namespace fs = std::filesystem;

	struct MoleculeGraph {
		struct Node {
			Node() {}
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

			int atomid{-1};

		private:
			std::string atomname{};
			int n_nonhydrogen_neighbors{};
			std::vector<Node*> neighbors;
		};


		template<typename NodePtr>
		class BFSRange : public std::ranges::view_interface<BFSRange<NodePtr>> {
			using NodeType = std::remove_pointer_t<NodePtr>;
		public:
			explicit BFSRange(NodePtr start_node) {
				if (start_node) {
					node_queue.push(start_node);
					visited.insert(start_node->atomid);
				}
			}

			// Iterator class for BFS traversal
			class Iterator {
			public:
				Iterator() = default;

				// Constructor that initializes from a BFSRange instance
				explicit Iterator(BFSRange* range) : range(range), current(range->next_node()) {}

				NodeType& operator*() const { return *current; }

				Iterator& operator++() {
					current = range->next_node();
					return *this;
				}

				bool operator==(std::default_sentinel_t) const { return !current; }

			private:
				BFSRange* range = nullptr;
				NodePtr current = nullptr;
			};

			// Begin and end for range-based for-loop support
			Iterator begin() { return Iterator(this); }
			std::default_sentinel_t end() const { return std::default_sentinel; }

		private:
			std::unordered_set<int> visited;
			std::queue<NodePtr> node_queue;

			// Generates the next node in BFS order
			NodePtr next_node() {
				if (node_queue.empty()) return nullptr;

				NodePtr current = node_queue.front();
				node_queue.pop();

				for (Node* neighbor : current->getNeighbors()) {
					NodePtr neighbor_ptr = neighbor; // Assign Node* to NodePtr
					if (visited.insert(neighbor->atomid).second) {
						node_queue.push(neighbor_ptr);
					}
				}
				return current;
			}
		};












		// Create a (possibly disconnected) graph from a MolType
		MoleculeGraph(const TopologyFile::Moleculetype&);

		// Create a connected graph from a node, and all nodes it is connected to
		MoleculeGraph(const Node* root);



		void addNode(int node_id, const std::string& atomname) {

			auto a = Node(node_id, atomname);
			nodes.emplace(node_id, Node(node_id, atomname) );
		}
		void connectNodes(int left_id, int right_id);

		std::unordered_map<int, Node> nodes;

		auto BFS(int start_node_id) const {
			return BFSRange(&nodes.at(start_node_id));
		}

		auto BFS(int start_node_id) {
			return BFSRange(&nodes.at(start_node_id));
		}


		bool GraphIsDisconnected() const;

		// A moleculeggraph may be disconnected. This function returns all subgraphs that are connected.
		// This is not cheap..
		std::vector<MoleculeGraph> GetListOfConnectedGraphs() const;


	};

	//MoleculeGraph createGraph(const TopologyFile::Moleculetype&);

	void reorderoleculeParticlesAccoringingToSubchains(GroFile&, TopologyFile::Moleculetype&);

};

