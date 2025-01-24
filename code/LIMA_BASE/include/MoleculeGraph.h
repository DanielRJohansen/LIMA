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
#include <map>


namespace LimaMoleculeGraph {
	namespace fs = std::filesystem;

	class MoleculeTree {
		std::unordered_map<int, std::vector<int>> tree;
	public:
		MoleculeTree(int rootId);
		void AddChild(int parentId, int childId) {
			tree.insert({ childId, {} });
			tree.at(parentId).emplace_back(childId);
		}

		// Get immediate children
		const std::vector<int>& GetChildIds(int parentId) const {
			return tree.at(parentId);
		}

		// Get all children recursively
		std::vector<int> GetAllChildIdsAndSelf(int parentId, std::unordered_map<int,int> nodeIdsNumDownstreamNodes) const;
	};



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
		MoleculeGraph(const TopologyFile::Moleculetype&, std::optional<const std::unordered_set<int>> allowedIds = std::nullopt);

		// Create a graph with 0-indexed consequtive nodes
		MoleculeGraph(const std::vector<std::pair<int, std::string>>& atoms, 
			const std::vector<std::array<int, 2>>& edges);
		MoleculeGraph() {};


		// Create a connected graph from a node, and all nodes it is connected to
		//MoleculeGraph(const Node* root);

		// Uses a BFS approach to construct a molecule without cycles
		MoleculeTree ConstructMoleculeTree() const;

		std::unordered_map<int, int> ComputeNumDownstreamNodes(const MoleculeTree& moleculeTree) const;

	/*	void addNode(int node_id, const std::string& atomname) {
			nodes.emplace(node_id, Node(node_id, atomname) );
		}*/
		void connectNodes(int left_id, int right_id);

		std::map<int, Node> nodes;
		Node* root = nullptr;

		auto BFS(int start_node_id) const {
			return BFSRange(&nodes.at(start_node_id));
		}

		auto BFS(int start_node_id) {
			return BFSRange(&nodes.at(start_node_id));
		}


		bool GraphIsDisconnected() const;

		// A moleculeggraph may be disconnected. This function returns all subgraphs that are connected.
		// This is not cheap..
		//std::vector<MoleculeGraph> GetListOfConnectedGraphs() const;
		std::vector<std::vector<int>> GetListOfListsofConnectedNodeids() const;


	};

	//MoleculeGraph createGraph(const TopologyFile::Moleculetype&);

	void reorderoleculeParticlesAccoringingToSubchains(GroFile&, TopologyFile::Moleculetype&);

};

