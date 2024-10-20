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

			const int atomid{-1};

		private:
			const std::string atomname{};
			int n_nonhydrogen_neighbors{};
			std::vector<Node*> neighbors;
		};


		void addNode(int node_id, const std::string& atomname) { 
			nodes.insert({ node_id, Node{node_id, atomname} }); 
			highestNodeId = std::max(highestNodeId, node_id);
		}
		void connectNodes(int left_id, int right_id);

		std::unordered_map<int, Node> nodes;
		int highestNodeId = -1;









        struct BFSRange {
            struct Iterator {
                using iterator_category = std::input_iterator_tag;
                using value_type = Node;
                using difference_type = std::ptrdiff_t;
                using pointer = Node*;
                using reference = Node&;

                Iterator(MoleculeGraph* graph, int start_node_id)
                    : graph(graph), current(nullptr) {
                    if (graph && graph->nodes.find(start_node_id) != graph->nodes.end()) {
                        visited.insert(start_node_id);
                        node_queue.push(&graph->nodes[start_node_id]);
                        ++(*this); // Initialize to first valid node
                    }
                }

                reference operator*() const { return *current; }
                pointer operator->() const { return current; }

                Iterator& operator++() {
                    if (!node_queue.empty()) {
                        current = node_queue.front();
                        node_queue.pop();
                        for (const auto neighbor : current->getNeighbors()) {
                            if (visited.insert(neighbor->atomid).second) {
                                node_queue.push(neighbor);
                            }
                        }
                    }
                    else {
                        current = nullptr;
                    }
                    return *this;
                }

                Iterator operator++(int) {
                    Iterator temp = *this;
                    ++(*this);
                    return temp;
                }

                bool operator==(const Iterator& other) const {
                    return current == other.current;
                }

                bool operator!=(const Iterator& other) const {
                    return !(*this == other);
                }

            private:
                MoleculeGraph* graph;
                Node* current;
                std::unordered_set<int> visited;
                std::queue<Node*> node_queue;
            };

            BFSRange(MoleculeGraph* graph, int start_node_id)
                : graph(graph), start_node_id(start_node_id) {}

            Iterator begin() { return Iterator(graph, start_node_id); }
            Iterator end() { return Iterator(nullptr, -1); }

        private:
            MoleculeGraph* graph;
            int start_node_id;
        };

        BFSRange BFS(int start_node_id) {
            return BFSRange(this, start_node_id);
        }
	};

	MoleculeGraph createGraph(const TopologyFile::Moleculetype1&);

    void reorderoleculeParticlesAccoringingToSubchains(GroFile&, TopologyFile::Moleculetype1&);

};

