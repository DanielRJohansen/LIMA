#pragma once

#include <array>
#include <filesystem>
#include <functional>
#include <iterator>
#include <map>
#include <optional>
#include <queue>
#include <ranges>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "Filehandling.h"
#include "MDFiles.h"

#include <cereal/access.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/vector.hpp>

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

		const std::vector<int>& GetChildIds(int parentId) const {
			return tree.at(parentId);
		}

		void ForSelfAndAllChildrenIds(int parentId, const std::function<void(int)>& visitor) const;
	};


	struct MoleculeGraph {
		struct Node {
			Node() {}
			Node(int id, const std::string& atomname) : atomid(id), atomname(atomname) {
				if (atomname == "") throw std::runtime_error("Cannot add noname atom to graph");
			};

			void addNeighbor(int neighborId, bool neighborIsHydrogen);

			bool isHydrogen() const { return atomname[0] == 'H'; }
			int getNNeighbors() const { return nNeighbors; }
			int getNNonhydrogenNeighbors() const { return n_nonhydrogen_neighbors; }
			const std::span<const int> getNeighbors() const { return std::span( neighborIds.data(), nNeighbors ); }

			bool isConnected(int id) const {
				for (const int neighborId : getNeighbors()) {
					if (neighborId == id)
						return true;
				}
				return false;
			}

			int atomid{ -1 };

		private:
			SmallString atomname{};
			
			static const int maxNeighbors = 8; // for easier serialization
			std::array<int, maxNeighbors> neighborIds;
			int nNeighbors=0;
			int n_nonhydrogen_neighbors{};


			// Serialization
			friend class cereal::access;
			template <class Archive>
			void serialize(Archive& archive) {
				archive(
					atomid,
					atomname,
					neighborIds,
					nNeighbors,
					n_nonhydrogen_neighbors
				);
			}

		};


		template<typename GraphType>
		class BFSRange : public std::ranges::view_interface<BFSRange<GraphType>> {
			using GraphNoRef = std::remove_reference_t<GraphType>;
			using NodeType = std::conditional_t<std::is_const_v<GraphNoRef>, const Node, Node>;

		public:
			BFSRange() = default;

			BFSRange(GraphType* graph, int startNodeId) : graph(graph) {
				if (graph && graph->nodes.contains(startNodeId)) {
					node_queue.push({ startNodeId, 0 });
					visited.insert(startNodeId);
				}
			}

			class Iterator {
			public:
				Iterator() = default;

				explicit Iterator(BFSRange* range) : range(range) {
					current = range->next_node();
				}

				NodeType& operator*() const {
					return range->graph->nodes.at(current->first);
				}

				int Depth() const {
					return current->second;
				}

				Iterator& operator++() {
					current = range->next_node();
					return *this;
				}

				bool operator==(std::default_sentinel_t) const {
					return !current.has_value();
				}

			private:
				BFSRange* range = nullptr;
				std::optional<std::pair<int, int>> current;
			};

			Iterator begin() { return Iterator(this); }
			std::default_sentinel_t end() const { return std::default_sentinel; }

		private:
			GraphType* graph = nullptr;
			std::unordered_set<int> visited;
			std::queue<std::pair<int, int>> node_queue;

			std::optional<std::pair<int, int>> next_node() {
				if (node_queue.empty())
					return std::nullopt;

				const auto [currentNodeId, currentDepth] = node_queue.front();
				node_queue.pop();

				const Node& currentNode = graph->nodes.at(currentNodeId);
				for (const int neighborId : currentNode.getNeighbors()) {
					if (!graph->nodes.contains(neighborId))
						continue;

					if (visited.insert(neighborId).second) {
						node_queue.push({ neighborId, currentDepth + 1 });
					}
				}

				return { { currentNodeId, currentDepth } };
			}
		};


		MoleculeGraph(const TopologyFile::Moleculetype&, std::optional<const std::unordered_set<int>> allowedIds = std::nullopt);

		MoleculeGraph(
			const std::vector<std::pair<int, std::string>>& atoms,
			const std::vector<std::array<int, 2>>& edges
		);

		MoleculeGraph() {};


		MoleculeTree ConstructMoleculeTree() const;

		std::unordered_map<int, int> ComputeNumDownstreamNodes(const MoleculeTree& moleculeTree) const;

		void connectNodes(int left_id, int right_id);

		std::map<int, Node> nodes;

		auto BFS(int start_node_id) const {
			return BFSRange<const MoleculeGraph>(this, start_node_id);
		}

		auto BFS(int start_node_id) {
			return BFSRange<MoleculeGraph>(this, start_node_id);
		}

		std::optional<int> DistanceBetweenNodes(int id0, int id1, int maxSearchDepth = 8) const;

		bool GraphIsDisconnected() const;

		std::vector<std::vector<int>> GetListOfListsofConnectedNodeids() const;
	};

	void reorderoleculeParticlesAccoringingToSubchains(GroFile&, TopologyFile::Moleculetype&);
};