//
// Created by Lei on 12/15/2023.
//

#ifndef LINEARSHARPOFFSET_GRAPHHANDLER_HPP
#define LINEARSHARPOFFSET_GRAPHHANDLER_HPP

#include "Graph.hpp"
#include "detail/Math.hpp"
#include <algorithm>

NAMESPACE_BEGIN(PCO)
namespace graph {

	class SparseGraphHandler : public SparseGraph {
	public:
		SparseGraphHandler() noexcept = default;

		explicit SparseGraphHandler(SparseGraph::SpGraph const& _G) : SparseGraph(_G) {}

		~SparseGraphHandler() noexcept override = default;

	protected:
		//////////////////////
		//!  Cliques Data  !//
		//////////////////////
		int _maxCliqueSize = 0;
		std::vector<std::list<int>> _maximalCliques;
		std::vector<std::list<int>> _maxCliques;

	public:
		///////////////////////
		//! Basic Operation !//
		///////////////////////
		size_t numNodes() { return _numNodes; }

		size_t numUEdges() { return _u_numEdges; }

		size_t numEdges() { return _numEdges; }

		size_t degeneracy() { return _degeneracy; }

		size_t maxDegree() { return _max_degree; }

		std::vector<size_t> nodeDegree() { return _nodeDegree; }

		// TODO
		/*void insertNode(int v) {
			if (v < 0) {
				fprintf(stderr, "[Graph] File %s at Line %d: Insert node %d is INVALID!\n", __FILE__, __LINE__, v);
				return;
			}
			this->_insertNode(v);
		}*/
		void insertNode(int v) = delete;

		// TODO
		void insertEdge(int v1, int v2) = delete;

		std::vector<int> softDeleteNode(int v) {
			if (v < 0 || v >= _numNodes) {
				fprintf(stderr, "[Graph] File %s at Line %d: Delete node %d is INVALID!\n", __FILE__, __LINE__, v);
				return std::vector<int>();
			}
			return this->_soft_deleteNode(v);
		}

		std::unordered_set<int>
			softDeleteMultiNodes(const std::unordered_set<int>& nodes) { return this->_soft_deleteMultiNodes(nodes); }

		// TODO
		void hardDeleteNode(int v) = delete;

		void deleteEdge(int v1, int v2) {
			if (v1 < 0 || v1 >= _numNodes ||
				v2 < 0 || v2 >= _numNodes) {
				fprintf(stderr, "[Graph] File %s at Line %d: Delete edge(%d, %d) is INVALID!\n", __FILE__, __LINE__,
					v1, v2);
				return;
			}
			this->_deleteEdge(v1, v2);
		}

		void updateEdges() { this->_updateEdges(); }

		void traverseGraph() { this->_traverseGraph(); }

		std::vector<int> getAdjNodes(int v) { return this->_getAdjNodes(v); }

		SpGraph
			getInducedSubgraph(int v, std::unordered_map<int, int>& subGToG);

	public:
		/**
		 * Compute degree for node v
		 */
		size_t degree(int v) { return this->_degree(v); }

		/**
		 * Compute degree for each node and update maximum degrees
		 */
		void computeAllDegree() { this->_computeAllDegree(); }

		/**
		 * Compute sorted degree for each node(ascending order by default)
		 */
		std::vector<IndexedValue<int>> computeSortedDegree(bool asc = true);

	public:
		/* Degeneracy */
		// TODO
		// void _computeDegeneracy() override;

		int getDegeneracy() {
			_computeDegeneracy();
			return this->_degeneracy;
		}

	public:
		/**
		 * Compute maximal clique(s)
		 * @param Algorithm name
		 */
		virtual std::vector<std::list<int>> run_maximalCliques(const std::string& name/* = "degeneracy"*/) = 0;

		/**
		 * Compute maximum clique(s)
		 * @param Algorithm name
		 */
		virtual std::vector<std::list<int>> run_maxCliques(const std::string& name/* = "degeneracy"*/) = 0;

		/**
		 * Get maximal clique(s)
		 * @return maximal clique(s)
		 */
		[[nodiscard]] std::vector<std::list<int>>
			maximalCliques() const { return _maximalCliques; }

		/**
		 * Get the size of largest(maximum) clique
		 * @return size of largest(maximum) clique
		 */
		[[nodiscard]] int maxCliqueSize() const { return _maxCliqueSize; }

		/**
		 * Get largest(maximum) clique(s)
		 * @return largest(maximum) clique(s)
		 */
		[[nodiscard]] std::vector<std::list<int>> maxCliques() const { return _maxCliques; }

		void outMaximalCliques(std::ostream& os = std::cout) const;

	public:
		void saveGraph(const std::string& save_file) { this->_saveGraph(save_file); }
	};

} // namespace graph
NAMESPACE_END(PCO)

#endif //LINEARSHARPOFFSET_GRAPHHANDLER_HPP
