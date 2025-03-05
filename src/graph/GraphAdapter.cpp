//
// Created by Lei on 12/15/2023.
//

#include "GraphAdaptor.hpp"
#include "utils/Log.hpp"

#include <ranges>

#include <quick-cliques/DegeneracyTools.h>
#include <quick-cliques/MemoryManager.h>
#include <quick-cliques/TomitaAlgorithm.h>
#include <quick-cliques/AdjacencyListAlgorithm.h>
#include <quick-cliques/HybridAlgorithm.h>
#include <quick-cliques/DegeneracyAlgorithm.h>

NAMESPACE_BEGIN(PCO)
namespace graph {

	void SparseGraphAdaptor::initCliqueAlgorithm(const std::string& name) {
		/*if (!isAdapted) */transformData();

		if (quickClique == nullptr) {
			if (name == "tomita") {
				for (int i = 0; i < _numNodes; i++) {
					adjacencyMatrix[i] = (char*)calloc(_numNodes == 0 ? 1 : _numNodes, sizeof(char));
					for (int const neighbor : adjacencyList[i]) {
						adjacencyMatrix[i][neighbor] = 1;
					}
				}
				quickClique = new TomitaAlgorithm(adjacencyMatrix, _numNodes);
			}
			else if (name == "adjlist") {
				std::vector<std::vector<int>> adjacencyArray;
				for (int i = 0; i < _numNodes; i++) {
					adjacencyArray[i].resize(adjacencyList[i].size());
					int j = 0;
					for (int const neighbor : adjacencyList[i]) {
						adjacencyArray[i][j++] = neighbor;
					}
				}
				quickClique = new AdjacencyListAlgorithm(adjacencyArray);
			}
			else if (name == "degeneracy") {
				quickClique = new DegeneracyAlgorithm(adjacencyList);
			}
			else if (name == "hybrid") {
				quickClique = new HybridAlgorithm(adjacencyList);
			}
			else {
				LOG::ERROR("ERROR: unrecognized algorithm name ", name);
				exit(1);
			}
		}
	}

	void SparseGraphAdaptor::transformData() {
		adjacencyList.clear();
		adjacencyList.resize(_numNodes);
		for (int k = 0; k < G.outerSize(); ++k) {
			for (SpGraph::InnerIterator it(G, k); it; ++it) {
				assert(it.row() < numNodes && it.row() > -1);
				assert(it.col() < numNodes && it.col() > -1);
				if (it.row() == it.col())
					fprintf(stderr, "[Graph] File %s at Line %d: Detected loop %td->%td\n", __FILE__, __LINE__,
						it.row(), it.col());
				//LOG::INFO("it.row: {}, it.col: {}", it.row(), it.col());
				adjacencyList[it.row()].emplace_back(it.col());
			}
		}
		isAdapted = true;
	}

	void SparseGraphAdaptor::_computeDegeneracy() {
		/*if (!isAdapted)*/ transformData();
		_degeneracy = computeDegeneracy(adjacencyList, _numNodes);
	}

	std::vector<std::list<int>> SparseGraphAdaptor::run_maximalCliques(const std::string& name) {
		initCliqueAlgorithm(name);
		std::vector<std::list<int>> maximalCliques;
		auto storeCliqueInList = [&maximalCliques](std::list<int> const& clique) {
			if (clique.size() >= 2) {
				maximalCliques.emplace_back(clique);
			}
			};
		quickClique->AddCallBack(storeCliqueInList);
		quickClique->Run();
		return maximalCliques;
	}

	std::vector<std::list<int>> SparseGraphAdaptor::run_maxCliques(const std::string& name) {
		initCliqueAlgorithm(name);

		int maxCliqueSize = 2;
		std::vector<std::list<int>> t_maxCliques;
		auto storeMaxCliqueInList = [&maxCliqueSize, &t_maxCliques](std::list<int> const& clique) {
			if (clique.size() >= maxCliqueSize) {
				maxCliqueSize = clique.size();
				t_maxCliques.emplace_back(clique);
			}
			};
		quickClique->AddCallBack(storeMaxCliqueInList);
		quickClique->Run();

		std::vector<std::list<int>> maxCliques;
		for (auto& clique : std::ranges::reverse_view(t_maxCliques)) {
			if (clique.size() < maxCliqueSize) break;
			maxCliques.emplace_back(clique);
		}
		return maxCliques;
	}

} // namespace graph
NAMESPACE_END(PCO)
