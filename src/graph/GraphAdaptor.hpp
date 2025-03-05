//
// Created by Lei on 12/15/2023.
//

#ifndef LINEARSHARPOFFSET_GRAPHADAPTOR_HPP
#define LINEARSHARPOFFSET_GRAPHADAPTOR_HPP

#include "GraphHandler.hpp"

#include <list>
#include <vector>

#include <quick-cliques/Algorithm.h>

NAMESPACE_BEGIN(PCO)
    namespace graph {
        using namespace QuickClique;

        class SparseGraphAdaptor : public SparseGraphHandler {
        private:
            bool isAdapted = false;

        private:
            /////////////////////////////////
            //!  Adapted data for Cliques
            /////////////////////////////////
            Algorithm *quickClique = nullptr;
            std::vector<std::list<int>> adjacencyList;
            char **adjacencyMatrix = nullptr;
            std::vector<std::vector<char>> vAdjacencyMatrix;

        public:
            SparseGraphAdaptor() noexcept = default;

            explicit SparseGraphAdaptor(const SparseGraph::SpGraph &_G) noexcept: SparseGraphHandler(_G) {}

            ~SparseGraphAdaptor() override {
                if (adjacencyMatrix != nullptr) {
                    int i = 0;
                    while (i < _numNodes) {
                        free(adjacencyMatrix[i]);
                        i++;
                    }
                    free(adjacencyMatrix);
                    adjacencyMatrix = nullptr;
                }
                if (quickClique != nullptr) {
                    delete quickClique;
                    quickClique = nullptr;
                }
            }

        private:
            void transformData();

        private:
            void initCliqueAlgorithm(const std::string &name);

        private:
            /////////////////////////////////
            //!  Degeneracy
            /////////////////////////////////
            void _computeDegeneracy() override;

        public:
            // 注意, 此处并没有拿到所有 maximal clique, 而是尽可能选择size比较大的clique
            std::vector<std::list<int>> run_maximalCliques(const std::string &name) override;

            std::vector<std::list<int>> run_maxCliques(const std::string &name) override;
        };

    } // namespace graph
NAMESPACE_END(PCO)

#endif //LINEARSHARPOFFSET_GRAPHADAPTOR_HPP
