//
// Created by Lei on 12/15/2023.
//

#include "GraphHandler.hpp"
#include <omp.h>

NAMESPACE_BEGIN(PCO)
    namespace graph {

        ///////////////////////
        //! Basic Operation !//
        ///////////////////////
        typename SparseGraph::SpGraph
        SparseGraphHandler::getInducedSubgraph(int v, std::unordered_map<int, int> &subGToG) {
            std::vector<SparseGraph::SpGraphTrip> subGTriplet;
            subGTriplet.reserve(G.rows());

            // 图 G 中节点索引到图 subG 节点索引的映射
            std::unordered_map<int, int> GToSubG;
            int subGIdx = 0;
            GToSubG[v] = subGIdx++;
            std::unordered_set<int> neighborSet;

            // 得到节点v的度数
            size_t v_degree = 0;
            if (!_nodeDegree.empty()) v_degree = _nodeDegree.at(v);
            else v_degree = this->degree(v);

            // 遍历节点v的邻居节点
            for (SpGraph::InnerIterator it(G, v); it; ++it) {
                int neighbor = it.row();
                GToSubG[neighbor] = subGIdx;
                subGToG[subGIdx] = neighbor;

                Scalar value = it.value();
                neighborSet.insert(neighbor);

                // 添加边(v, neighbor)和(neighbor, v)到诱导子图
                subGTriplet.emplace_back(0, subGIdx, value);
                subGTriplet.emplace_back(subGIdx, 0, value);
                ++subGIdx;
            }

            // v的诱导子图等于 节点v的度数 + 1
            int numNodes = v_degree + 1;
            // 创建一个稀疏矩阵表示诱导子图
            SpGraph subG(numNodes, numNodes);

            // 添加与节点v相邻的节点之间的边
            for (const int adj_v1: neighborSet) {
                for (const int adj_v2: neighborSet) {
                    if (adj_v1 != adj_v2 &&
                        G.coeff(adj_v1, adj_v2) != .0) {
                        subGTriplet.emplace_back(GToSubG.at(adj_v1), GToSubG.at(adj_v2),
                                                 G.coeff(adj_v1, adj_v2));
                    }
                }
            }

            subG.setFromTriplets(subGTriplet.begin(), subGTriplet.end());
            subG.makeCompressed();

            return subG;
        }

        ///////////////////////
        //!  Degree
        ///////////////////////
        std::vector<IndexedValue<int>>
        SparseGraphHandler::computeSortedDegree(bool asc) {
            if (_nodeDegree.empty()) computeAllDegree();
            std::vector<IndexedValue<int>> sortedDegree(_numNodes);
#pragma omp parallel for
            for (int i = 0; i < _numNodes; ++i) {
                sortedDegree[i] = IndexedValue<int>(i, _nodeDegree[i]);
            }
            if (asc)
                std::sort(sortedDegree.begin(), sortedDegree.end(), IndexedValue<int>::asc_cmp);
            else
                std::sort(sortedDegree.begin(), sortedDegree.end(), IndexedValue<int>::dsc_cmp);
            return sortedDegree;
        }

        ///////////////////////
        //!  Maximal Cliques
        ///////////////////////
        void SparseGraphHandler::outMaximalCliques(std::ostream &os) const {
            os << "[Graph] maxCliqueSize = " << _maxCliqueSize << std::endl;
            for (const auto &_maximalClique: _maximalCliques) {
                for (const auto &node: _maximalClique) {
                    os << node << " ";
                }
                os << "\n==========\n";
            }
        }

    } // namespace graph
NAMESPACE_END(PCO)