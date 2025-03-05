//
// Created by Lei on 12/15/2023.
//

#include "Graph.hpp"
#include "utils/LOG.hpp"
#include "utils/File.hpp"

#include <unsupported/Eigen/SparseExtra>

NAMESPACE_BEGIN(offset)
    namespace graph {

        std::vector<int> SparseGraph::_getAdjNodes(int v) {
            std::vector<int> adjNodes;
            // 使用列优先遍历 v 这一列
            for (SpGraph::InnerIterator it(G, v); it; ++it) {
                adjNodes.emplace_back(it.row());
            }
            return adjNodes;
        }

        std::vector<int> SparseGraph::_soft_deleteNode(int v) {
            std::vector<int> adjNodes = _getAdjNodes(v);

            // 先将与顶点 v 相关联的行和列对应值设为 .0
//            for (int k = G.outerIndexPtr()[v]; k < G.outerIndexPtr()[v + 1]; ++k) {
//                int colIndex = G.innerIndexPtr()[k];
//                G.coeffRef(v, colIndex) = .0;
//            }
            for (SpGraph::InnerIterator it(G, v); it; ++it) {
                G.coeffRef(it.row(), v) = .0;
                G.coeffRef(v, it.row()) = .0; // 对称矩阵
            }

            // 再 prune
            /// Be aware that this will trigger a costly memory copy of the remaining column entries.
            /// With sparse matrices, we usually never work in-place.
            G.prune(.0);
            G.makeCompressed();

            _updateEdges();

            // 更新与节点 v 相邻节点的度数
            if (!_nodeDegree.empty()) {
                _nodeDegree.at(v) = 0; // 并不直接删除
                for (int adjNode: adjNodes)
                    _nodeDegree.at(adjNode)--;
            }

            return adjNodes;
        }

        std::unordered_set<int> SparseGraph::_soft_deleteMultiNodes(const std::unordered_set<int> &nodes) {
            std::unordered_set<int> adjToNodes; // 与这些节点相邻的节点
            std::vector<int> _adjToNodes; // 与这些节点相邻的节点
            for (int v: nodes) {
                if (v < 0 || v >= _numNodes) {
                    fprintf(stderr, "[Graph] File %s at Line %d: Delete nodes %d is INVALID!\n", __FILE__, __LINE__, v);
                    exit(1);
                }
                for (SpGraph::InnerIterator it(G, v); it; ++it) {
                    if (it.row() >= G.rows()) break;
                    if (!nodes.count(it.row())) {
                        _adjToNodes.emplace_back(it.row());
                        adjToNodes.insert(it.row());
                    }
                }
                /*for (int k = 0; k < G.outerSize(); ++k)
                {
                    for (SpGraph::InnerIterator it(G, k); it; ++it) {
                        if (it.row() == v) {
                            _adjToNodes.emplace_back(it.row());
                            adjToNodes.insert(it.row());
                        }
                    }
                }*/
            }
            if (!_nodeDegree.empty()) {
                for (int v: _adjToNodes)
                    _nodeDegree.at(v)--;
            }


            for (int v: nodes) {
                //    // 先将与顶点 v 相关联的行和列对应值设为 .0
                //std::cout << "v = " << v;
                //std::cout << " G.rows = " << G.rows() << std::endl;
                //for (int k = G.outerIndexPtr()[v]; k < G.outerIndexPtr()[v + 1]; ++k) {
                //    int colIndex = G.innerIndexPtr()[k];
                //    //std::cout << "colIndex = " << colIndex << std::endl;
                //    G.coeffRef(v, colIndex) = .0;
                //}
                //    std::cout << "456!\n";
                for (SpGraph::InnerIterator it(G, v); it; ++it) {
                    if (it.row() >= G.rows()) break;
                    if (it.col() >= G.outerSize()) break;
                    G.coeffRef(it.row(), it.col()) = .0;
                    try {
                        G.coeffRef(it.col(), it.row()) = .0; // 对称矩阵
                    }
                    catch (const std::exception &e) {
                        // auto file_logger = spdlog::basic_logger_st("file_logger", LOG_FILE);
                        // file_logger->set_level(spdlog::level::err);
                        // file_logger->error("[GRAPH] Error in graph delete! Caught exception with message: {}",
                        // 	e.what());
                        // std::cerr << "error in graph: " << e.what() << std::endl;
                    }
                }

                //for (int k = 0; k < G.outerSize(); ++k)
                //{
                //	for (SpGraph::InnerIterator it(G, k); it; ++it) {
                //		if (it.row() == v || it.col() == v) {
                //			G.coeffRef(it.row(), it.col()) = .0;
                //			//G.coeffRef(it.col(), it.row()) = .0; // 对称矩阵
                //		}
                //	}
                //}

                // 再 prune
                /// Be aware that this will trigger a costly memory copy of the remaining column entries.
                /// With sparse matrices, we usually never work in-place.
                G.prune(.0);
                G.makeCompressed();

                // std::cout << "789!\n";
            }


            _updateEdges();

            return adjToNodes;
        }

        //        void SparseGraph::_soft_deleteSubGraph(const SpGraph &subG) {
        //
        //        }

        void SparseGraph::_deleteEdge(int v1, int v2) {
            G.coeffRef(v1, v2) = .0;
            G.coeffRef(v2, v1) = .0;

            G.prune(.0);
            G.makeCompressed();

            _updateEdges();
        }

        size_t SparseGraph::_degree(int v) {
            if (!_nodeDegree.empty()) return _nodeDegree.at(v);

            size_t degree = 0;
            for (SpGraph::InnerIterator it(G, v); it; ++it) {
                ++degree;
            }
            return degree;
        }

        void SparseGraph::_computeAllDegree() {
            _nodeDegree.resize(_numNodes, 0);
            for (int k = 0; k < G.outerSize(); ++k) {
                for (SpGraph::InnerIterator it(G, k); it; ++it) {
                    _nodeDegree[it.col()]++;
                }
            }
        }

        void SparseGraph::_saveGraph(const std::string &save_file) {
            utils::file::checkDir(save_file);
            Eigen::saveMarket(G, save_file);
        }
    }
NAMESPACE_END(offset)