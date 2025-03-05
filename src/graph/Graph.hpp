//
// Created by Lei on 12/15/2023.
//

#ifndef LINEARSHARPOFFSET_GRAPH_HPP
#define LINEARSHARPOFFSET_GRAPH_HPP

#include "Config.hpp"
#include "detail/Math.hpp"
#include "detail/BasicDataType.hpp"
#include <iostream>
#include <utility>

NAMESPACE_BEGIN(PCO)
    namespace graph {

        ////////////////////////
        //! Undirected Graph !//
        ////////////////////////
        // TODO: template
        class Graph {
        protected:
            /* Basic properties of a graph */
            size_t _numNodes; // 节点数
            size_t _u_numEdges; // 边的数量 *2
            size_t _numEdges; // 边的数量

        protected:
            /* Advanced properties of a graph */
            size_t _degeneracy;
            size_t _max_degree;
            std::vector<size_t> _nodeDegree;

        public:
            Graph() noexcept = default;

            virtual ~Graph() noexcept = default;

        protected:
            /**
             * 遍历整个图
             */
            virtual void _traverseGraph() = 0;

            /**
             * 获取一个节点的相邻节点
             * @param v
             * @return
             */
            virtual std::vector<int> _getAdjNodes(int v) = 0;

            /**
             * TODO: 在图中添加一个节点
             * @param v 待添加节点
             */
            // virtual void _insertNode(int v) = 0;

            /**
             * TODO: 在图中添加一条(无向)边
             * @param v1 待添加边的端点
             * @param v2 待添加边的端点
             */
            // virtual void _insertEdge(int v1, int v2) = 0;

            /**
             * 在图中删除一个节点及其相关联的边(但不真正在底层数据结构中删除, 即节点数并不会减少)
             * 并返回节点的邻居节点
             * @param v 待删除节点
             */
            virtual std::vector<int> _soft_deleteNode(int v) = 0;

            /**
             * 在图中删除多个节点及其相关联的边(但不真正在底层数据结构中删除, 即节点数并不会减少)
             * 并返回这些节点的邻居节点
             * @param v 待删除节点
             */
            virtual std::unordered_set<int> _soft_deleteMultiNodes(const std::unordered_set<int> &nodes) = 0;

            /**
             * TODO: 在图中删除一个节点及其相关联的边(但会在底层数据结构中删除, 即节点数会减少)
             * @param v 待删除节点
             */
            // virtual void _hard_deleteNode(int v) = 0;

            /**
             * 在图中删除一条(无向)边
             * @param v1 待删除边的端点
             * @param v2 待删除边的端点
             */
            virtual void _deleteEdge(int v1, int v2) = 0;

            virtual void _updateEdges() = 0;

        protected:
            virtual void _computeDegeneracy() = 0;

            // TODO
            // virtual void _computeDegeneracyOrder() = 0;

            virtual size_t _degree(int v) = 0;

            virtual void _computeAllDegree() = 0;

        protected:
            virtual void _saveGraph(const std::string &save_file) = 0;
        };

        class SparseGraph : public Graph {
        public:
            using SpGraph = SpMatrix;
            using SpGraphTrip = SpTriplet;

        public:
            SpGraph G;

        public:
            SparseGraph() noexcept = default;

            explicit SparseGraph(const SpGraph &_G) : G(_G) {
                _numNodes = G.rows();
                _u_numEdges = G.nonZeros();
                _numEdges = _u_numEdges / 2;
            }

            ~SparseGraph() noexcept override = default;

        protected:
            void _traverseGraph() override {
                traverseSparseMat(G);
            }

            std::vector<int> _getAdjNodes(int v) override;

            // TODO
            // void _insertNode(int v) override;

            // TODO
            // void _insertEdge(int v1, int v2) override;

            std::vector<int> _soft_deleteNode(int v) override;

            std::unordered_set<int> _soft_deleteMultiNodes(const std::unordered_set<int> &nodes) override;

            // TODO
            // void _hard_deleteNode(int v) override;

            void _deleteEdge(int v1, int v2) override;

            void _updateEdges() override {
//                _u_numEdges = G.nonZeros();
                _numEdges = G.nonZeros() / 2;
            }

        protected:
            size_t _degree(int v) override;

            void _computeAllDegree() override;

        protected:
            void _saveGraph(const std::string &save_file) override;

        public:
            /**
             * TODO: 在图中删除一个子图
             * @param v 待删除节点
             */
            std::vector<int> _soft_deleteSubGraph(const SpGraph &subG) = delete;
        };

        class DenseGraph : public Graph {
        public:
            using DeGraph = MatrixX;

        protected:
            DeGraph G;

        public:
            DenseGraph() noexcept = default;

            explicit DenseGraph(DeGraph _G) : G(std::move(_G)) {
                _numNodes = G.rows();
                _u_numEdges = G.nonZeros();
                _numEdges = _u_numEdges / 2;
            }

            ~DenseGraph() noexcept override = default;

        protected:

        };
    } // namespace graph

NAMESPACE_END(PCO)

#endif //LINEARSHARPOFFSET_GRAPH_HPP
