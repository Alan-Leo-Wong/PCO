#ifndef PCO_OCTREE_HPP
#define PCO_OCTREE_HPP

#include "OctreeTet.hpp"
#include "utils/Donut.hpp"
#include <memory>

NAMESPACE_BEGIN(PCO)
    namespace geometry {

        struct OctreeNode {
        public:
            OctreeNode() {
                for (int i = 0; i < 27; ++i) neighborsIdx[i] = -1;
            }

            OctreeNode(int _depth, const AABox<Vector3> &_bbox) : depth(_depth), bbox(_bbox), width(bbox.boxWidth) {
                center = bbox.center();
                circumRadius = 0.5 * width.norm();
                donut::Loop<int, 27>([&](int i) {
                    neighborsIdx[i] = -1;
                });
            }

        public:
            int depth = -1;
            int childIdx = 0;
            bool is_leaf = true;
            Vector3i xyzIdx = Vector3i(0, 0, 0);
            uint32_t mortonCode = -1;

            int globalIdx = -1;
            unsigned int neighborsIdx[27];
            unsigned int childsIdx[8];

            MatrixX corners;
            Vector3 center;
            double circumRadius;

            Vector3 width;
            AABox<Vector3> bbox;

            double centerDis;

            std::shared_ptr<OctreeNode> child[8];
            std::shared_ptr<OctreeNode> parent;

        public:
            std::array<OctreeTet, 5> tets;

        public:
            void setCorners() {
                corners.resize(8, 3);
                donut::Loop<int, 8>([&](int i) {
                    corners.row(i) = bbox.corner(i);
                });
            }

            void setLocalTets() {
                int type = tetsTypeInNodeLUT[childIdx];
                for (int i = 0; i < 5; ++i) {
                    tets[i].type = type;
                    for (int j = 0; j < 4; ++j) {
                        tets[i].coord.row(j) = bbox.corner(tetVertexInNodeLUT[type][i][j]);
                        tets[i].vertsInLocalNode[j] = tetVertexInNodeLUT[type][i][j];
                        tets[i].localFaces[j] = tetLocalFaceVertexLUT[type][i][j];
                        for (int k = 0; k < 3; ++k)
                            tets[i].facesInLocalNode[j][k] = tetTriInNodeLUT[type][i][j][k];

                        tets[i].localFaceNeighbors[j] = tetLocalFaceNeighborsLUT[type][i][j];
                    }

                    for (int j = 0; j < 6; ++j) {
                        const TetEdge &edge = OctreeTet::localEdges[j];
                        int x = edge[0];
                        int y = edge[1];
                        tets[i].edgesInLocalNode[j] = {tetVertexInNodeLUT[type][i][x], tetVertexInNodeLUT[type][i][y]};

                        for (int k = 0; k < 6; ++k)
                            tets[i].localEdgeNeighbors[j][k] = tetLocalEdgeNeighborsLUT[type][i][j][k];
                    }
                }
            }

            void setGlobalTets(const std::unordered_map<int, uint32_t> &localToGlobal) {
                int type = tetsTypeInNodeLUT[childIdx];
                for (int i = 0; i < 5; ++i) {
                    tets[i].type = type;
                    for (int j = 0; j < 4; ++j) {
                        tets[i].coord.row(j) = bbox.corner(tetVertexInNodeLUT[type][i][j]);
                        tets[i].globalVerts[j] = localToGlobal.at(tetVertexInNodeLUT[type][i][j]);
                        for (int k = 0; k < 3; ++k)
                            tets[i].globalFaces[j][k] = localToGlobal.at(tetTriInNodeLUT[type][i][j][k]);
                    }

                    for (int j = 0; j < 6; ++j) {
                        const TetEdge &edge = OctreeTet::localEdges[j];
                        int x = edge[0];
                        int y = edge[1];
                        tets[i].globalEdges[j] = {localToGlobal.at(tetVertexInNodeLUT[type][i][x]),
                                                  localToGlobal.at(tetVertexInNodeLUT[type][i][y])};
                    }
                }
            }

            void setLocalTetNeighbors() {
                int type = tetsTypeInNodeLUT[childIdx];
                for (int i = 0; i < 5; ++i) {
                    for (int j = 0; j < 6; ++j)
                        for (int k = 0; k < 6; ++k)
                            tets[i].localEdgeNeighbors[j][k] = tetLocalEdgeNeighborsLUT[type][i][j][k];

                    for (int j = 0; j < 4; ++j)
                        tets[i].localFaceNeighbors[j] = tetLocalFaceNeighborsLUT[type][i][j];
                }
            }

            OctreeTet &getTet(int i) {
                return tets[i];
            }

            OctreeTet getTet(int i) const {
                return tets[i];
            }

            TetCoord getTetCornersCoord(int i) const {
                TetCoord corners(4, 3);

                corners.row(0) = bbox.corner(tets[i].vertsInLocalNode[0]);
                corners.row(1) = bbox.corner(tets[i].vertsInLocalNode[1]);
                corners.row(2) = bbox.corner(tets[i].vertsInLocalNode[2]);
                corners.row(3) = bbox.corner(tets[i].vertsInLocalNode[3]);

                return corners;
            }

            Vector3 getTetCornerCoord(int i, int j) const {
                return bbox.corner(tets[i].vertsInLocalNode[j]);
            }

            TetEdgeNeighbor getTetEdgeLocalNeighbors(int i, int j) {
                assert(i >= 0 && i < 5 && j >= 0 && j < 6);
                return tets[i].localEdgeNeighbors[j];
            }

            TetEdgeGlobalNeighbor getTetEdgeGlobalNeighbors(int i, int j) {
                assert(i >= 0 && i < 5 && j >= 0 && j < 4);
                return tets[i].globalEdgeNeighbors[j];
            }

            TetFaceNeighbor getTetFaceLocalNeighbors(int i, int j) {
                assert(i >= 0 && i < 5 && j >= 0 && j < 6);
                return tets[i].localFaceNeighbors[j];
            }

            TetFaceGlobalNeighbor getTetFaceGlobalNeighbors(int i, int j) {
                assert(i >= 0 && i < 5 && j >= 0 && j < 4);
                return tets[i].globalFaceNeighbors[j];
            }

            void outputTet(std::ofstream &out, int tetIdx, int &vertBegIdx, int nodeIdx = -1) {
                if (tetIdx < 0 || tetIdx >= 5) return;
                if (vertBegIdx <= 0) vertBegIdx = 1;
                if (nodeIdx != -1) out << "# nodeIdx: " << nodeIdx << std::endl;

                int type = tetsTypeInNodeLUT[childIdx];
//                for (int j = 0; j < 8; ++j)
//                    out << "v " << bbox.corner(j).transpose() << std::endl;
//
//                {
//                    out << "# tetIdx: " << tetIdx << std::endl;
//                    for (int j = 0; j < 4; ++j)
//                        out << "f " << tetTriInNodeLUT[type][tetIdx][j][0] + vertBegIdx << " " <<
//                            tetTriInNodeLUT[type][tetIdx][j][1] + vertBegIdx << " " <<
//                            tetTriInNodeLUT[type][tetIdx][j][2] + vertBegIdx << std::endl;
//                }
//
//                vertBegIdx += 8;

                std::array<Vector3, 5> colors;
                colors[0] = Vector3(0.541, 0.725, 0.710);
                colors[1] = Vector3(0.956, 0.651, 0.596);
                colors[2] = Vector3(0.949, 0.776, 0.627);
                colors[3] = Vector3(0.722, 0.663, 0.788);
                colors[4] = Vector3(0.764, 0.741, 0.722);

                Eigen::RowVector3d color = colors[tetIdx];
                std::unordered_map<int, int> mp;
                for (int j = 0; j < 4; ++j) {
                    if (!mp.contains(tetTriInNodeLUT[type][tetIdx][j][0])) {
                        out << "v " << bbox.corner(tetTriInNodeLUT[type][tetIdx][j][0]).transpose() << " " << color
                            << std::endl;
                        mp[tetTriInNodeLUT[type][tetIdx][j][0]] = vertBegIdx++;
                    }
                    if (!mp.contains(tetTriInNodeLUT[type][tetIdx][j][1])) {
                        out << "v " << bbox.corner(tetTriInNodeLUT[type][tetIdx][j][1]).transpose() << " " << color
                            << std::endl;
                        mp[tetTriInNodeLUT[type][tetIdx][j][1]] = vertBegIdx++;
                    }
                    if (!mp.contains(tetTriInNodeLUT[type][tetIdx][j][2])) {
                        out << "v " << bbox.corner(tetTriInNodeLUT[type][tetIdx][j][2]).transpose() << " " << color
                            << std::endl;
                        mp[tetTriInNodeLUT[type][tetIdx][j][2]] = vertBegIdx++;
                    }
                    out << "f " << mp[tetTriInNodeLUT[type][tetIdx][j][0]] << " " <<
                        mp[tetTriInNodeLUT[type][tetIdx][j][1]] << " " <<
                        mp[tetTriInNodeLUT[type][tetIdx][j][2]] << std::endl;
                }
            }

            void outputTets(std::ofstream &out, int &vertBegIdx, int nodeIdx = -1) {
                if (vertBegIdx <= 0) vertBegIdx = 1;
                if (nodeIdx != -1) out << "# nodeIdx: " << nodeIdx << std::endl;

                int type = tetsTypeInNodeLUT[childIdx];
                //std::cout << corners << std::endl;
                out << "# tet type: " << type << std::endl;
//                for (int j = 0; j < 8; ++j)
//                    out << "v " << bbox.corner(j).transpose() << std::endl;

                std::array<Vector3, 5> colors;
                colors[0] = Vector3(0.541, 0.725, 0.710);
                colors[1] = Vector3(0.956, 0.651, 0.596);
                colors[2] = Vector3(0.949, 0.776, 0.627);
                colors[3] = Vector3(0.722, 0.663, 0.788);
                colors[4] = Vector3(0.764, 0.741, 0.722);
//                for (int i = 0; i < 5; ++i) {
//                    //if (i != 0 && i != 1) continue;
//                    out << "# tetIdx: " << i << std::endl;
//                    for (int j = 0; j < 4; ++j)
//                        out << "f " << tetTriInNodeLUT[type][i][j][0] + vertBegIdx << " " <<
//                            tetTriInNodeLUT[type][i][j][1] + vertBegIdx << " " <<
//                            tetTriInNodeLUT[type][i][j][2] + vertBegIdx << std::endl;
//                }
//                vertBegIdx += 8;

                for (int i = 0; i < 5; ++i) {
                    //if (i != 0 && i != 1) continue;
                    out << "# tetIdx: " << i << std::endl;
                    Eigen::RowVector3d color = colors[i];
                    std::unordered_map<int, int> mp;
                    for (int j = 0; j < 4; ++j) {
                        if (!mp.contains(tetTriInNodeLUT[type][i][j][0])) {
                            out << "v " << bbox.corner(tetTriInNodeLUT[type][i][j][0]).transpose() << " " << color
                                << std::endl;
                            mp[tetTriInNodeLUT[type][i][j][0]] = vertBegIdx++;
                        }
                        if (!mp.contains(tetTriInNodeLUT[type][i][j][1])) {
                            out << "v " << bbox.corner(tetTriInNodeLUT[type][i][j][1]).transpose() << " " << color
                                << std::endl;
                            mp[tetTriInNodeLUT[type][i][j][1]] = vertBegIdx++;
                        }
                        if (!mp.contains(tetTriInNodeLUT[type][i][j][2])) {
                            out << "v " << bbox.corner(tetTriInNodeLUT[type][i][j][2]).transpose() << " " << color
                                << std::endl;
                            mp[tetTriInNodeLUT[type][i][j][2]] = vertBegIdx++;
                        }
                        out << "f " << mp[tetTriInNodeLUT[type][i][j][0]] << " " <<
                            mp[tetTriInNodeLUT[type][i][j][1]] << " " <<
                            mp[tetTriInNodeLUT[type][i][j][2]] << std::endl;
                    }
                }
            }

        };

        struct OctreeVertex {
            uint32_t globalIdx = -1;
            uint32_t owner = -1;
            Vector3 coord;

            OctreeVertex() = default;

            OctreeVertex(uint32_t _globalIdx, uint32_t _owner, const Vector3 &_coord) : globalIdx(_globalIdx),
                                                                                        owner(_owner),
                                                                                        coord(_coord) {}
        };

        using EdgeOrder = std::pair<uint32_t, uint32_t>;

        struct OctreeEdge {
            uint32_t globalIdx = -1;
            uint32_t owner = -1;
            EdgeOrder vertGlobalIdx = {-1, -1}; // order isn't the matter
            OctreeEdge() = default;

            OctreeEdge(uint32_t _globalIdx, uint32_t _owner, const EdgeOrder &_vertGlobalIdx) : globalIdx(
                    _globalIdx), owner(_owner), vertGlobalIdx(_vertGlobalIdx) {}
        };

        using FaceOrder = std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>;

        struct OctreeFace {
            uint32_t globalIdx = -1;
            uint32_t owner = -1;
            FaceOrder vertGlobalIdx = {-1, -1, -1, -1}; // right order
            OctreeFace() = default;

            OctreeFace(uint32_t _globalIdx, uint32_t _owner, const FaceOrder &_vertGlobalIdx) : globalIdx(
                    _globalIdx), owner(_owner), vertGlobalIdx(_vertGlobalIdx) {}
        };

        struct OctreeTetVertex {
            uint32_t globalIdx = -1;
            uint32_t owner = -1; // ref to tet global idx
            Vector3 coord;

            OctreeTetVertex() = default;

            OctreeTetVertex(uint32_t _globalIdx, uint32_t _owner, const Vector3 &_coord) : globalIdx(_globalIdx),
                                                                                           owner(_owner),
                                                                                           coord(_coord) {}
        };

        using TetEdgeOrder = std::array<uint32_t, 2>;

        struct OctreeTetEdge {
            uint32_t globalIdx = -1;
            uint32_t owner = -1; // ref to tet global idx
            TetEdgeOrder edge = {-1u, -1u}; // order isn't the matter
            OctreeTetEdge() = default;

            OctreeTetEdge(uint32_t _globalIdx, uint32_t _owner, const TetEdgeOrder &_vertGlobalIdx) : globalIdx(
                    _globalIdx), owner(_owner), edge(_vertGlobalIdx) {}
        };

        using TetFaceOrder = std::array<uint32_t, 3>;

        struct OctreeTetFace {
            uint32_t globalIdx = -1;
            uint32_t owner = -1; // ref to tet global idx
            TetFaceOrder face = {-1u, -1u, -1u}; // right order
            OctreeTetFace() = default;

            OctreeTetFace(uint32_t _globalIdx, uint32_t _owner, const TetFaceOrder &_vertGlobalIdx) : globalIdx(
                    _globalIdx), owner(_owner), face(_vertGlobalIdx) {}
        };

        struct Octree {
        public:
            Octree() = default;

            Octree(int _maxOctreeDepth);

            void printValidNodesInfo();

            void constructAllNodes();

            void constructNodeNeighbors();

            void constructTetNeighbors();

            void setMortonCode(std::shared_ptr<OctreeNode> &node);

            static constexpr short int vertSharedLUT[64] = {
                    0, 1, 3, 4, 9, 10, 12, 13,

                    1, 2, 4, 5, 10, 11, 13, 14,

                    3, 4, 6, 7, 12, 13, 15, 16,

                    4, 5, 7, 8, 13, 14, 16, 17,

                    9, 10, 12, 13, 18, 19, 21, 22,

                    10, 11, 13, 14, 19, 20, 22, 23,

                    12, 13, 15, 16, 21, 22, 24, 25,

                    13, 14, 16, 17, 22, 23, 25, 26
            };

            // edge: 02 23 31 10   46 67 75 54   04 26 37 15
            static constexpr short int edgeLocalIdxPair[12][2] = {
                    {0, 2},
                    {2, 3},
                    {3, 1},
                    {1, 0},

                    {4, 6},
                    {6, 7},
                    {7, 5},
                    {5, 4},

                    {0, 4},
                    {2, 6},
                    {3, 7},
                    {1, 5},
            };

            static constexpr short int edgeSharedLUT[48] = {
                    3, 4, 12, 13,
                    4, 7, 13, 16,
                    4, 5, 13, 14,
                    1, 4, 10, 13,

                    12, 13, 21, 22,
                    13, 16, 22, 25,
                    13, 14, 22, 23,
                    10, 13, 19, 22,

                    9, 10, 12, 13,
                    12, 13, 15, 16,
                    13, 14, 16, 17,
                    10, 11, 13, 14
            };

            static constexpr short int faceLocalIdxPair[6][4] = {
                    {0, 2, 3, 1}, // bottom

                    {4, 5, 7, 6}, // up

                    {0, 1, 5, 4}, // left

                    {2, 6, 7, 3}, // right

                    {1, 3, 7, 5}, // front

                    {0, 4, 6, 2}, // behind
            };

            static constexpr short int symFaceLocalIdxPair[6][4] = {
                    {0, 1, 3, 2}, // bottom

                    {4, 6, 7, 5}, // up

                    {0, 4, 5, 1}, // left

                    {2, 3, 7, 6}, // right

                    {1, 5, 7, 3}, // front

                    {0, 2, 6, 4}, // behind
            };

            static constexpr short int faceSharedLUT[6] = {
                    4, 22, 10, 16, 14, 12
            };

            void constructNodePrimitives();

            void constructTetPrimitives();

            bool visOctree(const std::string &filename) const;

            bool visOctreeEdge(const std::string &filename) const;

            bool visOctreeFace(const std::string &filename) const;

            bool visOctreeTet(const std::string &filename) const;

            bool visOctreeValidLeafNodes(const std::string &filename) const;

        public:
            int maxOctreeDepth;

            std::shared_ptr<OctreeNode> root;

            std::vector<uint64_t> preDepthSumNodes;
            std::vector<std::vector<std::shared_ptr<OctreeNode>>> depthNodes;
            std::vector<std::shared_ptr<OctreeNode>> allNodes;
            std::vector<std::shared_ptr<OctreeNode>> allLeafNodes;
            std::vector<std::shared_ptr<OctreeNode>> validLeafNodes;

            std::vector<OctreeTet> allTets;

            std::vector<OctreeVertex> nodeVertexArray;
            std::vector<OctreeEdge> nodeEdgeArray;
            std::vector<OctreeFace> nodeFaceArray;
            std::vector<OctreeTetVertex> tetVertexArray;

            std::vector<OctreeTetEdge> tetEdgeArray;
            std::vector<OctreeTetFace> tetFaceArray;
        };
    } // namespace geometry
NAMESPACE_END(PCO)

#endif //PCO_NEWOCTREE_HPP
