//
// Created by Lei on 11/3/2024.
//

#include "Octree.hpp"
#include "MortonLUT.hpp"
#include "utils/File.hpp"
#include <set>

NAMESPACE_BEGIN(PCO)
    namespace geometry {
        Octree::Octree(int _maxOctreeDepth) : maxOctreeDepth(_maxOctreeDepth) {
            root = std::make_shared<OctreeNode>();
            depthNodes.resize(maxOctreeDepth);
        }

        void Octree::printValidNodesInfo() {
            std::ofstream out("xyzinfo.txt");
            Vector3 root_origin = root->bbox.boxOrigin;
            for (int dep = 0; dep < maxOctreeDepth; ++dep) {
                for (const auto &node: depthNodes[dep]) {
                    out << "Depth: " << node->depth << std::endl;
                    out << "Origin: " << node->bbox.boxOrigin.transpose() << std::endl;
                    out << "End: " << node->bbox.boxEnd.transpose() << std::endl;
                    out << "xyz: " << node->xyzIdx.transpose() << std::endl;
                    if (!(root_origin.array() +
                          node->bbox.boxWidth.array() * node->xyzIdx.cast<double>().array()).isApprox(
                            node->bbox.boxOrigin.array()))
                        std::cout << "no!" << std::endl;
                    if (node->parent != nullptr && node->mortonCode / 8 != node->parent->mortonCode)
                        std::cout << "no2!" << std::endl;
                    out << "mortonCode: " << node->mortonCode << std::endl;
                    out << "=============\n";
                }
            }
        }

        void Octree::setMortonCode(std::shared_ptr<OctreeNode> &node) {
            uint16_t x = node->xyzIdx(0);
            uint16_t y = node->xyzIdx(1);
            uint16_t z = node->xyzIdx(2);
            // order: zyx
            node->mortonCode = mortonEncode_LUT(x, y, z);
        }

        void Octree::constructAllNodes() {
            // 使用前缀和形式存储所有节点
            size_t preSumNodes = 0;
            preDepthSumNodes.resize(maxOctreeDepth + 1);
            preDepthSumNodes[0] = 0;
            //std::ofstream out("xyzinfo.txt");
            for (int dep = 0; dep < maxOctreeDepth; ++dep) {
                for (auto &node: depthNodes[dep]) {
                    if (node->parent != nullptr) node->parent->childsIdx[node->childIdx] = allNodes.size();
                    node->globalIdx = allNodes.size();

                    allNodes.push_back(node);
                }
                preSumNodes += depthNodes[dep].size();
                preDepthSumNodes[dep + 1] = preSumNodes;
            }
        }

        void Octree::constructNodePrimitives() {
            // 后面需要改成只存储valid的vertex、edge和face，
            // 而不是存储所有叶子节点的vertex、edge和face

            // vertex
            std::map<Vector3i, uint32_t> xyzToGlobalIdx;
            for (const auto &node: depthNodes[maxOctreeDepth - 1]) {
                Vector3i originXYZIdx = node->xyzIdx;
                std::unordered_map<int, uint32_t> localCornerToGlbal;
                for (int i = 0; i < 8; ++i) {
                    uint16_t x, y, z;
                    morton3D_32_decode(i, x, y, z);
                    Vector3i cornerXYZIdx = originXYZIdx + Vector3i(x, y, z);
                    if (xyzToGlobalIdx.count(cornerXYZIdx)) {
                        localCornerToGlbal.insert(std::make_pair(i, xyzToGlobalIdx.at(cornerXYZIdx)));
                        continue;
                    }
                    xyzToGlobalIdx.insert(std::make_pair(cornerXYZIdx, nodeVertexArray.size()));
                    localCornerToGlbal.insert(std::make_pair(i, nodeVertexArray.size()));

                    Vector3 corner = node->bbox.corner(i);
                    uint32_t minMorton = node->mortonCode;
                    uint32_t owner = node->globalIdx;
                    for (int j = 0; j < 8; ++j) {
                        unsigned int sharedNodeIdx = node->neighborsIdx[vertSharedLUT[i * 8 + j]];
                        if (sharedNodeIdx == -1) continue;
                        uint32_t morton = allNodes[sharedNodeIdx]->mortonCode;
                        if (morton != -1 && morton < minMorton) {
                            owner = sharedNodeIdx;
                            minMorton = morton;
                        }
                    }
                    nodeVertexArray.emplace_back(nodeVertexArray.size(), owner, corner);
                }
                node->setGlobalTets(localCornerToGlbal);
            }

            // edge
            std::map<EdgeOrder, uint32_t> edgeOrderToGlobalIdx; // 也可以用morton code来编码边
            for (const auto &node: depthNodes[maxOctreeDepth - 1]) {
                Vector3i originXYZIdx = node->xyzIdx;
                for (int i = 0; i < 12; ++i) {
                    short int localEdgeIdx[2] = {edgeLocalIdxPair[i][0], edgeLocalIdxPair[i][1]};

                    uint16_t x1, y1, z1;
                    morton3D_32_decode(localEdgeIdx[0], x1, y1, z1);
                    Vector3i cornerXYZIdx_1 = originXYZIdx + Vector3i(x1, y1, z1);

                    uint16_t x2, y2, z2;
                    morton3D_32_decode(localEdgeIdx[1], x2, y2, z2);
                    Vector3i cornerXYZIdx_2 = originXYZIdx + Vector3i(x2, y2, z2);

                    uint32_t vertexGlobalIdx_1 = xyzToGlobalIdx.at(cornerXYZIdx_1);
                    uint32_t vertexGlobalIdx_2 = xyzToGlobalIdx.at(cornerXYZIdx_2);
                    EdgeOrder edge_1 = std::make_pair(vertexGlobalIdx_1, vertexGlobalIdx_2);
                    EdgeOrder edge_2 = std::make_pair(vertexGlobalIdx_2, vertexGlobalIdx_1);

                    if (edgeOrderToGlobalIdx.count(edge_1) || edgeOrderToGlobalIdx.count(edge_2)) continue;
                    edgeOrderToGlobalIdx.insert(std::make_pair(edge_1, nodeEdgeArray.size()));

                    std::pair<uint32_t, uint32_t> edge = edge_1;

                    uint32_t minMorton = node->mortonCode;
                    uint32_t owner = node->globalIdx;
                    for (int j = 0; j < 4; ++j) {
                        unsigned int sharedNodeIdx = node->neighborsIdx[edgeSharedLUT[i * 4 + j]];
                        if (sharedNodeIdx == -1) continue;
                        uint32_t morton = allNodes[sharedNodeIdx]->mortonCode;
                        if (morton != -1 && morton < minMorton) {
                            owner = sharedNodeIdx;
                            minMorton = morton;
                        }
                    }
                    nodeEdgeArray.emplace_back(nodeEdgeArray.size(), owner, edge);
                }
            }

            // face
            std::map<FaceOrder, uint32_t> faceOrderToGlobalIdx; // 也可以用morton code来编码边
            for (const auto &node: depthNodes[maxOctreeDepth - 1]) {
                Vector3i originXYZIdx = node->xyzIdx;
                for (int i = 0; i < 6; ++i) {
                    short int localFaceIdx_1[4] = {faceLocalIdxPair[i][0], faceLocalIdxPair[i][1],
                                                   faceLocalIdxPair[i][2], faceLocalIdxPair[i][3]};

                    std::array<Vector3i, 4> cornerXYZIdx_1;
                    for (int m = 0; m < 4; ++m) {
                        uint16_t x, y, z;
                        morton3D_32_decode(localFaceIdx_1[m], x, y, z);
                        cornerXYZIdx_1[m] = originXYZIdx + Vector3i(x, y, z);
                    }
                    uint32_t vertexGlobalIdx_1 = xyzToGlobalIdx.at(cornerXYZIdx_1[0]);
                    uint32_t vertexGlobalIdx_2 = xyzToGlobalIdx.at(cornerXYZIdx_1[1]);
                    uint32_t vertexGlobalIdx_3 = xyzToGlobalIdx.at(cornerXYZIdx_1[2]);
                    uint32_t vertexGlobalIdx_4 = xyzToGlobalIdx.at(cornerXYZIdx_1[3]);
                    FaceOrder face_1 = std::make_tuple(vertexGlobalIdx_1, vertexGlobalIdx_2, vertexGlobalIdx_3,
                                                       vertexGlobalIdx_4);

                    short int localFaceIdx_2[4] = {symFaceLocalIdxPair[i][0], symFaceLocalIdxPair[i][1],
                                                   symFaceLocalIdxPair[i][2], symFaceLocalIdxPair[i][3]};
                    std::array<Vector3i, 4> cornerXYZIdx_2;
                    for (int m = 0; m < 4; ++m) {
                        uint16_t x, y, z;
                        morton3D_32_decode(localFaceIdx_2[m], x, y, z);
                        cornerXYZIdx_2[m] = originXYZIdx + Vector3i(x, y, z);
                    }
                    uint32_t vertexGlobalIdx_5 = xyzToGlobalIdx.at(cornerXYZIdx_2[0]);
                    uint32_t vertexGlobalIdx_6 = xyzToGlobalIdx.at(cornerXYZIdx_2[1]);
                    uint32_t vertexGlobalIdx_7 = xyzToGlobalIdx.at(cornerXYZIdx_2[2]);
                    uint32_t vertexGlobalIdx_8 = xyzToGlobalIdx.at(cornerXYZIdx_2[3]);
                    FaceOrder face_2 = std::make_tuple(vertexGlobalIdx_5, vertexGlobalIdx_6, vertexGlobalIdx_7,
                                                       vertexGlobalIdx_8);
                    if (faceOrderToGlobalIdx.count(face_1) || faceOrderToGlobalIdx.count(face_2)) continue;
                    faceOrderToGlobalIdx.insert(std::make_pair(face_1, nodeFaceArray.size()));

                    FaceOrder face = face_1;

                    uint32_t minMorton = node->mortonCode;
                    uint32_t owner = node->globalIdx;
                    {
                        unsigned int sharedNodeIdx = node->neighborsIdx[faceSharedLUT[i]];
                        if (sharedNodeIdx != -1) {
                            uint32_t morton = allNodes[sharedNodeIdx]->mortonCode;
                            if (morton != -1 && morton < minMorton) {
                                owner = sharedNodeIdx;
                                minMorton = morton;
                            }
                        }
                    }
                    nodeFaceArray.emplace_back(nodeFaceArray.size(), owner, face);
                }
            }
        }

        void Octree::constructTetPrimitives() {
            // 后面需要改成只存储valid的vertex、edge和face，
            // 而不是存储所有叶子节点的vertex、edge和face

            // tet and vertex
            std::map<Vector3i, uint32_t> xyzToGlobalIdx;
            std::map<TetEdgeOrder, uint32_t> edgeOrderToGlobalIdx; // 也可以用morton code来编码边
            std::map<TetFaceOrder, uint32_t> faceOrderToGlobalIdx;
            //for (auto& node : depthNodes[maxOctreeDepth - 1]) {
            for (auto &node: validLeafNodes) {
                Vector3i originXYZIdx = node->xyzIdx;
                //std::unordered_map<int, uint32_t> localCornerToGlbal;

                node->setLocalTets();
                for (int i = 0; i < 5; ++i) {
                    auto &tet = node->getTet(i);
                    tet.localIdx = i;
                    tet.globalIdx = allTets.size();
                    tet.globalNodeIdx = node->globalIdx;

                    // vertex
                    for (int j = 0; j < 4; ++j) {
                        int localCornerInNode = tet.vertsInLocalNode[j];
                        Vector3 corner = node->bbox.corner(localCornerInNode);
                        tet.coord.row(j) = corner;

                        uint16_t x, y, z;
                        morton3D_32_decode(localCornerInNode, x, y, z);
                        Vector3i cornerXYZIdx = originXYZIdx + Vector3i(x, y, z);

                        if (xyzToGlobalIdx.count(cornerXYZIdx)) {
                            tet.globalVerts[j] = xyzToGlobalIdx.at(cornerXYZIdx);
                            //localCornerToGlbal.insert(std::make_pair(i, xyzToGlobalIdx.at(cornerXYZIdx)));
                            tet.globalVertsId[j] = tet.globalVerts[j];
                            continue;
                        }
                        xyzToGlobalIdx.insert(std::make_pair(cornerXYZIdx, tetVertexArray.size()));
                        tet.globalVerts[j] = xyzToGlobalIdx.at(cornerXYZIdx);

                        tet.globalVertsId[j] = tetVertexArray.size();
                        // 直接将owner分配给第一个计算到该点的四面体
                        OctreeTetVertex octreeTetVertex;
                        octreeTetVertex.globalIdx = tetVertexArray.size();
                        octreeTetVertex.coord = tet.coord.row(j);
                        /*std::cout << "vertex: " << octreeTetVertex.coord.transpose() << std::endl;
                        std::cout << "cornerXYZIdx: " << cornerXYZIdx.transpose() << std::endl;*/
                        octreeTetVertex.owner = tet.globalIdx;
                        tetVertexArray.emplace_back(octreeTetVertex);
                    }

                    // edge
                    for (int j = 0; j < 6; ++j) {
                        int localEdgeIdx[2] = {tet.edgesInLocalNode[j][0],
                                               tet.edgesInLocalNode[j][1]}; // 此处edge存储的vertex id是一个node中的loca id，即0~7

                        std::array<Vector3i, 2> cornerXYZIdx;
                        for (int m = 0; m < 2; ++m) {
                            uint16_t x, y, z;
                            morton3D_32_decode(localEdgeIdx[m], x, y, z);
                            cornerXYZIdx[m] = originXYZIdx + Vector3i(x, y, z);
                        }

                        uint32_t vertexGlobalIdx_1 = xyzToGlobalIdx.at(cornerXYZIdx[0]);
                        uint32_t vertexGlobalIdx_2 = xyzToGlobalIdx.at(cornerXYZIdx[1]);
                        TetEdgeOrder edge_1 = {vertexGlobalIdx_1, vertexGlobalIdx_2};
                        tet.globalEdges[j] = edge_1;

                        TetEdgeOrder edge_2 = {vertexGlobalIdx_2, vertexGlobalIdx_1};
                        if (edgeOrderToGlobalIdx.count(edge_1) || edgeOrderToGlobalIdx.count(edge_2)) {
                            if (edgeOrderToGlobalIdx.count(edge_1))
                                tet.globalEdgesId[j] = edgeOrderToGlobalIdx.at(edge_1);
                            else tet.globalEdgesId[j] = edgeOrderToGlobalIdx.at(edge_2);

                            continue;
                        }
                        edgeOrderToGlobalIdx.insert(std::make_pair(edge_1, tetEdgeArray.size()));

                        tet.globalEdgesId[j] = tetEdgeArray.size();
                        OctreeTetEdge octreeTetEdge;
                        octreeTetEdge.globalIdx = tetEdgeArray.size();
                        octreeTetEdge.edge = edge_1;
                        octreeTetEdge.owner = tet.globalIdx;

                        tetEdgeArray.emplace_back(octreeTetEdge);
                    }

                    // face
                    for (int j = 0; j < 4; ++j) {
                        int localFaceIdx[3] = {tet.facesInLocalNode[j][0], tet.facesInLocalNode[j][1],
                                               tet.facesInLocalNode[j][2]};

                        std::array<Vector3i, 3> cornerXYZIdx;
                        for (int m = 0; m < 3; ++m) {
                            uint16_t x, y, z;
                            morton3D_32_decode(localFaceIdx[m], x, y, z);
                            cornerXYZIdx[m] = originXYZIdx + Vector3i(x, y, z);
                        }
                        uint32_t vertexGlobalIdx_1 = xyzToGlobalIdx.at(cornerXYZIdx[0]);
                        uint32_t vertexGlobalIdx_2 = xyzToGlobalIdx.at(cornerXYZIdx[1]);
                        uint32_t vertexGlobalIdx_3 = xyzToGlobalIdx.at(cornerXYZIdx[2]);
                        TetFaceOrder face = {vertexGlobalIdx_1, vertexGlobalIdx_2, vertexGlobalIdx_3};
                        tet.globalFaces[j] = face;

                        TetFaceOrder sortedFace = face;
                        std::sort(sortedFace.begin(), sortedFace.end());

                        if (faceOrderToGlobalIdx.count(sortedFace)) {
                            tet.globalFacesId[j] = faceOrderToGlobalIdx.at(sortedFace);
                            continue;
                        }
                        faceOrderToGlobalIdx.insert(std::make_pair(sortedFace, tetFaceArray.size()));

                        tet.globalFacesId[j] = tetFaceArray.size();
                        OctreeTetFace octreeTetFace;
                        octreeTetFace.globalIdx = tetFaceArray.size();
                        octreeTetFace.face = face;
                        octreeTetFace.owner = tet.globalIdx;

                        tetFaceArray.emplace_back(octreeTetFace);
                    }

                    allTets.emplace_back(tet);
                }
            }
        }

        void Octree::constructNodeNeighbors() {
            depthNodes[0][0]->neighborsIdx[13] = 0;
            for (int dep = 1; dep < maxOctreeDepth; ++dep) {
                //#pragma omp parallel for
                for (int i = 0; i < depthNodes[dep].size(); ++i) {
                    //for (auto& node : depthNodes[dep]) {
                    auto &node = depthNodes[dep][i];
                    for (int j = 0; j < 27; ++j) {
                        const uint8_t key = (node->mortonCode) & LOWER_3BIT_MASK;
                        const unsigned int neighborIdx = node->parent->neighborsIdx[neighbor_LUTparent[key][j]];
                        if (neighborIdx != -1) {
                            auto &h = allNodes[neighborIdx];
                            node->neighborsIdx[j] = h->childsIdx[neighbor_LUTchild[key][j]];
                        }
                    }
                }
            }
        }

        void Octree::constructTetNeighbors() {
            constructNodeNeighbors();

            int tetCnt = 0;
            // set global neighbors
//#pragma omp parallel for
            for (int i = 0; i < allTets.size(); ++i) {
                auto &tet = allTets[i];
                const auto &node = allNodes[tet.globalNodeIdx];

                //std::cout << "Check tet #" << tetCnt++ << " neighbors info...\n";

                //std::cout << "Check edge's neghbor...\n";
                for (int j = 0; j < 6; ++j) {
                    auto &globalEdgeNeighbor = tet.globalEdgeNeighbors[j];
                    const auto &localEdgeNeighbor = tet.localEdgeNeighbors[j];


                    for (int k = 0; k < 6; ++k) {

                        int nodeNeighbor = localEdgeNeighbor[k][0];
                        if (nodeNeighbor == -1 || node->neighborsIdx[nodeNeighbor] == -1) {
                            globalEdgeNeighbor[k] = -1;
                            continue;
                        }

                        const auto &neighborNode = allNodes[node->neighborsIdx[nodeNeighbor]];
                        int localTetIdx = localEdgeNeighbor[k][1];

                        globalEdgeNeighbor[k] = neighborNode->getTet(localTetIdx).globalIdx;

                        Eigen::RowVector3d ed_0 = tetVertexArray[tet.globalEdges[j][0]].coord.transpose();
                        Eigen::RowVector3d ed_1 = tetVertexArray[tet.globalEdges[j][1]].coord.transpose();
                        int correct_vertex = 0;
                        for (int m = 0; m < 4; ++m) {
                            if (allTets[globalEdgeNeighbor[k]].coord.row(m) == ed_0) ++correct_vertex;
                            else if (allTets[globalEdgeNeighbor[k]].coord.row(m) == ed_1) ++correct_vertex;
                        }
                    }
                }

                //std::cout << "Check face's neghbor...\n";
                for (int j = 0; j < 4; ++j) {
                    auto &globalFaceNeighbor = tet.globalFaceNeighbors[j];
                    const auto &localFaceNeighbor = tet.localFaceNeighbors[j];

                    int nodeNeighbor = localFaceNeighbor[0];
                    if (nodeNeighbor == -1 || node->neighborsIdx[nodeNeighbor] == -1) {
                        globalFaceNeighbor = -1;
                        continue;
                    }

                    const auto &neighborNode = allNodes[node->neighborsIdx[nodeNeighbor]];
                    int localTetIdx = localFaceNeighbor[1];
                    globalFaceNeighbor = neighborNode->getTet(localTetIdx).globalIdx;

                    Eigen::RowVector3d ed_0 = tetVertexArray[tet.globalFaces[j][0]].coord.transpose();
                    Eigen::RowVector3d ed_1 = tetVertexArray[tet.globalFaces[j][1]].coord.transpose();
                    Eigen::RowVector3d ed_2 = tetVertexArray[tet.globalFaces[j][2]].coord.transpose();
                    int correct_vertex = 0;
                    for (int m = 0; m < 4; ++m) {
                        if (allTets[globalFaceNeighbor].coord.row(m) == ed_0) ++correct_vertex;
                        else if (allTets[globalFaceNeighbor].coord.row(m) == ed_1) ++correct_vertex;
                        else if (allTets[globalFaceNeighbor].coord.row(m) == ed_2) ++correct_vertex;
                    }
                }
                ++tetCnt;
            }
        }

        bool
        Octree::visOctree(const std::string &filename) const {
            utils::file::checkDir(filename);
            std::ofstream out(filename);
            if (!out) {
                return false;
            }

            int vertBegIdx = 0;
            for (auto &leafNode: allLeafNodes) {
                leafNode->bbox.output(out, vertBegIdx);
            }

            out.close();
            return true;
        }

        bool Octree::visOctreeEdge(const std::string &filename) const {
            utils::file::checkDir(filename);
            std::ofstream out(filename);
            if (!out) return false;

            for (const auto &vertex: nodeVertexArray) {
                out << "v " << vertex.coord.transpose() << std::endl;
            }

            for (const auto &edge: nodeEdgeArray) {
                out << "l " << edge.vertGlobalIdx.first + 1 << " " << edge.vertGlobalIdx.second + 1 << std::endl;
            }

            return false;
        }

        bool Octree::visOctreeFace(const std::string &filename) const {
            utils::file::checkDir(filename);
            std::ofstream out(filename);
            if (!out) {
                return false;
            }

            for (const auto &vertex: nodeVertexArray) {
                out << "v " << vertex.coord.transpose() << std::endl;
            }

            for (const auto &face: nodeFaceArray) {
                out << "f "
                    << std::get<0>(face.vertGlobalIdx) + 1
                    << " " << std::get<1>(face.vertGlobalIdx) + 1
                    << " " << std::get<2>(face.vertGlobalIdx) + 1
                    << " " << std::get<3>(face.vertGlobalIdx) + 1
                    << std::endl;
            }

            return false;
        }

        bool Octree::visOctreeTet(const std::string &filename) const {
            utils::file::checkDir(filename);
            std::ofstream out(filename);
            if (!out) {
                return false;
            }

            for (const auto &tetVertex: tetVertexArray) {
                out << "v " << tetVertex.coord.transpose() << std::endl;
            }

            for (const auto &tetEdge: tetEdgeArray) {
                out << "l " << tetEdge.edge[0] + 1 << " " << tetEdge.edge[1] + 1 << std::endl;
            }

            for (const auto &tetFace: tetFaceArray) {
                out << "f " << tetFace.face[0] + 1 << " " << tetFace.face[1] + 1 << " " << tetFace.face[2] + 1
                    << std::endl;
            }

            return false;
        }


        bool
        Octree::visOctreeValidLeafNodes(const std::string &filename) const {
            utils::file::checkDir(filename);
            std::ofstream out(filename);
            if (!out) {
                return false;
            }

            int vertBegIdx = 0;
            for (size_t i = 0; i < validLeafNodes.size(); ++i) {
                out << "# i: " << i << std::endl;
                auto &leafNode = validLeafNodes[i];
                leafNode->bbox.output(out, vertBegIdx);
            }

            out.close();
            return true;
        }

    } // namespace geometry
NAMESPACE_END(PCO)