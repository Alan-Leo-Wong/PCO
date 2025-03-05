//
// Created by Lei on 11/3/2024.
//

#ifndef PCO_OCTREETET_HPP
#define PCO_OCTREETET_HPP

#include "OctreeTetLUT.hpp"
#include "Geometry.hpp"

NAMESPACE_BEGIN(PCO)
    namespace geometry {
        struct OctreeTet {
            int type;
            uint32_t localIdx = -1;
            uint32_t globalIdx = -1;
            uint32_t globalNodeIdx = -1;

            TetCoord coord;

            static constexpr std::array<TetVertex, 4> localVerts = {0, 1, 2, 3};
            static constexpr std::array<TetEdge, 6> localEdges = {
                    {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}}; // don't consider the orders
            std::array<TetFace, 4> localFaces;

            std::array<TetVertex, 4> vertsInLocalNode;
            std::array<TetEdge, 6> edgesInLocalNode;
            std::array<TetFace, 4> facesInLocalNode;

            std::array<TetGlobalVertex, 4> globalVerts;
            std::array<TetGlobalEdge, 6> globalEdges;
            std::array<TetGlobalFace, 4> globalFaces;

            std::array<TetEdgeNeighbor, 6> localEdgeNeighbors;
            std::array<TetFaceNeighbor, 4> localFaceNeighbors;

            std::array<TetEdgeGlobalNeighbor, 6> globalEdgeNeighbors;
            std::array<TetFaceGlobalNeighbor, 4> globalFaceNeighbors;

            std::array<size_t, 4> globalVertsId;
            std::array<size_t, 6> globalEdgesId;
            std::array<size_t, 4> globalFacesId;

            static std::array<int, 3> getPointAdjFace(int vertIdx) {
                int vert = localVerts[vertIdx];
                std::array<int, 3> res;
                int cnt = 0;
                for (int i = 0; i < 4; ++i) {
                    if (i != vert) res[cnt++] = i;
                }

                return res;
            }

            static std::array<int, 2> getEdgeAdjFace(int edgeIdx) {
                const auto &edge = localEdges[edgeIdx];
                std::array<int, 2> res;
                int cnt = 0;
                for (int i = 0; i < 4; ++i) {
                    if (i != edge[0] && i != edge[1]) res[cnt++] = i;
                }

                return res;
            }
        };

    } // namespace geometry
NAMESPACE_END(PCO)

#endif //PCO_OCTREETET_HPP
