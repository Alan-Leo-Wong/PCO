#pragma once

#include "detail/BasicDataType.hpp"

#include <array>
#include <vector>

NAMESPACE_BEGIN(PCO)
    namespace complex {
        static constexpr size_t INVALID = std::numeric_limits<size_t>::max();
        static constexpr int DIM = 3;

        /**
         * A point is represented as the intersection of planes. We store the index
         * of the plane here.
         */
        template<int DIM>
        using Point = std::array<size_t, DIM>;

        using Vertex = Point<3>;

        using VertexCoord = Eigen::Matrix<Scalar, 3, 1>;

        using CVertexCoord = Eigen::Matrix<Scalar, 4, 3>;

        struct Edge {
            size_t globalId = -1;
            std::array<size_t, 2> vertices = {INVALID, INVALID}; ///< unordered.
            std::array<size_t, 2> supporting_planes = {INVALID, INVALID};
        };

        struct Face {
            size_t globalId = -1;
            std::vector<size_t> edges; ///< ordered.
            size_t supporting_plane = INVALID;
            size_t positive_cell = INVALID;
            size_t negative_cell = INVALID;
        };

        struct Cell {
            std::vector<size_t> faces; ///< unordered.
            std::vector<bool> signs; ///< sign of each implicit function.
        };

        struct Complex {
            std::vector<Vertex> vertices;
            std::vector<Edge> edges;
            std::vector<Face> faces;
            std::vector<Cell> cells;

            std::array<size_t, 4> global_verticesId;
            std::array<size_t, 6> global_oriEdgesId;
            std::array<size_t, 4> global_oriFacesId;
            std::array<std::vector<size_t>, 4> oriFacesEdges;
        };

    } // namespace complex
NAMESPACE_END(PCO)
