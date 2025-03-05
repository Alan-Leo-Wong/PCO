#pragma once

#include "Complex.hpp"
#include "PlaneRepo.hpp"
#include "detail/BasicDataType.hpp"

#include <array>
#include <vector>

NAMESPACE_BEGIN(PCO)
    namespace complex {

        /**
         * A self-contained data structure for 2D or 3D arrangement representation.
         */
        struct OffsetSurface {
            static constexpr size_t None = std::numeric_limits<size_t>::max();

            struct Edge {
                size_t globalId = -1; // if the edge belongs the tetrahedron
                std::array<size_t, 2> vertices = { None, None };
                std::array<size_t, 2> supporting_planes = {None, None};
            };

            /**
             * A face represents a (DIM-1)-dimensional polytope.  For 2D and 3D, its
             * orientation is uniquely defined by ordering of its boundary vertices.
             */
            struct Face {
                size_t globalId = -1; // if the face belongs the tetrahedron
                /**
                 * An ordered list of boundary vertices.
                 *
                 * In 3D, the face is always oriented counterclockwise when viewed from
                 * the positive side of the supporting plane.
                 *
                 * In 2D, the face (aka edge) is oriented such that the positive side of
                 * the supporting plane (aka line) is on the right.
                 */
                std::vector<size_t> vertices;
                std::vector<size_t> global_verticesId;

                std::vector<Edge> bd_edges;

                /**
                 * Plane index of the supporting plane.
                 */
                size_t supporting_plane = None;

                /**
                 * The cell index on the positive and negative side of this face.
                 */
                size_t positive_cell = None;
                size_t negative_cell = None;
            };

            /**
             * A cell is a DIM-dimensional polytope.
             */
            struct Cell {
                /**
                 * A set of boundary face indices in no particular order.
                 */
                std::vector<size_t> faces;
            };

            void clear() {
                vertices.clear();
                vertices_coord.clear();
                faces.clear();
                cells.clear();
            }

            std::vector<Vertex> vertices;
            std::vector<VertexCoord> vertices_coord;
            std::vector<Face> faces;
            std::vector<Cell> cells;

            // Note: the following structure is only non-empty if input planes contain
            // duplicates.
            /*std::vector<size_t> unique_plane_indices;
            std::vector<std::vector<size_t>> unique_planes;
            std::vector<bool> unique_plane_orientations;*/
        };

    } // namespace complex
NAMESPACE_END(PCO)