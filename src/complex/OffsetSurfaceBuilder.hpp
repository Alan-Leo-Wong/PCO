#pragma once

#include "ComplexCut.hpp"
#include "PlaneRepo.hpp"
#include "OffsetSurface.hpp"
#include "utils/dsu/DSU.hpp"

#include <algorithm>
#include <array>
#include <type_traits>
#include <iostream>

NAMESPACE_BEGIN(PCO)
    namespace complex {

        class OffsetSurfaceBuilder {
        public:
            OffsetSurfaceBuilder(int tetLocalIdx,
                                 CVertexCoord complex_vertices_coord,
                                 bool _is_positive,
                                 const std::array<size_t, 4> &vertices_globalId,
                                 const std::array<size_t, 6> &edges_globalId,
                                 const std::array<size_t, 4> &faces_globalId);

            OffsetSurfaceBuilder(int tetLocalIdx,
                                 CVertexCoord complex_vertices_coord,
                                 const std::vector<Plane> &planes,
                                 bool _is_positive,
                                 const std::array<size_t, 4> &vertices_globalId,
                                 const std::array<size_t, 6> &edges_globalId,
                                 const std::array<size_t, 4> &faces_globalId);

        public:
            size_t add_new_plane(const Plane &plane, size_t global_planeId);

            const Complex &get_complex() { return complexCut.get_complex(); }

            const utils::dsu::DSU &get_coplanar_planes() { return m_coplanar_planes; }

            const std::vector<size_t> &get_global_extPlaneId() { return global_extPlaneId; }

            const OffsetSurface &get_offset_surface() const { return m_offset_surface; }

            bool is_empty() const;

            void extract_offset_surface();

            OffsetSurface &get_offset_surface() { return m_offset_surface; }

            OffsetSurface &&export_offset_surface() &&{ return std::move(m_offset_surface); }

            const PlaneRepo &getPlaneRepo() const { return m_planes; }

        private:
            void extract_offset_surface(const Complex &complex);

            VertexCoord compute_vertex_coord(const Vertex &p);

            Complex initialize_complex(int tetLocalIdx,
                                       const std::array<size_t, 4> &vertices_globalId,
                                       const std::array<size_t, 6> &edges_globalId,
                                       const std::array<size_t, 4> &faces_globalId);

            Complex initialize_complex(int tetLocalIdx, size_t num_planes,
                                       const std::array<size_t, 4> &vertices_globalId,
                                       const std::array<size_t, 6> &edges_globalId,
                                       const std::array<size_t, 4> &faces_globalId);

        private:
            bool is_positive = true;

            ComplexCut complexCut;
            CVertexCoord m_complex_vertices_coord;

            PlaneRepo m_planes;
            std::vector<size_t> global_extPlaneId;
            utils::dsu::DSU m_coplanar_planes;
            OffsetSurface m_offset_surface;
        };

    } // namespace complex
NAMESPACE_END(PCO)