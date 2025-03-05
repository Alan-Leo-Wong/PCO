//
// Created by Lei on 3/31/2024.
//

#include "OffsetSurfaceBuilder.hpp"

#include "ComplexCut.hpp"
#include "detail/Handler.hpp"
#include "utils/Log.hpp"
#include <utility>


NAMESPACE_BEGIN(PCO)
    namespace complex {

        OffsetSurfaceBuilder::OffsetSurfaceBuilder(int tetLocalIdx, CVertexCoord complex_vertices_coord,
                                                   bool _is_positive,
                                                   const std::array<size_t, 4> &vertices_globalId,
                                                   const std::array<size_t, 6> &edges_globalId,
                                                   const std::array<size_t, 4> &faces_globalId)
                : m_complex_vertices_coord(std::move(complex_vertices_coord)),
                  is_positive(_is_positive) {
            auto complex = initialize_complex(tetLocalIdx, vertices_globalId,
                                              edges_globalId, faces_globalId);
            m_coplanar_planes.init(DIM + 1);
            complexCut = ComplexCut(complex);
        }

        OffsetSurfaceBuilder::OffsetSurfaceBuilder(int tetLocalIdx, CVertexCoord complex_vertices_coord,
                                                   const std::vector<Plane> &planes, bool _is_positive,
                                                   const std::array<size_t, 4> &vertices_globalId,
                                                   const std::array<size_t, 6> &edges_globalId,
                                                   const std::array<size_t, 4> &faces_globalId)
                : m_complex_vertices_coord(std::move(complex_vertices_coord)),
                  m_planes(planes), is_positive(_is_positive) {
            const size_t num_planes = m_planes.get_num_planes();
            m_coplanar_planes.init(num_planes + DIM + 1);

            auto complex = initialize_complex(tetLocalIdx, num_planes + DIM + 1, vertices_globalId,
                                              edges_globalId, faces_globalId);
            complexCut = ComplexCut(m_planes, complex);

            for (size_t i = 0; i < num_planes; i++) {
                if (complexCut.get_complex().vertices.empty()) break;
                size_t plane_id = i + DIM + 1;
                size_t coplanar_plane_id = complexCut.add_plane(plane_id, is_positive);
                if (coplanar_plane_id != INVALID) {
                    LOG::DEBUG("coplanar plane id: {}", coplanar_plane_id);
                    m_coplanar_planes.merge(plane_id, coplanar_plane_id);
                }
            }

            extract_offset_surface(complexCut.get_complex());
        }

        size_t OffsetSurfaceBuilder::add_new_plane(const Plane &plane,
                                                   size_t global_planeId) {
            m_planes.add_new_plane(plane);
            size_t plane_id = m_planes.get_num_planes() + DIM;
            m_coplanar_planes.fa.push_back(plane_id);
            m_coplanar_planes.fa.back() = plane_id;
            size_t coplanar_plane_id = complexCut.add_plane(plane, is_positive);

            if (coplanar_plane_id != INVALID) {
                LOG::DEBUG("coplanar plane id: {}", coplanar_plane_id);
                if (coplanar_plane_id > 3) {
                    if (global_extPlaneId[coplanar_plane_id - 4] < global_planeId)
                        m_coplanar_planes.merge(plane_id, coplanar_plane_id);
                    else
                        m_coplanar_planes.merge(coplanar_plane_id, plane_id);
                } else {
                    m_coplanar_planes.merge(plane_id, coplanar_plane_id);
                }
            }

            global_extPlaneId.push_back(global_planeId);
            return plane_id;
        }

        bool OffsetSurfaceBuilder::is_empty() const {
            return complexCut.get_complex().faces.empty();
        }

        void OffsetSurfaceBuilder::extract_offset_surface() {
            extract_offset_surface(complexCut.get_complex());
        }

        void OffsetSurfaceBuilder::extract_offset_surface(const Complex &complex) {
            m_offset_surface.clear();
            m_offset_surface.vertices = complex.vertices;

            auto &edges = complex.edges;
            auto &faces = complex.faces;
            size_t num_faces = faces.size();
            m_offset_surface.faces.clear();
            m_offset_surface.faces.reserve(num_faces);

            for (size_t i = 0; i < num_faces; i++) {
                auto &cf = faces[i];
                size_t supporting_plane = m_coplanar_planes.find(cf.supporting_plane);
                if (supporting_plane <= DIM)
                    continue;

                m_offset_surface.faces.emplace_back();
                // auto &f = m_offset_surface.faces[i];
                auto &f = m_offset_surface.faces.back();
                if (cf.supporting_plane <= 3) {
                    if (cf.globalId == -1)
                        LOG::ERROR("Error global id in facet!");
                    f.globalId = cf.globalId;
                }
                const size_t num_bd_edges = cf.edges.size();
                PCO_ASSERT(num_bd_edges >= 3);
                f.bd_edges.reserve(num_bd_edges);
                f.vertices.reserve(num_bd_edges);

                size_t first_vertex = INVALID, last_vertex = INVALID;
                for (size_t j = 0; j < num_bd_edges; j++) {
                    auto &curr_e = edges[cf.edges[j]];
                    auto &next_e = edges[cf.edges[(j + 1) % num_bd_edges]];

                    if (curr_e.vertices[0] == next_e.vertices[0] ||
                        curr_e.vertices[0] == next_e.vertices[1]) {
                        f.vertices.push_back(curr_e.vertices[0]);

                        if (curr_e.vertices[0] <= 3)
                            f.global_verticesId.push_back(
                                    complex.global_verticesId[curr_e.vertices[0]]);
                        else
                            f.global_verticesId.push_back(-1);

                        if (last_vertex == INVALID) {
                            first_vertex = curr_e.vertices[0];
                            last_vertex = curr_e.vertices[0];
                        } else {
                            OffsetSurface::Edge bd_edge;
                            if (curr_e.supporting_planes[0] <= 3 &&
                                curr_e.supporting_planes[1] <= 3) {
                                if (curr_e.globalId == -1)
                                    LOG::ERROR("Error global id in edge!");
                                bd_edge.globalId = curr_e.globalId;
                            }
                            bd_edge.vertices = {last_vertex, curr_e.vertices[0]};
                            bd_edge.supporting_planes = curr_e.supporting_planes;
                            f.bd_edges.push_back(bd_edge);
                            last_vertex = curr_e.vertices[0];
                        }
                    } else {
                        PCO_ASSERT(curr_e.vertices[1] == next_e.vertices[0] ||
                                   curr_e.vertices[1] == next_e.vertices[1]);
                        f.vertices.push_back(curr_e.vertices[1]);

                        if (curr_e.vertices[1] <= 3)
                            f.global_verticesId.push_back(
                                    complex.global_verticesId[curr_e.vertices[1]]);
                        else
                            f.global_verticesId.push_back(-1);

                        if (last_vertex == INVALID) {
                            first_vertex = curr_e.vertices[1];
                            last_vertex = curr_e.vertices[1];
                        } else {
                            OffsetSurface::Edge bd_edge;
                            if (curr_e.supporting_planes[0] <= 3 &&
                                curr_e.supporting_planes[1] <= 3) {
                                if (curr_e.globalId == -1)
                                    std::cerr << "Error global id in edge!\n";
                                bd_edge.globalId = curr_e.globalId;
                            }
                            bd_edge.vertices = {last_vertex, curr_e.vertices[1]};
                            bd_edge.supporting_planes = curr_e.supporting_planes;
                            f.bd_edges.push_back(bd_edge);
                            last_vertex = curr_e.vertices[1];
                        }
                    }
                }
                OffsetSurface::Edge bd_edge;
                auto &curr_e = edges[cf.edges[0]];
                if (curr_e.supporting_planes[0] <= 3 && curr_e.supporting_planes[1] <= 3) {
                    if (curr_e.globalId == -1)
                        std::cerr << "Error global id in edge!\n";
                    bd_edge.globalId = curr_e.globalId;
                }
                bd_edge.vertices = {last_vertex, first_vertex};
                bd_edge.supporting_planes = curr_e.supporting_planes;
                f.bd_edges.push_back(bd_edge);

                f.positive_cell = cf.positive_cell;
                f.negative_cell = cf.negative_cell;
                f.supporting_plane = supporting_plane;
            }

            auto &cells = complex.cells;
            size_t num_cells = cells.size();
            m_offset_surface.cells.resize(num_cells);

            for (size_t i = 0; i < num_cells; i++) {
                auto &cc = cells[i];
                auto &c = m_offset_surface.cells[i];
                c.faces = std::move(cc.faces);
            }

            auto &vc = m_offset_surface.vertices_coord;
            for (const auto &vertex: m_offset_surface.vertices) {
                vc.push_back(compute_vertex_coord(vertex));
            }
        }

        VertexCoord OffsetSurfaceBuilder::compute_vertex_coord(const Vertex &p) {
            Eigen::Matrix<Scalar, DIM + 1, DIM + 1> A;
            for (int i = 0; i < DIM; ++i) {
                const auto &plane = m_planes.get_plane(p[i]);
                Eigen::Matrix<Scalar, DIM + 1, 1> vec;
                for (int j = 0; j < DIM + 1; ++j) {
                    vec(j) = plane[j];
                }
                A.row(i) = vec;
            }
            A.row(DIM) = Eigen::Matrix<Scalar, DIM + 1, 1>::Ones();

            Eigen::Matrix<Scalar, DIM + 1, 1> b;
            b.setZero();
            b(DIM) = 1;

#ifdef NON_ROBUST_COMPLEX_CUTTING
            auto x = A.inverse() * b;
#else
            auto x = A.fullPivHouseholderQr().solve(b);
#endif

            VertexCoord coord = m_complex_vertices_coord.transpose() * x;
            return coord;
        }

        Complex
        OffsetSurfaceBuilder::initialize_complex(int tetLocalIdx, const std::array<size_t, 4> &vertices_globalId,
                                                 const std::array<size_t, 6> &edges_globalId,
                                                 const std::array<size_t, 4> &faces_globalId) {
            Complex complex;
            complex.vertices.resize(DIM + 1);
            complex.global_verticesId = vertices_globalId;

            {
                complex.vertices[0] = {1, 2, 3};
                complex.vertices[1] = {2, 3, 0};
                complex.vertices[2] = {3, 0, 1};
                complex.vertices[3] = {0, 1, 2};

                complex.edges.resize(6);
                complex.edges[0].vertices = {0, 1};
                complex.edges[1].vertices = {0, 2};
                complex.edges[2].vertices = {0, 3};
                complex.edges[3].vertices = {1, 2};
                complex.edges[4].vertices = {1, 3};
                complex.edges[5].vertices = {2, 3};
                complex.edges[0].supporting_planes = {2, 3};
                complex.edges[1].supporting_planes = {1, 3};
                complex.edges[2].supporting_planes = {1, 2};
                complex.edges[3].supporting_planes = {0, 3};
                complex.edges[4].supporting_planes = {0, 2};
                complex.edges[5].supporting_planes = {0, 1};

                complex.global_oriEdgesId[0] = edges_globalId[0];
                complex.global_oriEdgesId[1] = edges_globalId[1];
                complex.global_oriEdgesId[2] = edges_globalId[2];
                complex.global_oriEdgesId[3] = edges_globalId[3];
                complex.global_oriEdgesId[4] = edges_globalId[4];
                complex.global_oriEdgesId[5] = edges_globalId[5];

                complex.faces.resize(4);
                complex.faces[0].edges = {5, 3, 4};
                complex.faces[1].edges = {2, 1, 5};
                complex.faces[2].edges = {4, 0, 2};
                complex.faces[3].edges = {1, 0, 3};
                complex.faces[0].supporting_plane = 0;
                complex.faces[1].supporting_plane = 1;
                complex.faces[2].supporting_plane = 2;
                complex.faces[3].supporting_plane = 3;
                complex.faces[0].positive_cell = 0;
                complex.faces[0].negative_cell = INVALID;
                complex.faces[1].positive_cell = 0;
                complex.faces[1].negative_cell = INVALID;
                complex.faces[2].positive_cell = 0;
                complex.faces[2].negative_cell = INVALID;
                complex.faces[3].positive_cell = 0;
                complex.faces[3].negative_cell = INVALID;

                complex.global_oriFacesId[0] = faces_globalId[0];
                complex.global_oriFacesId[1] = faces_globalId[1];
                complex.global_oriFacesId[2] = faces_globalId[2];
                complex.global_oriFacesId[3] = faces_globalId[3];
                complex.oriFacesEdges[0] = complex.faces[0].edges;
                complex.oriFacesEdges[1] = complex.faces[1].edges;
                complex.oriFacesEdges[2] = complex.faces[2].edges;
                complex.oriFacesEdges[3] = complex.faces[3].edges;

                complex.cells.resize(1);
                complex.cells[0].faces = {0, 1, 2, 3};
                {
                    if (tetLocalIdx == 0) {
                        complex.cells[0].signs = {false, false, false, false};
                    } else {
                        complex.cells[0].signs = {true, true, true, true};
                    }
                }
            }

            return complex;
        }

        Complex OffsetSurfaceBuilder::initialize_complex(
                int tetLocalIdx, size_t num_planes,
                const std::array<size_t, 4> &vertices_globalId,
                const std::array<size_t, 6> &edges_globalId,
                const std::array<size_t, 4> &faces_globalId) {
            Complex complex = initialize_complex(tetLocalIdx, vertices_globalId,
                                                 edges_globalId, faces_globalId);
            complex.cells[0].signs.resize(num_planes, false);

            return complex;
        }

    } // namespace complex
NAMESPACE_END(PCO)