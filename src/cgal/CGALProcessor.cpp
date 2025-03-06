//
// Created by Lei on 1/6/2024.
//

#include "CGALProcessor.hpp"
#include "utils/Log.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "detail/BasicDataType.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

NAMESPACE_BEGIN(PCO)
    namespace cgal {
        namespace PMP = CGAL::Polygon_mesh_processing;

        typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
        typedef CGAL::Surface_mesh<K::Point_3> CGALMesh;

        namespace internal {

            void to_cgal(const SurfaceMesh *input, CGALMesh &output) {
                output.clear();
                for (const SurfaceMesh::Vertex &v: input->getVertVec()) {
                    output.add_vertex(K::Point_3(v.x(), v.y(), v.z()));
                }

                for (const auto &f: input->getFaceVec()) {
                    std::vector<CGALMesh::Vertex_index> vts;
                    for (size_t v: f)
                        vts.emplace_back(v);
                    output.add_face(vts);
                }
            }

            void to_ours(const CGALMesh &input, SurfaceMesh *output) {
                output->clear();
                using ValType = SurfaceMesh::Vertex::Scalar;
                for (auto v: input.vertices()) {
                    const K::Point_3 &p = input.point(v);
                    output->addVertex(SurfaceMesh::Vertex(
                            static_cast<ValType>(CGAL::to_double(p.x())),
                            static_cast<ValType>(CGAL::to_double(p.y())),
                            static_cast<ValType>(CGAL::to_double(p.z()))
                    ));
                }

                for (const auto &f: input.faces()) {
                    SurfaceMesh::Polygon vts;
                    auto h = input.halfedge(f);
                    for (const auto &v: input.vertices_around_face(h))
                        vts.emplace_back(v.idx());
                    output->addFace(vts);
                }
            }

            void to_cgal(const SurfaceMesh *input,
                         std::vector<K::Point_3> &points,
                         std::vector<std::vector<std::size_t>> &polygons) {
                const std::vector<SurfaceMesh::Vertex> &input_points = input->getVertVec();
                const std::vector<SurfaceMesh::Polygon> &input_polygons = input->getFaceVec();
                points.resize(input_points.size());
                polygons.resize(input_polygons.size());

                for (size_t i = 0; i < input_points.size(); ++i) {
                    const SurfaceMesh::Vertex &p = input_points[i];
                    points[i] = K::Point_3(p.x(), p.y(), p.z());
                }

                for (size_t i = 0; i < input_polygons.size(); ++i) {
                    const auto &input_plg = input_polygons[i];
                    auto &plg = polygons[i];
                    plg.resize(input_plg.size());
                    for (size_t j = 0; j < input_plg.size(); ++j)
                        plg[j] = input_plg[j];
                }
            }

            void to_ours(const std::vector<K::Point_3> &input_points,
                         const std::vector<std::vector<std::size_t> > &input_polygons,
                         SurfaceMesh *output) {
                output->clear();
                std::vector<SurfaceMesh::Vertex> &points = output->getVertVec();
                std::vector<SurfaceMesh::Polygon> &polygons = output->getFaceVec();
                points.resize(input_points.size());
                polygons.resize(input_polygons.size());

                using ValType = SurfaceMesh::Vertex::Scalar;
                for (size_t i = 0; i < input_points.size(); ++i) {
                    const auto &p = input_points[i];
                    points[i] = SurfaceMesh::Vertex(
                            static_cast<ValType>(CGAL::to_double(p.x())),
                            static_cast<ValType>(CGAL::to_double(p.y())),
                            static_cast<ValType>(CGAL::to_double(p.z()))
                    );
                }

                for (size_t i = 0; i < input_polygons.size(); ++i) {
                    const auto &input_plg = input_polygons[i];
                    auto &plg = polygons[i];
                    plg.resize(input_plg.size());
                    for (size_t j = 0; j < input_plg.size(); ++j)
                        plg[j] = input_plg[j];
                }
            }

            void to_polygon_mesh(const std::vector<SurfaceMesh::Vertex> &points,
                                 const std::vector<SurfaceMesh::Polygon> &polygons,
                                 SurfaceMesh *mesh) {
                mesh->clear();
                for (const SurfaceMesh::Vertex &p: points)
                    mesh->addVertex(p);
                for (const SurfaceMesh::Polygon &plg: polygons) {
                    SurfaceMesh::Polygon vts;
                    for (size_t v: plg)
                        vts.emplace_back(v);
                    mesh->addFace(vts);
                }
            }

        } // namespace internal

        bool isVertexManifold(SurfaceMesh *surfaceMesh) {
            typedef boost::graph_traits<CGALMesh>::vertex_descriptor vertex_descriptor;
            CGALMesh cgalMesh;
            internal::to_cgal(surfaceMesh, cgalMesh);

            // Count non manifold vertices
            int counter = 0;
            for (vertex_descriptor v: vertices(cgalMesh)) {
                if (PMP::is_non_manifold_vertex(v, cgalMesh)) {
                    ++counter;
                }
            }
            return (counter != 0);
        }

        bool isSelfIntersect(SurfaceMesh *surfaceMesh) {
            typedef boost::graph_traits<CGALMesh>::face_descriptor face_descriptor;

            CGALMesh cgalMesh;
            internal::to_cgal(surfaceMesh, cgalMesh);

            bool flag = PMP::does_self_intersect<CGAL::Parallel_if_available_tag>(cgalMesh,
                                                                                  CGAL::parameters::vertex_point_map(
                                                                                          get(CGAL::vertex_point,
                                                                                              cgalMesh)));

            return flag;
        }

        bool orientPolygonSoup(SurfaceMesh *surfaceMesh) {
            // stitch
            std::vector<K::Point_3> points;
            std::vector<std::vector<std::size_t>> polygons(surfaceMesh->numFaces());
            internal::to_cgal(surfaceMesh, points, polygons);

            bool status = PMP::orient_polygon_soup(points, polygons);
            if (!PMP::is_polygon_soup_a_polygon_mesh(polygons)) {
                LOG::WARN("[CGAL] The polygons after orientation do not define a valid polygon mesh.");
                return false;
            }

            internal::to_ours(points, polygons, surfaceMesh);

            return status;
        }

        size_t removeDuplicateVertices(SurfaceMesh *surfaceMesh) {
            std::vector<K::Point_3> points;
            std::vector<std::vector<std::size_t> > polygons(surfaceMesh->numFaces());
            internal::to_cgal(surfaceMesh, points, polygons);

            size_t nb_v = PMP::merge_duplicate_points_in_polygon_soup(points, polygons);

            internal::to_ours(points, polygons, surfaceMesh);
            return nb_v;
        }

        size_t removeDuplicateFaces(SurfaceMesh *surfaceMesh) {
            std::vector<K::Point_3> points;
            std::vector<std::vector<std::size_t> > polygons(surfaceMesh->numFaces());
            internal::to_cgal(surfaceMesh, points, polygons);

            size_t nb_p = PMP::merge_duplicate_polygons_in_polygon_soup(points, polygons);

            internal::to_ours(points, polygons, surfaceMesh);
            return nb_p;
        }

        bool triangulatePolygon(SurfaceMesh *surfaceMesh) {
            CGALMesh cgalMesh;
            internal::to_cgal(surfaceMesh, cgalMesh);

            bool is_success = PMP::triangulate_faces(cgalMesh);
            internal::to_ours(cgalMesh, surfaceMesh);

            return is_success;
        }

        void nb_connectedComponents(SurfaceMesh *surfaceMesh) {
            namespace params = CGAL::parameters;

            CGALMesh cgalMesh;
            internal::to_cgal(surfaceMesh, cgalMesh);

            CGALMesh::Property_map <CGALMesh::Face_index, std::size_t> vol_id_map =
                    cgalMesh.add_property_map<CGALMesh::Face_index, std::size_t>().first;

            std::vector<PMP::Volume_error_code> err_codes;
            std::size_t nb_vol = PMP::volume_connected_components(cgalMesh, vol_id_map,
                                                                  params::error_codes(std::ref(err_codes)));

            // std::cout << "nb_vol: " << nb_vol << std::endl;
        }

        bool check(SurfaceMesh *surfaceMesh) {
            typedef boost::graph_traits<CGALMesh>::vertex_descriptor vertex_descriptor;
            typedef boost::graph_traits<CGALMesh>::face_descriptor face_descriptor;

            CGALMesh cgalMesh;
            internal::to_cgal(surfaceMesh, cgalMesh);

            // Count non manifold vertices
            int counter = 0;
            for (vertex_descriptor v: vertices(cgalMesh)) {
                if (PMP::is_non_manifold_vertex(v, cgalMesh)) {
                    ++counter;
                }
            }
            bool isVertexManifold = (counter == 0);
            bool isSelfIntersection = PMP::does_self_intersect<CGAL::Parallel_if_available_tag>(cgalMesh,
                                                                                                CGAL::parameters::vertex_point_map(
                                                                                                        get(CGAL::vertex_point,
                                                                                                            cgalMesh)));
            {
                LOG::SPLIT();
                LOG::INFO("Mesh correctness check:");
                LOG::INFO("-- Manifold? {}", isVertexManifold ? "yes" : "no");
                LOG::INFO("-- Self-intersecting? {}", isSelfIntersection ? "yes" : "no");
                if(isSelfIntersection){
                    std::vector<std::pair<face_descriptor, face_descriptor> > intersected_tris;
                    PMP::self_intersections<CGAL::Parallel_if_available_tag>(faces(cgalMesh), cgalMesh,
                                                                             std::back_inserter(intersected_tris));

                    CGALMesh intersectedMesh;
                    std::unordered_map<vertex_descriptor, vertex_descriptor> vertex_map;

                    // Copy only the intersected faces to the new mesh
                    for (const auto& face_pair : intersected_tris) {
                        const face_descriptor& face1 = face_pair.first;
                        const face_descriptor& face2 = face_pair.second;

                        // Extract the vertices of face1
                        std::vector<vertex_descriptor> face1_vertices;
                        for (auto v : CGAL::vertices_around_face(halfedge(face1, cgalMesh), cgalMesh)) {
                            face1_vertices.push_back(v);
                        }

                        // Extract the vertices of face2
                        std::vector<vertex_descriptor> face2_vertices;
                        for (auto v : CGAL::vertices_around_face(halfedge(face2, cgalMesh), cgalMesh)) {
                            face2_vertices.push_back(v);
                        }

                        // Add vertices and faces to the new mesh
                        std::vector<vertex_descriptor> new_face1_vertices;
                        for (auto v : face1_vertices) {
                            // Check if vertex already exists in the new mesh
                            if (vertex_map.find(v) == vertex_map.end()) {
                                // Add new vertex to the new mesh
                                new_face1_vertices.push_back(intersectedMesh.add_vertex(cgalMesh.point(v)));
                                vertex_map[v] = new_face1_vertices.back();
                            } else {
                                new_face1_vertices.push_back(vertex_map[v]);
                            }
                        }
                        intersectedMesh.add_face(new_face1_vertices);

                        std::vector<vertex_descriptor> new_face2_vertices;
                        for (auto v : face2_vertices) {
                            // Check if vertex already exists in the new mesh
                            if (vertex_map.find(v) == vertex_map.end()) {
                                // Add new vertex to the new mesh
                                new_face2_vertices.push_back(intersectedMesh.add_vertex(cgalMesh.point(v)));
                                vertex_map[v] = new_face2_vertices.back();
                            } else {
                                new_face2_vertices.push_back(vertex_map[v]);
                            }
                        }
                        intersectedMesh.add_face(new_face2_vertices);
                    }

                    // Write the intersected mesh to an OBJ file
                    if (!CGAL::IO::write_polygon_mesh(R"(F:\PCO\cube_2_intersected_faces.obj)", intersectedMesh)) {
                        std::cerr << "Failed to write OBJ file." << std::endl;
                        return false;
                    }
                }
                LOG::SPLIT();
            }

            return (isVertexManifold && isSelfIntersection);
        }

        void repair(SurfaceMesh *surfaceMesh) {
            std::vector<K::Point_3> points;
            std::vector<std::vector<std::size_t> > polygons(surfaceMesh->numFaces());
            internal::to_cgal(surfaceMesh, points, polygons);

            // std::cout << PMP::merge_duplicate_points_in_polygon_soup(points, polygons) << std::endl;
            PMP::merge_duplicate_polygons_in_polygon_soup(points, polygons);

            CGALMesh cgalMesh;
            PMP::polygon_soup_to_polygon_mesh(points, polygons, cgalMesh);
            //PMP::triangulate_faces(cgalMesh);

            internal::to_ours(cgalMesh, surfaceMesh);
        }

    } // namespace cgal
NAMESPACE_END(PCO)
