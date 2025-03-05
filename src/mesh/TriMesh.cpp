//
// Created by lei on 23-9-13.
//

#include "TriMesh.hpp"
#include "utils/Log.hpp"
#include "utils/Donut.hpp"

#include <omp.h>
#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>
#include <igl/per_edge_normals.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_vertex_normals.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/triangle_triangle_adjacency.h>

#ifdef ERROR
#   undef ERROR
#endif

NAMESPACE_BEGIN(PCO)
    namespace mesh {
        namespace internal {

            double signedTetVolumeOfTri(const Vector3 &triNormal,
                                        const Vector3 &triCorner,
                                        const Vector3 &p) {
                return triNormal.dot(p - triCorner);
            }

        } // namespace internal

        //////////////////////
        //   Constructors   //
        //////////////////////
        TriMesh::TriMesh(const std::string &file, bool lazyTag) {
            readMesh(file);

            if (!lazyTag) setInternalData();
        }

        TriMesh::TriMesh(MatrixX verts, MatrixXi faces) : vertMat(std::move(verts)),
                                                          faceMat(std::move(faces)) {
            setInternalData();
        }

        //////////////////////
        //   Internal data  //
        //////////////////////
        void TriMesh::setInternalData(bool vertChanged) {
            if (!vertChanged) {
                numMeshVerts = vertMat.rows();
                numMeshFaces = faceMat.rows();
            }

            // adjacency
            igl::triangle_triangle_adjacency(faceMat, TTAdj);
            igl::vertex_triangle_adjacency(vertMat, faceMat, VTAdj, VTiAdj);

            setBoundingBox<false>();

#ifndef NDEBUG
            modelBoundingBox.output(utils::file::concatFilePath(VIS_DIR, modelName,
                                                                "modelBoundingBox.obj"));
#endif
            LOG::INFO("[PCO] Diagonal length of the bounding-box = {}.", diagonalLengthOfBBox);

            // set faces' normal
            faceVertMats.resize(numMeshFaces);
            faceNormalMat.resize(numMeshFaces, 3);
            for (size_t i = 0; i < numMeshFaces; ++i) {
                MatrixX faceVertMat(3, 3);
                for (int j = 0; j < 3; ++j)
                    faceVertMat.row(j) = vertMat.row(faceMat(i, j));
                faceVertMats[i] = faceVertMat;

                Vector3 dir1 = faceVertMat.row(1) - faceVertMat.row(0);
                Vector3 dir2 = faceVertMat.row(2) - faceVertMat.row(0);
                Vector3 normal = dir1.cross(dir2).normalized();
                faceNormalMat.row(i) = normal;
            }

            // set vertices' normal
            igl::per_vertex_normals(vertMat, faceMat,
                                    igl::PerVertexNormalsWeightingType::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE,
                                    faceNormalMat, vertNormalMat);

            // set edges
            igl::per_edge_normals(vertMat, faceMat, igl::PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM, faceNormalMat,
                                  edgeNormalMat, meshEdges, meshEdgeMAP);

            // igl::edges(faceMat, meshEdges);
            numMeshEdges = meshEdges.rows();
            for (int i = 0; i < numMeshEdges; ++i) {
                Vector3 v0 = vertMat.row(meshEdges.row(i)(0));
                Vector3 v1 = vertMat.row(meshEdges.row(i)(1));
                double edge_length = (v1 - v0).norm();
                max_edge_length = std::max(max_edge_length, edge_length);
                avg_edge_length += edge_length;
            }
            avg_edge_length /= numMeshEdges;
            LOG::DEBUG("[PCO] Average edge length = ", avg_edge_length);
            LOG::DEBUG("[PCO] Maximum edge length = ", max_edge_length);

            meshAABB.init(vertMat, faceMat);

            cgal_triangles.resize(numMeshFaces);
#pragma omp parallel for
            for (int i = 0; i < numMeshFaces; ++i) {
                MatrixX faceVerts = this->faceVertMats[i];
                Point p0 = Point(faceVerts(0, 0), faceVerts(0, 1), faceVerts(0, 2));
                Point p1 = Point(faceVerts(1, 0), faceVerts(1, 1), faceVerts(1, 2));
                Point p2 = Point(faceVerts(2, 0), faceVerts(2, 1), faceVerts(2, 2));
                cgal_triangles[i] = Triangle(p0, p1, p2);
            }
            cgal_aabb = CGAL_AABB(cgal_triangles.begin(), cgal_triangles.end());
        }

        //////////////////////
        //     Utilities    //
        //////////////////////
        Scalar TriMesh::pointUDF(const Vector3 &p) const {
            int closestFace;
            Eigen::RowVector3d closestPoint;
            return std::sqrt(meshAABB.squared_distance(vertMat, faceMat, p, closestFace, closestPoint));
        }

        Scalar TriMesh::pointUDF(const Vector3 &p, int &closestFace, Eigen::RowVector3d &closestPoint) const {
            return std::sqrt(meshAABB.squared_distance(vertMat, faceMat, p, closestFace, closestPoint));
        }

        VectorX TriMesh::pointsUDF(const MatrixX &pts) const {
            VectorX res;
            VectorXi closestFace;
            MatrixX closestPoint;
            meshAABB.squared_distance(vertMat, faceMat, pts, res, closestFace, closestPoint);
            return res.array().sqrt();
        }

        VectorX TriMesh::pointsUDF(const MatrixX &pts, VectorXi &closestFace, MatrixX &closestPoint) const {
            VectorX res;
            meshAABB.squared_distance(vertMat, faceMat, pts, res, closestFace, closestPoint);
            return res.array().sqrt();
        }

        int TriMesh::pointSign(const Vector3 &p) {
            Point _p = Point(p.x(), p.y(), p.z());
            Ray ray(_p, K::Vector_3(1.0, 0.0, 0.0));

            bool is_inside = (cgal_aabb.number_of_intersected_primitives(ray) % 2 == 1);

            return (is_inside ? -1 : 1);
        }

        //////////////////////
        //  Low level APIs  //
        //////////////////////
        void TriMesh::printInfo() const {
            LOG::INFO("[PCO] Model: ", modelName);
            LOG::INFO("-- The number of vertices = {}.", numMeshVerts);
            LOG::INFO("-- The number of faces = {}.", numMeshFaces);
            LOG::INFO("-- The number of edges = {}.", numMeshEdges);
        }

        //////////////////////
        //       I / O      //
        //////////////////////
        /* Input */
        void TriMesh::readMesh(const std::string &in_file) {
            modelName = utils::file::getFileName(in_file);
            igl::read_triangle_mesh(in_file, vertMat, faceMat);
        }

        /* Output */
        void TriMesh::writeMesh(const std::string &filename) const {
            const std::string extension = utils::file::getFileExtension(filename);
            if (extension != ".obj" &&
                extension != ".off") {
                LOG::ERROR("[PCO] Unsupported file format \"{}\"!", extension);
                return;
            }
            if (extension == ".obj")
                igl::writeOBJ(filename, vertMat, faceMat);
            if (extension == ".off")
                igl::writeOFF(filename, vertMat, faceMat);
        }


    } // namespace mesh
NAMESPACE_END(PCO)