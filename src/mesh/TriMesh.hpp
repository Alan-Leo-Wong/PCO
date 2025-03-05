#pragma once

#include "detail/Geometry.hpp"

#include <string>
#include <vector>
#include <igl/AABB.h>
#include <igl/edges.h>
#include <igl/fast_winding_number.h>

#include <CGAL/Simple_cartesian.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Ray_3.h>

NAMESPACE_BEGIN(PCO)
    namespace mesh {
        using namespace geometry;

        /**
         * Triangle mesh/soup
         */
        class TriMesh {
        public:
            typedef TriMesh ThisClass;
            typedef CGAL::Simple_cartesian<double> K;
//            typedef CGAL::Exact_predicates_exact_constructions_kernel K;
            typedef K::FT FT;
            typedef K::Point_3 Point;
            typedef K::Ray_3 Ray;
            typedef K::Segment_3 Segment;
            typedef K::Triangle_3 Triangle;
            typedef K::Iso_cuboid_3 Iso_cuboid;
            typedef std::vector<Triangle>::iterator Iterator;
            typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
            typedef CGAL::AABB_traits_3<K, Primitive> AABB_triangle_traits;
            typedef CGAL::AABB_tree<AABB_triangle_traits> CGAL_AABB;

        public:
            /* Constructor and Destructor */
            TriMesh() noexcept = default;

            TriMesh(MatrixX verts, MatrixXi faces);

            explicit TriMesh(const std::string &file, bool lazyTag = false);

            virtual ~TriMesh() noexcept = default;

        public:
            /* Getter */
            [[nodiscard]] size_t numVertices() const noexcept { return vertMat.rows(); }

            [[nodiscard]] size_t numFaces() const noexcept { return faceMat.rows(); }

            [[nodiscard]] size_t numEdges() const noexcept { return meshEdges.rows(); }

            MatrixX &getVertMat() noexcept { return vertMat; }

            [[nodiscard]] const MatrixX &getVertMat() const noexcept { return vertMat; }

            std::vector<MatrixX> &getVertMats() noexcept { return faceVertMats; }

            [[nodiscard]] const std::vector<MatrixX> &getVertMats() const noexcept { return faceVertMats; }

            MatrixXi &getFaceMat() noexcept { return faceMat; }

            [[nodiscard]] const MatrixXi &getFaceMat() const noexcept { return faceMat; }

            MatrixX &getVertNormalMat() noexcept { return vertNormalMat; }

            [[nodiscard]] const MatrixX &getVertNormalMat() const noexcept { return vertNormalMat; }

            MatrixX &getFaceNormalMat() noexcept { return faceNormalMat; }

            [[nodiscard]] const MatrixX &getFaceNormalMat() const noexcept { return faceNormalMat; }

            AABox<Vector3> getModelBB() noexcept { return modelBoundingBox; }

            [[nodiscard]] const AABox<Vector3> &getModelBB() const noexcept { return modelBoundingBox; }

        private:
            void setInternalData(bool vertChanged = false);

            /* Set (uniform)bounding-box */
            template<bool is2Uniform>
            void setBoundingBox() {
                Vector3 minV = vertMat.colwise().minCoeff();
                Vector3 maxV = vertMat.colwise().maxCoeff();

                if constexpr (is2Uniform) {
                    modelBoundingBox = AABox<Vector3>(minV, maxV); // initialize answer
                    Vector3 lengths = maxV - minV; // check length of given bbox in every direction
                    const double max_length = fmaxf(lengths.x(), fmaxf(lengths.y(), lengths.z())); // find max length
                    for (unsigned int i = 0; i < 3; i++) { // for every direction (X,Y,Z)
                        if (max_length == lengths[i]) {
                            continue;
                        } else {
                            const double delta = max_length -
                                                 lengths[i]; // compute difference between largest length and current (X,Y or Z) length
                            modelBoundingBox.boxOrigin[i] =
                                    minV[i] - (delta / 2.0f); // pad with half the difference before current min
                            modelBoundingBox.boxEnd[i] =
                                    maxV[i] + (delta / 2.0f); // pad with half the difference behind current max
                        }
                    }

                    // Next snippet adresses the problem reported here: https://github.com/Forceflow/cuda_voxelizer/issues/7
                    // Suspected cause: If a triangle is axis-aligned and lies perfectly on a voxel edge, it sometimes gets counted / not counted
                    // Probably due to a numerical instability (division by zero?)
                    // Ugly fix: we pad the bounding box on all sides by 1/10001th of its total length, bringing all triangles ever so slightly off-grid
                    Vector3 epsilon = (modelBoundingBox.boxEnd - modelBoundingBox.boxOrigin) / 10001; // 之前是10001
                    modelBoundingBox.boxOrigin -= epsilon;
                    modelBoundingBox.boxEnd += epsilon;
                    modelBoundingBox.boxWidth = modelBoundingBox.boxEnd - modelBoundingBox.boxOrigin;
                } else {
                    modelBoundingBox = AABox(minV, maxV);
                }

                // diagonal length of the bounding-box
                diagonalLengthOfBBox = (modelBoundingBox.boxEnd - modelBoundingBox.boxOrigin).norm();
            }

        public:
            /* Query udf/sdf */
            Scalar pointUDF(const Vector3 &p) const;

            Scalar pointUDF(const Vector3 &p, int &closestFace, Eigen::RowVector3d &closestPoint) const;

            VectorX pointsUDF(const MatrixX &pts) const;

            VectorX pointsUDF(const MatrixX &pts, VectorXi &closestFace, MatrixX &closestPoint) const;

            int pointSign(const Vector3& p);

            Scalar pointPseudonormalSDF(const Vector3 &p);

            Scalar pointWnSDF(const Vector3 &p);

        public:
            /* Low level apis */
            void printInfo() const;

        public:
            /* I/O */
            /* Input */
            void readMesh(const std::string &in_file);

            /* Output */
            void writeMesh(const std::string &filename) const;

        protected:
            size_t numMeshVerts = 0;
            size_t numMeshEdges = 0;
            size_t numMeshFaces = 0;

            double avg_edge_length = 0.0;
            double max_edge_length = -1.0;

        protected:
            /* Basic data in matrix style(from libigl) */
            MatrixX vertMat;
            MatrixXi faceMat;
            std::vector<MatrixX> faceVertMats;

            MatrixXi meshEdges;
            VectorXi meshEdgeMAP;
            MatrixX vertNormalMat;
            MatrixX faceNormalMat;
            MatrixX edgeNormalMat;

        protected:
            MatrixXi TTAdj; // the element i,j is the id of the triangle adjacent to the j edge of triangle i
            // the first edge of a triangle is [0,1] the second [1,2] and the third [2,3].

            std::vector<std::vector<int>> VTAdj; // stores the index of the incident face with each point.
            std::vector<std::vector<int>> VTiAdj; // stores the index of vertices corresponding to each point within its associated face.

        protected:
            igl::AABB<MatrixX, 3> meshAABB;
            std::vector<Triangle> cgal_triangles;
            CGAL_AABB cgal_aabb;
            igl::FastWindingNumberBVH fwn_bvh;

        protected:
            /* Axis-aligned bounding box */
            double boxScaleFactor;
            double diagonalLengthOfBBox;
            AABox<Vector3> modelBoundingBox;

        public:
            std::string modelName;
            std::string textureName;
        };

    } // namespace mesh
NAMESPACE_END(PCO)