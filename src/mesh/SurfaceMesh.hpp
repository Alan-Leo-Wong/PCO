//
// Created by Lei on 4/3/2024.
//

#ifndef PCO_SURFACEMESH_HPP
#define PCO_SURFACEMESH_HPP

#include "Config.hpp"
#include "detail/BasicDataType.hpp"

#include <vector>

NAMESPACE_BEGIN(PCO)
    namespace mesh {

        /**
         * Polygonal mesh/soup
         */
        class SurfaceMesh {
        public:
            typedef Vector3f Vertex;
            typedef std::vector<size_t> Polygon;

        public:
            SurfaceMesh() noexcept = default;

            SurfaceMesh(const std::vector<Vertex> &_vertVec, const std::vector<Polygon> &_faceVec) noexcept;

            ~SurfaceMesh() noexcept = default;

        public:
            /* Getter */
            [[nodiscard]] size_t numVertices() const noexcept { return vertVec.size(); }

            [[nodiscard]] size_t numFaces() const noexcept { return faceVec.size(); }

            std::vector<Vertex> &getVertVec() noexcept { return vertVec; }

            [[nodiscard]] const std::vector<Vertex> &getVertVec() const noexcept { return vertVec; }

            std::vector<Polygon> &getFaceVec() noexcept { return faceVec; }

            [[nodiscard]] const std::vector<Polygon> &getFaceVec() const noexcept { return faceVec; }

            MatrixX &getFaceNormalMat() noexcept { return faceNormalMat; }

            [[nodiscard]] const MatrixX &getFaceNormalMat() const noexcept { return faceNormalMat; }

            MatrixX &getVertNormalMat() noexcept { return vertNormalMat; }

            [[nodiscard]] const MatrixX &getVertNormalMat() const noexcept { return vertNormalMat; }

            /* Setter */
            void setVertVec(const std::vector<Vertex> &_vertVec) noexcept { vertVec = _vertVec; }

            void setFaceVec(const std::vector<Polygon> &_faceVec) noexcept { faceVec = _faceVec; }

        private:
            void setInternalData();

        public:
            void clear();

            void addVertex(const Vertex &p);

            void addFace(const Polygon &vts);

        public:
            void removeDuplicatedVerts(int verbose = 0, bool lazyTag = false);

            void reverseOrientation(int verbose = 0);

            void triangulate(int verbose = 0);

            bool checkSelfIntersection();

            size_t removeDuplicateVertices();

            size_t removeDuplicateFaces();

            bool checkMesh();

            void repairMesh();

            void processOffsetMesh();

            void extractOffsetMesh();

            void extractTriMesh(MatrixX &V, MatrixXi &F, MatrixX &VN, MatrixX &FN);

        public:
            /* Low level apis */
            void printInfo() const;

        public:
            /* I/O */
            /* Output */
            void writeMesh(const std::string &outFile) const;

        protected:
            size_t numMeshVerts = 0;
            size_t numMeshFaces = 0;

        private:
            /* Pointer(point to address of data in MemoryPool) is stored in vector */
            std::vector<Vertex> vertVec;
            std::vector<Polygon> faceVec;

            MatrixX faceNormalMat;
            MatrixX vertNormalMat;
        };

    } // namespace mesh
NAMESPACE_END(PCO)

#endif //PCO_SURFACEMESH_HPP
