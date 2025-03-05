//
// Created by Lei on 4/3/2024.
//

#include "SurfaceMesh.hpp"
#include "utils/Log.hpp"
#include "utils/File.hpp"
#include "cgal/CGALProcessor.hpp"

#include <fstream>

#ifdef ERROR
#   undef ERROR
#endif

NAMESPACE_BEGIN(PCO)
    namespace mesh {

        SurfaceMesh::SurfaceMesh(const std::vector<Vertex> &_vertVec, const std::vector<Polygon> &_faceVec) noexcept:
                vertVec(_vertVec), faceVec(_faceVec) {
            setInternalData();
        }

        void SurfaceMesh::clear() {
            vertVec.clear();
            vertVec.shrink_to_fit();

            for (auto &plg: faceVec)
                plg.clear();
            faceVec.clear();
            faceVec.shrink_to_fit();
        }

        void SurfaceMesh::addVertex(const Vertex &p) {
            vertVec.emplace_back(p);
        }

        void SurfaceMesh::addFace(const Polygon &vts) {
            faceVec.emplace_back(vts);
        }

        void SurfaceMesh::setInternalData() {
            numMeshVerts = vertVec.size();
            numMeshFaces = faceVec.size();

            faceNormalMat.resize(numMeshFaces, 3);
            faceNormalMat.setZero();
#pragma omp parallel for
            for (int i = 0; i < numMeshFaces; ++i) {
                MatrixX faceVertMat(3, 3);
                for (int j = 0; j < 3; ++j)
                    faceVertMat.row(j) = vertVec[faceVec[i][j]].cast<double>();

                Vector3 dir1 = faceVertMat.row(1) - faceVertMat.row(0);
                Vector3 dir2 = faceVertMat.row(2) - faceVertMat.row(0);
                Vector3 normal = dir1.cross(dir2).normalized();

                faceNormalMat.row(i) = normal;
            }

            vertNormalMat.resize(numMeshVerts, 3);
            vertNormalMat.setZero();
            for (size_t i = 0; i < numMeshFaces; ++i) {
                size_t numFaceVerts = faceVec[i].size();
                for (size_t j = 0; j < numFaceVerts; ++j) {
                    Vertex v0 = vertVec[faceVec[i][j]];
                    Vertex v1 = vertVec[faceVec[i][(j + 1) % numFaceVerts]];
                    Vertex v2 = vertVec[faceVec[i][(j + numFaceVerts - 1) % numFaceVerts]];

                    Vertex v10 = (v1 - v0).normalized();
                    Vertex v20 = (v2 - v0).normalized();
                    Scalar angel = 2.0 * atan((v10 - v20).norm() / (v10 + v20).norm());

                    vertNormalMat.row(faceVec[i][j]) += angel * faceNormalMat.row(i);
                }
            }
            vertNormalMat.rowwise().normalize();
        }

        void SurfaceMesh::removeDuplicatedVerts(int verbose, bool lazyTag) {
            if (verbose == 1) {
                LOG::INFO("[PCO][Mesh] Remove duplicated vertices.");
                LOG::INFO("[PCO][Mesh] Before:");
                printInfo();
            }

            // helper class for identifying coincident vertices
            class CompareVec {
            public:
                explicit CompareVec(float _eps = FLT_MIN) : eps_(_eps) {}

                bool operator()(const Vertex &v0, const Vertex &v1) const {
                    if (fabs(v0[0] - v1[0]) <= eps_) {
                        if (fabs(v0[1] - v1[1]) <= eps_)
                            return (v0[2] < v1[2] - eps_);
                        else
                            return (v0[1] < v1[1] - eps_);
                    } else
                        return (v0[0] < v1[0] - eps_);
                }

            private:
                float eps_;
            };

            CompareVec comp(FLT_MIN);
            std::map<Vertex, size_t, CompareVec> vMap(comp);
            std::map<Vertex, size_t, CompareVec>::iterator vMapIt;
            std::map<size_t, size_t> index_map;

            std::vector<Vertex> result_points = vertVec;
            std::vector<Polygon> input_polygons = faceVec;
            clear();
            faceVec.resize(numMeshFaces);

            for (size_t i = 0; i < result_points.size(); ++i) {
                const Vertex &p = result_points.at(i);

                // has the point been referenced before?
                vMapIt = vMap.find(p);
                if (vMapIt == vMap.end()) {
                    // No : add vertex and remember idx/vector mapping
                    vMap[p] = vertVec.size();
                    index_map[i] = vertVec.size();
                    vertVec.emplace_back(p);
                } else {
                    // Yes : get index from map
                    index_map[i] = vMapIt->second;
                }
            }
            for (size_t i = 0; i < input_polygons.size(); ++i) {
                for (auto id: input_polygons[i])
                    faceVec[i].push_back(index_map[id]);
            }

            if (!lazyTag) {
                setInternalData();
                if (verbose == 1) {
                    LOG::INFO("[PCO][Mesh] After:");
                    printInfo();
                }
            } else {
                if (verbose == 1) {
                    LOG::INFO("[PCO][Mesh] After(lazy):");
                    printInfo();
                }
            }
        }

        void SurfaceMesh::reverseOrientation(int verbose) {
            if (verbose == 1) {
                LOG::INFO("[PCO][Mesh] Reverse orientation.");
                LOG::INFO("[PCO][Mesh] Before:");
                printInfo();
            }

#pragma omp parallel for
            for (int i = 0; i < numMeshFaces; ++i) {
                std::reverse(faceVec[i].begin(), faceVec[i].end());
                faceNormalMat.row(i) *= -1.0;
            }
#pragma omp parallel for
            for (int i = 0; i < numMeshVerts; ++i) {
                vertNormalMat.row(i) *= -1.0;
            }

            if (verbose == 1) {
                LOG::INFO("[PCO][Mesh] After:");
                printInfo();
            }
        }

        void SurfaceMesh::triangulate(int verbose) {
            if (verbose == 1) {
                LOG::INFO("[PCO][Mesh] Triangulate polygon mesh.");
                LOG::INFO("[PCO][Mesh] Before:");
                printInfo();
            }

            ::PCO::cgal::triangulatePolygon(this);
            setInternalData();

            if (verbose == 1) {
                LOG::INFO("[PCO][Mesh] After:");
                printInfo();
            }
        }

        bool SurfaceMesh::checkSelfIntersection() {
            return ::PCO::cgal::isSelfIntersect(this);
        }

        size_t SurfaceMesh::removeDuplicateVertices() {
            return ::PCO::cgal::removeDuplicateVertices(this);
        }

        size_t SurfaceMesh::removeDuplicateFaces() {
            return ::PCO::cgal::removeDuplicateFaces(this);
        }

        bool SurfaceMesh::checkMesh() {
            ::PCO::cgal::nb_connectedComponents(this);
            return ::PCO::cgal::check(this);
        }

        void SurfaceMesh::repairMesh() {
            ::PCO::cgal::repair(this);
        }

        void SurfaceMesh::processOffsetMesh() {
            repairMesh();
            // checkMesh();
        }

        void SurfaceMesh::extractOffsetMesh() {

        }

        void SurfaceMesh::extractTriMesh(MatrixX &V, MatrixXi &F, MatrixX &VN, MatrixX &FN) {
            V.resize(vertVec.size(), 3);
            F.resize(faceVec.size(), 3);
            VN = vertNormalMat;
            FN = faceNormalMat;

#pragma omp parallel for
            for (int i = 0; i < vertVec.size(); ++i)
                V.row(i) = vertVec[i].cast<Scalar>();

#pragma omp parallel for
            for (int i = 0; i < faceVec.size(); ++i) {
                if (faceVec[i].size() < 3) continue;
                F.row(i) = Vector3i(faceVec[i][0], faceVec[i][1], faceVec[i][2]);
            }
        }

        void SurfaceMesh::printInfo() const {
            LOG::INFO("[PCO][Mesh] Surface Mesh Information");
            LOG::INFO("-- The number of vertices = {}.", numMeshVerts);
            LOG::INFO("-- The number of faces = {}.", numMeshFaces);
        }

        void SurfaceMesh::writeMesh(const std::string &outFile) const {
            utils::file::checkDir(outFile);
            std::ofstream out(outFile);
            if (!out) {
                LOG::ERROR("[PCO][I/O] File \"{}\" could not be opened!", outFile);
                return;
            }

            const std::string &fileExt = utils::file::getFileExtension(outFile);
            if (fileExt != ".obj") {
                LOG::ERROR("[PCO][I/O] Unsupported file format \"{}\"!", fileExt);
                out.close();
                return;
            }

            if (fileExt == ".obj") {
                // Write vertices
                for (const auto &vertex: vertVec) {
                    out << "v " << vertex.x() << " " << vertex.y() << " " << vertex.z() << std::endl;
                }

                // Write faces
                for (const auto &face: faceVec) {
                    if (face.empty()) continue;
                    std::stringstream facess;
                    facess << "f";
                    for (const auto &v: face) {
                        facess << " " << v + 1;
                    }
                    out << facess.str() << std::endl;
                }
            }

            out.close();
        }

    } // namespace mesh
NAMESPACE_END(PCO)
