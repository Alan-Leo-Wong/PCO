//
// Created by Lei on 2/9/2024.
//
#ifndef PCO_MESHVIEWER_HPP
#define PCO_MESHVIEWER_HPP

#include "Config.hpp"
#include "mesh/TriMesh.hpp"
#include "mesh/SurfaceMesh.hpp"

#include <polyscope/polyscope.h>
#include <polyscope/camera_view.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>
#include <polyscope/simple_triangle_mesh.h>
#include <polyscope/combining_hash_functions.h>

NAMESPACE_BEGIN(PCO)
    namespace viewer {
        using namespace mesh;
        using namespace polyscope;

        struct MeshData {
            std::string meshName, curveName;

            size_t nVertices = 0;
            size_t nFaces = 0;
            size_t nEdges = 0;

            std::vector<glm::vec3> vertexPositionsGLM;
            std::vector<std::vector<size_t>> faceIndices;

            std::unordered_set<std::pair<size_t, size_t>,
                    polyscope::hash_combine::hash<std::pair<size_t, size_t>>> seenEdges;
            std::vector<std::array<size_t, 2>> edges;

            std::vector<glm::vec3> fNormals;
            std::vector<glm::vec3> fCenters;
            std::vector<glm::vec3> vNormals;

            /// surface mesh quantities
            std::map<std::string, SurfaceVertexScalarQuantity *> vertexScalarQuantities;
            std::map<std::string, SurfaceVertexColorQuantity *> vertexColorQuantities;

            std::map<std::string, SurfaceFaceScalarQuantity *> faceScalarQuantities;
            std::map<std::string, SurfaceFaceColorQuantity *> faceColorQuantities;

            std::map<std::string, SurfaceEdgeScalarQuantity *> edgeScalarQuantities;
            std::map<std::string, SurfaceVertexVectorQuantity *> vertexVectorQuantities;
            std::map<std::string, SurfaceVertexTangentVectorQuantity *> vertexTangentVectorQuantities;
            std::map<std::string, SurfaceVertexParameterizationQuantity *> vertexParameterizationQuantities;

            std::map<std::string, SurfaceFaceVectorQuantity *> faceVectorQuantities;
            std::map<std::string, SurfaceFaceTangentVectorQuantity *> faceTangentVectorQuantities;

            std::map<std::string, SurfaceHalfedgeScalarQuantity *> halfEdgeScalarQuantities;
            std::map<std::string, SurfaceCornerScalarQuantity *> cornerScalarQuantities;

            std::map<std::string, SurfaceOneFormTangentVectorQuantity *> oneFormTangentVectorQuantities;

            /// curve network quantities
            std::map<std::string, CurveNetworkNodeScalarQuantity *> curveNetworkNodeScalarQuantities;
            std::map<std::string, CurveNetworkNodeColorQuantity *> curveNetworkNodeColorQuantities;
            std::map<std::string, CurveNetworkNodeVectorQuantity *> curveNetworkNodeVectorQuantities;

            std::map<std::string, CurveNetworkEdgeScalarQuantity *> curveNetworkEdgeScalarQuantities;
            std::map<std::string, CurveNetworkEdgeColorQuantity *> curveNetworkEdgeColorQuantities;
            std::map<std::string, CurveNetworkEdgeVectorQuantity *> curveNetworkEdgeVectorQuantities;
        };

        class MeshViewer {
        public:
            typedef polyscope::SurfaceMesh pSurfaceMesh;
            typedef polyscope::SimpleTriangleMesh pSimpleTriMesh;
            typedef polyscope::CurveNetwork pCurveNetwork;

        public:
            MeshViewer(const std::string &_meshName, SurfaceMesh *_surfaceMesh);

            MeshViewer(const std::string &_triMeshName, TriMesh *_triMesh,
                       const std::string &_surfaceMeshName, SurfaceMesh *_surfaceMesh);

        private:
            void render();
        };

    } // namespace viewer
NAMESPACE_END(PCO)

#endif // PCO_MESHVIEWER_HPP
