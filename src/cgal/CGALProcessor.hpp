//
// Created by Lei on 1/6/2024.
//

#ifndef LINEARSHARPOFFSET_CGALPROCESSOR_HPP
#define LINEARSHARPOFFSET_CGALPROCESSOR_HPP

#include "Config.hpp"

NAMESPACE_BEGIN(PCO)
    namespace mesh {
        // forward declarations
        class SurfaceMesh;
    }
    namespace cgal {
        using ::PCO::mesh::SurfaceMesh;

        bool isSelfIntersect(SurfaceMesh *surfaceMesh);

        bool isVertexManifold(SurfaceMesh *surfaceMesh);

        bool orientPolygonSoup(SurfaceMesh *surfaceMesh);

        size_t removeDuplicateVertices(SurfaceMesh *surfaceMesh);

        size_t removeDuplicateFaces(SurfaceMesh *surfaceMesh);

        bool triangulatePolygon(SurfaceMesh *surfaceMesh);

        void nb_connectedComponents(SurfaceMesh *surfaceMesh);

        bool check(SurfaceMesh *surfaceMesh);

        void repair(SurfaceMesh *surfaceMesh);

    } // namespace cgal
NAMESPACE_END(PCO)

#endif //LINEARSHARPOFFSET_CGALPROCESSOR_HPP
