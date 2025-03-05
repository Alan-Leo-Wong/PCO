//
// Created by Lei on 11/3/2024.
//

#ifndef PCO_OCTREETETLUT_HPP
#define PCO_OCTREETETLUT_HPP

#include "BasicDataType.hpp"

NAMESPACE_BEGIN(PCO)
    namespace geometry {

        using TetVertex = int;
        using TetGlobalVertex = uint32_t;

        using TetEdge = std::array<int, 2>;
        using TetGlobalEdge = std::array<uint32_t, 2>;

        using TetFace = std::array<int, 3>;
        using TetGlobalFace = std::array<uint32_t, 3>;

        using TetNeighbor = std::array<int, 2>; // {node's neighbor, local tet index in a node}
        using TetGlobalNeighbor = std::array<uint32_t, 2>; // {node's neighbor, local tet index in a node}

        using TetEdgeNeighbor = std::array<TetNeighbor, 6>; // {node's neighbor, local tet index in a node}
        using TetEdgeGlobalNeighbor = std::array<uint32_t, 6>; // {node's neighbor, local tet index in a node}

        using TetFaceNeighbor = TetNeighbor; // {node's neighbor, local tet index in a node}
        using TetFaceGlobalNeighbor = uint32_t; // {node's neighbor, local tet index in a node}

        ////// The following look-up tables are used to construct tetrahedra important data within a node
        // child idx -> tet type
        constexpr static int tetsTypeInNodeLUT[8] =
                {
                        0, 1,
                        1, 0,
                        1, 0,
                        0, 1
                };

        constexpr static int tetLocalIdxLUT[2][5] =
                {
                        {1, 1, 1, 0, 0},
                        {0, 1, 1, 0, 0}
                };

        // each element stores the indices of three local tetrahedron vertices.
        // type - tet - face
        constexpr static TetFace tetLocalFaceVertexLUT[2][5][4] =
                {
                        // type 0
                        {
                                {{2, 3, 1}, {2, 0, 3}, {0, 1, 3}, {2, 1, 0}}, // tet 0
                                {{3, 1, 2}, {3, 2, 0}, {3, 0, 1}, {0, 2, 1}}, // tet 1
                                {{1, 2, 3}, {2, 0, 3}, {1, 3, 0}, {1, 0, 2}}, // tet 2
                                {{2, 3, 1}, {2, 0, 3}, {0, 1, 3}, {2, 1, 0}}, // tet 3
                                {{2, 1, 3}, {0, 2, 3}, {0, 3, 1}, {0, 1, 2}}  // tet 4
                        },

                        // type 1
                        {
                                {{2, 1, 3}, {0, 2, 3}, {0, 3, 1}, {0, 1, 2}}, // tet 0
                                {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}}, // tet 1
                                {{3, 1, 2}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}}, // tet 2
                                {{3, 2, 1}, {0, 2, 3}, {0, 3, 1}, {0, 1, 2}}, // tet 3
                                {{1, 3, 2}, {2, 3, 0}, {0, 3, 1}, {1, 2, 0}}  // tet 4
                        },
                };

        // each element stores the neighbor information for each tet face.
        // type - tet - face
        constexpr static TetNeighbor tetLocalFaceNeighborsLUT[2][5][4] =
                {
                        // type 0
                        {
                                {{13, 3}, {13, 2}, {13, 4}, {13, 1}}, // tet 0
                                {{13, 0}, {12, 2}, {10, 3}, {4,  1}}, // tet 1
                                {{22, 2}, {14, 1}, {13, 0}, {10, 4}}, // tet 2
                                {{22, 3}, {16, 1}, {13, 0}, {12, 4}}, // tet 3
                                {{16, 2}, {14, 3}, {13, 0}, {4,  4}}, // tet 4
                        },

                        // type 1
                        {
                                {{13, 4}, {13, 1}, {13, 3}, {13, 2}}, // tet 0
                                {{22, 1}, {13, 0}, {12, 2}, {10, 3}}, // tet 1
                                {{14, 1}, {13, 0}, {10, 4}, {4,  2}}, // tet 2
                                {{16, 1}, {13, 0}, {12, 4}, {4,  3}}, // tet 3
                                {{22, 4}, {16, 2}, {14, 3}, {13, 0}}, // tet 4
                        },
                };

        // each element stores the neighbor information for each tet edge.
        // type - tet - edge - neighbor
        constexpr static TetNeighbor tetLocalEdgeNeighborsLUT[2][5][6][6] =
                {
                        // type 0
                        {
                                // tet 0
                                {
                                        {{4, 0}, {4, 1}, {4, 4}, {13, 0}, {13, 1}, {13, 4}}, // edge 0

                                        {{10, 0}, {10, 3}, {10, 4}, {13, 0}, {13, 1}, {13, 2}}, // edge 1

                                        {{13, 0}, {13, 2}, {13, 4}, {14, 0}, {14, 1}, {14, 3}}, // edge 2

                                        {{12, 0}, {12, 2}, {12, 4}, {13, 0}, {13, 1}, {13, 3}}, // edge 3

                                        {{13, 0}, {13, 3}, {13, 4}, {16, 0}, {16, 1}, {16, 2}}, // edge 4

                                        {{13, 0}, {13, 2}, {13, 3}, {22, 0}, {22, 2}, {22, 3}}  // edge 5
                                },
                                // tet 1
                                {
                                        {{1, 3}, {4,  1}, {10, 3}, {13, 1}, {-1, -1}, {-1, -1}}, // edge 0

                                        {{3,  2}, {4,  1}, {12, 2}, {13, 1}, {-1, -1}, {-1, -1}}, // edge 1

                                        {{9,  4}, {10, 3}, {12, 2}, {13, 1}, {-1, -1}, {-1, -1}}, // edge 2

                                        {{4,  0}, {4,  1}, {4,  4}, {13, 0}, {13, 1},  {13, 4}}, // edge 3

                                        {{10, 0}, {10, 3}, {10, 4}, {13, 0}, {13, 1},  {13, 2}}, // edge 4

                                        {{12, 0}, {12, 2}, {12, 4}, {13, 0}, {13, 1}, {13, 3}}  // edge 5
                                },
                                // tet 2
                                {
                                        {{10, 0}, {10, 3}, {10, 4}, {13, 0}, {13, 1},  {13, 2}}, // edge 0

                                        {{10, 4}, {11, 3}, {13, 2}, {14, 1}, {-1, -1}, {-1, -1}}, // edge 1

                                        {{13, 0}, {13, 2}, {13, 4}, {14, 0}, {14, 1}, {14, 3}}, // edge 2

                                        {{10, 4}, {13, 2}, {19, 4}, {22, 2}, {-1, -1}, {-1, -1}}, // edge 3

                                        {{13, 0}, {13, 2}, {13, 3}, {22, 0}, {22, 2},  {22, 3}}, // edge 4

                                        {{13, 2}, {14, 1}, {22, 2}, {23, 1}, {-1, -1}, {-1, -1}}  // edge 5
                                },
                                // tet 3
                                {
                                        {{12, 0}, {12, 2}, {12, 4}, {13, 0}, {13, 1},  {13, 3}}, // edge 0

                                        {{12, 4}, {13, 3}, {15, 2}, {16, 1}, {-1, -1}, {-1, -1}}, // edge 1

                                        {{13, 0}, {13, 3}, {13, 4}, {16, 0}, {16, 1}, {16, 2}}, // edge 2

                                        {{12, 4}, {13, 3}, {21, 4}, {22, 3}, {-1, -1}, {-1, -1}}, // edge 3

                                        {{13, 0}, {13, 2}, {13, 3}, {22, 0}, {22, 2},  {22, 3}}, // edge 4

                                        {{13, 3}, {16, 1}, {22, 3}, {25, 1}, {-1, -1}, {-1, -1}}  // edge 5
                                },
                                // tet 4
                                {
                                        {{4,  0}, {4,  1}, {4,  4}, {13, 0}, {13, 1}, {13, 4}}, // edge 0

                                        {{4,  4}, {5,  3}, {13, 4}, {14, 3}, {-1, -1}, {-1, -1}}, // edge 1

                                        {{13, 0}, {13, 2}, {13, 4}, {14, 0}, {14, 1},  {14, 3}}, // edge 2

                                        {{4,  4}, {7,  2}, {13, 4}, {16, 2}, {-1, -1}, {-1, -1}}, // edge 3

                                        {{13, 0}, {13, 3}, {13, 4}, {16, 0}, {16, 1},  {16, 2}}, // edge 4

                                        {{13, 4}, {14, 3}, {16, 2}, {17, 1}, {-1, -1}, {-1, -1}}  // edge 5
                                }
                        },

                        // type 1
                        {
                                // tet 0
                                {
                                        {{4, 0}, {4, 2}, {4, 3}, {13, 0}, {13, 2}, {13, 3}}, // edge 0

                                        {{10, 0}, {10, 3}, {10, 4}, {13, 0}, {13, 1}, {13, 2}}, // edge 1

                                        {{12, 0}, {12, 2}, {12, 4}, {13, 0}, {13, 1}, {13, 3}}, // edge 2

                                        {{13, 0}, {13, 2}, {13, 4}, {14, 0}, {14, 1}, {14, 3}}, // edge 3

                                        {{13, 0}, {13, 3}, {13, 4}, {16, 0}, {16, 1}, {16, 2}}, // edge 4

                                        {{13, 0}, {13, 1}, {13, 4}, {22, 0}, {22, 1}, {22, 4}}  // edge 5
                                },
                                // tet 1
                                {
                                        {{9, 4}, {10, 3}, {12, 2}, {13, 1}, {-1, -1}, {-1, -1}}, // edge 0

                                        {{10, 0}, {10, 3}, {10, 4}, {13, 0}, {13, 1},  {13, 2}}, // edge 1

                                        {{12, 0}, {12, 2}, {12, 4}, {13, 0}, {13, 1},  {13, 3}}, // edge 2

                                        {{10, 3}, {13, 1}, {19, 3}, {22, 1}, {-1, -1}, {-1, -1}}, // edge 3

                                        {{12, 2}, {13, 1}, {21, 2}, {22, 1}, {-1, -1}, {-1, -1}}, // edge 4

                                        {{13, 0}, {13, 1}, {13, 4}, {22, 0}, {22, 1}, {22, 4}}  // edge 5
                                },
                                // tet 2
                                {
                                        {{1,  4}, {4,  2}, {10, 4}, {13, 2}, {-1, -1}, {-1, -1}}, // edge 0

                                        {{4,  0}, {4,  2}, {4,  3}, {13, 0}, {13, 2},  {13, 3}}, // edge 1

                                        {{10, 0}, {10, 3}, {10, 4}, {13, 0}, {13, 1}, {13, 2}}, // edge 2

                                        {{4,  2}, {5,  1}, {13, 2}, {14, 1}, {-1, -1}, {-1, -1}}, // edge 3

                                        {{10, 4}, {11, 3}, {13, 2}, {14, 1}, {-1, -1}, {-1, -1}}, // edge 4

                                        {{13, 0}, {13, 2}, {13, 4}, {14, 0}, {14, 1},  {14, 3}}  // edge 5
                                },
                                // tet 3
                                {
                                        {{3,  4}, {4,  3}, {12, 4}, {13, 3}, {-1, -1}, {-1, -1}}, // edge 0

                                        {{4,  0}, {4,  2}, {4,  3}, {13, 0}, {13, 2},  {13, 3}}, // edge 1

                                        {{12, 0}, {12, 2}, {12, 4}, {13, 0}, {13, 1}, {13, 3}}, // edge 2

                                        {{4,  3}, {7,  1}, {13, 3}, {16, 1}, {-1, -1}, {-1, -1}}, // edge 3

                                        {{12, 4}, {13, 3}, {15, 2}, {16, 1}, {-1, -1}, {-1, -1}}, // edge 4

                                        {{13, 0}, {13, 3}, {13, 4}, {16, 0}, {16, 1},  {16, 2}}  // edge 5
                                },
                                // tet 4
                                {
                                        {{13, 0}, {13, 2}, {13, 4}, {14, 0}, {14, 1}, {14, 3}}, // edge 0

                                        {{13, 0}, {13, 3}, {13, 4}, {16, 0}, {16, 1},  {16, 2}}, // edge 1

                                        {{13, 4}, {14, 3}, {16, 2}, {17, 1}, {-1, -1}, {-1, -1}}, // edge 2

                                        {{13, 0}, {13, 1}, {13, 4}, {22, 0}, {22, 1},  {22, 4}}, // edge 3

                                        {{13, 4}, {14, 3}, {22, 4}, {23, 3}, {-1, -1}, {-1, -1}}, // edge 4

                                        {{13, 4}, {16, 2}, {22, 4}, {25, 2}, {-1, -1}, {-1, -1}}  // edge 5
                                }
                        }
                };

        // each element stores the mapping from the (local) vertex indices of tet to the (local) indices of the node corners.
        // local tet vertex idx -> local node corner idx
        // type - tet - vertex
        constexpr static int tetVertexInNodeLUT[2][5][4] =
                {
                        {
                                {1, 2, 4, 7},
                                {0, 1, 2, 4},
                                {1, 4, 5, 7},
                                {2, 4, 6, 7},
                                {1, 2, 3, 7}
                        },
                        {
                                {0, 3, 5, 6},
                                {0, 4, 5, 6},
                                {0, 1, 3, 5},
                                {0, 2, 3, 6},
                                {3, 5, 6, 7}
                        }
                };

        // each element stores the mapping from the (local) vertex indices of each tet face to the (local) indices of the node corners.
        // local tet tris vertex idx -> local node corner idx
        // type - tet - facet - vertex
        constexpr static int tetTriInNodeLUT[2][5][4][3] =
                {
                        {
                                {{4, 7, 2}, {4, 1, 7}, {1, 2, 7}, {4, 2, 1}},

                                {{4, 1, 2}, {4, 2, 0}, {4, 0, 1}, {0, 2, 1}},

                                {{4, 5, 7}, {5, 1, 7}, {4, 7, 1}, {4, 1, 5}},

                                {{4, 7, 6}, {6, 7, 2}, {4, 2, 7}, {4, 6, 2}},

                                {{3, 2, 7}, {1, 3, 7}, {1, 7, 2}, {1, 2, 3}}
                        },
                        {
                                {{5, 3, 6}, {0, 5, 6}, {0, 6, 3}, {0, 3, 5}},

                                {{4, 5, 6}, {0, 6, 5}, {0, 4, 6}, {0, 5, 4}},

                                {{5, 1, 3}, {0, 5, 3}, {0, 1, 5}, {0, 3, 1}},

                                {{6, 3, 2}, {0, 3, 6}, {0, 6, 2}, {0, 2, 3}},

                                {{5, 7, 6}, {6, 7, 3}, {3, 7, 5}, {5, 6, 3}}
                        }
                };

    }
NAMESPACE_END(PCO)

#endif //PCO_OCTREETETLUT_HPP
