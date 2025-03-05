if (NOT TARGET Eigen3::Eigen)
    include(eigen)
endif ()

if (NOT TARGET igl::core)
    include(libigl)
endif ()

if (NOT TARGET CLI11::CLI11)
    include(CLI)
    target_compile_definitions(CLI11 INTERFACE -DCLI11_STD_OPTIONAL=0)
    target_compile_definitions(CLI11 INTERFACE -DCLI11_EXPERIMENTAL_OPTIONAL=0)
endif ()

if (NOT TARGET spdlog::spdlog)
    include(spdlog)
endif ()

if (NOT TARGET indirect_predicates::indirect_predicates)
    include(indirect_predicates)
endif ()

if (NOT TARGET implicit_predicates::implicit_predicates)
    include(implicit_predicates)
endif ()

if (NOT TARGET CGAL::CGAL)
    include(cgal)
endif ()