project(small_cplex_soc
    LANGUAGES C CXX)

if (NOT TARGET CoMISo::CoMISo)
    find_package(CoMISo MODULE REQUIRED)
endif()

add_executable(small_cplex_soc "main.cc")
target_link_libraries (small_cplex_soc CoMISo::CoMISo)
