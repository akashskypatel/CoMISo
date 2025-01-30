project(small_symmetric_dirichlet
    LANGUAGES C CXX)

if (NOT TARGET CoMISo::CoMISo)
    find_package(CoMISo MODULE REQUIRED)
endif()

add_executable(small_symmetric_dirichlet "main.cc")
target_link_libraries (small_symmetric_dirichlet CoMISo::CoMISo)
