project(small_quadratic_resolve_example
    LANGUAGES C CXX)

if (NOT TARGET CoMISo::CoMISo)
    find_package(CoMISo MODULE REQUIRED)
endif()

add_executable(small_quadratic_resolve_example "main.cc")
target_link_libraries (small_quadratic_resolve_example CoMISo::CoMISo)
