project(small_quadratic_example
    LANGUAGES C CXX)

if (NOT TARGET CoMISo::CoMISo)
    find_package(CoMISo MODULE REQUIRED)
endif()

add_executable(small_quadratic_example "main.cc")
target_link_libraries (small_quadratic_example CoMISo::CoMISo)
