project(small_linear_problem
    LANGUAGES C CXX)

if (NOT TARGET CoMISo::CoMISo)
    find_package(CoMISo MODULE REQUIRED)
endif()

add_executable(small_linear_problem "main.cc")
target_link_libraries (small_linear_problem CoMISo::CoMISo)
