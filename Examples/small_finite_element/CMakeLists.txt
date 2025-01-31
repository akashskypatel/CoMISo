project(small_finite_element
    LANGUAGES C CXX)

if (NOT TARGET CoMISo::CoMISo)
    find_package(CoMISo MODULE REQUIRED)
endif()

add_executable(small_finite_element "main.cc")
target_link_libraries (small_finite_element CoMISo::CoMISo)
