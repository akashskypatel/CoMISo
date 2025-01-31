project(small_nsolver
    LANGUAGES C CXX)

if (NOT TARGET CoMISo::CoMISo)
    find_package(CoMISo MODULE REQUIRED)
endif()

add_executable(small_nsolver "main.cc")
target_link_libraries (small_nsolver CoMISo::CoMISo)
