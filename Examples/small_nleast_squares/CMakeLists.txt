project(small_nleast_squares
    LANGUAGES C CXX)

if (NOT TARGET CoMISo::CoMISo)
    find_package(CoMISo MODULE REQUIRED)
endif()

add_executable(small_nleast_squares "main.cc")
target_link_libraries (small_nleast_squares CoMISo::CoMISo)
