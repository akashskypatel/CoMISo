project(small_sparseqr
    LANGUAGES C CXX)

if (NOT TARGET CoMISo::CoMISo)
    find_package(CoMISo MODULE REQUIRED)
endif()

add_executable(small_sparseqr "main.cc")
target_link_libraries (small_sparseqr CoMISo::CoMISo)
