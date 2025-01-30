project(small_adolc
    LANGUAGES C CXX)

if (NOT TARGET CoMISo::CoMISo)
    find_package(CoMISo MODULE REQUIRED)
endif()

add_executable(small_adolc "main.cc")
target_link_libraries (small_adolc CoMISo::CoMISo)
