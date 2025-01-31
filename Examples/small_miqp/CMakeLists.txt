project(small_miqp
    LANGUAGES C CXX)

if (NOT TARGET CoMISo::CoMISo)
    find_package(CoMISo MODULE REQUIRED)
endif()

add_executable(small_miqp "main.cc")
target_link_libraries (small_miqp CoMISo::CoMISo)
