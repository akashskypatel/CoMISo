project(small_tinyad
    LANGUAGES C CXX)

if (NOT TARGET CoMISo::CoMISo)
    find_package(CoMISo MODULE REQUIRED)
endif()

add_executable(small_tinyad "main.cc")
target_link_libraries (small_tinyad CoMISo::CoMISo)
