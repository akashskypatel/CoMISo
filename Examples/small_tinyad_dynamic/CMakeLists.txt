project(small_tinyad_dynamic
    LANGUAGES C CXX)

if (NOT TARGET CoMISo::CoMISo)
    find_package(CoMISo MODULE REQUIRED)
endif()

add_executable(small_tinyad_dynamic "main.cc")
target_link_libraries (small_tinyad_dynamic CoMISo::CoMISo)
