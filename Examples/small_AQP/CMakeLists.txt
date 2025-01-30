project(small_AQP
    LANGUAGES C CXX)

if (NOT TARGET CoMISo::CoMISo)
    find_package(CoMISo MODULE REQUIRED)
endif()

add_executable(small_AQP "main.cc")
target_link_libraries (small_AQP CoMISo::CoMISo)
