project(small_factored_example
    LANGUAGES C CXX)

if (NOT TARGET CoMISo::CoMISo)
    find_package(CoMISo MODULE REQUIRED)
endif()

add_executable(small_factored_example "main.cc")
target_link_libraries (small_factored_example CoMISo::CoMISo)
