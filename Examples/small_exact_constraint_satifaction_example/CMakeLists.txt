project(small_exact_constraint_satifaction_example
    LANGUAGES C CXX)

if (NOT TARGET CoMISo::CoMISo)
    find_package(CoMISo MODULE REQUIRED)
endif()

add_executable(small_exact_constraint_satifaction_example "main.cc")
target_link_libraries (small_exact_constraint_satifaction_example CoMISo::CoMISo)
