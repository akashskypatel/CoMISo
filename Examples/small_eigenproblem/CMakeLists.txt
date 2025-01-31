project(small_eigenproblem
    LANGUAGES C CXX)

if (NOT TARGET CoMISo::CoMISo)
    find_package(CoMISo MODULE REQUIRED)
endif()

add_executable(small_eigenproblem "main.cc")
target_link_libraries (small_eigenproblem CoMISo::CoMISo)
