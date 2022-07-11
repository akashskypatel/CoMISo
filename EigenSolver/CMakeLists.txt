SET(my_headers
   ${CMAKE_CURRENT_SOURCE_DIR}/ArpackSolver.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/EigenArpackMatrixT.hh
   PARENT_SCOPE
)

SET(my_t_impls
   ${CMAKE_CURRENT_SOURCE_DIR}/EigenArpackMatrixT_impl.hh
   PARENT_SCOPE
)

SET(my_sources
   ${CMAKE_CURRENT_SOURCE_DIR}/ArpackSolver.cc
   PARENT_SCOPE
)

