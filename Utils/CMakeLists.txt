set(my_headers
   ${CMAKE_CURRENT_SOURCE_DIR}/CoMISoError.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/CoMISoErrorInc.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/ExactConstraintSatisfaction.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/MutablePriorityQueueT.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/SmartPointer.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/Tools.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/VSToolsT.hh
   PARENT_SCOPE
)

set(my_sources
   ${CMAKE_CURRENT_SOURCE_DIR}/CoMISoError.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/ExactConstraintSatisfaction.cc
   PARENT_SCOPE
)
