SET(my_headers
   ${CMAKE_CURRENT_SOURCE_DIR}/BoundConstraint.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/COMISOSolver.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/ConeConstraint.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/CBCSolver.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/CPLEXSolver.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/DOCloudSolver.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/GurobiHelper.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/GUROBISolver.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/IPOPTSolver.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/LeastSquaresProblem.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/LinearConstraint.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/LinearConstraintHandlerElimination.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/LinearConstraintHandlerPenalty.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/LinearProblem.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/NConstraintInterface.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/NewtonSolver.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/NPDerivativeChecker.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/NPLinearConstraints.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/NProblemGmmInterface.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/NProblemInterface.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/NProblemInterfaceADOLC.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/NPTiming.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/QuadraticProblem.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/SuperSparseMatrixT.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/TAOSolver.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/TapeIDSingleton.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/TestInterface.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/VariableType.hh
   PARENT_SCOPE
)

set(my_t_impls
   ${CMAKE_CURRENT_SOURCE_DIR}/LinearConstraintHandlerPenaltyT.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/LinearConstraintHandlerEliminationT.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/NPLinearConstraintsT.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/SuperSparseMatrixT.cc
   PARENT_SCOPE
)

SET(my_sources
   ${CMAKE_CURRENT_SOURCE_DIR}/BoundConstraint.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/COMISOSolver.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/ConeConstraint.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/CBCSolver.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/CPLEXSolver.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/DOCloudSolver.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/GurobiHelper.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/GUROBISolver.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/IPOPTSolver.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/LeastSquaresProblem.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/LinearConstraint.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/LinearConstraintHandlerElimination.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/LinearConstraintHandlerPenalty.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/NewtonSolver.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/NPDerivativeChecker.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/NPLinearConstraints.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/NProblemInterface.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/NPTiming.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/TAOSolver.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/TapeIDSingleton.cc
   PARENT_SCOPE
)



   