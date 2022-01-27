// Copyright 2022 Autodesk, Inc. All rights reserved.


//=============================================================================
//
//  CLASS MultiDimConstrainedSolver
//
//=============================================================================


#ifndef COMISO_SUBSYSTEMSOLVERT_HH
#define COMISO_SUBSYSTEMSOLVERT_HH


//== INCLUDES =================================================================
#include <CoMISo/Config/CoMISoDefines.hh>
#include <CoMISo/Config/StdTypes.hh>
#include <CoMISo/Solver/SolverBaseT.hh>
#include <CoMISo/Solver/MultiDimConstrainedSolverT.hh>
#include <CoMISo/Utils/ProblemSubsetMapT.hh>

#include <vector>

//== NAMESPACES ===============================================================

namespace COMISO
{


// Class to solve linear systems of the form:
// Minimize ||Ax - b||^2, subject to linear constraints Cx=d,
// fix point constraints, and integer constraints.
// This class handles efficiently the case of solving a system with many unused
// variables. For solving many such systems its worth keeping the solver object
// alive as its construction is in O(n) where n is largest variable index.
template <int DIM>
class COMISODLLEXPORT SubProblemSolverT
{
public:
  using Point            = typename SolverBaseT<DIM>::Point;
  using LinearEquation   = typename SolverBaseT<DIM>::LinearEquation;
  using ValueVector      = typename SolverBaseT<DIM>::ValueVector;
  using IndexVector      = std::vector<int>;
  using Result           = ValueVector;


  SubProblemSolverT(
      size_t _max_n_vars, const ValueVector& _fixed_values = ValueVector());

  // Add an equation to the system. solve() will minimize sum of the quadratic
  // errors of all equations
  void add_equation(LinearEquation _eq);

  // Add a linear constraint to the system
  void add_constraint(LinearEquation _eq);

  // Add multiple linear constraints to the system
  void add_constraints(std::vector<LinearEquation> _eqs);

  // Set the fix point constraints which are special linear constraints for
  // fixing individual variables to given values. They are handled more
  // efficiently than the more general constraints specified with
  // add_constraints() by being removed from the system immediatly.
  // Also clears all equations, linear constraints, and integer constraints.
  void reset(ValueVector _fixed_values);

  // Set the integer constraints
  void set_integers(IndexVector _int_var_indcs);

  // Solve the system that has been setup with the calls above.
  // Clears all equations, constraints, and integer constraints that
  // have been setup using the functions above. Keeps fixed values.
  void solve(Result& _result);

private:

  ProblemSubsetMapT<DIM> sbst_map_;
  MultiDimConstrainedSolverT<DIM> solver_; // Solver used to solve subproblem
};


//=============================================================================
} // namespace COMISO
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(COMISO_SUBSYSTEMSOLVERT_C)
#define COMISO_SUBSYSTEMSOLVERT_TEMPLATES
#include "SubProblemSolverT.cc"
#endif
//=============================================================================
#endif // COMISO_CONSTRAINEDSOLVER_HH defined
//=============================================================================

