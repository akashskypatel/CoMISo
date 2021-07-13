// Copyright 2021 Autodesk, Inc. All rights reserved.

//=============================================================================
//
//  CLASS EigenLSQSolverT
//
//=============================================================================

#ifndef COMISO_EIGEN_LSQC_SOLVER_HH
#define COMISO_EIGEN_LSQC_SOLVER_HH

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if (COMISO_EIGEN3_AVAILABLE)
//== INCLUDES =================================================================

#include <CoMISo/Config/CoMISoDefines.hh>
#include <CoMISo/Config/StdTypes.hh>

//== NAMESPACES ===============================================================

namespace COMISO
{

//== CLASS DEFINITION =========================================================

/** \class EigenLSQSolver EigenLSQSolver.hh

    Least Square Minimization Solver.
    Problem
      find x that minimize             Sum_i(Sum_j((a_ij*x_jh - b_ih)^2))
    i is in [0, number of equations]
    j is in [0, number of variables]
    h is in [0,DIM), where DIM is the dimension of the point we want to solve

    x is an array of solution (var_name, Point)
*/

template <size_t DIM>
class COMISODLLEXPORT EigenLSQSolverT : public COMISO_STD::SolverBaseT<DIM>
{
public:
  using ValueVector = PointVector;
  void add_equation(LinearEquation&& _lin_eq)
  {
    lin_eqs_.emplace_back(std::move(_lin_eq));
  }

  const ValueVector& solve();

private:
  std::vector<LinearEquation> lin_eqs_;
  // System solution
  ValueVector result_;
};

//=============================================================================
} // namespace COMISO
//=============================================================================

//=============================================================================
#endif // COMISO_EIGEN3_AVAILABLE
//=============================================================================
#endif // COMISO_EIGEN_LSQC_SOLVER_HH defined
//=============================================================================
