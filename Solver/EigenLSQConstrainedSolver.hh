// Copyright 2021 Autodesk, Inc. All rights reserved.

//=============================================================================
//
//  CLASS EigenLSQConstrainedSolverT
//
//=============================================================================

#ifndef COMISO_EIGEN_LSQC_CONSTRAINED_SOLVER_HH
#define COMISO_EIGEN_LSQC_CONSTRAINED_SOLVER_HH

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if (COMISO_EIGEN3_AVAILABLE)
//== INCLUDES =================================================================

#include <CoMISo/Config/CoMISoDefines.hh>

#include <array>
#include <memory>
#include <vector>

//== NAMESPACES ===============================================================

namespace COMISO
{

//== CLASS DEFINITION =========================================================

/** \class EigenLSQConstrainedSolver EigenLSQConstrainedSolver.hh

    Constrained Least Square Minimization Solver.
    Problem
      find x that minimize             Sum_i(Sum_j((a_ij*x_jh - b_ih)^2))
      with a set of fixed constraints  x_i = b_ih
    and a set of linear constraints    Sum (c_ij * x_j) = d_jh
    i is in [0, number of equations]
    j is in a set variable names
    h is in [0,DIM), where DIM is the dimension of the point we want to solve

    x is an array of solution (var_name, Point)
*/

template <size_t DIM>
class COMISODLLEXPORT EigenLSQConstrainedSolverT
{
public:
  using Point = std::array<double, DIM>; // Object to minimize
  struct Value // It means that X_name_ = val_. 
               // Used to get the results and to set the fixed variables
  {
    size_t var_name; // Variable name.
    Point point; // Value of the variable
    bool operator<(const Value& _vl) const { return var_name < _vl.var_name; }
    bool operator==(const Value& _vl) const { return var_name == _vl.var_name; }
  };
  using ValueVector = std::vector<Value>;
  struct LinearTerm // It means that coeff_ * X_name_
  {
    size_t var_name; // Variable name.
    double coeff;
    bool operator<(const LinearTerm& _lt) const
    {
      return var_name < _lt.var_name;
    }
  };
  using LinearTermVector = std::vector<LinearTerm>;
  struct LinearEquation // a linear equation in the form Sum(c_i * x_i) =
                        // constant_term
  {
    LinearTermVector linear_terms;
    Point const_term;
  };

  /** \constructor
   * _var_nmbr: the number of points to optimize.
                The variable names are [0, 1, .._vars]
   * _fixed:    vector of fixed positions in the form (var_name, point)
   *            These are fixed points, that will not be minimized.
   */
  EigenLSQConstrainedSolverT(size_t _var_nmbr, ValueVector&& _fixed);

  /** \constructor
   * _var_names: a vector of positive integers representing the names of the
                 variable to optimize. Can be not sequential.
   * _fixed: vector of fixed positions in the form (var_name, point)
   */
  EigenLSQConstrainedSolverT(
      std::vector<size_t>&& _var_names, ValueVector&& _fixed);

  ~EigenLSQConstrainedSolverT();

  /** \Add a linear equation in the form Sum_i(a_i * x_i) = b_i
  * 
  * The solution will minimize the sum of the square of all the added
   * equations: min Sum((Sum_i(a_i * x_i) - b_i)^2)
   */
  void add_equation(const LinearEquation& _lin_eq);

  /** \Add a linear constraint
  * The solution must satisfy it.
  * The input data is consumed by the function.
  */
  void add_linear_constraint(LinearEquation&& _lin_cnstr);

  /** \Return the list of fixed points, without duplications
   */
  const ValueVector& fixed_points() const;

  /** \Solve the problem and return the solution.
   * Note: do not use the output after the class is destroyed
   * The function can throw if it find a reference to a variable outside
   * the valid range or if the system is under-constraint, i.e., it does not
   * have a unique minimum.
   */
  const ValueVector& solve();

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

//=============================================================================
} // namespace COMISO
//=============================================================================

//=============================================================================
#endif // COMISO_EIGEN3_AVAILABLE
//=============================================================================
#endif // COMISO_EIGEN_LSQC_CONSTRAINED_SOLVER_HH defined
//=============================================================================
