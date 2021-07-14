// Copyright 2021 Autodesk, Inc. All rights reserved.

//=============================================================================
//
//  CLASS SolverBaseT
//
//=============================================================================

#ifndef COMISO_SOLVERBASET_HH
#define COMISO_SOLVERBASET_HH

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if (COMISO_EIGEN3_AVAILABLE)

//== INCLUDES =================================================================

#include <CoMISo/Config/CoMISoDefines.hh>

#include <array>
#include <vector>

//== NAMESPACES ===============================================================

namespace COMISO
{

template <size_t DIM> struct SolverBaseT
{
  using Point = std::array<double, DIM>; // Object to minimize

  using PointVector = std::vector<Point>;

  struct LinearTerm // It means that coeff_ * X_name_
  {
    size_t var_name; // Variable name.
    double coeff;
    bool operator<(const LinearTerm& _lt) const
    {
      return var_name < _lt.var_name;
    }
  }; // struct LinearTerm

  using LinearTermVector = std::vector<LinearTerm>;

  struct LinearEquation // a linear equation in the form Sum(c_i * x_i) =
                        // constant_term
  {
    LinearTermVector linear_terms;
    // Construction must create an equation 'nothing' = 0, so const_term is
    // initialized with zeros.
    Point const_term{};
  }; // struct LinearEquation

}; // struct Types

//=============================================================================
} // namespace COMISO
//=============================================================================

//=============================================================================
#endif // COMISO_EIGEN3_AVAILABLE
//=============================================================================
#endif // COMISO_SOLVERBASET_HH defined
//=============================================================================
