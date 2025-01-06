//=============================================================================
//
//  CLASS ConstraintTools
//
//=============================================================================


#ifndef COMISO_CONSTRAINTTOOLS_HH
#define COMISO_CONSTRAINTTOOLS_HH


//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_EIGEN3_AVAILABLE

//== INCLUDES =================================================================

#include <stdio.h>
#include <iostream>
#include <vector>

#include <CoMISo/Utils/gmm.hh>

#include <CoMISo/Config/CoMISoDefines.hh>
#include <CoMISo/NSolver/NConstraintInterface.hh>
#include <CoMISo/Solver/Eigen_Tools.hh>
#include <CoMISo/NSolver/LinearConstraintConverter.hh>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO 
{

namespace ConstraintTools
{
using ConstraintVector = std::vector<NConstraintInterface*>;
const double DEFAULT_EPS = 1e-8; // TODO: document

// struct to return the result of the constraint elimination
struct ConstraintRemovalResult {
  size_t n_constraints_eliminated = 0;
  size_t n_infeasible_detected = 0;
};

// Remove all linear dependent linear equality constraints. The remaining
// constraints are a subset of the original ones. Non-linear or equality
// constraints are preserved.
COMISODLLEXPORT ConstraintRemovalResult remove_dependent_linear_constraints(
    ConstraintVector& _constraints, const double _eps = DEFAULT_EPS);

// As above but assumes that all constraints are linear equality constraints
COMISODLLEXPORT ConstraintRemovalResult remove_dependent_linear_constraints_only_linear_equality(
    ConstraintVector& _constraints, const double _eps = DEFAULT_EPS);

// same as above but designed for Eigen::SparseMatrix
template< typename SMatrixT, typename VectorT>
COMISODLLEXPORT ConstraintRemovalResult remove_dependent_linear_constraints(SMatrixT& _A, VectorT& _b, const double _eps = DEFAULT_EPS)
{
  // convert into NConstraints
  LinearConstraintConverter lcc(_A,_b);
  std::vector<NConstraintInterface*> constraints = lcc.constraints_nsolver();

  // process
  auto result = remove_dependent_linear_constraints_only_linear_equality(constraints, _eps);

  // convert back
  LinearConstraintConverter::nsolver_to_eigen(constraints, _A, _b);
  return result;
}


using HalfSparseRowMatrix = COMISO_EIGEN::HalfSparseRowMatrix<double>;
using HalfSparseColMatrix = COMISO_EIGEN::HalfSparseColMatrix<double>;
using SparseVector = Eigen::SparseVector<double>;
using uint = unsigned int;
using IntVector = std::vector<int>;

enum Flags // TODO: document flags
{
  FL_NONE = 0,
  FL_DO_GCD,
  FL_REORDER = FL_DO_GCD << 1,
  FL_DEFAULT = FL_DO_GCD | FL_REORDER
};


// struct to return the result of the constraint elimination
struct GaussEliminationResult {
  size_t n_rows_linearly_dependent = 0;
  size_t n_rows_contradicting = 0;
};


/*!
Perform Gauss elimination on the constraint matrix to facilitate constraint
elimination downstream.
 
\note Contradicting constraints are ignored.
  
\warning Care must be taken downstream when non-trivial constraints occur
where some of the variables contain integer-variables (to be rounded) as
the optimal result might not always occur.
*/
COMISODLLEXPORT GaussEliminationResult gauss_elimination(
    HalfSparseRowMatrix& _constraints, // constraint matrix
    IntVector& _elmn_clmn_indcs, // return the variable indices and the order in
                                 // which they can be eliminated
    const IntVector& _indcs_to_round = IntVector(), // variables to be rounded
    HalfSparseRowMatrix* _update_D = nullptr,       // TODO: document
    const double _eps = DEFAULT_EPS,                // TODO: document
    const uint _flags = Flags::FL_DEFAULT           // control execution flags
);

} // namespace ConstraintTools

//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // COMISO_GMM_AVAILABLE
//=============================================================================
#endif // COMISO_CONSTRAINTTOOLS_HH defined
//=============================================================================

