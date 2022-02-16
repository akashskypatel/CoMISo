//=============================================================================
//
//  CLASS CoonstraintTools
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

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================


/** \class ConstraintTools ConstraintTools.hh <CoMISo/NSolver/ConstraintTools.hh>

    A more elaborate description follows.
*/
class COMISODLLEXPORT ConstraintTools
{
public:

  using HalfSparseRowMatrix = COMISO_EIGEN::HalfSparseRowMatrix<double>;
  using HalfSparseColMatrix = COMISO_EIGEN::HalfSparseColMatrix<double>;
  typedef Eigen::SparseVector<double> SparseVector;

  ConstraintTools(double _epsilon = 1e-8, bool _do_gcd = true,
      bool _use_reordering = true, HalfSparseRowMatrix* _update_D = nullptr)
      : epsilon_(_epsilon), do_gcd_(_do_gcd), use_reordering_(_use_reordering),
        update_D_(_update_D)
  {
  }


  // remove all linear dependent linear equality constraints. the remaining
  // constraints are a subset of the original ones nonlinear or equality
  // constraints are preserved.
  static void remove_dependent_linear_constraints(
      std::vector<NConstraintInterface*>& _constraints,
      const double _eps = 1e-8);

  // same as above but assumes already that all constraints are linear equality
  // constraints
  static void remove_dependent_linear_constraints_only_linear_equality(
      std::vector<NConstraintInterface*>& _constraints,
      const double _eps = 1e-8);


  /// Make constraints independent
  /**
   *  This function performs a Gauss elimination on the constraint matrix making
   * the constraints easier to eliminate. \note A certain amount of independence
   * of the constraints is assumed. \note contradicting constraints will be
   * ignored. \warning care must be taken when non-trivial constraints occur
   * where some of the variables contain integer-variables (to be rounded) as
   * the optimal result might not always occur.
   *  @param _constraints  row matrix with constraints
   *  @param _idx_to_round indices of variables to be rounded (these must be
   * considered.)
   *  @param _c_elim the "returned" vector of variable indices and the order in
   * which the can be eliminated.
   */
  void make_constraints_independent(
          HalfSparseRowMatrix& _constraints,
          const std::vector<int>& _idx_to_round,
                std::vector<int>& _c_elim);

  // Same as above but without rounded indices
  void make_constraints_independent(
      HalfSparseRowMatrix& _constraints, std::vector<int>& _c_elim);

private:

  // adjust _constraints such that col _col is zero except in row _row and
  // and those rows that should be ignored.
  // col will be chosen automatically.
  // _changed_rows returns the rows that were changed.
  void make_constraint_independent(
          HalfSparseRowMatrix& _constraints,
          HalfSparseColMatrix& _constraints_c,
          int                  _row,
          int&                 _col,
    const std::vector<bool>&   _integer,
    const std::vector<bool>&   _ignore,
          std::vector<int>&    _changed_rows);

  // same as above but without returning changed rows
  void make_constraint_independent(
            HalfSparseRowMatrix& _constraints,
            HalfSparseColMatrix& _constraints_c,
      const int                  _row,
            int&                 _col,
      const std::vector<bool>&   _integer,
      const std::vector<bool>&   _ignore);

  void make_constraints_independent_reordering(
          HalfSparseRowMatrix& _constraints,
          const std::vector<int>& _idx_to_round,
                std::vector<int>& _c_elim);

  void make_constraints_independent_no_reordering(
          HalfSparseRowMatrix& _constraints,
          const std::vector<int>& _idx_to_round,
                std::vector<int>& _c_elim);

  // add _coeff * (_source_row of _source_mat)  to _target_row of _target_rmat
  // and target_cmat. set element in _zero_col to 0 if it exists
  void add_row_simultaneously(
    const Eigen::Index         _target_row,
    const double               _coeff,
    const HalfSparseRowMatrix& _source_mat,
    const Eigen::Index         _source_row,
          HalfSparseRowMatrix& _target_rmat,
          HalfSparseColMatrix& _target_cmat,
    const Eigen::Index         _zero_col = -1);


  // TODO if no gcd correction was possible, at least use a variable divisible
  // by 2 as new elim_j (to avoid in-exactness e.g. 1/3)
  bool update_constraint_gcd(
          SparseVector&     _row,
    const int               _elim_j,
          std::vector<int>& _v_gcd,
          int&              _n_ints);

  static inline int gcd(int _a, int _b)
  {
    while (_b != 0)
    {
      int t(_b);
      _b = _a % _b;
      _a = t;
    }
    return _a;
  }

  static int find_gcd(std::vector<int>& _v_gcd, int& _n_ints);

private:
  double epsilon_;
  bool do_gcd_;
  bool use_reordering_;
  HalfSparseRowMatrix* update_D_;

};


//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // COMISO_GMM_AVAILABLE
//=============================================================================
#endif // COMISO_CONSTRAINTTOOLS_HH defined
//=============================================================================

