//=============================================================================
//
//  CLASS IterativeSolverT
//
//=============================================================================

#ifndef COMISO_ITERATIVESOLVERT_HH
#define COMISO_ITERATIVESOLVERT_HH

//== INCLUDES =================================================================

#include <CoMISo/Utils/gmm.hh>
#include <deque>
#include <set>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO
{

//== CLASS DEFINITION =========================================================

/** \class IterativeSolverT IterativeSolverT.hh <COMISO/.../IterativeSolverT.hh>

    Brief Description.

    A more elaborate description follows.
*/

template <class RealT> class IterativeSolverT
{
public:
  typedef RealT Real;
  typedef std::vector<Real> Vector;
  typedef gmm::csc_matrix<Real> Matrix;

  // local gauss_seidel
  bool gauss_seidel_local(const Matrix& _A, Vector& _x, const Vector& _rhs,
      const std::vector<unsigned int>& _idxs, const int _max_iter,
      const Real& _tolerance);

  // local gauss_seidel
  bool gauss_seidel_local2(const Matrix& _A, Vector& _x, const Vector& _rhs,
      const std::vector<unsigned int>& _idxs, const int _max_iter,
      const Real& _tolerance);

  // conjugate gradient
  bool conjugate_gradient(const Matrix& _A, Vector& _x, const Vector& _rhs,
      int& _max_iter, Real& _tolerance);

private:
  // compute relative norm
  Real vect_norm_rel(const Vector& _v, const Vector& _diag) const;

private:
  // helper for conjugate gradient
  Vector p_;
  Vector q_;
  Vector r_;
  Vector d_;

  //  helper for local gauss seidel
  std::vector<unsigned int> i_temp;
  std::deque<unsigned int> q;
  std::set<int> s;
};

//=============================================================================
} // namespace COMISO
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(COMISO_ITERATIVESOLVERT_C)
#define COMISO_ITERATIVESOLVERT_TEMPLATES
#include "IterativeSolverT.cc"
#endif
//=============================================================================
#endif // COMISO_ITERATIVESOLVERT_HH defined
//=============================================================================
