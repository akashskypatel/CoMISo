//=============================================================================
//
//  CLASS IterativeSolverT
//
//=============================================================================

#ifndef COMISO_ITERATIVESOLVERT_HH
#define COMISO_ITERATIVESOLVERT_HH

//== INCLUDES =================================================================

#include <CoMISo/Utils/gmm.hh>
#include <Eigen/Sparse>
#include <deque>

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
  typedef unsigned int uint;
  typedef RealT Real;
  typedef std::vector<Real> Vector;
  typedef std::vector<uint> IndexVector;
  typedef gmm::csc_matrix<Real> Matrix;
  typedef Eigen::SparseMatrix<Real, Eigen::ColMajor> EigenMatrix;
  typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> EigenVector;

  // local Gauss-Seidel
  bool gauss_seidel_local(const Matrix& _A, Vector& _x, const Vector& _rhs,
      const IndexVector& _idxs, const int _max_iter, const Real& _tolerance);

  bool gauss_seidel_local_eigen(const Matrix& _A, Vector& _x, const Vector& _rhs,
      const IndexVector& _idxs, const int _max_iter, const Real& _tolerance);

  bool gauss_seidel_local(const EigenMatrix& _A, EigenVector& _x, const EigenVector& _rhs,
      const IndexVector& _idxs, const int _max_iter, const Real& _tolerance);

  // get the indices of any variables updated during the last local Gauss-Seidel
  const IndexVector& updated_variable_indices() const
  {
    return updt_vrbl_indcs_;
  }

  // conjugate gradient
  bool conjugate_gradient(const Matrix& _A, Vector& _x, const Vector& _rhs,
      int& _max_iter, Real& _tolerance);

  // conjugate gradient
  bool conjugate_gradient_eigen(const Matrix& _A, Vector& _x, const Vector& _rhs,
      int& _max_iter, Real& _tolerance);

  // conjugate gradient
  bool conjugate_gradient(const EigenMatrix& _A, EigenVector& _x, const EigenVector& _rhs,
      int& _max_iter, Real& _tolerance);

private:
  // compute relative norm
  Real vect_norm_rel(const Vector& _v, const Vector& _diag) const;

  // compute relative norm
  Real vect_norm_rel(const EigenVector& _v, const EigenVector& _diag) const;

private:
  // context  for Conjugate Gradient
  Vector p_;
  Vector q_;
  Vector r_;
  Vector d_;

  EigenVector ep_;
  EigenVector eq_;
  EigenVector er_;
  EigenVector ed_;

  // context for local Gauss-Seidel
  IndexVector indx_temp_;
  std::deque<uint> indx_queue_;
  IndexVector updt_vrbl_indcs_; // updated variable indices
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
