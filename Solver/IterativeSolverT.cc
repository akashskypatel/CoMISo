//=============================================================================
//
//  CLASS IterativeSolverT - IMPLEMENTATION
//
//=============================================================================

#define COMISO_ITERATIVESOLVERT_C

//== INCLUDES =================================================================

#include "IterativeSolverT.hh"
#include <Base/Debug/DebOut.hh>
#include <CoMISo/Solver/Eigen_Tools.hh>
#include <CoMISo/Solver/GMM_Tools.hh>

//== NAMESPACES ===============================================================

namespace COMISO
{

//== IMPLEMENTATION ==========================================================

template <class RealT>
bool IterativeSolverT<RealT>::gauss_seidel_local(const Matrix& _A, Vector& _x,
    const Vector& _rhs, const IndexVector& _idxs, const int _max_iter,
    const Real& _tolerance)
{
  if (COMISO_GMM::use_eigen())
    return gauss_seidel_local_eigen(_A, _x, _rhs, _idxs, _max_iter, _tolerance);

  if (_max_iter == 0)
    return false;

  typedef typename gmm::linalg_traits<Matrix>::const_sub_col_type ColT;
  typedef typename gmm::linalg_traits<ColT>::const_iterator CIter;

  updt_vrbl_indcs_.clear();
  indx_queue_.clear();

  for (unsigned int i = 0; i < _idxs.size(); ++i)
    indx_queue_.push_back(_idxs[i]);

  int it_count = 0;

  while (!indx_queue_.empty() && it_count < _max_iter)
  {
    ++it_count;
    const auto i = indx_queue_.front();
    indx_queue_.pop_front();
    indx_temp_.clear();

    double res_i = -_rhs[i];
    double x_i_new = _rhs[i];
    double diag = 1.0;

    const ColT col = mat_const_col(_A, i);
    for (auto it = gmm::vect_const_begin(col), ite = gmm::vect_const_end(col);
         it != ite; ++it)
    {
      const auto j = static_cast<unsigned>(it.index());
      res_i += (*it) * _x[j];
      x_i_new -= (*it) * _x[j];
      if (j != i)
        indx_temp_.push_back(j);
      else
        diag = *it;
    }

    // take inverse of diag
    diag = 1.0 / diag;

    // compare relative residuum normalized by diagonal entry
    if (fabs(res_i * diag) > _tolerance)
    {
      _x[i] += x_i_new * diag;
      updt_vrbl_indcs_.push_back(i);
      for (unsigned int j = 0; j < indx_temp_.size(); ++j)
        indx_queue_.push_back(indx_temp_[j]);
    }
  }

  return indx_queue_.empty(); // converged?
}

//-----------------------------------------------------------------------------

template <class RealT>
bool IterativeSolverT<RealT>::gauss_seidel_local_eigen(const Matrix& _A,
    Vector& _x, const Vector& _rhs, const IndexVector& _idxs,
    const int _max_iter, const Real& _tolerance)
{
  EigenMatrix A;
  EigenVector x;
  EigenVector rhs;
  COMISO_EIGEN::gmm_to_eigen(_A, A);
  COMISO_EIGEN::to_eigen_vec(_x, x);
  COMISO_EIGEN::to_eigen_vec(_rhs, rhs);

  auto res = gauss_seidel_local(A, x, rhs, _idxs, _max_iter, _tolerance);

  COMISO_EIGEN::from_eigen_vec(x, _x);

  return res;
}

//-----------------------------------------------------------------------------

template <class RealT>
bool IterativeSolverT<RealT>::gauss_seidel_local(const EigenMatrix& _A,
    EigenVector& _x, const EigenVector& _rhs, const IndexVector& _idxs,
    const int _max_iter, const Real& _tolerance)
{
  if (_max_iter == 0)
    return false;

  updt_vrbl_indcs_.clear();
  indx_queue_.clear();

  for (size_t i = 0; i < _idxs.size(); ++i)
    indx_queue_.push_back(_idxs[i]);

  int it_count = 0;

  while (!indx_queue_.empty() && it_count < _max_iter)
  {
    ++it_count;
    const auto i = indx_queue_.front();
    indx_queue_.pop_front();
    indx_temp_.clear();

    double res_i = -_rhs[i];
    double x_i_new = _rhs[i];
    double diag = 1.0;

    for (typename EigenMatrix::InnerIterator it(_A, i); it; ++it)
    {
      const auto j = static_cast<unsigned>(it.row());
      res_i += it.value() * _x[j];
      x_i_new -= it.value() * _x[j];
      if (j != i)
        indx_temp_.push_back(j);
      else
        diag = it.value();
    }

    // take inverse of diag
    diag = 1.0 / diag;

    // compare relative residuum normalized by diagonal entry
    if (std::abs(res_i * diag) > _tolerance)
    {
      _x[i] += x_i_new * diag;
      updt_vrbl_indcs_.push_back(i);
      for (size_t j = 0; j < indx_temp_.size(); ++j)
        indx_queue_.push_back(indx_temp_[j]);
    }
  }

  return indx_queue_.empty(); // converged?
}

//-----------------------------------------------------------------------------

template <class RealT>
bool IterativeSolverT<RealT>::conjugate_gradient(const Matrix& _A, Vector& _x,
    const Vector& _rhs, int& _max_iter, Real& _tolerance)
{
  if (COMISO_GMM::use_eigen())
    return conjugate_gradient_eigen(_A, _x, _rhs, _max_iter, _tolerance);

  DEB_enter_func;

  Real rho, rho_1(0), a;

  // initialize vectors
  p_.resize(_x.size());
  q_.resize(_x.size());
  r_.resize(_x.size());
  d_.resize(_x.size());
  gmm::copy(_x, p_); // gets overwritten before being used?

  // initialize diagonal (for relative norm)
  for (unsigned int i = 0; i < _x.size(); ++i)
    d_[i] = 1.0 / _A(i, i);

  // start with iteration 0
  int cur_iter(0);

  gmm::mult(_A, gmm::scaled(_x, Real(-1)), _rhs, r_); // r_ = _A * -_x + _rhs
  rho = gmm::vect_sp(r_, r_);
  gmm::copy(r_, p_);

  bool not_converged = true;
  Real res_norm(0);

  // while not converged
  while ((not_converged = ((res_norm = vect_norm_rel(r_, d_)) > _tolerance)) &&
         cur_iter < _max_iter)
  {
    DEB_line(11, "iter " << cur_iter << "  res " << res_norm);

    if (cur_iter != 0)
    {
      rho = gmm::vect_sp(r_, r_);
      gmm::add(r_, gmm::scaled(p_, rho / rho_1), p_);
    }

    gmm::mult(_A, p_, q_);

    a = rho / gmm::vect_sp(q_, p_);
    gmm::add(gmm::scaled(p_, a), _x);
    gmm::add(gmm::scaled(q_, -a), r_);
    rho_1 = rho;

    ++cur_iter;
  }

  _max_iter = cur_iter;
  _tolerance = res_norm;

  return !not_converged;
}

//-----------------------------------------------------------------------------


template <class RealT>
bool IterativeSolverT<RealT>::conjugate_gradient_eigen(const Matrix& _A,
    Vector& _x, const Vector& _rhs, int& _max_iter, Real& _tolerance)
{

  EigenMatrix A;
  EigenVector x;
  EigenVector rhs;
  COMISO_EIGEN::gmm_to_eigen(_A, A);
  COMISO_EIGEN::to_eigen_vec(_x, x);
  COMISO_EIGEN::to_eigen_vec(_rhs, rhs);

  bool res = conjugate_gradient(A, x, rhs, _max_iter, _tolerance);

  COMISO_EIGEN::from_eigen_vec(x, _x);

  return res;
}

//-----------------------------------------------------------------------------


template <class RealT>
bool IterativeSolverT<RealT>::conjugate_gradient(const EigenMatrix& _A,
    EigenVector& _x, const EigenVector& _rhs, int& _max_iter, Real& _tolerance)
{
  DEB_enter_func;
  Real rho, rho_1(0), a;

  // initialize vectors
  ep_.resize(_x.size());
  eq_.resize(_x.size());
  er_.resize(_x.size());
  ed_.resize(_x.size());
  ep_ = _x; // gets overwritten before being used?

  // initialize diagonal (for relative norm)
  ed_ = _A.diagonal().cwiseInverse();

  // start with iteration 0
  int cur_iter(0);

  er_ = _A * -_x + _rhs;
  rho = er_.dot(er_);
  ep_ = er_;

  bool not_converged = true;
  Real res_norm(0);

  // while not converged
  while ((not_converged = ((res_norm = vect_norm_rel(er_, ed_)) > _tolerance)) &&
         cur_iter < _max_iter)
  {
    DEB_line(11, "iter " << cur_iter << "  res " << res_norm);

    if (cur_iter != 0)
    {
      rho = er_.dot(er_);
      ep_ = er_ + rho / rho_1 * ep_;
    }

    eq_ = _A * ep_;

    a = rho / eq_.dot(ep_);
    _x += a * ep_;
    er_ -= a * eq_;
    rho_1 = rho;

    ++cur_iter;
  }

  _max_iter = cur_iter;
  _tolerance = res_norm;

  return !not_converged;
}

//-----------------------------------------------------------------------------

template <class RealT>
typename IterativeSolverT<RealT>::Real IterativeSolverT<RealT>::vect_norm_rel(
    const Vector& _v, const Vector& _diag) const
{
  Real res = 0.0;
  for (unsigned int i = 0; i < _v.size(); ++i)
    res = std::max(fabs(_v[i] * _diag[i]), res);
  return res;
}

//-----------------------------------------------------------------------------

template <class RealT>
typename IterativeSolverT<RealT>::Real IterativeSolverT<RealT>::vect_norm_rel(
    const EigenVector& _v, const EigenVector& _diag) const
{
  // compute component wise product
  const auto cwise_product = _v.array() * _diag.array();
  // return largest coefficient
  return cwise_product.abs().maxCoeff();
}


//=============================================================================
} // namespace COMISO
//=============================================================================
