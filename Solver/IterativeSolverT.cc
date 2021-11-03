//=============================================================================
//
//  CLASS IterativeSolverT - IMPLEMENTATION
//
//=============================================================================

#define COMISO_ITERATIVESOLVERT_C

//== INCLUDES =================================================================

#include "IterativeSolverT.hh"
#include <Base/Debug/DebOut.hh>

//== NAMESPACES ===============================================================

namespace COMISO
{

//== IMPLEMENTATION ==========================================================

template <class RealT>
bool IterativeSolverT<RealT>::gauss_seidel_local(const Matrix& _A, Vector& _x,
    const Vector& _rhs, const IndexVector& _idxs, const int _max_iter,
    const Real& _tolerance)
{
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
bool IterativeSolverT<RealT>::conjugate_gradient(const Matrix& _A, Vector& _x,
    const Vector& _rhs, int& _max_iter, Real& _tolerance)
{
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
typename IterativeSolverT<RealT>::Real IterativeSolverT<RealT>::vect_norm_rel(
    const Vector& _v, const Vector& _diag) const
{
  Real res = 0.0;
  for (unsigned int i = 0; i < _v.size(); ++i)
    res = std::max(fabs(_v[i] * _diag[i]), res);
  return res;
}


//=============================================================================
} // namespace COMISO
//=============================================================================
