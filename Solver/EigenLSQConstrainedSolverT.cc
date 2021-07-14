// Copyright 2021 Autodesk, Inc. All rights reserved.

#include "EigenLSQConstrainedSolver.hh"

#include <CoMISo/Utils/CoMISoError.hh>

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#if (COMISO_EIGEN3_AVAILABLE)
//== INCLUDES =================================================================

#include <Base/Code/Quality.hh>
LOW_CODE_QUALITY_SECTION_BEGIN
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SparseCholesky>
LOW_CODE_QUALITY_SECTION_END

#include <fstream>
#include <numeric>

namespace COMISO
{

namespace
{
static std::vector<size_t> make_sequence_vector(size_t _n)
{
  std::vector<size_t> v(_n);
  std::iota(v.begin(), v.end(), 0);
  return v;
}

template <class vectT> static void sort_and_compact(vectT& _v)
{
  std::sort(_v.begin(), _v.end());
  _v.erase(std::unique(_v.begin(), _v.end()), _v.end());
}

} // namespace

template <size_t DIM> class EigenLSQConstrainedSolverT<DIM>::Impl
{
public:
  Impl(std::vector<size_t>&& _var_names, ValueVector&& _fixed);
  Impl(size_t _var_nmbr, ValueVector&& _fixed)
      : Impl(make_sequence_vector(_var_nmbr), std::move(_fixed))
  {
  }

  void add_equation(
      const LinearEquation& _lin_eq);

  void add_linear_constraint(
      LinearEquation&& lin_condtr)
  {
    cnstrs_.emplace_back(lin_condtr);
  }

  const EigenLSQConstrainedSolverT<DIM>::ValueVector& fixed() const
  {
    return fixed_;
  }

  const ValueVector& solve();

private:

  bool fixed_variable(size_t _var_name,
      Point* _pt = nullptr) const
  {
    auto it = std::lower_bound(fixed_.begin(), fixed_.end(), _var_name,
        [](const Value& _val, size_t _var) { return _val.var_name < _var; });
    auto res = it != fixed_.end() && it->var_name == _var_name;
    if (_pt != nullptr && res)
      *_pt = it->point;
    return res;
  }

  bool indexed_variable(size_t _var_name, size_t& _var_idx) const
  {
    auto it = std::lower_bound(var_names_.begin(), var_names_.end(), _var_name);
    auto res = it != var_names_.end() && *it == _var_name;
    if (res)
      _var_idx = it - var_names_.begin();
    return res;
  }

  void remove_linearly_dependent_constraints();

private:
  // Fixed position points [x_i, position]
  ValueVector fixed_;

  // Map x_i (i can be any sequence of unsigned) to [0 .. n]
  std::vector<size_t> var_names_;

  // Minimization problem: (Ax-B)'*(Ax-B)
  std::vector<Eigen::Triplet<double, size_t>> A_coeff_;
  std::vector<Point> B_coeff_;
  size_t row_nmbr_ = 0;

  // Non trivial linear constraints: Sum(a_i * x_i = b_i)
  std::vector<LinearEquation> cnstrs_;

  // System solution
  ValueVector result_;
};

template <size_t DIM>
EigenLSQConstrainedSolverT<DIM>::Impl::Impl(
    std::vector<size_t>&& _var_names, ValueVector&& _fixed)
    : fixed_(std::move(_fixed)), var_names_(std::move(_var_names))
{
  sort_and_compact(fixed_);
  sort_and_compact(var_names_);
  auto new_end = std::remove_if(var_names_.begin(), var_names_.end(),
      [this](size_t _var_name) { return fixed_variable(_var_name); });
  var_names_.erase(new_end, var_names_.end());
}

template <size_t DIM>
void EigenLSQConstrainedSolverT<DIM>::Impl::add_equation(
    const LinearEquation& _lin_eq)
{
  B_coeff_.resize(row_nmbr_ + 1);
  auto& b = B_coeff_[row_nmbr_];
  b = _lin_eq.const_term;
  bool relevant_equation = false;
  for (auto& lnr_trm : _lin_eq.linear_terms)
  {
    size_t var_ind;
    Point pt;
    if (indexed_variable(lnr_trm.var_name, var_ind))
    {
      A_coeff_.push_back({row_nmbr_, var_ind, lnr_trm.coeff});
      relevant_equation = true;
    }
    else if (fixed_variable(lnr_trm.var_name, &pt))
    {
      for (auto i = 0; i < pt.size(); ++i)
        b[i] -= lnr_trm.coeff * pt[i];
    }
    else
      COMISO_THROW(LSQC_UNEXPECTED_VARIABLE);
  }
  if (relevant_equation)
    ++row_nmbr_;
}

// Remove linearly dependent constraints.
// Method: create a matrix with all constraints in a column all equation
// constraints. The constant terms are in the last DIM rows of the equation. In
// order to remove the dependent constraints, we factorize the matrix using
// FullPivLU, then, if the detected rank is less than the number of original
// constraints (so some of them are linearly dependent), we use the image method
// to get a set of linearly independent columns and we convert back into a
// reduced set of linear constraints.
template <size_t DIM>
void EigenLSQConstrainedSolverT<
    DIM>::Impl::remove_linearly_dependent_constraints()
{
  if (cnstrs_.empty())
    return;
  std::vector<size_t> var_name_to_idx;
  for (const auto& cnstr : cnstrs_)
  {
    for (const auto& lnr_trm : cnstr.linear_terms)
      var_name_to_idx.push_back(lnr_trm.var_name);
  }
  sort_and_compact(var_name_to_idx);
  auto get_index = [&var_name_to_idx](size_t var_name)
  {
    return std::lower_bound(
               var_name_to_idx.begin(), var_name_to_idx.end(), var_name) -
           var_name_to_idx.begin();
  };
  Eigen::MatrixXd M(var_name_to_idx.size() + DIM, cnstrs_.size());
  M.setZero();
  for (size_t i = 0; i < cnstrs_.size(); ++i)
  {
    for (const auto& monm : cnstrs_[i].linear_terms)
      M(get_index(monm.var_name), i) = monm.coeff;
    const auto& pt = cnstrs_[i].const_term;
    auto j = var_name_to_idx.size();
    for (const auto coord : pt)
      M(j++, i) = coord;
  }
  Eigen::FullPivLU<Eigen::MatrixXd> lu(M);
  auto rank = lu.rank();
  if (static_cast<size_t>(rank) >= cnstrs_.size())
    return;
  Eigen::MatrixXd CxD = lu.image(M);
  cnstrs_.resize(CxD.cols());
  for (size_t j = 0; j < static_cast<size_t>(CxD.cols()); ++j)
  {
    cnstrs_[j].linear_terms.clear();
    for (size_t i = 0; i < var_name_to_idx.size(); ++i)
    {
      if (CxD(i, j) != 0)
        cnstrs_[j].linear_terms.push_back({var_name_to_idx[i], CxD(i, j)});
    }
    auto k = var_name_to_idx.size();
    for (auto& coord : cnstrs_[j].const_term)
      coord = CxD(k++, j);
  }
}

template <size_t DIM>
const typename EigenLSQConstrainedSolverT<DIM>::ValueVector&
EigenLSQConstrainedSolverT<DIM>::Impl::solve()
{
  remove_linearly_dependent_constraints();
  using SparseMatrix = Eigen::SparseMatrix<double>;
  using ColumnMatrix = Eigen::Matrix<double, Eigen::Dynamic, DIM>;
  const auto size_M = var_names_.size() + cnstrs_.size();
  SparseMatrix A(row_nmbr_, size_M);
  A.setFromTriplets(A_coeff_.begin(), A_coeff_.end());
  // M is symmetric. We need to compute only the Lower part for SimplicialLDLT
  // factorization
  SparseMatrix M = (A.transpose() * A).triangularView<Eigen::Lower>();
  ColumnMatrix B(row_nmbr_, DIM);
  for (size_t i = 0; i < row_nmbr_; ++i)
  {
    for (size_t j = 0; j < DIM; ++j)
      B(i, j) = B_coeff_[i][j];
  }
  ColumnMatrix N = A.transpose() * B;
  for (auto i = 0; i != cnstrs_.size(); ++i)
  {
    auto b_var_idx = i + var_names_.size();
    const auto& pt = cnstrs_[i].const_term;
    for (int j = 0; j < DIM; ++j)
      N(b_var_idx, j) = pt[j];
    for (const auto& lin_trm : cnstrs_[i].linear_terms)
    {
      size_t idx;
      Point pos;
      if (indexed_variable(lin_trm.var_name, idx))
        M.coeffRef(i + var_names_.size(), idx) += lin_trm.coeff;
      else if (fixed_variable(lin_trm.var_name, &pos))
      {
        for (int j = 0; j < DIM; ++j)
          N(b_var_idx, j) -= lin_trm.coeff * pos[j];
      }
      else
        COMISO_THROW(LSQC_UNEXPECTED_VARIABLE);
    }
  }
  Eigen::SimplicialLDLT<SparseMatrix> solver(M);
  ColumnMatrix X = solver.solve(N);
  if (solver.info() != Eigen::Success)
  {
    // We may want to try a different factorization here, for example SparseLU
    // is to be able to find a result. Nevertheless if we are here the
    // minimization problem has multiple solutions and randomly picking one of
    // them is dangerous, for example the solution may have unpleasant
    // oscillations.
    COMISO_THROW(LSQC_SINGULAR);
  }
  result_.resize(var_names_.size());
  for (auto i = 0; i < result_.size(); ++i)
  {
    result_[i].var_name = var_names_[i];
    for (auto j = 0; j < DIM; ++j)
      result_[i].point[j] = X(i, j);
  }

  return result_;
}

template <size_t DIM>
EigenLSQConstrainedSolverT<DIM>::EigenLSQConstrainedSolverT(
    size_t _var_nmbr, ValueVector&& _fixed)
    : impl_(std::make_unique<Impl>(_var_nmbr, std::move(_fixed)))
{
}

template <size_t DIM>
EigenLSQConstrainedSolverT<DIM>::EigenLSQConstrainedSolverT(
    std::vector<size_t>&& _var_names, ValueVector&& _fixed)
    : impl_(std::make_unique<Impl>(std::move(_var_names), std::move(_fixed)))
{
}

template <size_t DIM>
EigenLSQConstrainedSolverT<DIM>::~EigenLSQConstrainedSolverT() = default;

template <size_t DIM>
void EigenLSQConstrainedSolverT<DIM>::add_equation(
    const LinearEquation& _lin_eq)
{
  impl_->add_equation(_lin_eq);
}

template <size_t DIM>
void EigenLSQConstrainedSolverT<DIM>::add_linear_constraint(
    LinearEquation&& _lin_cnstr)
{
  impl_->add_linear_constraint(std::move(_lin_cnstr));
}

template <size_t DIM>
const typename EigenLSQConstrainedSolverT<DIM>::ValueVector&
EigenLSQConstrainedSolverT<DIM>::fixed_points() const
{
  return impl_->fixed();
}

template <size_t DIM>
const typename EigenLSQConstrainedSolverT<DIM>::Result&
EigenLSQConstrainedSolverT<DIM>::solve()
{
  return impl_->solve();
}

} // namespace COMISO

//=============================================================================
#endif // COMISO_EIGEN3_AVAILABLE
//=============================================================================
