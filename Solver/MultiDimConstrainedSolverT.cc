// Copyright 2022 Autodesk, Inc. All rights reserved.

#define COMISO_MULTIDIMCONSTRAINEDSOLVERT_C

#include "MultiDimConstrainedSolverT.hh"

#include <CoMISo/Solver/ConstrainedSolver.hh>

#include <Eigen/Sparse>

namespace COMISO
{

template <int DIM>
class MultiDimConstrainedSolverT<DIM>::Impl
{
  using RowMatrix = ConstrainedSolver::RowMatrix;
  using Vector = ConstrainedSolver::Vector;

  using TripletVector = std::vector<Eigen::Triplet<double>>;

public:
  Impl()
      : var_nmbr_(0)
  {
  }

  void add_equation(const LinearEquation& _eq) { add(A_, b_, _eq); }
  void add_constraint(const LinearEquation& _eq) { add(C_, d_, _eq); }

  void set_integers(IndexVector _int_var_indcs)
  {
    int_var_indcs_ = std::move(_int_var_indcs);
  }

  // Solve problem. Resets all equations, constraints and integer constraints
  void solve(PointVector& _result);

private:

  // Create matrix from triplets with first dim of rhs in last column
  static void to_matrix(TripletVector& _triplets, size_t _col_num,
      const PointVector& _rhs, RowMatrix& _mat)
  {
    const auto row_num = _rhs.size();
    for (size_t i = 0; i < row_num; ++i)
      _triplets.emplace_back((int)i, (int)_col_num, -_rhs[i][0]);
    _mat.resize(row_num, _col_num + 1);
    _mat.reserve(_triplets.size());
    _mat.setFromTriplets(_triplets.begin(), _triplets.end());
  }

  // Extract simple one dimensional vector for column _col of a point vector
  static void to_vector(const PointVector& _point_vec, int _col, Vector& _vec)
  {
    _vec.resize(_point_vec.size());
    for (int i = 0; i < _point_vec.size(); ++i)
      _vec[i] = _point_vec[i][_col];
  }

  // Exchange rhs (last column) of _mat with column _col of _rhs
  static void update_rhs(RowMatrix& _mat, const PointVector& _rhs, int _col)
  {
    for (int i = 0; i < _mat.rows(); ++i)
      _mat.coeffRef(i, _mat.cols() - 1) = -_rhs[i][_col];
  }

  // Add a linear equation into triplets and rhs vector
  void add(TripletVector& _M, PointVector& _rhs, const LinearEquation _eq)
  {
    if (_eq.linear_terms.empty())
      return; // don't add empty equations
    const auto row_idx = _rhs.size();
    for (const auto& term : _eq.linear_terms)
    {
      _M.emplace_back((int)row_idx, (int)term.var_name, term.coeff);
      var_nmbr_ = std::max(var_nmbr_, term.var_name+1);
    }
    _rhs.push_back(_eq.const_term);
  }

  // Reset all data accumulated by calling add_equation add_constraint and
  // set_integer
  void reset()
  {
    A_.clear();
    b_.clear();
    C_.clear();
    d_.clear();
    int_var_indcs_.clear();
    var_nmbr_ = 0;
  }

  TripletVector A_; // Matrix A defining the minimization objective ||Ax-b||^2
  PointVector   b_; // b of the minimization objective ||Ax-b||^2
  TripletVector C_; // Matrix C defining the linear equality constraints Cx=d
  PointVector   d_; // rhs of equality constraints Cx = d
  IndexVector   int_var_indcs_; // List of variables which should be rounded to
                                // integers.

  size_t var_nmbr_; // Number of variables in the problem
};

template <int DIM>
void
MultiDimConstrainedSolverT<DIM>::Impl::solve(PointVector& _result)
{
  ConstrainedSolver solver;

  // create systems with first dimension of _b,_d as rhs in last column of A,C
  RowMatrix A;
  RowMatrix C;
  to_matrix(A_, var_nmbr_, b_, A);
  to_matrix(C_, var_nmbr_, d_, C);

  ConstrainedSolver::Vector solution(var_nmbr_);

  // solve system for first dimension
  solver.solve(C, A, solution, int_var_indcs_, 0.0, false);
  _result.resize(var_nmbr_);

  // store result
  for (int i = 0; i < var_nmbr_; ++i)
    _result[i][0] = solution[i];

  // solve system for remaining dimensions
  for (int k = 1; k < DIM; ++k)
  {
    // exchange right hand sides
    update_rhs(A, b_, k);
    Vector rhs;
    to_vector(d_, k, rhs);

    // resolve for updated rhs
    solver.resolve(A, solution, &rhs);

    // store result
    for (int j = 0; j < var_nmbr_; ++j)
      _result[j][k] = solution[j];
  }

  reset();
}


template <int DIM>
MultiDimConstrainedSolverT<DIM>::MultiDimConstrainedSolverT()
    : impl_(new Impl())
{
}

template <int DIM>
MultiDimConstrainedSolverT<DIM>::~MultiDimConstrainedSolverT()
{
  if (impl_ != nullptr)
    delete impl_;
}

template <int DIM>
void
MultiDimConstrainedSolverT<DIM>::add_equation(const LinearEquation& _eq)
{
  impl_-> add_equation(_eq);
}

template <int DIM>
void
MultiDimConstrainedSolverT<DIM>::add_constraint(const LinearEquation& _eq)
{
  impl_->add_constraint(_eq);
}

template <int DIM>
void
MultiDimConstrainedSolverT<DIM>::set_integers(IndexVector _int_var_indcs)
{
  impl_->set_integers(std::move(_int_var_indcs));
}

template <int DIM>
void
MultiDimConstrainedSolverT<DIM>::solve(PointVector& _result)
{
  return impl_->solve(_result);
}

}//namespace COMISO


