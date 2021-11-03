// Copyright 2021 Autodesk, Inc. All rights reserved.

#include "EigenLSQConstrainedSolver.hh"

#include <CoMISo/Utils/CoMISoError.hh>

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#if (COMISO_EIGEN3_AVAILABLE)
//== INCLUDES =================================================================

#include <Base/Code/Quality.hh>

#include <CoMISo/Solver/Eigen_Tools.hh>

LOW_CODE_QUALITY_SECTION_BEGIN
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/SparseCholesky>
LOW_CODE_QUALITY_SECTION_END

#include <fstream>
#include <numeric>
#include <set>

namespace COMISO
{

//#define COMISO_EIGENLSQCONSTRAINEDSOLVER_DUMP_SYSTEMS

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

  void add_equation(const LinearEquation& _lnr_eq);

  void add_linear_constraints(std::vector<LinearEquation>& _lnr_cnstrs);

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

#ifdef COMISO_EIGENLSQCONSTRAINEDSOLVER_DUMP_SYSTEMS
  bool fixed_variable_index(size_t _var_name, size_t& _var_idx) const
  {
    auto it = std::lower_bound(fixed_.begin(), fixed_.end(), _var_name,
        [](const Value& _val, size_t _var) { return _val.var_name < _var; });
    auto res = it != fixed_.end() && it->var_name == _var_name;
    if (res)
      _var_idx = var_names_.size() + (it - fixed_.begin());
    return res;
  }
#endif

  bool indexed_variable(size_t _var_name, size_t& _var_idx) const
  {
    auto it = std::lower_bound(var_names_.begin(), var_names_.end(), _var_name);
    auto res = it != var_names_.end() && *it == _var_name;
    if (res)
      _var_idx = it - var_names_.begin();
    return res;
  }

  using PermutationMatrix =
      Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, size_t>;

  PermutationMatrix compute_permutation_matrix();

private:
  // Fixed position points [x_i, position]
  ValueVector fixed_;

  // Constraints in explicit form, so a variable is a linear function of other
  // variables. Any entry of the map means that
  // variable[entry.first] = entry.second, with entry.second = Sum(c_i *
  // variable(x_i)) + d_i. Any variable that is made explicit can appear only
  // as a key in the map. Variables that appear on the right cannot be keys.
  std::map<size_t, LinearEquation> subst_vars_;

  // Map x_i (i can be any sequence of unsigned) to [0 .. n]
  std::vector<size_t> var_names_;

  // Minimization problem: (Ax-B)'*(Ax-B)
  std::vector<Eigen::Triplet<double, size_t>> A_coeff_;
  std::vector<Point> B_coeff_;
  size_t row_nmbr_ = 0;


#ifdef COMISO_EIGENLSQCONSTRAINEDSOLVER_DUMP_SYSTEMS
  // count number of constraint sets (used for output filenames)
  int n_constraint_sets_ = 0;
  // id for writing out multiple systems to individual files
  int solver_id_;

  // add equation reshapes the input to eliminate fixed variables.
  // store original coefficients for debug output
  std::vector<Eigen::Triplet<double, size_t>> A_coeff_with_fixed_;
  std::vector<Point> B_coeff_with_fixed_;
#endif

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

#ifdef COMISO_EIGENLSQCONSTRAINEDSOLVER_DUMP_SYSTEMS
  static int solver_id = 0;
  solver_id_ = solver_id++;

  // write vector of fixed variables
  {
    std::vector<Eigen::Triplet<double>> trips;
    for (int i = 0; i < fixed_.size(); ++i)
    {
      for (int j = 0; j < DIM; ++j)
        trips.emplace_back(i, j, fixed_[i].point[j]);
    }
    Eigen::SparseMatrix<double> fixed(fixed_.size(), 3);
    fixed.setFromTriplets(trips.begin(), trips.end());
    std::string filename = "EigenLSQConstrainedSolver_" +
                           std::to_string(solver_id_) + "_fixed.mtx";
    COMISO_EIGEN::write_matrix(filename, fixed);
  }
#endif
}

template <size_t DIM>
void EigenLSQConstrainedSolverT<DIM>::Impl::add_equation(
    const LinearEquation& _lnr_eq)
{
  B_coeff_.resize(row_nmbr_ + 1);
  auto& b = B_coeff_[row_nmbr_];
  b = _lnr_eq.const_term;
  bool relevant_equation = false;

#ifdef COMISO_EIGENLSQCONSTRAINEDSOLVER_DUMP_SYSTEMS
  B_coeff_with_fixed_.push_back(_lnr_eq.const_term);
  size_t A_coeff_with_fixed_size_before = A_coeff_with_fixed_.size();
#endif
  for (auto& lnr_trm : _lnr_eq.linear_terms)
  {
    size_t var_ind;
    Point pt;
    if (indexed_variable(lnr_trm.var_name, var_ind))
    {
      A_coeff_.push_back({row_nmbr_, var_ind, lnr_trm.coeff});

#ifdef COMISO_EIGENLSQCONSTRAINEDSOLVER_DUMP_SYSTEMS
      A_coeff_with_fixed_.push_back({row_nmbr_, var_ind, lnr_trm.coeff});
#endif
      relevant_equation = true;
    }
    else if (fixed_variable(lnr_trm.var_name, &pt))
    {
      for (auto i = 0; i < pt.size(); ++i)
        b[i] -= lnr_trm.coeff * pt[i];
#ifdef COMISO_EIGENLSQCONSTRAINEDSOLVER_DUMP_SYSTEMS
      fixed_variable_index(lnr_trm.var_name, var_ind);
      A_coeff_with_fixed_.push_back({row_nmbr_, var_ind, lnr_trm.coeff});
#endif
    }
    else
      COMISO_THROW(LSQC_UNEXPECTED_VARIABLE);
  }
  if (relevant_equation)
    ++row_nmbr_;
  else
  {
    B_coeff_.resize(row_nmbr_);
#ifdef COMISO_EIGENLSQCONSTRAINEDSOLVER_DUMP_SYSTEMS
    B_coeff_with_fixed_.resize(B_coeff_with_fixed_.size() - 1);
    A_coeff_with_fixed_.resize(A_coeff_with_fixed_size_before);
#endif
  }
}

template <size_t DIM>
void EigenLSQConstrainedSolverT<DIM>::Impl::add_linear_constraints(
    std::vector<LinearEquation>& _lnr_cnstrs)
{
  if (_lnr_cnstrs.empty())
    return;

#ifdef COMISO_EIGENLSQCONSTRAINEDSOLVER_DUMP_SYSTEMS
  // write out data
  {
    int n = static_cast<int>(var_names_.size() + fixed_.size());
    Eigen::SparseMatrix<double, Eigen::RowMajor> C(_lnr_cnstrs.size(), n);
    Eigen::SparseMatrix<double, Eigen::ColMajor> rhs(n, DIM);
    std::vector<Eigen::Triplet<double>> trips_C;
    std::vector<Eigen::Triplet<double>> trips_rhs;
    for (int i = 0; i < _lnr_cnstrs.size(); ++i)
    {
      const auto& eq = _lnr_cnstrs[i];
      for (const auto& term : eq.linear_terms)
      {
        size_t var_idx = 0;
        if (indexed_variable(term.var_name, var_idx) ||
            fixed_variable_index(term.var_name, var_idx))
        {
          trips_C.emplace_back(i, static_cast<int>(var_idx), term.coeff);
        }
      }
      for (int j = 0; j < DIM; ++j)
        trips_rhs.emplace_back(i, j, eq.const_term[j]);
    }
    C.setFromTriplets(trips_C.begin(), trips_C.end());
    rhs.setFromTriplets(trips_rhs.begin(), trips_rhs.end());

    std::string filename = "EigenLSQConstrainedSolver_" +
                           std::to_string(solver_id_) +
                           "_constraint_set_";
    filename += std::to_string(n_constraint_sets_++);
    COMISO_EIGEN::write_matrix(filename + "_C.mtx", C);
    COMISO_EIGEN::write_matrix(filename + "_rhs.mtx", rhs);
  }
#endif


  // Collect the variables used by the set of constraints and substitute fixed
  // variables with their value.
  std::vector<size_t> idx_to_var_name;
  for (auto& cnstr : _lnr_cnstrs)
  {
    for (auto& lnr_trm : cnstr.linear_terms)
    {
      Point pos;
      if (!fixed_variable(lnr_trm.var_name, &pos))
        idx_to_var_name.push_back(lnr_trm.var_name);
      else
      {
        for (size_t i = 0; i < DIM; ++i)
          cnstr.const_term[i] -= lnr_trm.coeff * pos[i];
        lnr_trm.coeff = 0;
      }
    }
    cnstr.linear_terms.erase(
        std::remove_if(cnstr.linear_terms.begin(), cnstr.linear_terms.end(),
            [](const LinearTerm& _term) { return _term.coeff == 0; }),
        cnstr.linear_terms.end());
    COMISO_THROW_if(cnstr.infeasible(1.e-6), LSQC_INFEASIBLE);
  }
  // Remove empty constraints
  _lnr_cnstrs.erase(std::remove_if(_lnr_cnstrs.begin(), _lnr_cnstrs.end(),
                        [](const LinearEquation& _lnr_eq)
                        { return _lnr_eq.linear_terms.empty(); }),
      _lnr_cnstrs.end());

  sort_and_compact(idx_to_var_name);
  auto get_index = [&idx_to_var_name](size_t _var_name)
  {
    return std::lower_bound(
               idx_to_var_name.begin(), idx_to_var_name.end(), _var_name) -
           idx_to_var_name.begin();
  };

  // Compute 2 matrices C and D that can express the linear constraints as
  // C * X = D (dimensions C(n, m), X(n, 1), D(n, DIM)
  // Here X is the unknown vector. X(i) is the variable with name idx_to_var_name(i).
  Eigen::MatrixXd C(_lnr_cnstrs.size(), idx_to_var_name.size()),
      D(_lnr_cnstrs.size(), DIM);
  C.setZero();
  D.setZero();
  for (size_t i = 0; i < _lnr_cnstrs.size(); ++i)
  {
    for (const auto& monm : _lnr_cnstrs[i].linear_terms)
      C(i, get_index(monm.var_name)) = monm.coeff;
    const auto& pt = _lnr_cnstrs[i].const_term;
    for (auto j = 0; j < DIM; ++j)
      D(i, j) = pt[j];
  }

  // Compute the QR factorization of C (C = Q * R * P_inv), with
  // Q a unit matrix (Q' = inverse(Q))
  // P_inv is the inverse of a permutation matrix P provided by the
  // factorization. R an upper trapezoidal matrix. R is 0 in last (n - rank)
  // rows if rank < n. R can be split into 4 matrices:
  // R = [R1 R2]
  //     [ 0  0]
  // with R1(rank, rank) square and upper triangular.
  // Q can be split in Q = [Q1 Q2], with Q1(n, rank), Q2(n, n - rank)
  // Note: Q1 * Q1' = I, Q1' * Q1 = I, ...
  // Note: Q2 is multiplied by 0, so disappear (is the null space of C * P)
  // So: C * X = D  ==>  Q * R * P_inv * X = D ==> ...
  // Define:
  //   Y = P_inv * X = |Y1|  (Y1(rank, DIM), Y2(n - rank, DIM))
  //                   |Y2|
  // ... ==> Q1 * (R1 * Y1 + R2 * Y2) = D  (Remember Q1' = Q1_inverse)
  // R1 * Y1 + R2 * Y2 = Q1' * D
  // R1 * Y1 = -R2 * Y2 + Q1' * D
  // Y1 = -R1_inverse * R2 * Y2 + R1_inverse * Q1' * D
  // Y1 = Y2_coeff * Y2 + Y_cnst_term
  // So we can find the substitution for the variables in Y1 in the LSQ problem.

  Eigen::ColPivHouseholderQR<Eigen::Ref<Eigen::MatrixXd>> qr_dec(C);
  Eigen::MatrixXd Q = qr_dec.matrixQ();
  const auto rank = qr_dec.rank();
  const auto tol = qr_dec.threshold();
  if (rank < C.rows())
  {
    // Q.rightCols(C.rows() - rank) is the null space  (or kernel) of C * P.
    // If the projection of D on this null space is not zero, the system is
    // infeasible.
    // The projection can be computed as a dot product because the columns of Q
    // are orthonormal vectors.
    COMISO_THROW_if(
        !(Q.rightCols(C.rows() - rank).transpose() * D).isZero(tol),
        LSQC_INFEASIBLE);
    Q = Q.leftCols(rank).eval();
  }
  auto P = qr_dec.colsPermutation();
  Eigen::MatrixXd R1 = qr_dec.matrixR()
                            .topLeftCorner(rank, rank)
                            .template triangularView<Eigen::Upper>();
  Eigen::MatrixXd R1_inv = R1.inverse();
  Eigen::MatrixXd R2 = qr_dec.matrixR()
                           .topRightCorner(rank, idx_to_var_name.size() - rank);
  Eigen::MatrixXd Y_cnst_term = R1_inv * Q.transpose() * D;
  Eigen::MatrixXd Y2_coeff = -R1_inv * R2;

  // Move the equation into a map of linear equations that uses the original
  // variable names.
  for (auto i = 0; i < rank; ++i)
  {
    auto var_name = idx_to_var_name[P.indices()[i]];
    auto& lnr_eq = subst_vars_[var_name];
    for (auto j = 0; j < Y_cnst_term.cols(); ++j)
      lnr_eq.const_term[j] = Y_cnst_term(i, j);
    for (auto j = 0; j < Y2_coeff.cols(); ++j)
    {
      // If the coefficient is smaller than the threshold used by the QR
      // decomposition to find the rank, we ignore the term. This will reduce
      // the number of coefficients in the sparse matrix inside the
      // EigenLSQConstrainedSolverT::solve
      if (std::fabs(Y2_coeff(i, j)) > tol)
      {
        lnr_eq.linear_terms.push_back(LinearTerm(
            {idx_to_var_name[P.indices()[j + rank]], Y2_coeff(i, j)}));
      }
    }
  }
}

template <size_t DIM>
const typename EigenLSQConstrainedSolverT<DIM>::ValueVector&
EigenLSQConstrainedSolverT<DIM>::Impl::solve()
{
  // We solve min||A * X - B||^2, where a subset Xc of the variables X = {x_i}
  // are eliminated by a set of equality constraints CX = D
  // (see \ref set_linear_constraints) and substituted with linear combinations
  // of the other {x_i} not in Xc
  //
  // Method:
  // Find a permutation matrix such that Y = P * X = [Y1 Y2]'. Here we assume Y2
  // contains all the substituted variables (in subst_vars_) ==> Y2 = F * y1 + G.
  //
  // The algorithm compute Y1 as LSQ unconstrained problem and then Y2 by
  // back substitution.
  //
  // The expression A * X - B can be simplified as follow:
  // A * X - B = A * P_inv * P * X - B = A * P_inv * Y - B =
  // = A * P_inv * [Y1 Y2]' - B = ...
  // ( Let's call AP = A * P_inv = [AP1 AP2] )
  // ... = AP * [Y1 Y2]' - B = [AP1 AP2] * [Y1 Y2]' = AP1 * Y1 + AP2 * Y2 - B =
  // = AP1 * Y1 + AP2 * (F * y1 + G) - B =
  // = (AP1 + AP2 * F) * Y1 - (B - AP2 * G)
  // In the end we can find Y2 solving the unconstrained LSQ problem:
  // min ||A_reduced * Y1 - B_reduced||^2, with
  //                      A_reduced = AP1 + AP2 * F and
  //                      B_reduced = B0 - AP2 * G
  // Then Y2 = F * Y1 + G;    X = P_inv * [Y1 Y2]'
  //
  // Matrices size:
  //   A(m, n), B(m, DIM), subst_vars_(k, n-k)
  //   Y1(n-k, DIM), Y2(k, DIM), AP1(m, n-k), AP2(m, k), F(k, n-k), G(k, DIM)
  //   A_reduced(m, n-k), B_reduced(m, DIM)


#ifdef COMISO_EIGENLSQCONSTRAINEDSOLVER_DUMP_SYSTEMS
  // write data for tests
  {

    int n = static_cast<int>(var_names_.size() + fixed_.size());
    Eigen::SparseMatrix<double> A(row_nmbr_, n);
    A.setFromTriplets(A_coeff_with_fixed_.begin(), A_coeff_with_fixed_.end());
    auto flnm_prfx = "EigenLSQConstrainedSolver_" + std::to_string(solver_id_);
    COMISO_EIGEN::write_matrix(flnm_prfx + "_A.mtx", A);

    Eigen::MatrixXd B(row_nmbr_, DIM);
    for (int i = 0; i < DIM; ++i)
    {
      for (int j = 0; j < B_coeff_with_fixed_.size(); ++j)
        B(j, i) = B_coeff_with_fixed_[j][i];
    }
    COMISO_EIGEN::write_matrix(flnm_prfx + "_b.mtx", B);
  }
#endif

  using PermutationMatrix =
      Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, size_t>;
  PermutationMatrix P = compute_permutation_matrix();
  PermutationMatrix P_inv = P.inverse();

  auto variable_name_to_index = [this, &P_inv](size_t var_name)
  {
    size_t idx;
    indexed_variable(var_name, idx);
    return P_inv.indices()[idx];
  };

  using SparseMatrix = Eigen::SparseMatrix<double>;
  using ColumnMatrix = Eigen::Matrix<double, Eigen::Dynamic, DIM>;
  size_t unconstr_size = var_names_.size() - subst_vars_.size();

  // Parse the solved constraints to build the matrices F and G
  SparseMatrix F(subst_vars_.size(), unconstr_size);
  ColumnMatrix G(subst_vars_.size(), DIM);
  {
    std::vector<Eigen::Triplet<double, size_t>> F_coeff;
    for (const auto& slvd_cnstr : subst_vars_)
    {
      auto row = variable_name_to_index(slvd_cnstr.first) - unconstr_size;
      for (const auto& lnr_prt : slvd_cnstr.second.linear_terms)
      {
        auto col = variable_name_to_index(lnr_prt.var_name);
        F_coeff.push_back({row, col, lnr_prt.coeff});
      }
      for (size_t j = 0; j < DIM; ++j)
        G(row, j) = slvd_cnstr.second.const_term[j];
    }
    F.setFromTriplets(F_coeff.begin(), F_coeff.end());
  }

  // Use the substitution to compute matrices A_reduced and B_reduced
  SparseMatrix A_reduced;
  ColumnMatrix B_reduced;
  {
    // Update the triplet data including the column permutation. This is
    // equivalent to multiply A * P_inv;
    for (auto& tri : A_coeff_)
    {
      tri = Eigen::Triplet<double, size_t>(
          tri.row(), P_inv.indices()[tri.col()], tri.value());
    }

    const auto size_M = var_names_.size();
    SparseMatrix AP(row_nmbr_, size_M);
    AP.setFromTriplets(A_coeff_.begin(), A_coeff_.end());
    auto AP2 = AP.rightCols(subst_vars_.size());
    A_reduced = AP.leftCols(unconstr_size) + AP2 * F;

    ColumnMatrix B(row_nmbr_, DIM);
    for (size_t i = 0; i < row_nmbr_; ++i)
    {
      for (size_t j = 0; j < DIM; ++j)
        B(i, j) = B_coeff_[i][j];
    }
    B_reduced = B - AP2 * G;
  }

  // Solve the unconstrained LSQ problem min ||A_reduced * Y1 - B_reduced||^2
  ColumnMatrix Y1;
  {
    // M is symmetric. We need to compute only the Lower part for SimplicialLDLT
    // factorization
    SparseMatrix M =
        (A_reduced.transpose() * A_reduced).triangularView<Eigen::Lower>();

    ColumnMatrix N = A_reduced.transpose() * B_reduced;
    Eigen::SimplicialLDLT<SparseMatrix> solver(M);
    Y1 = solver.solve(N);
    if (solver.info() != Eigen::Success)
    {
      // We may want to try a different factorization here, for example SparseLU
      // is to be able to find a result. Nevertheless if we are here the
      // minimization problem has multiple solutions and randomly picking one of
      // them is dangerous, for example the solution may have unpleasant
      // oscillations.
      COMISO_THROW(LSQC_SINGULAR);
    }
  }

  ColumnMatrix Y2 = F * Y1 + G;

  // Fill the result vector
  result_.resize(var_names_.size());
  for (auto i = 0; i < result_.size(); ++i)
  {
    result_[i].var_name = var_names_[i];
    size_t idx = P_inv.indices()[i];
    auto Y_ptr = &Y1;
    if (idx >= static_cast<size_t>(Y1.rows()))
    {
      idx -= Y1.rows();
      Y_ptr = &Y2;
    }
    for (auto j = 0; j < DIM; ++j)
      result_[i].point[j] = (*Y_ptr)(idx, j);
  }

#ifdef COMISO_EIGENLSQCONSTRAINEDSOLVER_DUMP_SYSTEMS
  // write solution to file
  {
    int n = static_cast<int>(var_names_.size() + fixed_.size());
    Eigen::MatrixXd sol(n, DIM);
    size_t var_idx = 0;
    for (const auto& val : result_)
    {
      if (indexed_variable(val.var_name, var_idx))
      {
        for (int i = 0; i < DIM; ++i)
          sol(static_cast<int>(var_idx), i) = val.point[i];
      }
      else
      {
        COMISO_THROW_TODO("Shouldn't all solution values be in the set of"
                          "indexed variables?");
      }
    }

    for (const auto& val : fixed_)
    {
      if (fixed_variable_index(val.var_name, var_idx))
      {
        for (int i = 0; i < DIM; ++i)
          sol(static_cast<int>(var_idx), i) = val.point[i];
      }
      else
      {
        COMISO_THROW_TODO(
            "Shouldn't all fixed values be in the set of fixed variables?");
      }
    }
    auto filename =
        "EigenLSQConstrainedSolver_" + std::to_string(solver_id_) + "_sol.mtx";
    COMISO_EIGEN::write_matrix(filename, sol);
  }
#endif

  return result_;
}

template <size_t DIM>
typename EigenLSQConstrainedSolverT<DIM>::Impl::PermutationMatrix
EigenLSQConstrainedSolverT<DIM>::Impl::compute_permutation_matrix()
{
  // Compute a permutation matrix P in such a way that
  // Y = P * X = [Y1 Y2]'
  // Where Y2 contains all variables that can be substituted using the
  // constrains, i.e., from the constrains equations:
  // Y2 = F * Y1 + G
  // P will be an array of indices that index(i) is the index j such that
  // x_i <==> y_j
  // All the variable solved in subst_vars_ must have the biggest indices.

  // Find the indices of the solved variables
  std::set<size_t> solved_vars;
  for (const auto& slv_cnstr : subst_vars_)
  {
    size_t idx;
    indexed_variable(slv_cnstr.first, idx);
    solved_vars.insert(idx);
  }

  // Compute the index map putting at the beginning all the indices not in
  // solved_vars, and then the others.
  Eigen::Matrix<size_t, 1, Eigen::Dynamic> indices(var_names_.size());
  size_t j_ins = 0;
  auto it = solved_vars.begin();
  for (size_t i = 0; i < var_names_.size(); ++i)
  {
    if (it == solved_vars.end() || i < *it)
      indices(0, j_ins++) = i;
    else if (it != solved_vars.end())
      ++it;
  }
  for (const auto& sv : solved_vars)
    indices(0, j_ins++) = sv;

  // Create the inverse of the permutation matrix.
  return PermutationMatrix(indices);
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
    const LinearEquation& _lnr_eq)
{
  impl_->add_equation(_lnr_eq);
}


template <size_t DIM>
void EigenLSQConstrainedSolverT<DIM>::add_linear_constraints(
    std::vector<LinearEquation>& _lnr_cnstrs)
{
  impl_->add_linear_constraints(_lnr_cnstrs);
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
