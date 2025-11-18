//=============================================================================
//
//  CLASS OSQPSolver - IMPLEMENTATION
//
//=============================================================================

//== INCLUDES =================================================================

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_OSQP_AVAILABLE
//=============================================================================
#include "OSQPSolver.hh"
#include "OSQPEigen.hh"

//#include <osqp/osqp.h>

#include <CoMISo/Utils/CoMISoError.hh>
#include <CoMISo/Utils/StopWatch.hh>
#include <CoMISo/NSolver/LazyConstraintSolver.hh>

#include <Base/Debug/DebTime.hh>

#include <Eigen/Sparse>

//== NAMESPACES ===============================================================

namespace COMISO {

//== IMPLEMENTATION ========================================================== 

namespace
{

using ContraintVector = OSQPSolver::ContraintVector;
using SMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, c_int>;

void regularize_hessian(NProblemInterface::SMatrixNP& _H, const double _reg_factor = 1e-8)
{
  NProblemInterface::SMatrixNP id;
  id.resize(_H.rows(), _H.cols());
  id.setIdentity();
  auto diag = _H.diagonal();
  _H = _H + _reg_factor * diag.sum() / diag.rows() * id; // perturbation?!
}

NProblemInterface::SMatrixNP get_hessian(NProblemInterface* _problem)
{
  std::vector<double> zero(_problem->n_unknowns(), 0);
  NProblemInterface::SMatrixNP H;
  _problem->eval_hessian(zero.data(), H);
  regularize_hessian(H, 0.0); // make sure that diagonal entries are created!
  return H;
}

Eigen::VectorXd get_linear_energy_coefficients(NProblemInterface* _problem)
{
  std::vector<double> zero(_problem->n_unknowns(), 0);
  Eigen::VectorXd q;
  q.resize(_problem->n_unknowns());
  _problem->eval_gradient(zero.data(), q.data());
  return q;
}
double get_constant_coefficient(NProblemInterface* _problem)
{
  std::vector<double> zero(_problem->n_unknowns(), 0);
  return _problem->eval_f(zero.data());
}

void get_constraints(int _n_cols, const ContraintVector& _constraints,
    SMatrix& _C, Eigen::VectorXd& _lower_bounds,
    Eigen::VectorXd& _upper_bounds)
{
  size_t n_rows = _constraints.size();
  _C.resize(n_rows, _n_cols);
  _lower_bounds.resize(n_rows);
  _upper_bounds.resize(n_rows);

  std::vector<double> x(_n_cols, 0.0);
  std::vector<Eigen::Triplet<double>> triplets;
  int current_row = 0;
  for (const auto& c : _constraints)
  {
    if (!c->is_linear())
    {
      DEB_error(
          "OSQP non-linear constraints are not supported and thus ignored.");
      continue;
    }

    NConstraintInterface::SVectorNC gc;
    c->eval_gradient(x.data(), gc);
    for (NConstraintInterface::SVectorNC::InnerIterator v_it(gc); v_it; ++v_it)
      triplets.emplace_back(current_row, v_it.index(), v_it.value());

    const auto b = c->eval_constraint(x.data());
    switch (c->constraint_type())
    {
    case NConstraintInterface::NC_EQUAL:
      _lower_bounds[current_row] = -b;
      _upper_bounds[current_row] = -b;
      break;

    case NConstraintInterface::NC_LESS_EQUAL:
      _lower_bounds[current_row] = -std::numeric_limits<double>::max();
      _upper_bounds[current_row] = -b;
      break;

    case NConstraintInterface::NC_GREATER_EQUAL:
      _lower_bounds[current_row] = -b;
      _upper_bounds[current_row] = std::numeric_limits<double>::max();
      break;
    }

    ++current_row;
  }

  _C.setFromTriplets(triplets.begin(), triplets.end());
}

void throw_solve_failure(const c_int _status)
{
  DEB_enter_func;
  DEB_warning(1, " OSQP solve failure code is " << int(_status));
  switch (_status)
  {
  case OSQP_MAX_ITER_REACHED:
    COMISO_THROW(QP_MAXIMUM_ITERATIONS_EXCEEDED);
  // case Ipopt::NonIpopt_Exception_Thrown: // TODO: handle interrupts
  //  // this could be due to a thrown PROGRESS_ABORTED exception, ...
  //  PROGRESS_RESUME_ABORT; // ... so check if we need to resume it
  default:
    COMISO_THROW(QP_OPTIMIZATION_FAILED);
  }
}

void check_solve_status(const c_int _status)
{
  if (_status != OSQP_SOLVED)
    throw_solve_failure(_status);
}


class Impl // manage OSQP objects
{
public:

  Impl()
  {
    osqp_eigen_.settings().alpha = 1.0; // this value works better than the default
    osqp_eigen_.settings().max_iter = 10000;
#if COMISO_OSQP_NEW_API
    osqp_eigen_.settings().warm_starting = true;
    osqp_eigen_.settings().polishing = 1;
#else
    osqp_eigen_.settings().warm_start = true;
    osqp_eigen_.settings().polish = 1;
#endif
    osqp_eigen_.settings().polish_refine_iter = 5;
    osqp_eigen_.settings().eps_abs = 1e-5;      // absolute convergence tolerance
    osqp_eigen_.settings().eps_rel = 1e-5;      // relative convergence tolerance
    osqp_eigen_.settings().eps_prim_inf = 1e-6; // primal infeasibility tolerance
    osqp_eigen_.settings().eps_dual_inf = 1e-8;   // dual infeasibility tolerance
    // osqp_eigen_.settings().linsys_solver = MKL_PARDISO_SOLVER;
  }

  ~Impl() {}

  void solve(NProblemInterface* _problem, const ContraintVector& _constraints);
  const double* get_solution() const { return osqp_eigen_.get_x(); }
        double  objective_value() const{ return const_coeff_ + osqp_eigen_.objective_value();}

private:
  OSQPEigen osqp_eigen_;
  double const_coeff_ = 0.;
};

void Impl::solve(
    NProblemInterface* _problem, const ContraintVector& _constraints)
{
  const auto H = get_hessian(_problem);
  const auto lin_q = get_linear_energy_coefficients(_problem);
  const_coeff_ = get_constant_coefficient(_problem);



  SMatrix A; // inequality constraints
  Eigen::VectorXd lower;                  // lower bounds
  Eigen::VectorXd upper;                  // upper bounds
  get_constraints(_problem->n_unknowns(), _constraints, A, lower, upper);

  auto exitflag = osqp_eigen_.setup(H, lin_q, A, lower, upper);

  DEB_error_if(exitflag != 0, "OSQP Setup failed with exit flag " << int(exitflag));
  COMISO_THROW_if(exitflag != 0, QP_INITIALIZATION_FAILED);

  // Solve Problem
  exitflag = osqp_eigen_.solve();
  DEB_error_if(exitflag != 0, "OSQP solve() failed with exit flag " << int(exitflag));
  COMISO_THROW_if(exitflag != 0, QP_OPTIMIZATION_FAILED);


  _problem->store_result(osqp_eigen_.get_x());
}

} // namespace

//class Impl // manage OSQP objects
//{
//public:
//
//  Impl()
//  {
//    osqp_set_default_settings(&settings);
//    settings.alpha = 1.0; // this value works better than the default
//    settings.max_iter = 10000;
//    settings.warm_start = true;
//    settings.polish = 1;
//    settings.polish_refine_iter = 5;
//    settings.eps_abs = 1e-5;      // absolute convergence tolerance
//    settings.eps_rel = 1e-5;      // relative convergence tolerance
//    settings.eps_prim_inf = 1e-6; // primal infeasibility tolerance
//    settings.eps_dual_inf = 1.;   // dual infeasibility tolerance
//    // settings.linsys_solver = MKL_PARDISO_SOLVER;
//
//    data.n = 0;
//    data.m = 0;
//    data.P = nullptr;
//    data.A = nullptr;
//    data.q = nullptr;
//    data.l = nullptr;
//    data.u = nullptr;
//  }
//
//  ~Impl()
//  {
//    delete data.P;
//    delete data.A;
//    delete work;
//  }
//
//  void solve(NProblemInterface* _problem, const ContraintVector& _constraints);
//  const double* get_solution() const { return work->solution->x; }
//
//private:
//  OSQPSettings settings;
//  OSQPData data;
//  OSQPWorkspace* work = nullptr;
//};
//
//void Impl::solve(
//    NProblemInterface* _problem, const ContraintVector& _constraints)
//{
//  const auto H = get_hessian(_problem);
//  const auto lin_q = get_linear_energy_coefficients(_problem);
//
//  COMISO::NProblemInterface::SMatrixNP A; // inequality constraints
//  Eigen::VectorXd lower;                  // lower bounds
//  Eigen::VectorXd upper;                  // upper bounds
//  get_constraints(_problem->n_unknowns(), _constraints, A, lower, upper);
//
//  COMISO::NProblemInterface::SMatrixNP HupperTriangle =
//      H.triangularView<Eigen::Upper>();
//  HupperTriangle.makeCompressed();
//
//  data.n = static_cast<int>(HupperTriangle.cols()); // number of variables n
//  data.m = static_cast<int>(A.rows());              // number of constraints m
//
//  c_float* P_x =                 // the upper triangular part of the quadratic
//     HupperTriangle.valuePtr(); // cost matrix P in csc format (size n x n).
//  c_int P_nnz =
//  static_cast<int>(HupperTriangle.nonZeros()); // number of non zeros
//  std::vector<c_int> P_i(P_nnz), P_p(P_nnz);
//  for(c_int i=0; i<P_nnz; ++i)
//  {
//    P_i[i] = HupperTriangle.innerIndexPtr()[i];
//    P_p[i] = HupperTriangle.outerIndexPtr()[i];
//  }
//
////   c_int == long long int but Eigen uses int
////    c_int* P_i = (c_int*) HupperTriangle.innerIndexPtr();     // row indices
////    c_int* P_p = (c_int*) HupperTriangle.outerIndexPtr();     // column pointers
//
//    data.P = csc_matrix(data.n, data.n, P_nnz, P_x, P_i.data(), P_p.data());
//
//  c_float* A_x =
//      A.valuePtr(); // linear constraints matrix A in csc format (size m x n)
//  c_int A_nnz = static_cast<int>(A.nonZeros()); // number of non zeros
//
//  std::vector<c_int> A_i(A_nnz), A_p(A_nnz);
//  for(c_int i=0; i<P_nnz; ++i)
//  {
//    A_i[i] = A.innerIndexPtr()[i];
//    A_p[i] = A.outerIndexPtr()[i];
//  }
//
////   c_int == long long int but Eigen uses int
////  c_int* A_i = (c_int*) A.innerIndexPtr();               // number of non z
////  c_int* A_p = (c_int*) A.outerIndexPtr();               // row indices
//  data.A = csc_matrix(data.m, data.n, A_nnz, A_x, A_i.data(), A_p.data());
//
//
//  data.q = const_cast<c_float*>(
//      lin_q.data()); // dense array for linear part of cost function (size n)
//  data.l = lower.data(); // dense array for lower bound (size m)
//  data.u = upper.data(); // dense array for upper bound (size m)
//
//  auto exitflag = osqp_setup(&work, &data, &settings); // Setup workspace
//  DEB_error_if(exitflag != 0, "OSQP Setup failed with exit flag " << int(exitflag));
//  COMISO_THROW_if(exitflag != 0, QP_INITIALIZATION_FAILED);
//
//  // Solve Problem
//  exitflag = osqp_solve(work);
//  DEB_error_if(exitflag != 0, "OSQP Setup failed with exit flag " << int(exitflag));
//  COMISO_THROW_if(exitflag != 0, QP_OPTIMIZATION_FAILED);
//  check_solve_status(work->info->status_val);
//
//  _problem->store_result(work->solution->x);
//}
//
//} // namespace

void OSQPSolver::solve(
    NProblemInterface* _problem, const ContraintVector& _constraints)
{
  Impl impl;
  impl.solve(_problem, _constraints);

  obj_val_ = impl.objective_value();
}

void OSQPSolver::solve(NProblemInterface* _problem,
    const ContraintVector& _constraints,
    const ContraintVector& _lazy_constraints, double _acceptable_tolerance,
    double _almost_infeasible_threshold, int _max_passes,
    bool _final_step_with_all_constraints)
{
  Impl impl;
  const auto solve_function = [&](
      NProblemInterface* _problem, const ContraintVector _constraints)
  {
    impl.solve(_problem, _constraints);
    obj_val_ = impl.objective_value();
  };
  const auto result_function = [&impl]() { return impl.get_solution(); };

  return solve_with_lazy_constraints(solve_function, result_function, _problem,
      _constraints, _lazy_constraints, _acceptable_tolerance,
      _almost_infeasible_threshold, _max_passes,
      _final_step_with_all_constraints);

}


//-----------------------------------------------------------------------------

//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // COMISO_OSQP_AVAILABLE
//=============================================================================
