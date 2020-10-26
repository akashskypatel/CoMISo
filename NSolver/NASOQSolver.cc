//=============================================================================
//
//  CLASS NASOQSolver - IMPLEMENTATION
//
//=============================================================================

//== INCLUDES =================================================================

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
//#if COMISO_NASOQ_AVAILABLE // TODO
//=============================================================================
#include "NASOQSolver.hh"
#include <nasoq_eigen.h>

#include <CoMISo/Utils/CoMISoError.hh>
#include <CoMISo/Utils/StopWatch.hh>

#include <Base/Debug/DebTime.hh>

//== NAMESPACES ===============================================================

namespace COMISO {

//== IMPLEMENTATION ========================================================== 

void NASOQSolver::regularize_hessian(NProblemInterface::SMatrixNP& _H)
{
  NProblemInterface::SMatrixNP id;
  id.resize(_H.rows(), _H.cols());
  id.setIdentity();
  double reg_factor = 1e-8;
  auto diag = _H.diagonal();
  _H = _H + reg_factor * diag.sum()/diag.rows() * id;
}

NProblemInterface::SMatrixNP NASOQSolver::get_hessian(NProblemInterface* _problem)
{
  std::vector<double> zero(_problem->n_unknowns(), 0);
  NProblemInterface::SMatrixNP H;
  _problem->eval_hessian(zero.data(), H);
  regularize_hessian(H);
  return H;
}

Eigen::VectorXd NASOQSolver::get_q(NProblemInterface* _problem)
{
  std::vector<double> zero(_problem->n_unknowns(), 0);
  Eigen::VectorXd q;
  q.resize(_problem->n_unknowns());
  _problem->eval_gradient(zero.data(), q.data());
  return q;
}

void NASOQSolver::get_equality_constraints(int _n_cols, const std::vector<NConstraintInterface*>& _constraints, COMISO::NProblemInterface::SMatrixNP& _A, Eigen::VectorXd& _rhs)
{
  int n_rows = 0;
  for (auto c : _constraints)
    if (c->linear_equality())
      ++n_rows;
  _A.resize(n_rows, _n_cols);
  _rhs.resize(n_rows);

  std::vector<double> x(_n_cols, 0.0);
  std::vector<Eigen::Triplet<double>> triplets;
  int current_row = 0;
  for (auto c : _constraints)
    if (c->linear_equality())
    {
      NConstraintInterface::SVectorNC gc;
      c->eval_gradient(x.data(), gc);
      for(NConstraintInterface::SVectorNC::InnerIterator v_it(gc); v_it; ++v_it)
        triplets.emplace_back(current_row, v_it.index(), v_it.value());

      double b = c->eval_constraint(x.data());
      _rhs(current_row) = -b;

      ++current_row;
    }

  _A.setFromTriplets(triplets.begin(), triplets.end());
}

void NASOQSolver::get_inequality_constraints(int _n_cols, const std::vector<NConstraintInterface*>& _constraints, COMISO::NProblemInterface::SMatrixNP& _C, Eigen::VectorXd& _rhs)
{
  int n_rows = 0;
  for (auto c : _constraints)
    if (c->is_linear() && !c->linear_equality())
      ++n_rows;
  _C.resize(n_rows, _n_cols);
  _rhs.resize(n_rows);

  std::vector<double> x(_n_cols, 0.0);
  std::vector<Eigen::Triplet<double>> triplets;
  int current_row = 0;
  for (auto c : _constraints)
    if (c->is_linear() && !c->linear_equality())
    {
      // NASOQ expects all constraints to be c^T x <= rhs, multiply row with -1 for to turn >= into <= constraints
      int sign = c->constraint_type() == NConstraintInterface::NC_LESS_EQUAL ? 1 : -1;
      NConstraintInterface::SVectorNC gc;
      c->eval_gradient(x.data(), gc);
      for(NConstraintInterface::SVectorNC::InnerIterator v_it(gc); v_it; ++v_it)
        triplets.emplace_back(current_row, v_it.index(), sign * v_it.value());

      double b = sign * c->eval_constraint(x.data());
      _rhs(current_row) = -b;

      ++current_row;
    }

  _C.setFromTriplets(triplets.begin(), triplets.end());
}

bool NASOQSolver::solve(NProblemInterface* _problem, const std::vector<NConstraintInterface*>& _constraints)
{
  DEB_enter_func;

  auto H = get_hessian(_problem);
  auto q = get_q(_problem);

  COMISO::NProblemInterface::SMatrixNP A; // equality constraints
  Eigen::VectorXd b; // equality right hand side
  get_equality_constraints(_problem->n_unknowns(), _constraints, A, b);

  COMISO::NProblemInterface::SMatrixNP C; // inequality constraints
  Eigen::VectorXd d; // inequality constraints rhs
  get_inequality_constraints(_problem->n_unknowns(), _constraints, C, d);

  Eigen::VectorXd x; // solution vector of primal vars
  Eigen::VectorXd y; // solution vector of dual vars for equality    // probably, TODO: verify
  Eigen::VectorXd z; // solution vector of dual vars for inequality  // probably, TODO: verify

  COMISO::NProblemInterface::SMatrixNP HlowerTriangle = H.triangularView<Eigen::Lower>();

  auto res = nasoq::quadprog(HlowerTriangle, q, A, b, C, d, x, y, z, nullptr);

  DEB_line_if(res == nasoq::nasoq_status::Optimal, 2, "NASOQ Solver converged.");
  DEB_warning_if(res == nasoq::nasoq_status::Inaccurate, 2, "NASOQ Solver converged to an inaccurate solution."); // TODO: what does this mean?
  DEB_error_if(res == nasoq::nasoq_status::NotConverged, "NASOQ Solver did not converge!");
  DEB_error_if(res == nasoq::nasoq_status::Infeasible, "NASOQ Solver reports infeasible problem!");


  if (res == nasoq::nasoq_status::Optimal || res == nasoq::nasoq_status::Inaccurate)
  {
    _problem->store_result(x.data());
  };

  return res;
}


//-----------------------------------------------------------------------------

//=============================================================================
} // namespace COMISO
//=============================================================================
//#endif // COMISO_NASOQ_AVAILABLE
//=============================================================================
