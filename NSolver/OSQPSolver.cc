//=============================================================================
//
//  CLASS OSQPSolver - IMPLEMENTATION
//
//=============================================================================

//== INCLUDES =================================================================

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
//#if COMISO_OSQP_AVAILABLE // TODO
//=============================================================================
#include "OSQPSolver.hh"

#include <osqp.h>

#include <CoMISo/Utils/CoMISoError.hh>
#include <CoMISo/Utils/StopWatch.hh>
#include <CoMISo/NSolver/LazyConstraintSolver.hh>

#include <Base/Debug/DebTime.hh>

#include <Eigen/Sparse>

//== NAMESPACES ===============================================================

namespace COMISO {

//== IMPLEMENTATION ========================================================== 


// helper struct to create and destruct OSQP objects
struct OSQPSolver::OSQPStructures
{
  OSQPStructures()
    :
      workspace(nullptr), // will be created by OSQP
      settings(new OSQPSettings()),
      data(new OSQPData)
  {
    osqp_set_default_settings(settings);
    data->n = 0;
    data->m = 0;
    data->P = nullptr;
    data->A = nullptr;
    data->q = nullptr;
    data->l = nullptr;
    data->u = nullptr;
  }

  ~OSQPStructures()
  {
    // Cleanup
    if (workspace != nullptr)
      delete  workspace;
    if (data != nullptr)
    {
      if (data->A)
        delete data->A;
      if (data->P)
        delete data->P;
      delete data;
    }
    if (settings != nullptr)
      delete settings;
  }

  OSQPWorkspace* workspace;
  OSQPSettings*  settings;
  OSQPData*      data;
};


void OSQPSolver::regularize_hessian(NProblemInterface::SMatrixNP& _H)
{
  NProblemInterface::SMatrixNP id;
  id.resize(_H.rows(), _H.cols());
  id.setIdentity();
  double reg_factor = 1e-8;
  auto diag = _H.diagonal();
  _H = _H + reg_factor * diag.sum()/diag.rows() * id;
}

NProblemInterface::SMatrixNP OSQPSolver::get_hessian(NProblemInterface* _problem)
{
  std::vector<double> zero(_problem->n_unknowns(), 0);
  NProblemInterface::SMatrixNP H;
  _problem->eval_hessian(zero.data(), H);
  regularize_hessian(H);
  return H;
}

Eigen::VectorXd OSQPSolver::get_linear_energy_coefficients(NProblemInterface* _problem)
{
  std::vector<double> zero(_problem->n_unknowns(), 0);
  Eigen::VectorXd q;
  q.resize(_problem->n_unknowns());
  _problem->eval_gradient(zero.data(), q.data());
  return q;
}

void OSQPSolver::get_constraints(int _n_cols, const std::vector<NConstraintInterface*>& _constraints, COMISO::NProblemInterface::SMatrixNP& _C, Eigen::VectorXd& _lower_bounds, Eigen::VectorXd& _upper_bounds)
{
  size_t n_rows = _constraints.size();
  _C.resize(n_rows, _n_cols);
  _lower_bounds.resize(n_rows);
  _upper_bounds.resize(n_rows);

  std::vector<double> x(_n_cols, 0.0);
  std::vector<Eigen::Triplet<double>> triplets;
  int current_row = 0;
  for (auto c : _constraints)
    if (c->is_linear())
    {
      NConstraintInterface::SVectorNC gc;
      c->eval_gradient(x.data(), gc);
      for(NConstraintInterface::SVectorNC::InnerIterator v_it(gc); v_it; ++v_it)
        triplets.emplace_back(current_row, v_it.index(), v_it.value());

      double b = c->eval_constraint(x.data());

      if (c->constraint_type() == NConstraintInterface::NC_EQUAL)
      {
        _lower_bounds[current_row] = -b;
        _upper_bounds[current_row] = -b;
      }
      else if (c->constraint_type() == NConstraintInterface::NC_LESS_EQUAL)
      {
        _lower_bounds[current_row] = -std::numeric_limits<double>::max();
        _upper_bounds[current_row] = -b;
      }
      else if (c->constraint_type() == NConstraintInterface::NC_GREATER_EQUAL)
      {
        _lower_bounds[current_row] = -b;
        _upper_bounds[current_row] = std::numeric_limits<double>::max();
      }

      ++current_row;
    }
    else
    {
      DEB_error("OSQP received non linear constraints which is not supported and thus ignored.");
    }

  _C.setFromTriplets(triplets.begin(), triplets.end());
}



bool OSQPSolver::solve(NProblemInterface* _problem, const std::vector<NConstraintInterface*>& _constraints)
{
  OSQPStructures osqp_structures;
  return solve(_problem, _constraints, osqp_structures);
}


bool OSQPSolver::solve(NProblemInterface* _problem, const std::vector<NConstraintInterface*>& _constraints, const std::vector<NConstraintInterface*>& _lazy_constraints, double _acceptable_tolerance, double _almost_infeasible_threshold)
{

  OSQPStructures osqp_structures;
  auto solve_function = [this, &osqp_structures](NProblemInterface* _problem, const std::vector<NConstraintInterface*> _constraints)
  {
    return solve(_problem, _constraints, osqp_structures);
  };

  auto get_res_function = [&osqp_structures]()
  {
    return osqp_structures.workspace->solution->x;
  };

  return solve_with_lazy_constraints(solve_function, get_res_function, _problem, _constraints, _lazy_constraints, _acceptable_tolerance, _almost_infeasible_threshold);
}

bool OSQPSolver::solve(NProblemInterface* _problem, const std::vector<NConstraintInterface*>& _constraints, OSQPSolver::OSQPStructures& _osqp_structures)
{
  // Load problem data
//  c_float P_x[3] = {4.0, 1.0, 2.0, };       // the upper triangular part of the quadratic cost matrix P in csc format (size n x n).
//  c_int P_nnz = 3;                          // number of non zeros
//  c_int P_i[3] = {0, 0, 1, };               // row indices
//  c_int P_p[3] = {0, 1, 3, };               // column pointers
//  c_float q[2] = {1.0, 1.0, };              // dense array for linear part of cost function (size n)
//  c_float A_x[4] = {1.0, 1.0, 1.0, 1.0, };  // linear constraints matrix A in csc format (size m x n)
//  c_int A_nnz = 4;                          // number of non zeros
//  c_int A_i[4] = {0, 1, 0, 2, };            // number of non z
//  c_int A_p[3] = {0, 2, 4, };               // row indices
//  c_float l[3] = {1.0, 0.0, 0.0, };         // dense array for lower bound (size m)
//  c_float u[3] = {1.0, 0.7, 0.7, };         // dense array for upper bound (size m)
//  c_int n = 2;                              // number of variables n
//  c_int m = 3;                              // number of constraints m

  auto H = get_hessian(_problem);
  auto lin_q = get_linear_energy_coefficients(_problem);

  COMISO::NProblemInterface::SMatrixNP A; // inequality constraints
  Eigen::VectorXd lower;                  // lower bounds
  Eigen::VectorXd upper;                  // upper bounds
  get_constraints(_problem->n_unknowns(), _constraints, A, lower, upper);

  COMISO::NProblemInterface::SMatrixNP HupperTriangle = H.triangularView<Eigen::Upper>();
  HupperTriangle.makeCompressed();


  c_float* P_x   = HupperTriangle.valuePtr();                    // the upper triangular part of the quadratic cost matrix P in csc format (size n x n).
  c_int    P_nnz = static_cast<int>(HupperTriangle.nonZeros());  // number of non zeros
  c_int*   P_i   = HupperTriangle.innerIndexPtr();               // row indices
  c_int*   P_p   = HupperTriangle.outerIndexPtr();               // column pointers
  c_float* q     = lin_q.data();                                 // dense array for linear part of cost function (size n)
  c_float* A_x   = A.valuePtr();                                 // linear constraints matrix A in csc format (size m x n)
  c_int    A_nnz = static_cast<int>(A.nonZeros());               // number of non zeros
  c_int*   A_i   = A.innerIndexPtr();                            // number of non z
  c_int*   A_p   = A.outerIndexPtr();                            // row indices
  c_float* l     = lower.data();                                 // dense array for lower bound (size m)
  c_float* u     = upper.data();                                 // dense array for upper bound (size m)
  c_int    n     = static_cast<int>(HupperTriangle.cols());      // number of variables n
  c_int    m     = static_cast<int>(A.rows());                   // number of constraints m

  // Exitflag
  c_int exitflag = 0;

  // Workspace structures
  OSQPWorkspace *& work     = _osqp_structures.workspace;
  OSQPSettings  *& settings = _osqp_structures.settings;
  OSQPData      *& data     = _osqp_structures.data;

  // Populate data
  if (data) {
    data->n = n;
    data->m = m;
    if (data->P != nullptr)
      delete data->P;
    data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
    data->q = q;
    if (data->A != nullptr)
      delete data->A;
    data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
    data->l = l;
    data->u = u;
  }

  // Define solver settings
  if (settings) {
    osqp_set_default_settings(settings);
    settings->alpha = 1.0; // Change alpha parameter
    settings->max_iter = 20000;
    settings->warm_start = true;
  }

  // Setup workspace
  exitflag = osqp_setup(&work, data, settings);

  if (exitflag != 0)
  {
    DEB_error("OSQP Setup failed with exit flag " << exitflag);
    return false;
  }

  // Solve Problem
  exitflag = osqp_solve(work);

  _problem->store_result(work->solution->x);

  if (exitflag != 0)
  {
    DEB_error("OSQP failed solving with exit flag " << exitflag);
    return false;
  }

  return exitflag == 0;
}



//-----------------------------------------------------------------------------

//=============================================================================
} // namespace COMISO
//=============================================================================
//#endif // COMISO_OSQP_AVAILABLE
//=============================================================================
