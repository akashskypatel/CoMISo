//=============================================================================
//
//  CLASS NewtonSolver - IMPLEMENTATION
//
//=============================================================================

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>

//== INCLUDES =================================================================

#include "NewtonSolver.hh"
#include <CoMISo/Solver/CholmodSolver.hh>
#include <CoMISo/NSolver/ConstraintTools.hh>
#include <Base/Debug/DebTime.hh>
//== NAMESPACES ===============================================================

namespace COMISO {

//== IMPLEMENTATION ========================================================== 

#if COMISO_GMM_AVAILABLE
// solve
int
NewtonSolver::
solve(NProblemGmmInterface* _problem)
{
  DEB_enter_func;
#if COMISO_SUITESPARSE_CHOLMOD_AVAILABLE
  converged_ = true;

  // get problem size
  int n = _problem->n_unknowns();

  // hesse matrix
  NProblemGmmInterface::SMatrixNP H;
  // gradient
  std::vector<double> x(n), x_new(n), dx(n), g(n);

  // get initial x, initial grad and initial f
  _problem->initial_x(P(x));
  double f = _problem->eval_f(P(x));

  double reg = 1e-3;
  COMISO::CholmodSolver chol;

  for(int i=0; i<max_iters_; ++i)
  {
    _problem->eval_gradient(P(x), P(g));
    // check for convergence
    if( gmm::vect_norm2(g) < eps_)
    {
      DEB_line(2, "Newton Solver converged after " << i << " iterations");
      _problem->store_result(P(x));
      return true;
    }

    // get current hessian
    _problem->eval_hessian(P(x), H);

    // regularize
    double reg_comp = reg*gmm::mat_trace(H)/double(n);
    for(int j=0; j<n; ++j)
      H(j,j) += reg_comp;

    // solve linear system
    bool factorization_ok = false;
    if(constant_hessian_structure_ && i != 0)
      factorization_ok = chol.update_system_gmm(H);
    else
      factorization_ok = chol.calc_system_gmm(H);

    bool improvement = false;
    if(factorization_ok)
      if(chol.solve( dx, g))
      {
        gmm::add(x, gmm::scaled(dx,-1.0),x_new);
        double f_new = _problem->eval_f(P(x_new));

        if( f_new < f)
        {
          // swap x and x_new (and f and f_new)
          x_new.swap(x);
          f = f_new;
          improvement = true;

          DEB_line(6, "energy improved to " << f);
        }
      }

    // adapt regularization
    if(improvement)
    {
      if(reg > 1e-9)
        reg *= 0.1;
    }
    else
    {
      if(reg < 1e4)
        reg *= 10.0;
      else
      {
        _problem->store_result(P(x));
        DEB_warning(2, "Newton solver reached max regularization but did not "
          "converge");
        converged_ = false;
        return false;
      }
    }
  }
  _problem->store_result(P(x));
  DEB_warning(2, "Newton Solver did not converge!!! after iterations.");
  converged_ = false;
  return false;

#else
  DEB_warning(1,"NewtonSolver requires not-available CholmodSolver");
  converged_ = false;
  return false;
#endif	    
}
#endif // COMISO_GMM_AVAILABLE

//-----------------------------------------------------------------------------


int NewtonSolver::solve(NProblemInterface* _problem, const SMatrixD& _A, 
  const VectorD& _b)
{
  DEB_time_func_def;
  converged_ = false;

  double KKT_res_eps = 1e-6;
  int    max_KKT_regularization_iters = 40;
  double regularize_constraints_limit = 1e-6;
  double max_allowed_constraint_violation2 = 1e-12;

  // require lower accuracy for iterative solvers
  if(solver_type_ == LS_EigenCG || solver_type_ == LS_EigenBiCGSTAB)
  {
    cg_solver_.setTolerance(1e-8);
    bicgstab_solver_.setTolerance(1e-8);

    KKT_res_eps = 1e-3;
    max_allowed_constraint_violation2 = 1e-6;
  }

  // number of unknowns
  size_t n = _problem->n_unknowns();
  // number of constraints
  size_t m = _b.size();

  DEB_line(2, "optimize via Newton with " << n << " unknowns and " << m <<
    " linear constraints");
  DEB_line(2,"using linear solver " << name_of(solver_type_));

  // initialize vectors of unknowns
  VectorD x(n);
  _problem->initial_x(x.data());

  double initial_constraint_violation2 = (_A*x-_b).squaredNorm();

  // storage of update vector dx and rhs of KKT system
  VectorD dx(n+m), rhs(n+m), g(n);
  rhs.setZero();

  // resize temp vector for line search (and set to x1 to approx Hessian correctly if problem is non-quadratic!)
  x_ls_ = x;

  // indicate that system matrix is symmetric
  lu_solver_.isSymmetric(true);

  // start with no regularization
  double regularize_hessian(0.0);
  double regularize_constraints(0.0);
  int iter=0;
  bool first_factorization = true;
  while( iter < max_iters_)
  {
    double kkt_res2(0.0);
    double constraint_res2(0.0);
    int    reg_iters(0);
    bool fact_ok = true;
    do
    {
      // get Newton search direction by solving LSE
      fact_ok = factorize(_problem, _A, _b, x, regularize_hessian, regularize_constraints, first_factorization);
      first_factorization = false;

      if(fact_ok)
      {
        // get rhs
        _problem->eval_gradient(x.data(), g.data());
        rhs.head(n) = -g;
        rhs.tail(m) = _b - _A*x;

        // solve KKT system
        solve_kkt_system(rhs, dx);

        // check numerical stability of KKT system and regularize if necessary
        kkt_res2 = (KKT_*dx-rhs).squaredNorm();
        constraint_res2 = (_A*dx.head(n)-rhs.tail(m)).squaredNorm();
      }

      if(!fact_ok || kkt_res2 > KKT_res_eps || constraint_res2 > max_allowed_constraint_violation2)
      {
        DEB_warning(2, "Numerical issues in KKT system");
        DEB_warning_if(!fact_ok, 2, "Factorization not ok");
        DEB_line_if(
            kkt_res2 > KKT_res_eps, 3, "KKT Residuum too high: " << kkt_res2);
        DEB_line_if(constraint_res2 > max_allowed_constraint_violation2, 3,
            "Constraint violation too high: " << constraint_res2);
        // alternate hessian and constraints regularization
        if(reg_iters % 2 == 0 || regularize_constraints >= regularize_constraints_limit)
        {
          DEB_line(2, "residual ^ 2 " << kkt_res2 << "->regularize hessian");
          if(regularize_hessian == 0.0)
            regularize_hessian = 1e-6;
          else
            regularize_hessian *= 2.0;
        }
        else
        {
          DEB_line(2, "residual^2 " << kkt_res2 << " -> regularize constraints");
          if(regularize_constraints == 0.0)
            regularize_constraints = 1e-8;
          else
            regularize_constraints *= 2.0;
        }
      }
      ++reg_iters;
    }
    while( (!fact_ok || kkt_res2 > KKT_res_eps || constraint_res2 > max_allowed_constraint_violation2) && reg_iters < max_KKT_regularization_iters);

    // no valid step could be found?
    if(kkt_res2 > KKT_res_eps || constraint_res2 > max_allowed_constraint_violation2 || reg_iters >= max_KKT_regularization_iters)
    {
      DEB_error("numerical issues in KKT system could not be resolved "
        "-> terminating NewtonSolver with current solution");
      _problem->store_result(x.data());
      return 0;
    }

    // get maximal reasonable step
    double t_max  = std::min(1.0, 
			     max_feasible_step_safety_factor_ * _problem->max_feasible_step(x.data(), dx.data()));

    // perform line-search
    double newton_decrement(0.0);
    double fx(0.0);
    double t = backtracking_line_search(_problem, x, g, dx, newton_decrement, fx, t_max);

    // perform update
    x += dx.head(n)*t;

    double constraint_violation2 = (_A*x-_b).squaredNorm();

    if(constraint_violation2 > 2*initial_constraint_violation2 && constraint_violation2 > max_allowed_constraint_violation2)
    {
      DEB_warning(2, "Numerical issues in KKT system lead to "
        "constraint violation -> recovery phase");
      // restore old solution
      x -= dx.head(n)*t;

      regularize_constraints *= 0.5;
      regularize_constraints_limit = regularize_constraints;
    }

    DEB_line(4,
        "iter: " << iter << ", f(x) = " << fx << ", t = " << t
                 << " (tmax=" << t_max << ")" << (t < t_max ? " _clamped_" : "")
                 << ", eps = [Newton decrement] = " << newton_decrement
                 << ", constraint violation prior = " << rhs.tail(m).norm()
                 << ", after = " << (_b - _A * x).norm()
                 << ", KKT residual^2 = " << kkt_res2);

    // converged?
    if (newton_decrement < eps_ || std::abs(t) < eps_ls_)
    {
      converged_ = true;
      break;
    }

    ++iter;
  }

  // store result
  _problem->store_result(x.data());

  // return success
  return 1;
}


//-----------------------------------------------------------------------------


int NewtonSolver::solve_infeasible_start(NProblemInterface* _problem, const SMatrixD& _A,
                        const VectorD& _b)
{
  DEB_time_func_def;
  converged_ = false;
  mu_merit_  = 0.0; // reset

  int    max_regularization_iters          = 10;
  double max_regularization_value          = 1e6;//1e6
  double max_allowed_constraint_violation2 = 1e-12;
  double KKT_res_eps                       = 1e-6;

  // require lower accuracy for iterative solvers
  if(solver_type_ == LS_EigenCG || solver_type_ == LS_EigenBiCGSTAB)
  {
    cg_solver_.setTolerance(1e-8);
    bicgstab_solver_.setTolerance(1e-8);
    max_allowed_constraint_violation2 = 1e-6;
  }

  // copy constraints (since they might be made linear independent if necessary)
  SMatrixD A = _A;
  VectorD  b = _b;

  // number of unknowns
  size_t n = _problem->n_unknowns();
  // number of constraints
  size_t m = b.size();

  // initialize vector of unknowns
  VectorD x(n);
  _problem->initial_x(x.data());
  // initialize vector of dual variables
  VectorD nue(m);
  nue.setZero();

  // resize temp vector for line search (and set to x0 to approx Hessian correctly if problem is non-quadratic!)
  x_ls_ = x;
  nue_ls_ = nue;

  // storage of update vector dz=(dx,dnue) and rhs of KKT system
  VectorD dz(n+m), rhs(n+m), g(n);
  rhs.setZero();
  g_ls_ = g;

  // update gradient
  _problem->eval_gradient(x.data(), g.data());

  double res_primal2 = (A*x-b).squaredNorm();
  double res_dual2   = (g + A.transpose()*nue).squaredNorm();

  DEB_line(2, "***** optimize via Newton (infeasible start version) with " << n << " unknowns and " << m <<
                                                                           " linear constraints (initial residuals r_primal = ||Ax-b||^2 = " << res_primal2 << ", r_dual = ||g+A^T nue||^2 = " << res_dual2 << ")" );
  DEB_line(2,"using linear solver " << name_of(solver_type_));
  // indicate that system matrix is symmetric
  lu_solver_.isSymmetric(true);

  // start with no regularization
  double regularize_constraints = prescribed_constraint_regularization_absolute_;
  double regularize_hessian     = prescribed_hessian_regularization_relative_;
  int reg_iters=0;
  iters_ = 0;
  double kkt_res2=0.0;
  bool first_factorization = true;
  bool constraints_checked = false;
  while( iters_ < max_iters_)
  {
    // step length (default -1 indicates that line-search was not performed)
    double t = -1.0;

    // set to infinity
    kkt_res2 = std::numeric_limits<double>::infinity();

    // get hessian
    SMatrixD H(n, n);
    _problem->eval_hessian(x.data(), H);

    // get rhs
    rhs.head(n) = -g - A.transpose() * nue;
    rhs.tail(m) = b - A * x;

    // get Newton search direction by solving LSE
    bool fact_ok = factorize(H, rhs.head(n), A, rhs.tail(m), regularize_hessian, regularize_constraints, first_factorization);
    first_factorization = false;

    // factorization succeeded?
    if (fact_ok)
    {
      // solve KKT system
      solve_kkt_system(rhs, dz);
      // check numerical stability of KKT system and regularize if necessary
      kkt_res2 = (KKT_ * dz - rhs).squaredNorm();

      if(kkt_res2 < KKT_res_eps)
      {
        // get maximal reasonable step
        double t_max = std::min(1.0,
                                max_feasible_step_safety_factor_ * _problem->max_feasible_step(x.data(), dz.data()));
        t = t_max;

        // decide on line search strategy:
        // ||Ax-b|| too large --> line search on residual
        // else               --> feasible Newton line search
        bool feasible_mode = (res_primal2 < max_allowed_constraint_violation2);
        double newton_decrement(DBL_MAX);
        double fx(0.0);
        if (feasible_mode)
          t = backtracking_line_search(_problem, x, g, dz, newton_decrement, fx, t_max);
        else
          t = backtracking_line_search_infeasible_merit_l1(_problem, H, A, b, x, g, dz, fx, t_max,
                                                  100);

        // old infeasible line search works well only for convex problems
        //        t = backtracking_line_search_infeasible(_problem, A, b, x, nue, g, dz, fx, res_primal2, res_dual2, t_max,
        //                                                100);

        // perform update
//        if (feasible_mode) // infeasible mode updates already during line search
        {
          // update primal variables
          x += dz.head(n) * t;
          // update dual variables
          nue += dz.tail(m) * t;
          // update gradient
          _problem->eval_gradient(x.data(), g.data());
          // update residuals
          res_primal2 = (A * x - b).squaredNorm();
          res_dual2 = (g + A.transpose() * nue).squaredNorm();
        }

        if (feasible_mode)
        {
          DEB_line(4,
                   "iter: " << iters_ << ", feasible, f(x) = " << fx << ", t = " << t
                            << " (tmax=" << t_max << ")" << (t < t_max ? " _clamped_" : "")
                            << ", res_primal = " << res_primal2
                            << ", res_dual   = " << res_dual2
                            << ", regularization = " << regularize_hessian
                            << ", eps = [Newton decrement] = " << newton_decrement
                            << ", KKT residual^2 = " << kkt_res2);
        }
        else
          {
          DEB_line(4,
                   "iter: " << iters_ << ", infeasible, f(x) = " << fx << ", t = " << t
                            << " (tmax=" << t_max << ")" << (t < t_max ? " _clamped_" : "")
                            << ", res_primal = " << res_primal2
                            << ", res_dual   = " << res_dual2
                            << ", regularization = " << regularize_hessian
                            << ", KKT residual^2 = " << kkt_res2);

          }

        // converged?
        if (feasible_mode && newton_decrement < eps_) {
          converged_ = true;
          break;
        }
//        else if (t <= eps_ls_) // line search was not able to determine a valid step
//        {
//          DEB_warning(1, "NEWTON TERMINATE because line search could not find a valid nonzero step");
//          converged_ = false;
//          break;
//        }

        if(t > eps_ls_ && use_trust_region_regularization_)
        {
          // adapt hessian regularization if necessary (similar to trust-region method)
          if (t < 0.1) // increase regularization if step length small
            regularize_hessian = std::min(std::max(1e-3, 10.0 * regularize_hessian), max_regularization_value);
          if (t == 1.0) // decrease regularization after full step
          {
            regularize_hessian = std::max(prescribed_hessian_regularization_relative_, 0.5 * regularize_hessian);
            if(regularize_hessian < 1e-3) // set to lowest possible value when below 1e-3
              regularize_hessian = prescribed_hessian_regularization_relative_;
          }
        }
        ++iters_;
      }
    }

    // solution of KKT system was not successfull, or step length vanishing? --> try to increase hessian regularization
    if(!fact_ok || kkt_res2 > KKT_res_eps || t<=eps_ls_)
    {
      // check for linear dependency of constraints
      bool constraint_update_successfull = false;
      if (!constraints_checked)
      {
        // make sure that constraints are linearly independent
        COMISO::ConstraintTools::remove_dependent_linear_constraints(A, b);

        // if number of constraitns m changed, adapt everything related to m
        if ((size_t)b.size() < m)
        {
          DEB_line(2, "--> constraints are linearly dependent reduce from " << m << " to " << (int)b.size());

          // adapt size
          m = b.size();
          nue.resize(m);
          nue.setZero();
          nue_ls_ = nue;

          dz.resize(n + m);
          rhs.resize(n + m);

          // try again to factorize
          constraint_update_successfull = factorize(H, rhs.head(n), A, rhs.tail(m), regularize_hessian, regularize_constraints, true);
        }
        constraints_checked = true;
      }

      // still not ok?
      if(!constraint_update_successfull)
      {
        regularize_hessian = std::max(1e-3, 10.0 * regularize_hessian);
        DEB_line(2, "--> valid step could not be found: factorization_ok =  " << int(fact_ok) << ", KKT residual^2 = " << kkt_res2 << ", t = " << t << ", ||grad||^2 = " <<  g.squaredNorm());
        DEB_line(2, "--> regularize Hessian with value " << regularize_hessian);

        ++reg_iters;

        if (reg_iters >
            max_regularization_iters || regularize_hessian > max_regularization_value) // regularization failed after max allowed number of regularization steps
        {
          DEB_warning(1, "NEWTON TERMINATE because of reaching max of allowed regularization");
          converged_ = false;
          break;
        }
      }
    }
  }

  // store result
  _problem->store_result(x.data());

  // return success
  return converged_;
}


//-----------------------------------------------------------------------------


bool NewtonSolver::factorize(NProblemInterface* _problem,
  const SMatrixD& _A, const VectorD& _b, const VectorD& _x, double& _regularize_hessian_relative, double& _regularize_constraints_absolute,
  const bool _first_factorization)
{
  DEB_enter_func;

  const int n  = _problem->n_unknowns();
  const int m  = static_cast<int>(_A.rows());
  const int nf = n+m;

  // get hessian of quadratic problem
  SMatrixD H(n,n);
  _problem->eval_hessian(_x.data(), H);

  // set up KKT matrix
  // create sparse matrix
  trips_.clear();
  trips_.reserve(H.nonZeros() + 2*_A.nonZeros());

  // add elements of H
  for (int k=0; k<H.outerSize(); ++k)
    for (SMatrixD::InnerIterator it(H,k); it; ++it)
      trips_.push_back(Triplet(static_cast<int>(it.row()),static_cast<int>(it.col()),it.value()));

  // add elements of _A
  for (int k=0; k<_A.outerSize(); ++k)
    for (SMatrixD::InnerIterator it(_A,k); it; ++it)
    {
      // insert _A block below
      trips_.push_back(Triplet(static_cast<int>(it.row())+n,static_cast<int>(it.col()),it.value()));

      // insert _A^T block right
      trips_.push_back(Triplet(static_cast<int>(it.col()),static_cast<int>(it.row())+n,it.value()));
    }

  // regularize constraints
//  if(_regularize_constraints != 0.0)
  for( int i=0; i<m; ++i)
    trips_.push_back(Triplet(n+i,n+i,_regularize_constraints_absolute));

  // regularize Hessian
//  if(_regularize_hessian != 0.0)
  {
    double ad(0.0);
    for( int i=0; i<n; ++i)
      ad += H.coeffRef(i,i);
    ad *= _regularize_hessian_relative/double(n);
    for( int i=0; i<n; ++i)
      trips_.push_back(Triplet(i,i,ad));
  }

  // create KKT matrix
  KKT_.resize(nf,nf);
  KKT_.setFromTriplets(trips_.begin(), trips_.end());

  // compute LU factorization
  if(_first_factorization)
    analyze_pattern(KKT_);

  return numerical_factorization(KKT_);
}


//-----------------------------------------------------------------------------


bool NewtonSolver::factorize(const SMatrixD& _H, const VectorD& _g, const SMatrixD& _A, const VectorD& _b,
                             const double _regularize_hessian, const double _regularize_constraints, const bool _update_pattern)
{
  DEB_enter_func;

  const int n  = static_cast<int>(_H.cols());
  const int m  = static_cast<int>(_A.rows());
  const int nf = n+m;

  // set up KKT matrix
  // create sparse matrix
  trips_.clear();
  trips_.reserve(_H.nonZeros() + 2*_A.nonZeros());

  // add elements of H
  for (int k=0; k<_H.outerSize(); ++k)
    for (SMatrixD::InnerIterator it(_H,k); it; ++it)
      trips_.push_back(Triplet(static_cast<int>(it.row()),static_cast<int>(it.col()),it.value()));

  // add elements of _A
  for (int k=0; k<_A.outerSize(); ++k)
    for (SMatrixD::InnerIterator it(_A,k); it; ++it)
    {
      // insert _A block below
      trips_.push_back(Triplet(static_cast<int>(it.row())+n,static_cast<int>(it.col()),it.value()));

      // insert _A^T block right
      trips_.push_back(Triplet(static_cast<int>(it.col()),static_cast<int>(it.row())+n,it.value()));
    }

  // regularize constraints (add also if zero to obtain constant pattern)
  for( int i=0; i<m; ++i)
    trips_.push_back(Triplet(n+i,n+i,_regularize_constraints));

  // regularize Hessian (add also if zero to obtain constant pattern)
  for( int i=0; i<n; ++i)
    trips_.push_back(Triplet(i,i,_regularize_hessian));

  // create KKT matrix
  KKT_.resize(nf,nf);
  KKT_.setFromTriplets(trips_.begin(), trips_.end());

#if COMISO_OSQP_AVAILABLE
  if(solver_type_ == LS_OSQP_DIRECT || solver_type_ == LS_OSQP_CG ||  solver_type_ == LS_OSQP_MKL)
  {
    if(_update_pattern)
      return setup_osqp(_H, -_g, _A, _b);
    else
      return update_osqp(_H, -_g, _b);
  }
#else
  if(solver_type_ == LS_OSQP_DIRECT || solver_type_ == LS_OSQP_CG ||  solver_type_ == LS_OSQP_MKL)
  {
    DEB_warning(1, "OSQP solver selected but not available!");
  }
#endif



  // compute LU factorization
  if(_update_pattern)
    analyze_pattern(KKT_);

  return numerical_factorization(KKT_);
}


//-----------------------------------------------------------------------------


double NewtonSolver::backtracking_line_search(NProblemInterface* _problem, 
  VectorD& _x, VectorD& _g, VectorD& _dx, double& _newton_decrement, 
  double& _fx, const double _t_start)
{
  DEB_enter_func;
  size_t n = _x.size();

  // pre-compute objective
  double fx = _problem->eval_f(_x.data());

  // pre-compute dot product
  double gtdx = _g.transpose()*_dx.head(n);
  _newton_decrement = std::abs(gtdx);

  // current step size
  double t = _t_start;

  // backtracking (stable in case of NAN and with max 100 iterations)
  for(int i=0; i<100; ++i)
  {
    // current update
    x_ls_ = _x + _dx.head(n)*t;
    double fx_ls = _problem->eval_f(x_ls_.data());

    if( fx_ls <= fx + alpha_ls_*t*gtdx )
    {
      _fx = fx_ls;
      return t;
    }
    else
      t *= beta_ls_;
  }

  DEB_warning(1, "line search could not find a valid step within 100 "
    "iterations");
  _fx = fx;
  return 0.0;
}


//-----------------------------------------------------------------------------


double NewtonSolver::backtracking_line_search_infeasible(NProblemInterface* _problem, const SMatrixD& _A, const VectorD& _b,
                                              VectorD& _x, VectorD& _nue, VectorD& _g, VectorD& _dz, double& _fx, double& _res_primal2, double& _res_dual2,
                                              const double _t_start, const int _max_ls_iters)
{
  DEB_enter_func;
  size_t n = _x.size();
  size_t m = _b.size();

  // current step size
  double t = _t_start;
  // current residual
  double r0  = _res_primal2 + _res_dual2;
  double rp2 = _res_primal2;
  double rd2 = _res_dual2;

  double fx(0.0);

  // backtracking (stable in case of NAN and with max iterations)
  for(int i=0; i<_max_ls_iters; ++i)
  {
    // current update of x, nue and g
    x_ls_   = _x + _dz.head(n)*t;
    fx = _problem->eval_f(x_ls_.data());
    // check if update is inside domain
    if(std::isfinite(fx))
    {
      nue_ls_ = _nue + _dz.tail(m) * t;
      _problem->eval_gradient(x_ls_.data(), g_ls_.data());

      // get current residual
      rp2 = (_A * x_ls_ - _b).squaredNorm();
      rd2 = (g_ls_ + _A.transpose() * nue_ls_).squaredNorm();
      double r = rp2 + rd2;

 //     std::cerr << " line search r=" << r << " vs r0=" << r0 << " r0-r=" << r0-r << std::endl;

      // sufficient decrease in residual?
      if (r < std::pow(1.0 - alpha_ls_*t,2) * r0)
      {
        // succesfull line search
        _res_primal2 = rp2;
        _res_dual2 = rd2;
        _fx = fx;
        _g.swap(g_ls_);
        _x.swap(x_ls_);
        _nue.swap(nue_ls_);
        return t;
      }
    }

    // shrink with factor beta
    t *= beta_ls_;
  }

  DEB_warning(1, "line search could not find a valid step");
  _fx = _problem->eval_f(_x.data());
  _problem->eval_gradient(_x.data(), _g.data());
  return 0.0;
}


//-----------------------------------------------------------------------------


double NewtonSolver::backtracking_line_search_infeasible_merit_l1(NProblemInterface* _problem, const SMatrixD& _H,
                                                                  const SMatrixD& _A, const VectorD& _b,
                                                                  const VectorD& _x, const VectorD& _g, VectorD& _dz, double& _fx,
                                                                  const double _t_start, const int _max_ls_iters)
{
  DEB_enter_func;
  size_t n = _x.size();

  // update mu
  double res_primal_1 = (_A*_x-_b).template lpNorm<1>();
  double gdx          = _g.transpose()*_dz.head(n);
  double dxHdx        = _dz.head(n).transpose()*_H*_dz.head(n);
  double mu_new = 1.2*(gdx+0.5*std::max(0.0,dxHdx))/((1.0-rho_merit_)*res_primal_1);
  mu_merit_ = std::max(mu_new, mu_merit_);

  // current step size
  double t = _t_start;

  // merit function and directional derivative for t=0
//  double merit_0  = _fx + mu_merit_*res_primal_1;
  double merit_0  = _problem->eval_f(_x.data()) + mu_merit_*res_primal_1;
  double D_merit_0 = gdx - mu_merit_*res_primal_1;

//  std::cerr << "mu_merit=" << mu_merit_ << std::endl;
//  std::cerr << "D_merit_0=" << D_merit_0 << std::endl;

  double fx(0.0);

  // backtracking (stable in case of NAN and with max iterations)
  for(int i=0; i<_max_ls_iters; ++i)
  {
    // current update of x, nue and g
    x_ls_   = _x + _dz.head(n)*t;
    fx = _problem->eval_f(x_ls_.data());
    // check if update is inside domain
    if(std::isfinite(fx))
    {
      double merit_t = fx + mu_merit_*(_A*x_ls_-_b).template lpNorm<1>();

 //     std::cerr << "t=" << t << ", merit(t)=" << merit_t << ", merit(0)+alpha*t*Dmerit(0)=" << merit_0 + alpha_ls_*t*D_merit_0 << std::endl;

      // sufficient decrease in residual?
      if (merit_t <= merit_0 + alpha_ls_*t*D_merit_0)
      {
        // succesfull line search
        _fx = fx;
        return t;
      }
    }
 //   else std::cerr << "t=" << t << ", merit(t)=" << std::numeric_limits<double>::infinity() << std::endl;

    // shrink with factor beta
    t *= beta_ls_;
  }

  DEB_warning(1, "line search could not find a valid step");
  return 0.0;
}


//-----------------------------------------------------------------------------


void NewtonSolver::analyze_pattern(SMatrixD& _KKT)
{
  DEB_enter_func;
  switch(solver_type_)
  {
    case LS_EigenLU:      lu_solver_.analyzePattern(_KKT); break;
#if COMISO_SUITESPARSE_UMFPACK_AVAILABLE
    case LS_Umfpack: umfpack_solver_.analyzePattern(_KKT); break;
#endif
    case LS_SPQR:
    case LS_EigenCG:
    case LS_EigenBiCGSTAB: break; // nothing to do

    default: DEB_warning(1, "selected linear solver not availble");
  }
}


//-----------------------------------------------------------------------------


bool NewtonSolver::numerical_factorization(SMatrixD& _KKT)
{
  DEB_enter_func;
  switch(solver_type_)
  {
    case LS_EigenLU:
      lu_solver_.factorize(_KKT); 
      return (lu_solver_.info() == Eigen::Success);
#if COMISO_SUITESPARSE_UMFPACK_AVAILABLE
    case LS_Umfpack:
      umfpack_solver_.factorize(_KKT); 
      return (umfpack_solver_.info() == Eigen::Success);
#endif
#if COMISO_SUITESPARSE_SPQR_AVAILABLE &&  COMISO_SUITESPARSE_CHOLMOD_AVAILABLE && EIGEN_VERSION_AT_LEAST(3,4,90)
      case LS_SPQR:
//      spqr_solver_.factorize(_KKT);
      spqr_solver_.compute(_KKT);
      return (spqr_solver_.info() == Eigen::Success);
#endif
    case LS_EigenCG:
      cg_solver_.compute(_KKT);
      return (cg_solver_.info() == Eigen::Success);
    case LS_EigenBiCGSTAB:
      bicgstab_solver_.compute(_KKT);
      return (bicgstab_solver_.info() == Eigen::Success);
    default: 
      DEB_warning(1, "selected linear solver not availble!"); 
      return false;
  }
}


//-----------------------------------------------------------------------------


void NewtonSolver::solve_kkt_system(const VectorD& _rhs, VectorD& _dx)
{
  DEB_enter_func;
  switch(solver_type_)
  {
    case LS_EigenLU: _dx =      lu_solver_.solve(_rhs); break;
#if COMISO_SUITESPARSE_UMFPACK_AVAILABLE
    case LS_Umfpack: _dx = umfpack_solver_.solve(_rhs); break;
#endif
#if COMISO_SUITESPARSE_SPQR_AVAILABLE && COMISO_SUITESPARSE_CHOLMOD_AVAILABLE && EIGEN_VERSION_AT_LEAST(3,4,90)
    case LS_SPQR: _dx = spqr_solver_.solve(_rhs); break;
#endif
    case LS_EigenCG: _dx = cg_solver_.solve(_rhs); break;
    case LS_EigenBiCGSTAB: _dx = bicgstab_solver_.solve(_rhs); break;

    case LS_OSQP_DIRECT:
    case LS_OSQP_CG:
    case LS_OSQP_MKL:
#if COMISO_OSQP_AVAILABLE
      solve_kkt_system_osqp(_dx); break;
#else
      DEB_warning(1, "OSQP Solver selected but not available"); break;
#endif
    default: DEB_warning(1, "selected linear solver not available"); break;
  }
}


//-----------------------------------------------------------------------------


double NewtonSolver::determine_spd_regularization(const SMatrixD& _A, const double _r0, const double _rscale, const int _max_iters)
{
  DEB_enter_func;

  size_t n = _A.rows();

  if(_A.rows() != _A.cols())
  {
    DEB_warning(1, "determine_spd_regularization received non-square matrix: " << (int)_A.rows() << " times " << (int)_A.cols());
    return 0.0;
  }

  Eigen::SimplicialLLT<SMatrixD> chol;

  SMatrixD A = _A;

  // determine initial regularization
  double r0(0.0);
  double trace(0.0);
  for(size_t i=0; i<n; ++i)
    trace += A.coeffRef(i,i);
  trace /= double(n);
  if(trace < 0.0)       r0 = -2.0*trace;
  else if(trace == 0.0) r0 = _r0;
  else                  r0 = _r0*trace;

  int iter = 0;
  double r =0.0;

  while(1)
  {
    // factorize and break if successfull
    chol.compute(A);
    if(chol.info() == Eigen::Success)
    {
      break;
    }

    if(r == 0.0)
    {
      r = r0;
      A.diagonal().array() += r0;
    }
    else
      {
        r *= _rscale;
        A.diagonal().array() += _rscale-1.0;
      }

    ++iter;
    if(iter >= _max_iters)
    {
      DEB_warning(1, "determine_spd_regularization failed after " << iter << " iterations ");
      break;
    }
  }

  return r;
}

//-----------------------------------------------------------------------------

#if COMISO_OSQP_AVAILABLE

int NewtonSolver::setup_osqp(const SMatrixD& _H, const VectorD& _g, const SMatrixD& _A, const VectorD& _b)
{
#if COMISO_OSQP_NEW_API
    // we now have these options:
    //osqp_eigen_.settings().linsys_solver = OSQP_DIRECT_SOLVER;
    //osqp_eigen_.settings().linsys_solver = OSQP_INDIRECT_SOLVER;
#else
  if(solver_type_ == LS_OSQP_DIRECT)
    osqp_eigen_.settings().linsys_solver = QDLDL_SOLVER;
  if(solver_type_ == LS_OSQP_MKL)
    osqp_eigen_.settings().linsys_solver = MKL_PARDISO_SOLVER;
#endif

  return (osqp_eigen_.setup(_H, _g, _A, _b, _b) == 0);
}


//-----------------------------------------------------------------------------


int NewtonSolver::update_osqp(const SMatrixD& _H, const VectorD& _g, const VectorD& _b)
{
  auto upd_obj_exit   = osqp_eigen_.update_objective(_H, _g);
  auto upd_const_exit = osqp_eigen_.update_bounds(_b, _b);

  return ( upd_obj_exit == 0 && upd_const_exit == 0);
}


//-----------------------------------------------------------------------------


void NewtonSolver::solve_kkt_system_osqp(VectorD& _dx)
{
  osqp_eigen_.solve();
  osqp_eigen_.get_x( &(_dx[0]));
  osqp_eigen_.get_y( &(_dx[x_ls_.size()]));
}

#endif // COMISO_OSQP_AVAILABLE

//=============================================================================
} // namespace COMISO
//=============================================================================
