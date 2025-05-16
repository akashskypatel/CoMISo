/*===========================================================================*\
 *                                                                           *
 *                            TruncatedNewtonPCG                             *
 *      Copyright (C) 2024 by Computer Graphics Group, University of Bern    *
 *                           http://cgg.unibe.ch                             *
 *                                                                           *
 *      Author: David Bommes                                                 *
 *                                                                           *
\*===========================================================================*/


#include "TruncatedNewtonPCG.hh"

#include <ios>
#include <iomanip>
#include <cassert>
#include <sstream>

#include <CoMISo/Utils/StopWatch.hh>
#include <CoMISo/NSolver/LinearConstraintConverter.hh>

#include <Base/Debug/DebTime.hh>
#include <Eigen/IterativeLinearSolvers>


#if COMISO_SUITESPARSE_CHOLMOD_AVAILABLE
#include <Eigen/CholmodSupport>
#endif

#if COMISO_METIS_AVAILABLE
#include <Eigen/MetisSupport>
#endif


namespace COMISO
{

int
TruncatedNewtonPCG::
solve(NProblemInterface* _problem)
{
  DEB_enter_func;

  // reset status
  status_ = OptimizerStatus();
  status_.feasible = true; // there are no constraints

  // number of unknowns
  Eigen::Index n = _problem->n_unknowns();
  DEB_line_if(!config_.silent, 2, "optimize via TruncatedNewtonPCG with " << (int)n << " unknowns");

  // Newton parameters
  const int    max_iters  = config_.max_iters;
  const double newton_tol = config_.eps;
//    const double relative_improvement_thres = 1e-3;
  // CG parameters
//    const int max_pcg_iters = std::max(100, _problem->n_unknowns()/10);
  const int max_pcg_iters = config_.max_pcg_iters;

  // allow minimal steps of 1e-6
  const int max_iter_ls   = int ( std::log(1e-6)/std::log(config_.beta_ls));

  // initialize vectors of unknowns
  VectorD x(n);
  _problem->initial_x(x.data());

  // storage of update vector dx and gradient
  VectorD dx(n), g(n), xt(n);
  dx.setZero();

  // hessian matrix
  SMatrixD H(n,n);

  // Prconditioned Conjugate Gradient Solver
//    Eigen::ConjugateGradient<SMatrixD, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double> > pcg;
  Eigen::ConjugateGradient<SMatrixD, Eigen::Lower|Eigen::Upper> pcg;
  pcg.setMaxIterations(max_pcg_iters);
  pcg.setTolerance(config_.pcg_tolerance);

  // statistics
  int n_pcg_iters(0);

  bool ls_succeeded = true;

  // get function value at current point
  double fx = _problem->eval_f(x.data());
  double fx_old = fx;

  DEB_line_if(!config_.silent, 2, "intial objective value =" << fx);
  if(!std::isfinite(fx))
  {
    DEB_line_if(!config_.silent, 2, "ERROR: intial objective value is not finite ---> abort");
    return false;
  }

  // data for adaptive hessian update
  int count_hessian_skip = 0;
  bool hessian_updated = false;
  double alpha_line_search = 0.0;
  double rel_objective_decrease = 0.0;

  for(int iter=0; iter<max_iters; ++iter)
  {
    // get gradient and Hessian
    _problem->eval_gradient(x.data(), g.data());

    if(iter == 0 || // first iteration
       count_hessian_skip >= config_.hessian_max_skips || // maximum number of Hessians have been skipped
       alpha_line_search < config_.hessian_min_acceptable_alpha || // line search truncated too much
       rel_objective_decrease < config_.hessian_min_acceptable_rel_objective_decrease // objective did not decrease sufficiently
            ) // update Hessian?
    {
      _problem->eval_hessian(x.data(), H);
      hessian_updated = true;
      count_hessian_skip = 0;
    }
    else
    {
      ++count_hessian_skip;
      hessian_updated = false;
    }

    // stopping criterion: gradient norm
    double gn = g.norm();
    status_.projected_gradient_norm = gn;
    if(status_.converged_to_local_optimum())
    {
      DEB_line_if(!config_.silent, 4,
                  "iter = " << iter << ", f(x) = " << fx
                            << ", |grad| = " << gn
                            << ", #PCG_iters = " << n_pcg_iters
                            << " -----> converged to desired accuracy!");
      break;
    }

//      // stopping criterion: relative progress
//      if( iter > 0 && (fx_old-fx) < relative_improvement_thres * std::abs(fx_old))
//      {
//        DEB_line(4,
//                 "iter = " << iter << ", f(x) = " << fx
//                           << ", |grad| = " << gn
//                           << ", #PCG_iters = " << n_pcg_iters
//                           << ", f(x_old) - f(x) = " << (fx_old-fx)
//                           << " -----> converged because of small relative progress!");
//        converged_ = true;
//        break;
//      }

    // stop if line search could not make any progress
    if(status_.line_search_t == 0.0)
    {
      DEB_line_if(!config_.silent, 4,
                  "iter = " << iter << ", f(x) = " << fx
                            << ", |grad| = " << gn
                            << ", #PCG_iters = " << n_pcg_iters
                            << " -----> stopped because line search could not make any progress!");
      break;
    }

    // set new CG tolerance
//      double eta = std::min(0.1,gn);
    double eta = std::min(0.1,std::sqrt(gn)); // ToDo: is this one suggested in Nocedal better?
//      double eta = std::min(0.1, std::sqrt(gn))*gn;  // eigen already uses the *gn factor internally!
    pcg.setTolerance(eta);

    pcg.compute(H);
    if(pcg.info() != Eigen::Success)
    {
      std::cerr << "Warning: pcg.compute(H) did not succeed" << std::endl;
    }

    // is previous step a descent direction?
    bool warm_start = false;
    double gdx = g.dot(dx);
    status_.newton_decrement = gdx;
    if(gdx < 0.0 && config_.allow_warmstart)
    {
      warm_start = true;
      double s = -gdx/(dx.transpose()*H*dx);
      xt = s*dx;
      dx = pcg.solveWithGuess(-g, xt);

      // make sure that at least 1 iteration is done
      if(pcg.iterations() == 0)
      {
        // temporarily set other parameters
        pcg.setTolerance(1e-12);
        pcg.setMaxIterations(1);

        dx = pcg.solveWithGuess(-g, xt);

        // reset parameters
        pcg.setTolerance(eta);
        pcg.setMaxIterations(max_pcg_iters);
      }
    }
    else
    {
      dx = pcg.solve(-g);

      // make sure that at least 1 iteration is done
      if(pcg.iterations() == 0)
      {
        // temporarily set other parameters
        pcg.setTolerance(1e-12);
        pcg.setMaxIterations(1);

        dx = pcg.solve(-g);

        // reset parameters
        pcg.setTolerance(eta);
        pcg.setMaxIterations(max_pcg_iters);
      }
    }

    status_.cg_converged = (pcg.info() == Eigen::Success);

    // update statistics
    status_.cg_iterations_total += pcg.iterations();

    // get maximal reasonable step
    double t_max  = std::min(1.0,
                             config_.max_feasible_step_safety_factor * _problem->max_feasible_step(x.data(), dx.data()));
    status_.line_search_t_max_feasible = t_max;

    // backtracking line search
    gdx = g.dot(dx);
    double t=t_max;
    xt = x + t*dx;
    double fxt = _problem->eval_f(xt.data());
    int iter_ls = 0;
    ls_succeeded = false;
    while(! (fxt <= fx + config_.alpha_ls*gdx*t) )
    {
      t *= config_.beta_ls;
      xt = x + t*dx;
      fxt = _problem->eval_f(xt.data());
      ++iter_ls;
      if(iter_ls >= max_iter_ls) // maximum number of allowed iterations reached?
      {
        t = 0.0;
        break;
      }
    }

    // store line search truncation
    status_.line_search_t = t;
    status_.line_search_iterations = iter_ls;
    // store relative reduction
    rel_objective_decrease = std::abs((fx-fxt)/fx);

    // update x, fx and fx_old
    if( status_.line_search_t > 0.0)
    {
      x.swap(xt);
      fx_old = fx;
      fx = fxt;
      status_.fx = fx;
      ls_succeeded = true;
    }

    if(status_.converged_to_local_optimum())
    {
      DEB_line_if(!config_.silent, 4,
                  "iter = " << iter << ", f(x) = " << fx
                            << ", |grad| = " << gn
                            << ", gdx = " << gdx
                            << ", #PCG_iters = " << n_pcg_iters
                            << " -----> converged to desired accuracy!");
      break;
    }


    DEB_line_if(!config_.silent, 4,
                "iter = " << iter << ", f(x) = " << fx << ", t = " << t
                          << " (tmax=" << status_.line_search_t_max_feasible << "), " << "#ls = "
                          << status_.line_search_iterations
                          << ", |grad| = " << status_.projected_gradient_norm
                          << ", " << "PCG_tol = " << eta
                          << ", " << "PCG_iters = " << int(pcg.iterations())
                          << " (total = " << n_pcg_iters << ") "
                          << ", " << "PCG_converged = " << int(status_.cg_converged)
                          << ", " << "PCG_warmstart = " << int(warm_start)
                          << ", " << "gdx = " << status_.newton_decrement
                          << ", " << "hessian_update = " << int(hessian_updated)
                          << ", " << "rel_obj_decrease = " << rel_objective_decrease
    );
  }

  // store result
  _problem->store_result(x.data());

  // return success
  return status_.converged_to_local_optimum();
}


//-----------------------------------------------------------------------------

int
TruncatedNewtonPCG::
solve(NProblemInterface *_problem, const SMatrixD &_A, const VectorD &_b)
{
  // partially based on ideas of
  // Gould, Nicholas IM, Mary E. Hribar, and Jorge Nocedal.
  // "On the solution of equality constrained quadratic programming problems arising in optimization."
  //  SIAM Journal on Scientific Computing 23.4 (2001): 1376-1395.

//    DEB_time_func_def;
  DEB_enter_func;

  // reset status
  status_ = OptimizerStatus();
  status_.cg_iterations_total = 0;
  status_.n_newton_iters = 0;
  status_.refinement_iters_total = 0;

//    converged_ = false;
//    feasible_solution_found_ = false;

  // number of unknowns
  Eigen::Index n = _problem->n_unknowns();
  Eigen::Index m = _A.rows();

  DEB_line_if(!config_.silent, 2,
              "optimize via TruncatedNewton with " << (int) n << " unknowns and " << (int) m
                                                                              << " linear constraints");

  // allow only steps larger than eps_ls_ in line search
  const int max_iter_ls = int(std::log(config_.eps_ls) / std::log(config_.beta_ls));
  const int max_iter_ls_negative_curvature = 20;

  // initialize mu for merit-function based line search
  double mu_merit = 0.0;

  // initialize vectors of unknowns
  VectorD x(n);
  _problem->initial_x(x.data());

  status_.constraint_violation_inf_norm = (_A * x - _b).lpNorm<Eigen::Infinity>();
  if (status_.constraint_violation_inf_norm < config_.eps_constraints_violation)
  {
    feasible_solution_found_ = true;
    status_.feasible = true;
  }
//  DEB_line_if(!silent_, 2, "initial inf-norm constraint violation = " << status_.constraint_violation_inf_norm
//                                                                      << " numerically feasible = "
//                                                                      << int(status_.feasible));

  // storage of update vector dx and gradient
  VectorD g(n), dx(n), xn(n), gz(n), dz(n);
  dx.setZero(); // required for warmstart

  // storage CG
  VectorD r(n), v(n), q(n), p(n), r2(n), q2(n), Hp(n);

  // storage for iterative refinement
  // ToDo: can we reduce number of temporaries by using them in several non-overlapping places?
  VectorD rho_q, rho_v, dq, dv, A_inv_row_norm;
  if (config_.enable_iterative_refinement)
  {
    rho_q.resize(n);
    rho_v.resize(n);
    dq.resize(n);
    dv.resize(n);

    // initialize row norms
    A_inv_row_norm.resize(_A.rows());
    A_inv_row_norm.setZero();
    for (int k = 0; k < _A.outerSize(); ++k)
      for (SMatrixD::InnerIterator it(_A, k); it; ++it)
        A_inv_row_norm[it.row()] += std::pow(it.value(), 2);
    for (int k = 0; k < A_inv_row_norm.size(); ++k)
    {
      double d = 1.0 / std::sqrt(A_inv_row_norm[k]);
      if (std::isfinite(d))
        A_inv_row_norm[k] = d;
      else
      {
        std::cerr << "Warning: row " << k << " of _A has degenerate norm of " << A_inv_row_norm[k] << std::endl;
        A_inv_row_norm[k] = 1.0;
      }
    }
  }

  // hessian matrix
  SMatrixD H(n, n);

  // Diagonal Preconditioner
  VectorD W(n), Wi(n);

  std::unique_ptr<MatrixDecomposition<double>> decomposed_projection;

  // get function value at current point
  status_.fx = _problem->eval_f(x.data());
//  DEB_line_if(!silent_, 2, "initial objective value = " << status_.fx);
  if (!std::isfinite(status_.fx))
  {
    DEB_line_if(!config_.silent, 2, "ERROR: intial objective value is not finite ---> abort");
    return false;
  }

  // data for adaptive hessian update
  int count_hessian_skip = 0;
  status_.hessian_updated = false;
  status_.line_search_t = 0.0;
  double rel_objective_decrease = 0.0;

  // print initial data
  print_iteration_data(status_);

  // link iter to status_.n_newton_iters
  int& iter = status_.n_newton_iters;
  for (iter=1; iter <= config_.max_iters; ++iter)
  {
    // initialize
    status_.refinement_iters = 0;

    // get gradient and Hessian
    _problem->eval_gradient(x.data(), g.data());

    const bool update_hessian =
            (iter == 0 || // first iteration
             count_hessian_skip >= config_.hessian_max_skips || // maximum number of Hessians have been skipped
             status_.line_search_t < config_.hessian_min_acceptable_alpha || // line search truncated too much
             rel_objective_decrease < config_.hessian_min_acceptable_rel_objective_decrease); // objective did not decrease sufficiently


    if(update_hessian) // update Hessian?
    {
      _problem->eval_hessian(x.data(), H);
      status_.hessian_updated = true;
      count_hessian_skip = 0;
    }
    else
    {
      ++count_hessian_skip;
      status_.hessian_updated = false;
    }

    // pre-factor projection matrix and update preconditioner
    if (iter == 0 ||
        (config_.always_update_preconditioner && status_.hessian_updated))
    {
      // setup preconditioner W = (|diag(H)|)
      // and Wi = (|diag(H)|)^-1

      for (int j = 0; j < n; ++j)
      {
        double d = std::abs(H.coeffRef(j, j));

        if (!std::isfinite(d))
        {
          d = 1.0;
          DEB_line_if(!config_.silent, 2, "ERROR: diagonal entry of Hessian is not finite ---> ignore for preconditioner");
        }

        // clamp preconditioner range
        if (d < config_.min_preconditioner_value)
        {
//          DEB_line_if(!silent_, 2, "iter = " << iter << ", increase preconditioner value from " << d << " to " << min_preconditioner_value_);
          d = config_.min_preconditioner_value;
        }

        if (d > config_.max_preconditioner_value)
        {
//          DEB_line_if(!silent_, 2, "iter = " << iter << ", decrease preconditioner value from " << d << " to " << max_preconditioner_value_);
          d = config_.max_preconditioner_value;
        }

        W[j] = d;
        Wi[j] = 1.0 / W[j];

//        minW = std::min(minW,W[j]);
//        maxW = std::max(maxW,W[j]);

//          // no preconditioner
//          W[j]  = 1.0;
//          Wi[j] = 1.0;
      }

      // prepare constraint projection
//      if (iter == 0) DEB_line_if(!silent_, 2, "prepare constraint projection ...");

      SMatrixD AWiAt = _A * Wi.asDiagonal() * _A.transpose();

      // update factorization (if valid matrix)
      if (_A.rows() > 0 && _A.cols() > 0)
      {
        COMISO::StopWatch swf;
        swf.start();
        if (!decomposed_projection)
        {
          decomposed_projection = make_decomposition<double>(config_.matrix_decomposition_algo);
          decomposed_projection->analyzePattern(AWiAt);
        }
        decomposed_projection->factorize(AWiAt);
//        std::cerr << "factorization took = " << swf.stop() / 1000.0 << " seconds" << std::endl;

        // ldlt.compute(AWiAt); // old update

        if (decomposed_projection->info() != Eigen::Success)
        {
          for (unsigned int j = 0; j < 10; ++j)
          {
            double reg = 1e-8 * AWiAt.diagonal().sum() / double(m);
            DEB_line_if(!config_.silent, 2,
                        "Warning: LDLT factorization failed --> regularize ADAt (could lead to reduced accuracy), reg="
                                << reg);
            for (Eigen::Index j = 0; j < AWiAt.rows(); ++j)
              AWiAt.coeffRef(j, j) += reg;  // operation is safe since all diagonal entries are nonzero!!!

            decomposed_projection->compute(AWiAt);
            if (decomposed_projection->info() == Eigen::Success)
              break;
          }
        }
      }
//      if (iter == 0) DEB_line_if(!silent_, 2, "done!");
    }

    // perform feasiblity step?
    status_.feasibility_step_productive = false;
    if (status_.constraint_violation_inf_norm >= config_.eps_constraints_violation_desirable)
    {
//        // optimize constraint violation in Hessian-diagonal norm
//        v  = ldlt.solve(_b-_A*x);
//        dz = Wi.asDiagonal()*_A.transpose()*v;

      // optimize full KKT-system but approximate Hessian only through diagonal
      if (_A.rows() > 0 && _A.cols() > 0)
        v = decomposed_projection->solve(_b - _A * x + _A * Wi.asDiagonal() * g);
      else
        v.setZero();

      dz = Wi.asDiagonal() * (_A.transpose() * v - g);

      // debug
//        std::cerr << "debug: constraint violation after full projection = " << (_A*(x+dz)-_b).norm() << std::endl;

      // truncate feasiblity step
      double t_max = std::min(1.0,
                              config_.max_infeasibility_step_safety_factor * _problem->max_feasible_step(x.data(), dz.data()));

      status_.line_search_t_inf_max_feasible = t_max;

      //      std::cerr << "debug: t_max = " << t_max << std::endl;

      if (config_.line_search_feasibility_step)
      {
        double fxn = status_.fx;
        double t = backtracking_line_search_infeasible_merit_l1(_problem, H, _A, _b, x, status_.fx,
                                                                g, dz,
                                                                xn, fxn, mu_merit,
                                                                t_max, max_iter_ls, status_.line_search_inf_iterations);
        status_.line_search_t_inf = t;

        if (t > 0.0)
        {

          // line search alrady computes xn
          // xn = x + t*dz;

          double constraint_violation_new = (_A * xn - _b).lpNorm<Eigen::Infinity>();

//          DEB_line_if(!silent_, 2,
//                      "iter = " << iter << ", feasibility step (merit-based) reducing constraint violation "
//                                << status_.constraint_violation_inf_norm << " ---> "
//                                << constraint_violation_new << "  (t=" << t_max << ")");

          // update objective
          double fx_new = _problem->eval_f(xn.data());
          if (std::isfinite(fx_new)) // valid objective function?
          {
            // update
            xn.swap(x);
            status_.feasibility_step_productive = true;
            status_.constraint_violation_inf_norm = constraint_violation_new;
            status_.fx = fx_new;

            // feasible solution found?
            if (status_.constraint_violation_inf_norm < config_.eps_constraints_violation)
              status_.feasible = true;

            // update gradient
            _problem->eval_gradient(x.data(), g.data());
          }
          else
          {
            DEB_line_if(!config_.silent, 2, "Warning: feasibility step line search truncated to 0.0");
            status_.feasibility_step_productive = false;
          }
        }
      }
      else // do not perform line search ---> old option
      {
        // new x
        xn = x + t_max * dz;

        double constraint_violation_new = (_A * xn - _b).lpNorm<Eigen::Infinity>();

        // is step productive? (can be numerically non-productive)
        if (constraint_violation_new < status_.constraint_violation_inf_norm)
        {
//          DEB_line_if(!silent_, 2, "iter = " << iter << ", feasibility step reducing constraint violation "
//                                             << status_.constraint_violation_inf_norm << " ---> "
//                                             << constraint_violation_new << "  (t=" << t_max << ")");

          // update objective
          double fx_new = _problem->eval_f(xn.data());
          if (std::isfinite(fx_new)) // valid objective function?
          {
            // update
            xn.swap(x);
            status_.feasibility_step_productive = true;
            status_.constraint_violation_inf_norm = constraint_violation_new;
            status_.fx = fx_new;

            // feasible solution found?
            if (status_.constraint_violation_inf_norm < config_.eps_constraints_violation)
              status_.feasible = true;

            // update gradient
            _problem->eval_gradient(x.data(), g.data());
          }
          else
          {
            DEB_line_if(!config_.silent, 2, "Warning: feasibility step resulted in non-finite objective value --> revert to previous x");
          }
        }
        else
        {
          DEB_line_if(!config_.silent, 2, "Warning: feasibility step was not productive --> revert to previous x");
        }
      }
    }


    // ---------------------------
    // Projected-Preconditioned-CG

    // choose starting point
    bool warmstart = false;
    if(config_.allow_warmstart && g.dot(dx) < 0.0)
      warmstart = true;
    else
      dx.setZero();

    r = g;

    auto preconditioned_projection = [&](const VectorD &_v, VectorD &_v_proj, VectorD &_v_help) {
      // preconditioned-projection _v --> _v_proj
      // note: _v_help is used later for constraint refinement

      if (_A.rows() == 0 || _A.cols() == 0)
      {
        _v_help.setZero();
        _v_proj = Wi.asDiagonal() * _v;
        return;
      }

      _v_help = decomposed_projection->solve(_A * Wi.asDiagonal() * _v);
      _v_proj = Wi.asDiagonal() * (_v - _A.transpose() * _v_help);

      // perform iterative refinement
      if (config_.enable_iterative_refinement)
      {
        int n_iters = 0;
        for (; n_iters < config_.max_iterative_refinement_iters; ++n_iters)
        {
          // desired accuracy reached?
          if (max_abs_cos_angle(_A, A_inv_row_norm, _v_proj) < config_.iterative_refinement_cos_angle_threshold)
            break;

          rho_q = _v - W.asDiagonal() * _v_proj - _A.transpose() * _v_help;
          rho_v = -_A * _v_proj;

          dv = decomposed_projection->solve(_A * Wi.asDiagonal() * rho_q - rho_v);
          dq = Wi.asDiagonal() * (rho_q - _A.transpose() * dv);

          _v_proj += dq;
          _v_help += dv;
        }
        // store number of refinement iterations
        status_.refinement_iters += n_iters;
        status_.refinement_iters_total += n_iters;
      }
    };

    // project residual/gradient
    preconditioned_projection(r, q, v);
    // constraint refinement
    gz = r - _A.transpose() * v;
    r = gz;
    // initialize p
    p = -q;

    double rtq = r.dot(q);

    // check convergence
    status_.projected_gradient_norm = gz.norm();
    status_.projected_gradient_norm_within_tolerance = (status_.projected_gradient_norm <= config_.eps);
    if (status_.converged_to_local_optimum())
      break;

    // choose CG tolerance
    double eta;
    if(config_.adaptive_tolerance)
      eta = config_.adaptive_tolerance_modifier*std::min(0.1, std::sqrt(status_.projected_gradient_norm))*status_.projected_gradient_norm;
    else
      eta = config_.pcg_tolerance*status_.projected_gradient_norm;

//      if(gzn >= newton_tol) // avoid numerical issues when still infeasible but projected gradient is vanishing
    {
      int n_pcg_iters = 0;
      for (int pcg_iter = 0; pcg_iter < config_.max_pcg_iters; ++pcg_iter)
      {
        // compute norm of projected residual (in original norm)
        double rpn = (W.asDiagonal() * q).norm();

        // stop if accuracy sufficient
        if (rpn < eta)
        {
          status_.cg_converged = true;
          break;
        }

        // cache re-used quantities
        Hp = H * p;
        double pHp = p.dot(Hp);

        // handle directions of negative curvature (important for Hessians of non-convex problems)
        status_.negative_curvature_step = false;
        if (pHp <= 0.0)
        {
          status_.negative_curvature_step = true;
          if (pcg_iter == 0)
          {
            // use p = proj(-grad f)
            // step is not well scaled but line search will take care of this
            dx = p;
          }

          // quit PCG iteration
          break;
        }

        // optimal step length along p (minimizer of quadratic)
        double alpha = rtq / pHp;

        // update dx
        dx += alpha * p;
        r2 = r + alpha * Hp;

        // preconditioned-projection of residual/gradient r2 --> q2
        // note: v is used later for constraint refinement
        preconditioned_projection(r2, q2, v);

        double rtq2 = r2.dot(q2);
        double beta = rtq2 / rtq;
        p = -q2 + beta * p;

        // swap vectors
        q.swap(q2);
        // r.swap(r2);

        // constraint refinement
        r = r2 - _A.transpose() * v;

        // update rtq
        rtq = rtq2;

        // count number of pcg iters
        ++n_pcg_iters;
      }

      if (config_.project_dx_before_update)
      {
        // project dx
        if (_A.rows() > 0 && _A.cols() > 0)
          v = decomposed_projection->solve(_A * dx);
        else
          v.setZero();
        gz = dx - Wi.asDiagonal() * _A.transpose() * v;
        dx = gz;
      }

      // store number of pcg iterations
      status_.cg_iterations = n_pcg_iters;
      status_.cg_iterations_total += n_pcg_iters;

      // negative curvature steps
      if(status_.negative_curvature_step)
        status_.n_negative_curvature_iters += 1;

      // backtracking line search
      status_.line_search_t = 0.0;
      status_.line_search_t_max_feasible = 0.0;
      status_.line_search_iterations = 0;
      double gdx = g.dot(dx);
      // update Newton decrement
      status_.newton_decrement = gdx;
      status_.newton_decrement_within_tolerance = (std::abs(status_.newton_decrement) < config_.eps_gdx);

      if (gdx < 0.0) // require descent direction!!!
      {
        // get maximal reasonable step
        status_.line_search_t_max_feasible =
                config_.max_feasible_step_safety_factor * _problem->max_feasible_step(x.data(), dx.data());
        double t_max = status_.line_search_t_max_feasible;

        double fxn = status_.fx;
        int    iter_ls = 0;
        double t = 0.0;

        if(status_.negative_curvature_step)
        {
          t = line_search_negative_curvature(_problem, x, status_.fx,
                                             dx, t_max, max_iter_ls_negative_curvature,
                                             xn, fxn, iter_ls);
        }
        else // regular step
        {
          // for regular case target full Newton steps
          if(t_max > 1.0)
            t_max = 1.0;

          t = backtracking_line_search(_problem, x, status_.fx,
                                       dx, gdx, t_max, max_iter_ls,
                                       xn, fxn, iter_ls);
        }

        // store line search truncation
        status_.line_search_t = t;
        status_.line_search_iterations = iter_ls;
        // store relative reduction
        rel_objective_decrease = std::abs((status_.fx - fxn) / status_.fx);

        // update x, fx, xz (if line search productive
        if (status_.line_search_t > 0.0)
        {
          // TODO: original dx is needed to compute dual variables such that update here is harmful
          if (config_.allow_warmstart) // update dx if needed afterwards
            dx = xn - x;
          x.swap(xn);
          status_.fx = fxn;
          // update constraint violation
          status_.constraint_violation_inf_norm = (_A * x - _b).lpNorm<Eigen::Infinity>();
          status_.feasible = (status_.constraint_violation_inf_norm < config_.eps_constraints_violation);
        }
      }
      else
      {
        DEB_line_if(!config_.silent, 2, "Warning: dx is not a descent direction, gdx = " << gdx);
      }

      // line-search not productive?
      if (status_.line_search_t == 0.0)
      {
        if (!status_.feasible && status_.feasibility_step_productive)
        {
          DEB_line_if(!config_.silent, 4,
                      "Info: line search failed but infeasible ---> continue since feasibility_step_productive");
        }
        else
        {
          DEB_line_if(!config_.silent, 4,
                      "Warning: line search failed at infeasible point and feasibility step non-productive ---> terminate");

          break;
        }
      }

      // output iteration data
      print_iteration_data(status_);

      // converged?
      if (status_.converged_to_local_optimum())
        break;
    }
  }

  n_iterations_used_ = iter; // TODO: or increase? may make more sense for incremental solves.
  status_.n_newton_iters = iter;

  // store result
  _problem->store_result(x.data());

  // compute dual variables
  if (config_.compute_dual_variables)
  {
    if (!decomposed_projection)
    {
      std::cerr << "ERROR: dual variables cannot be computed since ldlt is not initialized!!!" << std::endl;
      return status_.converged_to_local_optimum();
    }

    // Variant I ---> compute dual variables with preconditioner (does not require additional factorization)
    // ToDo: is dx uptodate here?
    if (_A.rows() > 0 && _A.cols() > 0)
      nue_ = decomposed_projection->solve(-_A * Wi.asDiagonal() * (g + H * dx));
    else
      nue_.setZero();

//      // Variant II ---> no preconditioner to compute dual variables (requires additional factorization)
//      ldlt.compute(_A*_A.transpose());
//      nue_ = ldlt.solve(-_A*(g+H*dx));

//      std::cerr << "Truncated Newton dual variables residual ||A^T nue + g + H dx|| = " << (_A.transpose()*nue_+g+H*dx).norm() << std::endl;
//      std::cerr << "||dx|| = " << dx.norm() << std::endl;
  }

  // print summary
  print_summary(status_);

  // return success
  return status_.converged_to_local_optimum();
}


//-----------------------------------------------------------------------------


int
TruncatedNewtonPCG::
solve_projected_normal_equation(NProblemInterface *_problem, const SMatrixD &_A, const VectorD &_b)
{
  return solve(_problem, _A, _b);
}


//-----------------------------------------------------------------------------


int
TruncatedNewtonPCG::
solve(NProblemInterface *_problem, std::vector<LinearConstraint> &_constraints)
{
  // convert constraints
  SMatrixD A;
  VectorD b;
  LinearConstraintConverter::nsolver_to_eigen(_constraints, A, b, _problem->n_unknowns());

  return solve(_problem, A, b);
}


//-----------------------------------------------------------------------------

int
TruncatedNewtonPCG::
solve(NProblemInterface *_problem, std::vector<NConstraintInterface *> &_constraints)
{
  // convert constraints
  SMatrixD A;
  VectorD b;
  LinearConstraintConverter::nsolver_to_eigen(_constraints, A, b, _problem->n_unknowns());

  return solve(_problem, A, b);
}


//-----------------------------------------------------------------------------


int
TruncatedNewtonPCG::
solve_projected_normal_equation(NProblemInterface *_problem, std::vector<LinearConstraint> &_constraints)
{
  // convert constraints
  SMatrixD A;
  VectorD b;
  LinearConstraintConverter::nsolver_to_eigen(_constraints, A, b, _problem->n_unknowns());

  return solve_projected_normal_equation(_problem, A, b);
}


//-----------------------------------------------------------------------------

int
TruncatedNewtonPCG::
solve_projected_normal_equation(NProblemInterface *_problem, std::vector<NConstraintInterface *> &_constraints)
{
  // convert constraints
  SMatrixD A;
  VectorD b;
  LinearConstraintConverter::nsolver_to_eigen(_constraints, A, b, _problem->n_unknowns());

  return solve_projected_normal_equation(_problem, A, b);
}


//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------


double
TruncatedNewtonPCG::
backtracking_line_search_infeasible_merit_l1(NProblemInterface *_problem, const SMatrixD &_H,
                                             const SMatrixD &_A, const VectorD &_b,
                                             const VectorD &_x, const double &_fx,
                                             const VectorD &_g, VectorD &_dx,
                                             VectorD &_x_new, double &_fx_new,
                                             double &_mu_merit,
                                             const double _t_start, const int _max_ls_iters, int& _n_iters) const
{
  DEB_enter_func;
  size_t n = _x.size();

  // update mu
  double res_primal_1 = (_A * _x - _b).template lpNorm<1>();
  double gdx = _g.transpose() * _dx;
  double dxHdx = _dx.transpose() * _H * _dx;
  double mu_new = 1.2 * (gdx + 0.5 * std::max(0.0, dxHdx)) / ((1.0 - config_.rho_merit) * res_primal_1);
  _mu_merit = std::max(mu_new, _mu_merit);

  // current step size
  double t = _t_start;

  // merit function and directional derivative for t=0
  double merit_0 = _fx + _mu_merit * res_primal_1;
//  double merit_0  = _problem->eval_f(_x.data()) + _mu_merit*res_primal_1;
  double D_merit_0 = gdx - _mu_merit * res_primal_1;

//  std::cerr << "mu_merit=" << mu_merit_ << std::endl;
//  std::cerr << "D_merit_0=" << D_merit_0 << std::endl;

  double fx(0.0);

  _n_iters = 0;

  // backtracking (stable in case of NAN and with max iterations)
  for (int i = 0; i < _max_ls_iters; ++i)
  {
    // current update of x, nue and g
    _x_new = _x + _dx * t;
    fx = _problem->eval_f(_x_new.data());
    // check if update is inside domain
    if (std::isfinite(fx))
    {
      double merit_t = fx + _mu_merit * (_A * _x_new - _b).template lpNorm<1>();

      //     std::cerr << "t=" << t << ", merit(t)=" << merit_t << ", merit(0)+alpha*t*Dmerit(0)=" << merit_0 + alpha_ls_*t*D_merit_0 << std::endl;

      // sufficient decrease in residual?
      if (merit_t <= merit_0 + config_.alpha_ls * t * D_merit_0)
      {
        // succesfull line search
        _fx_new = fx;
        return t;
      }
    }
    //   else std::cerr << "t=" << t << ", merit(t)=" << std::numeric_limits<double>::infinity() << std::endl;

    // shrink with factor beta
    t *= config_.beta_ls;
    ++_n_iters;
  }

  // restore
  _x_new = _x;

  DEB_warning(1, "line search could not find a valid step");
  return 0.0;
}

//-----------------------------------------------------------------------------


double
TruncatedNewtonPCG::
backtracking_line_search(NProblemInterface* _problem,
                         const VectorD& _x, const double _fx,
                         const VectorD& _dx,
                         const double _gdx, const double _t_max, const int _max_iter_ls,
                         VectorD& _x_new, double& _fx_new, int& _iter_ls) const
{
  double t = _t_max;

  _x_new = _x + t * _dx;
  _fx_new = _problem->eval_f(_x_new.data());

  _iter_ls = 0;
  while (!(_fx_new <= _fx + config_.alpha_ls * _gdx * t))
  {
    t *= config_.beta_ls;
    _x_new = _x + t * _dx;
    _fx_new = _problem->eval_f(_x_new.data());
    ++_iter_ls;
    // maximum number of steps reached without finding a valid step?
    if (_iter_ls >= _max_iter_ls)
    {
      _x_new  = _x;
      _fx_new = _fx;
      return 0.0;
    }
  }

  return t;
}


//-----------------------------------------------------------------------------


double
TruncatedNewtonPCG::
line_search_negative_curvature( NProblemInterface* _problem,
                                const VectorD& _x, const double _fx,
                                const VectorD& _dx,
                                const double _t_max, const int _max_iter_ls,
                                VectorD& _x_new, double& _fx_new, int& _iter_ls) const
{
  const double beta_inc_negative_curvature = 2.0;
  const double beta_dec_negative_curvature = 0.5;
  const double expected_decrease = 0.0;

  double t = std::min(1.0,_t_max);

  _x_new = _x + t * _dx;
  _fx_new = _problem->eval_f(_x_new.data());

  double  t_test, fx_test;
  VectorD x_test;

  _iter_ls = 0;

  // full step already sufficiently better --> try to grow
  if(_fx_new < (1.0-t*expected_decrease)*_fx)
  {
    while (_iter_ls < _max_iter_ls)
    {
      // cannot grow further?
      if(t == _t_max)
        return t;

      t_test  = std::min(beta_inc_negative_curvature * t, _t_max);
      x_test  = _x + t_test * _dx;
      fx_test = _problem->eval_f(x_test.data());
      ++_iter_ls;

      // better than before?
      if(fx_test < _fx_new)
      {
        _fx_new = fx_test;
        _x_new.swap(x_test);
        t = t_test;
      }
      else
        return t; // no improvement --> terminate
    }
    return t;
  }

  // backtrack while improving
  while (_iter_ls < _max_iter_ls)
  {
    // shrink t
    t *= beta_dec_negative_curvature;

    t_test  = beta_dec_negative_curvature * t;
    x_test  = _x + t_test * _dx;
    fx_test = _problem->eval_f(x_test.data());
    ++_iter_ls;

    // better than before?
    if(fx_test < _fx_new)
    {
      _fx_new = fx_test;
      _x_new.swap(x_test);
      t = t_test;
    }
    else
    {
      // sufficiently good --> done
      if (fx_test < (1.0-t*expected_decrease)*_fx)
        return t;
    }
  }

  // maximum number of steps reached without finding a valid step?
  if (_iter_ls >= _max_iter_ls)
  {
    _x_new  = _x;
    _fx_new = _fx;
    return 0.0;
  }

  return t;
}

//-----------------------------------------------------------------------------


double
TruncatedNewtonPCG::
max_abs_cos_angle(const SMatrixD& _A, const VectorD& _A_inv_row_norm, const VectorD& _v) const
{
  assert(_A_inv_row_norm.size() == _A.rows());

  double s = 1.0/_v.norm();
  if(!std::isfinite(s))
  {
    std::cerr << "Warning: max_abs_cos_angle called with degenerate vector" << std::endl;
    return 0.0;
  }

  // compute dot-products with rows
  VectorD dp = _A*_v;

  double max_abs_cos = 0.0;
  for(int i=0; i<dp.size(); ++i)
    max_abs_cos = std::max(max_abs_cos, std::abs(dp[i]*s*_A_inv_row_norm[i]));

  return max_abs_cos;
}

//-----------------------------------------------------------------------------


void
TruncatedNewtonPCG::
print_iteration_data(const OptimizerStatus& _status) const
{
  DEB_enter_func;

  std::stringstream ss_out;

  if(!config_.silent && _status.n_newton_iters % 10 == 0)
  {
    ss_out    << std::left << std::setw(6) << "iter" << " | "
              << std::left << std::setw(11) << "f(x)" << " | "
              << std::left << std::setw(11) << "||Ax-b||" << " | "
              << std::left << std::setw(8) << "||g_p||" << " | "
              << std::left << std::setw(9) << "gdx" << " | "
              << std::left << std::setw(8) << "t" << " | "
              << std::left << std::setw(8) << "t_max" << " | "
              << std::left << std::setw(3) << "#ls" << " | "
              << std::left << std::setw(8) << "ti" << " | "
              << std::left << std::setw(8) << "ti_max" << " | "
              << std::left << std::setw(4) << "#lsi" << " | "
              << std::left << std::setw(3) << "PCG" << " | "
              << std::left << std::setw(3) << "REF" << " | "
              << std::left << std::setw(1) << "C" << " | "
              << std::left << std::setw(1) << "N" << " | "
              << std::left << std::setw(1) << "H" << " | "
              << std::endl;
  }

  if(!config_.silent)
  {
    char feasible = ' ';
    if(_status.feasible)
      feasible = '*';

    if(_status.n_newton_iters == 0)
    {
      ss_out    << std::setprecision(4) << std::scientific
                << std::left << std::setw(2) << feasible
                << std::left << std::setw(4) << _status.n_newton_iters << " | "
                << std::left << std::setw(11) << _status.fx << " | "
                << std::left << std::setw(11) << _status.constraint_violation_inf_norm << " | "
                << std::left << std::setw(8) << "-" << " | "
                << std::left << std::setw(9) << "-" << " | "
                << std::left << std::setw(8) << "-" << " | "
                << std::left << std::setw(8) << "-" << " | "
                << std::left << std::setw(3) << "-" << " | "
                << std::left << std::setw(8) << "-" << " | "
                << std::left << std::setw(8) << "-" << " | "
                << std::left << std::setw(4) << "-" << " | "
                << std::left << std::setw(3) << "-" << " | "
                << std::left << std::setw(3) << "-" << " | "
                << std::left << std::setw(1) << "-" << " | "
                << std::left << std::setw(1) << "-" << " | "
                << std::left << std::setw(1) << "-" << " | ";
    }
    else
    {
      ss_out    << std::setprecision(4) << std::scientific
                << std::left << std::setw(2) << feasible
                << std::left << std::setw(4) << _status.n_newton_iters << " | "
                << std::left << std::setw(11) << _status.fx << " | "
                << std::left << std::setw(11) << _status.constraint_violation_inf_norm << " | "
                << std::setprecision(2)
                << std::left << std::setw(8) << _status.projected_gradient_norm << " | "
                << std::left << std::setw(9) << _status.newton_decrement << " | "
                << std::left << std::setw(8) << _status.line_search_t << " | "
                << std::left << std::setw(8) << _status.line_search_t_max_feasible << " | "
                << std::left << std::setw(3) << _status.line_search_iterations << " | "
                << std::left << std::setw(8) << _status.line_search_t_inf << " | "
                << std::left << std::setw(8) << _status.line_search_t_inf_max_feasible << " | "
                << std::left << std::setw(4) << _status.line_search_inf_iterations << " | "
                << std::left << std::setw(3) << _status.cg_iterations << " | "
                << std::left << std::setw(3) << _status.refinement_iters << " | "
                << std::left << std::setw(1) << int(_status.cg_converged) << " | "
                << std::left << std::setw(1) << int(_status.negative_curvature_step) << " | "
                << std::left << std::setw(1) << int(_status.hessian_updated) << " | ";
    }
  }

  DEB_line_if(!config_.silent, 2, ss_out.str());
}


//-----------------------------------------------------------------------------


void
TruncatedNewtonPCG::
print_summary(const OptimizerStatus& _status) const
{
  DEB_enter_func;

  DEB_line_if(!config_.silent, 2, "------------------ TruncatedNewtonSummary BEGIN ------------------");
  DEB_line_if(!config_.silent, 2, "Converged to local optimum = " << int(_status.converged_to_local_optimum()));
  DEB_line_if(!config_.silent, 2, "objective f(x)   = " << _status.fx);
  DEB_line_if(!config_.silent, 2, "||Ax-b||_inf     = " << _status.constraint_violation_inf_norm);
  DEB_line_if(!config_.silent, 2, "||proj(grad)||   = " << status_.projected_gradient_norm);
  DEB_line_if(!config_.silent, 2, "Newton decrement = " << status_.newton_decrement);
  DEB_line_if(!config_.silent, 2, "#iters Newton    = " << status_.n_newton_iters);
  DEB_line_if(!config_.silent, 2, "#iters PCG       = " << status_.cg_iterations_total);
  DEB_line_if(!config_.silent, 2, "#iters Neg. curv.= " << status_.n_negative_curvature_iters);
  DEB_line_if(!config_.silent, 2, "#iters refinement= " << status_.refinement_iters_total);
  DEB_line_if(!config_.silent, 2, "------------------ TruncatedNewtonSummary END ------------------");
}


// below are experimental versions, which are currently not used anymore

//  // solve with linear constraints
//  // Warning: so far only feasible starting points with (_A*_problem->initial_x() == b) are supported!
//  // It is also required that the constraints are linearly independent
//  int solve_experimental(NProblemInterface* _problem, const SMatrixD& _A, const VectorD& _b)
//  {
//    DEB_time_func_def;
//    converged_ = false;
//
//    {
//      std::cerr << "compute QR decomposition of A^T...\n";
//      SMatrixD At = _A.transpose();
//      At.makeCompressed();
//      StopWatch swq; swq.start();
//      Eigen::SparseQR<SMatrixD,Eigen::COLAMDOrdering<int> > sqr(At);
//      std::cerr << "done! time = " << swq.stop()/1000.0 << "s\n";
//    }
//
//    // number of unknowns
//    Eigen::Index n = _problem->n_unknowns();
//    Eigen::Index m = _A.rows();
//    DEB_line(2, "optimize via TruncatedNewtonPCG with " << (int)n << " unknowns and " << (int)m << " linear constraints");
//
//    // Newton parameters
//    const int    max_iters     = max_iters_;
//    const double newton_tol = eps_;
//    //    const double relative_improvement_thres = 1e-3;
//    // CG parameters
////    const int max_pcg_iters = 2*(n-m);
//    const int max_pcg_iters = 500;
//
//    // allow minimal steps of eps_ls_
//    const int max_iter_ls   = int ( std::log(eps_ls_)/std::log(beta_ls_));
//
//    // initialize vectors of unknowns
//    VectorD x(n);
//    _problem->initial_x(x.data());
//
//    // storage of update vector dx and gradient
//    VectorD dx(n), g(n), g2(n), xt(n);
//    dx.setZero();
//
//    // hessian matrix
//    SMatrixD H(n,n), H2(n,n), A2(m,n);
//
//    // preconditioner
//    VectorD D(n), D2(n);
//
//    // storage for CG
//    VectorD r(n), q(n), v(m), d(n), Hd(n), r1(n), q1(n);
//
//    // cholesky for constraint projection
//    Eigen::SimplicialLDLT<SMatrixD> ldlt;
//
//    double constraint_violation(0.0);
//
//    // statistics
//    int n_pcg_iters(0);
//
//    bool ls_succeeded = true;
//
//    // get function value at current point
//    double fx = _problem->eval_f(x.data());
//    double fx_old = fx;
//
//    for(int iter=0; iter<max_iters; ++iter)
//    {
//      // get gradient and Hessian
//      _problem->eval_gradient(x.data(), g.data());
//      _problem->eval_hessian(x.data(), H);
//
//      if(iter == 0)
//      {
//        // setup preconditioner D2 = (|diag(H)|)^(-1/2)
//        // and D = (|diag(H)|)^(1/2)
//        for(int j=0; j<n ;++j)
//        {
//          D[j]  = std::sqrt(std::abs(H.coeffRef(j,j)));
//          D2[j] = 1.0/D[j];
//        }
//
//        // check
//        if(0)
//        {
//          H2 = D2.asDiagonal() * H * D2.asDiagonal();
//          std::cerr << "************ Hessian diagonal: " << H.diagonal() << std::endl;
//          std::cerr << "************ Preconditioned Hessian diagonal: " << H2.diagonal() << std::endl;
//        }
//
//        // prepare constraint projection
//        DEB_line(2, "prepare constraint projection ..." );
//        constraint_violation = (_A*x-_b).norm();
//        DEB_line(2, "initial constraint violation = " << constraint_violation);
//        A2 = _A*D2.asDiagonal();
//        SMatrixD AAt = A2*A2.transpose();
//        ldlt.compute(AAt);
//        if(ldlt.info() != Eigen::Success)
//        {
//          for(unsigned int j=0; j<10; ++j)
//          {
//            double reg = 1e-8 * AAt.diagonal().sum() / double(m);
//            DEB_line(2,
//                     "Warning: LDLT factorization failed --> regularize AAt (could lead to reduced accuracy), reg=" << reg);
//            for (Eigen::Index j = 0; j < AAt.rows(); ++j)
//              AAt.coeffRef(j, j) += reg;
//
//            ldlt.compute(AAt);
//            if(ldlt.info() == Eigen::Success)
//              break;
//          }
//
//        }
//        DEB_line(2, "done!" );
//
//      }
//
//      // stopping criterion: gradient norm of eliminated problem (or equivalently norm of projected gradient)
//      // project g  --> g2
//      g2 = D2.asDiagonal()*g;
//
//      v = ldlt.solve(A2*g2);
//      g2 -= A2.transpose()*v;
//      g2 = D.asDiagonal()*g2;
//
//      double gn = g2.norm();
//      if(gn < newton_tol)
//      {
//        DEB_line(4,
//                 "iter = " << iter << ", f(x) = " << fx
//                           << ", |projected grad| = " << gn
//                           << ", #PCG_iters = " << n_pcg_iters
//                           << " -----> converged to desired accuracy!");
//        converged_ = true;
//        break;
//      }
//
//      // stop if line search did not make any progress
//      if(!ls_succeeded)
//      {
//        DEB_line(4,
//                 "iter = " << iter << ", f(x) = " << fx
//                           << ", |grad| = " << gn
//                           << ", #PCG_iters = " << n_pcg_iters
//                           << " -----> stopped because line search could not make any progress!");
//        converged_ = false;
//        break;
//      }
//
//      // set new CG tolerance
//      //      double eta = std::min(0.1,gn);
//      double eta = std::min(0.1,std::sqrt(gn)); // ToDo: is this one suggested in Nocedal better?
//      //      pcg.setTolerance(eta);
//
//      //      pcg.compute(H);
//      //      if(pcg.info() != Eigen::Success)
//      //      {
//      //        std::cerr << "Warning: pcg.compute(H) did not succeed" << std::endl;
//      //      }
//
//
//      // is previous step a descent direction?
//      //      double gdx = g.dot(dx);
//      //      if(gdx < 0.0)
//      //      {
//      //        double s = -gdx/(dx.transpose()*H*dx);
//      //        xt = s*dx;
//      //        dx = pcg.solveWithGuess(-g, xt);
//      //      }
//      //      else
//      //      {
//      //        dx = pcg.solve(-g);
//      //      }
//      //
//      bool cg_converged = false;
//
//      // projeced CG
//      {
//        // precondition QP
//        H2 = D2.asDiagonal() * H * D2.asDiagonal();
//        g2 = D2.asDiagonal()*g;
//
//        dx.setZero();
//        r = H2*dx+g2;
//
//        // norm of pre-conditioned gradient of QP
//        double gn2 = g2.norm();
//
//        double eta = std::min(0.1,std::sqrt(gn2));
//        eta = 1e-12;
//
//        // project r  --> q
//        v = ldlt.solve(A2*r);
//        q = r - A2.transpose()*v;
//        d = -q;
//
//        double rtq = r.transpose()*q;
//
//        for (int i = 0; i < max_pcg_iters; ++i)
//        {
//          Hd = H2*d;
//
//          double alpha = rtq/(d.transpose()*Hd);
//          dx += alpha*d;
//          r1 = r + alpha*Hd;
//
//          // project r1  --> q1
//          v = ldlt.solve(A2*r1);
//          q1 = r1 - A2.transpose()*v;
//
//          //          double beta = (r1.transpose()*q1).eval()/rtq;
//          double beta = r1.transpose()*q1;
//          beta /= rtq;
//          d = -q1 + beta*d;
//
//          q.swap(q1);
//          r.swap(r1);
//
//          // check convergence
//          rtq = r.transpose()*q;
//
//          std::cerr << "CG iter: rtq = " << rtq << ", eps = " << eta*gn2 << std::endl;
//
//          if(rtq < eta*gn2 || i == max_pcg_iters-1)
//          {
//            n_pcg_iters += i+1;
//            // transfrom dx to original coordinates
//            dx = (D2.asDiagonal()*dx).eval();
//            constraint_violation = (_A*(x+dx) -_b).norm();
//            std::cerr << "projected CG terminated with constraint violation " << constraint_violation << " and r^Tq = " << rtq << std::endl;
//            cg_converged = !(i == max_pcg_iters-1);
//            break;
//          }
//        }
//      }
//
//      // get maximal reasonable step
//      double t_max  = std::min(1.0,
//                               max_feasible_step_safety_factor_ * _problem->max_feasible_step(x.data(), dx.data()));
//
//      // backtracking line search
//      double  gdx = g.dot(dx);
//      //      std::cerr << "gdx = " << gdx << std::endl;
//      double t=t_max;
//      xt = x + t*dx;
//      double fxt = _problem->eval_f(xt.data());
//      int iter_ls = 0;
//      ls_succeeded = false;
//      while(! (fxt <= fx + alpha_ls_*gdx*t) && iter_ls < max_iter_ls)
//      {
//        t *= beta_ls_;
//        xt = x + t*dx;
//        fxt = _problem->eval_f(xt.data());
//        ++iter_ls;
//      }
//
//      // update x, fx and fx_old
//      if( iter_ls < max_iter_ls)
//      {
//        x.swap(xt);
//        fx_old = fx;
//        fx = fxt;
//        ls_succeeded = true;
//      }
//
//      DEB_line(4,
//               "iter = " << iter << ", f(x) = " << fx << ", t = " << t
//                         << " (tmax=" << t_max << "), " << "#ls = " << iter_ls
//                         << ", |proj grad| = " << gn
//                         << ", " << "PCG_tol = " << eta
//                         << ", " << "PCG_iters = " << n_pcg_iters
//                         << ", " << "PCG_converged = " << int(cg_converged) );
//    }
//
//    // store result
//    _problem->store_result(x.data());
//
//    // return success
//    return converged_;
//  }


//  // solve with linear constraints
//  // Warning: so far only feasible starting points with (_A*_problem->initial_x() == b) are supported!
//  // It is also required that the constraints are linearly independent
//  int solve_reduced_system( NProblemInterface* _problem, const SMatrixD& _A, const VectorD& _b )
//  {
//    DEB_time_func_def;
//    converged_ = false;
//
//    // number of unknowns
//    Eigen::Index n = _problem->n_unknowns();
//    Eigen::Index m = _A.rows();
//
//    std::cerr << "compute QR decomposition of A^T...\n";
//    SMatrixD At = _A.transpose();
//    At.makeCompressed();
//    StopWatch swq; swq.start();
//    Eigen::SparseQR<SMatrixD,Eigen::COLAMDOrdering<int> > sqr(At);
//    int mq = sqr.rank();
//    // #variables of reduced problem
//    int nz = n - mq;
//    std::cerr << "done! time = " << swq.stop()/1000.0 << "s\n";
//    // get null-space basis Z
//    SMatrixD Q;
//    Q = sqr.matrixQ();
//    SMatrixD Z(n, nz);
//    Z = Q.middleCols(mq,nz); //ToDo: Is the performance of this operation ok?
//    std::cerr << "||A*Z|| = " << (_A*Z).norm() << std::endl;
//
//    DEB_line(2, "optimize via ReducedSystemNewtonPCG with " << (int)n << " unknowns and " << (int)m << " linear constraints (" << (int)mq << " linearly independent)");
//    DEB_line(2, "----> reduced (unconstrained) problem has " << (int)nz << " unknowns");
//
//    // Newton parameters
//    const int    max_iters  = max_iters_;
//    const double newton_tol = eps_; // norm of gradient of reduced (unconstrained) problem
//    // CG parameters
////    const int    max_pcg_iters = n-m;
//    const int    max_pcg_iters = n-m;
//    const double pcg_tol = 1e-4;
//    // allow only steps larger than 1e-8 in line search
//    const int max_iter_ls   = int ( std::log(eps_)/std::log(beta_ls_));
//
//    // initialize vectors of unknowns
//    VectorD x(n);
//    _problem->initial_x(x.data());
//
//    double constraint_violation = (_A*x-_b).norm();
//    DEB_line(2, "initial constraint violation = " << constraint_violation);
//    if(constraint_violation > 1e-6)
//      std::cerr << "Warning: this method expects a feasible starting point!" << std::endl;
//
//    // storage of update vector dx and gradient
//    VectorD g(n), dx(n), xn(n);
//
//    // storage for reduced problem and CG
//    VectorD xz(nz), cz(nz), rz(nz), gz(nz), dz(nz), ZtHzdz(nz), rz2(nz), gz2(nz),cz0(nz);
//    xz.setZero();
//
//    // hessian matrix
//    SMatrixD H(n,n);
//
//    // preconditioner
//    VectorD W(n), Wi(n);
//
//    // statistics
//    int n_pcg_iters(0);
//
//    double eta = pcg_tol;
//
//    bool ls_succeeded = true;
//
//    bool cg_converged = false;
//
//    // get function value at current point
//    double fx = _problem->eval_f(x.data());
//
//    for(int iter=0; iter<max_iters; ++iter)
//    {
//      // get gradient and Hessian
//      _problem->eval_gradient(x.data(), g.data());
//      _problem->eval_hessian(x.data(), H);
//
//      // diagonal preconditioner
//      for (int j = 0; j < n; ++j)
//      {
//        W[j] = std::abs(H.coeffRef(j, j));
//        Wi[j] = 1.0 / W[j];
////
////        // identity preconditioner
////        W[j] = 1.0;
////        Wi[j] = 1.0;
//      }
//
//      //--------------------------
//
//      // solve Hz xz = -gz via CG
//      cg_converged = false;
//      xz.setZero();
//      // project gradient
//      cz = Z.transpose() * g;
//      rz = cz; // Note: rz = Z.transpose()*H*Z*xz if xz \neq 0
//      // apply preconditioner
//      gz = Z.transpose() * (Wi.asDiagonal() * (Z * rz));
//      dz = -gz;
//
//      // check convergence
//      double gzn = cz.norm();
//      if(gzn < newton_tol)
//      {
//        constraint_violation = (_A*x-_b).norm();
//
//        DEB_line(4,"converged" <<  ", f(x) = " << fx
//                               << ", ||gz|| = " << gzn
//                               << ", constraint_violation = " << constraint_violation
//                               << ", " << "PCG_iters = " << n_pcg_iters );
//        break;
//      }
//
//      // choose CG tolerance
//      if(adaptive_tolerance_)
//        eta = std::min(0.1, std::sqrt(gzn))*gzn;
//      else
//        eta = pcg_tol*gzn;
//
//      cz0 = cz; // debug
//
//      for (int iter_cg = 0; iter_cg < max_pcg_iters; ++iter_cg)
//      {
//        double rztgz = rz.transpose() * gz;
//        ZtHzdz = Z.transpose()*(H * (Z * dz));
//        double alpha = rztgz / (dz.transpose() * ZtHzdz);
//
//        xz += alpha * dz;
//
//        rz2 = rz + alpha * ZtHzdz;
//
//        gz2 = Z.transpose() * (Wi.asDiagonal() * (Z * rz2));
//
//        double rz2tgz2 = rz2.transpose() * gz2;
//        double beta = rz2tgz2 / rztgz;
//
//        dz = (-gz2 + beta * dz).eval();
//
//        // update current gz and rz
//        gz.swap(gz2);
//        rz.swap(rz2);
//
//        ++n_pcg_iters;
//
//        // is this correct?
//        if(rz.norm() < eta) break;
//
//        if(0)
//        {
//          // debug output
//          VectorD rr = Z.transpose()*H*Z*xz + cz0;
//          std::cerr << "CG iter=" << iter_cg
//                    << ", ||xz||=" << xz.norm()
//                    << ", ||rz||=" << rz.norm()
//                    << ", ||rr||=" << rr.norm()
//                    << ", ||gz||=" << gz.norm()
//                    << std::endl;
//        }
//      }
//
//      // compute dx of full problem
//      //      dx = (x0 + Z * xz) - x;
//      dx = Z * xz;
//
//      // get maximal reasonable step
//      double t_max = std::min(1.0,
//                              max_feasible_step_safety_factor_ * _problem->max_feasible_step(x.data(), dx.data()));
//
//      // backtracking line search
//      double gdx = g.dot(dx);
//      double t = t_max;
//      xn = x + t*dx;
//      double fxn = _problem->eval_f(xn.data());
//      int iter_ls = 0;
//      ls_succeeded = false;
//      while (!(fxn <= fx + alpha_ls_ * gdx * t) && iter_ls < max_iter_ls) {
//        t *= beta_ls_;
//        xn = x + t * dx;
//        fxn = _problem->eval_f(xn.data());
//        ++iter_ls;
//      }
//
//      // update x, fx, xz
//      if (iter_ls < max_iter_ls)
//      {
//        x.swap(xn);
//        fx = fxn;
//        ls_succeeded = true;
//        // update reduced variables
//        //        xz = Z.transpose()*(x-x0);
//      }
//      else
//      {
//        std::cerr << "Warning: line search failed ---> skip update" << std::endl;
//        break;
//      }
//
//
//      DEB_line(4,
//               "iter = " << iter << ", f(x) = " << fx << ", t = " << t
//                         << " (tmax=" << t_max << "), " << "#ls = " << iter_ls
//                         << ", |reduced grad| = " << gzn
//                         << ", " << "PCG_tol = " << eta
//                         << ", " << "PCG_iters = " << n_pcg_iters
//                         << ", " << "PCG_converged = " << int(cg_converged));
//
//
//    }
//
//    // store result
//    _problem->store_result(x.data());
//
//    // return success
//    return converged_;
//  }


//  // solve with linear constraints
//  // Warning: so far only feasible starting points with (_A*_problem->initial_x() == b) are supported!
//  // It is also required that the constraints are linearly independent
//  int solve_reduced_system_EigenCG( NProblemInterface* _problem, const SMatrixD& _A, const VectorD& _b )
//  {
//    DEB_time_func_def;
//    converged_ = false;
//
//    // number of unknowns
//    Eigen::Index n = _problem->n_unknowns();
//    Eigen::Index m = _A.rows();
//
//    std::cerr << "compute QR decomposition of A^T...\n";
//    StopWatch swq; swq.start();
//    SMatrixD At = _A.transpose();
//    At.makeCompressed();
//    Eigen::SparseQR<SMatrixD,Eigen::COLAMDOrdering<int> > sqr(At);
//    int mq = sqr.rank();
//    // #variables of reduced problem
//    int nz = n - mq;
//    std::cerr << "done! time = " << swq.stop()/1000.0 << "s\n";
//    // get null-space basis Z
//    SMatrixD Q;
//    Q = sqr.matrixQ();
//    SMatrixD Z(n, nz);
//    Z = Q.middleCols(mq,nz); //ToDo: Is the performance of this operation ok?
////    std::cerr << "||A*Z|| = " << (_A*Z).norm() << std::endl;
////    std::cerr << "||Z^T Z||^2 = " << std::pow((Z.transpose()*Z).norm(), 2) << std::endl;
////    std::cerr << "||Z Z^T||^2 = " << std::pow((Z*Z.transpose()).norm(), 2) << std::endl;
//
//    DEB_line(2, "optimize via ReducedSystemNewtonPCGSimple with " << (int)n << " unknowns and " << (int)m << " linear constraints (" << (int)mq << " linearly independent)");
//    DEB_line(2, "----> reduced (unconstrained) problem has " << (int)nz << " unknowns");
//
//    // Newton parameters
//    const int    max_iters  = max_iters_;
//    const double newton_tol = eps_; // norm of gradient of reduced (unconstrained) problem
//    // CG parameters
//    const int    max_pcg_iters = n-m;
//    const double pcg_tol = 1e-4;
//    // allow only steps larger than eps_ls_ in line search
//    const int max_iter_ls   = int ( std::log(eps_ls_)/std::log(beta_ls_));
//
//    // initialize vectors of unknowns
//    VectorD x(n);
//    _problem->initial_x(x.data());
//
//    double constraint_violation = (_A*x-_b).norm();
//    DEB_line(2, "initial constraint violation = " << constraint_violation);
//    if(constraint_violation > 1e-6)
//      std::cerr << "Warning: this method expects a feasible starting point!" << std::endl;
//
//    // storage of update vector dx and gradient
//    VectorD g(n), dx(n), xn(n);
//
//    // storage of update vector dx and gradient
//    VectorD  dxz(nz), gz(nz);
//
//    // hessian matrix
//    SMatrixD H(n,n), Hz(nz,nz);
//
//    // statistics
//    int n_pcg_iters(0);
//
//    // choose precision of CG solve
//    double eta = pcg_tol;
//
//    bool ls_succeeded = true;
//
//    bool cg_converged = false;
//
//    // get function value at current point
//    double fx = _problem->eval_f(x.data());
//
//    for(int iter=0; iter<max_iters; ++iter)
//    {
//      // get gradient and Hessian
//      _problem->eval_gradient(x.data(), g.data());
//      _problem->eval_hessian(x.data(), H);
//
//      // get gradient and Hessian of reduced problem
//      gz = Z.transpose() * g;
//      Hz = Z.transpose()*H*Z;
//
//      // check convergence
//      double gzn = gz.norm();
//      if(gzn < newton_tol)
//      {
//        constraint_violation = (_A*x-_b).norm();
//
//        DEB_line(4,"converged" <<  ", f(x) = " << fx
//                                        << ", ||gz|| = " << gzn
//                                        << ", constraint_violation = " << constraint_violation
//                                        << ", " << "PCG_iters = " << n_pcg_iters );
//        break;
//      }
//
//      // choose CG tolerance
//      if(adaptive_tolerance_)
//        eta = std::min(0.1, std::sqrt(gzn));
//      else
//        eta = pcg_tol;
//
//      Eigen::ConjugateGradient<SMatrixD, Eigen::Lower|Eigen::Upper> pcg;
////      Eigen::ConjugateGradient<SMatrixD, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double> > pcg;
//      pcg.setMaxIterations(max_pcg_iters);
//      pcg.setTolerance(eta);
//      pcg.compute(Hz);
//      dxz=pcg.solve(-gz);
//      n_pcg_iters += pcg.iterations();
//      cg_converged = (pcg.info() == Eigen::Success);
//
//      if(0)
//      { // debug
//        std::cerr << "pcg residual=" << (Z.transpose() * H * Z * dxz + gz).norm() / gz.norm()
//                  << std::endl;
//      }
//
//      // compute dx of full problem
//      dx = Z * dxz;
//
//      // get maximal reasonable step
//      double t_max = std::min(1.0,
//                              max_feasible_step_safety_factor_ * _problem->max_feasible_step(x.data(), dx.data()));
//
//      // backtracking line search
//      double gdx = g.dot(dx);
//      double t = t_max;
//      xn = x + t*dx;
//      double fxn = _problem->eval_f(xn.data());
//      int iter_ls = 0;
//      ls_succeeded = false;
//      while (!(fxn <= fx + alpha_ls_ * gdx * t) && iter_ls < max_iter_ls) {
//        t *= beta_ls_;
//        xn = x + t * dx;
//        fxn = _problem->eval_f(xn.data());
//        ++iter_ls;
//      }
//
//      // update x, fx, xz
//      if (iter_ls < max_iter_ls)
//      {
//        x.swap(xn);
//        fx = fxn;
//        ls_succeeded = true;
//        // update reduced variables
//        //        xz = Z.transpose()*(x-x0);
//      }
//      else
//      {
//        std::cerr << "Warning: line search failed ---> terminate" << std::endl;
//        break;
//      }
//
//
//      DEB_line(4,
//               "iter = " << iter << ", f(x) = " << fx << ", t = " << t
//                         << " (tmax=" << t_max << "), " << "#ls = " << iter_ls
//                         << ", |reduced grad| = " << gzn
//                         << ", " << "PCG_tol = " << eta
//                         << ", " << "PCG_iters = " << n_pcg_iters
//                         << ", " << "PCG_converged = " << int(cg_converged));
//
//
//    }
//
//    // store result
//    _problem->store_result(x.data());
//
//    // return success
//    return converged_;
//  }


//  int solve_experimental(NProblemInterface* _problem, std::vector<LinearConstraint>& _constraints)
//  {
//    // convert constraints
//    SMatrixD A;
//    VectorD b;
//    LinearConstraintConverter::nsolver_to_eigen(_constraints, A, b, _problem->n_unknowns());
//
//    return solve_experimental(_problem, A, b);
//  }
//
//  int solve_reduced_system(NProblemInterface* _problem, std::vector<LinearConstraint>& _constraints)
//  {
//    // convert constraints
//    SMatrixD A;
//    VectorD b;
//    LinearConstraintConverter::nsolver_to_eigen(_constraints, A, b, _problem->n_unknowns());
//
//    return solve_reduced_system(_problem, A, b);
//  }
//
//  int solve_reduced_system_EigenCG(NProblemInterface* _problem, std::vector<LinearConstraint>& _constraints)
//  {
//    // convert constraints
//    SMatrixD A;
//    VectorD b;
//    LinearConstraintConverter::nsolver_to_eigen(_constraints, A, b, _problem->n_unknowns());
//
//    return solve_reduced_system_EigenCG(_problem, A, b);
//  }


}
