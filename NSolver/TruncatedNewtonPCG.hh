
/*===========================================================================*\
 *                                                                           *
 *                            TruncatedNewtonPCG                             *
 *      Copyright (C) 2024 by Computer Graphics Group, University of Bern    *
 *                           http://cgg.unibe.ch                             *
 *                                                                           *
 *      Author: David Bommes                                                 *
 *                                                                           *
\*===========================================================================*/



#ifndef COMISO_TRUNCATEDNEWTONPCG_HH
#define COMISO_TRUNCATEDNEWTONPCG_HH

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>

//== INCLUDES =================================================================

#include <iomanip>

#include <CoMISo/Config/CoMISoDefines.hh>
#include <CoMISo/Utils/StopWatch.hh>
#include <CoMISo/NSolver/NProblemInterface.hh>
#include <CoMISo/NSolver/NConstraintInterface.hh>
#include <CoMISo/NSolver/LinearConstraintConverter.hh>

#include <Base/Debug/DebTime.hh>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

#include <Eigen/CholmodSupport>
#include <Eigen/MetisSupport>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================

	      

/** \class NewtonSolver NewtonSolver.hh <CoMISo/.../TruncatedNewtonPCG.hh>

    Brief Description.
  
    A more elaborate description follows.
*/
class COMISODLLEXPORT TruncatedNewtonPCG
{
public:

  typedef Eigen::VectorXd             VectorD;
  typedef Eigen::SparseMatrix<double> SMatrixD;
  typedef Eigen::Triplet<double>      Triplet;

  /// Default constructor
  TruncatedNewtonPCG(const double _eps = 1e-3, const double _eps_line_search = 1e-8,
      const int _max_iters = 500, const double _alpha_ls = 0.2,
      const double _beta_ls = 0.6)
      : eps_(_eps), eps_ls_(_eps_line_search), eps_gdx_(0.1), max_iters_(_max_iters), max_pcg_iters_(500),
        alpha_ls_(_alpha_ls), beta_ls_(_beta_ls), max_feasible_step_safety_factor_(0.5),
        adaptive_tolerance_(true), always_update_preconditioner_(true), allow_warmstart_(false), adaptive_tolerance_modifier_(1.0), converged_(false), silent_(false), reduced_gradient_norm_(0.0)
  {
  }

  int solve(NProblemInterface* _problem)
  {
    DEB_enter_func;
//    DEB_time_func_def;
    converged_ = false;

    // number of unknowns
    Eigen::Index n = _problem->n_unknowns();
    DEB_line_if(!silent_, 2, "optimize via TruncatedNewtonPCG with " << (int)n << " unknowns");

    // Newton parameters
    const int    max_iters     = max_iters_;
    const double newton_tol = eps_;
//    const double relative_improvement_thres = 1e-3;
    // CG parameters
//    const int max_pcg_iters = std::max(100, _problem->n_unknowns()/10);
    const int max_pcg_iters = max_pcg_iters_;

    // allow minimal steps of 1e-6
    const int max_iter_ls   = int ( std::log(1e-6)/std::log(beta_ls_));

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
    pcg.setTolerance(1e-6);

    // statistics
    int n_pcg_iters(0);

    bool ls_succeeded = true;

    // get function value at current point
    double fx = _problem->eval_f(x.data());
    double fx_old = fx;

    DEB_line_if(!silent_, 2, "intial objective value =" << fx);
    if(!std::isfinite(fx))
    {
      DEB_line_if(!silent_, 2, "ERROR: intial objective value is not finite ---> abort");
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
         count_hessian_skip >= hessian_max_skips_ || // maximum number of Hessians have been skipped
         alpha_line_search < hessian_min_acceptable_alpha_ || // line search truncated too much
         rel_objective_decrease < hessian_min_acceptable_rel_objective_decrease_ // objective did not decrease sufficiently
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
      if(gn < newton_tol)
      {
        DEB_line_if(!silent_, 4,
                 "iter = " << iter << ", f(x) = " << fx
                           << ", |grad| = " << gn
                           << ", #PCG_iters = " << n_pcg_iters
                           << " -----> converged to desired accuracy!");
        converged_ = true;
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

      // stop if line search did not make any progress
      if(!ls_succeeded)
      {
        DEB_line_if(!silent_, 4,
                 "iter = " << iter << ", f(x) = " << fx
                           << ", |grad| = " << gn
                           << ", #PCG_iters = " << n_pcg_iters
                           << " -----> stopped because line search could not make any progress!");
        converged_ = false;
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
      if(gdx < 0.0 && allow_warmstart_)
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

      bool cg_converged = (pcg.info() == Eigen::Success);

      // update statistics
      n_pcg_iters += pcg.iterations();

      // get maximal reasonable step
      double t_max  = std::min(1.0,
                               max_feasible_step_safety_factor_ * _problem->max_feasible_step(x.data(), dx.data()));

      // backtracking line search
      gdx = g.dot(dx);
      double t=t_max;
      xt = x + t*dx;
      double fxt = _problem->eval_f(xt.data());
      int iter_ls = 0;
      ls_succeeded = false;
      while(! (fxt <= fx + alpha_ls_*gdx*t) && iter_ls < max_iter_ls)
      {
        t *= beta_ls_;
        xt = x + t*dx;
        fxt = _problem->eval_f(xt.data());
        ++iter_ls;
      }

      // store line search truncation
      alpha_line_search = t;
      // store relative reduction
      rel_objective_decrease = std::abs((fx-fxt)/fx);

      // update x, fx and fx_old
      if( iter_ls < max_iter_ls)
      {
        x.swap(xt);
        fx_old = fx;
        fx = fxt;
        ls_succeeded = true;
      }

      if(std::abs(gdx) < eps_gdx_)
      {
        DEB_line_if(!silent_, 4,
                 "iter = " << iter << ", f(x) = " << fx
                           << ", |grad| = " << gn
                           << ", gdx = " << gdx
                           << ", #PCG_iters = " << n_pcg_iters
                           << " -----> converged to desired accuracy!");
        converged_ = true;
        break;
      }


      DEB_line_if(!silent_, 4,
               "iter = " << iter << ", f(x) = " << fx << ", t = " << t
                         << " (tmax=" << t_max << "), " << "#ls = " << iter_ls
                         << ", |grad| = " << gn
                         << ", " << "PCG_tol = " << eta
                         << ", " << "PCG_iters = " << int(pcg.iterations())
                         << " (total = " << n_pcg_iters << ") "
                         << ", " << "PCG_converged = " << int(cg_converged)
                         << ", " << "PCG_warmstart = " << int(warm_start)
                         << ", " << "gdx = " << gdx
                         << ", " << "hessian_update = " << int(hessian_updated)
                         << ", " << "rel_obj_decrease = " << rel_objective_decrease
                         );
    }

    // store result
    _problem->store_result(x.data());

    // return success
    return converged_;
  }

  // solve with linear constraints
  // Warning: so far only feasible starting points with (_A*_problem->initial_x() == b) are supported!
  // It is also required that the constraints are linearly independent
  int solve(NProblemInterface* _problem, const SMatrixD& _A, const VectorD& _b)
  {
    DEB_time_func_def;
    converged_ = false;

    // number of unknowns
    Eigen::Index n = _problem->n_unknowns();
    Eigen::Index m = _A.rows();
    DEB_line(2, "optimize via TruncatedNewtonPCG with " << (int)n << " unknowns and " << (int)m << " linear constraints");

    // Newton parameters
    const int    max_iters     = max_iters_;
    const double newton_tol = eps_;
    //    const double relative_improvement_thres = 1e-3;
    // CG parameters
    const int max_pcg_iters = 1000;

    // allow minimal steps of 1e-6
    const int max_iter_ls   = int ( std::log(eps_ls_)/std::log(beta_ls_));

    // initialize vectors of unknowns
    VectorD x(n);
    _problem->initial_x(x.data());

    // storage of update vector dx and gradient
    VectorD dx(n), g(n), xt(n);
    dx.setZero();

    // hessian matrix
    SMatrixD H(n,n);

    // storage for CG
    VectorD r(n), q(n), v(m), d(n), Hd(n), r1(n), q1(n);

    // prepare constraint projection
    DEB_line(2, "prepare constraint projection ..." );
    double constraint_violation = (_A*x-_b).norm();
    DEB_line(2, "initial constraint violation = " << constraint_violation);
    SMatrixD AAt = _A*_A.transpose();
    Eigen::SimplicialLDLT<SMatrixD> ldlt(AAt);
    if(ldlt.info() != Eigen::Success)
    {
      for(unsigned int j=0; j<10; ++j)
      {
        double reg = 1e-8 * AAt.diagonal().sum() / double(m);
        DEB_line(2,
                 "Warning: LDLT factorization failed --> regularize AAt (could lead to reduced accuracy), reg=" << reg);
        for (Eigen::Index j = 0; j < AAt.rows(); ++j)
          AAt.coeffRef(j, j) += reg;

        ldlt.compute(AAt);
        if(ldlt.info() == Eigen::Success)
          break;
      }

    }
    DEB_line(2, "done!" );

    // statistics
    int n_pcg_iters(0);

    bool ls_succeeded = true;

    // get function value at current point
    double fx = _problem->eval_f(x.data());
    double fx_old = fx;

    for(int iter=0; iter<max_iters; ++iter)
    {
      // get gradient and Hessian
      _problem->eval_gradient(x.data(), g.data());
      _problem->eval_hessian(x.data(), H);

      // stopping criterion: gradient norm
      double gn = g.norm();
      if(gn < newton_tol)
      {
        DEB_line(4,
                 "iter = " << iter << ", f(x) = " << fx
                           << ", |grad| = " << gn
                           << ", #PCG_iters = " << n_pcg_iters
                           << " -----> converged to desired accuracy!");
        converged_ = true;
        break;
      }

      // stop if line search did not make any progress
      if(!ls_succeeded)
      {
        DEB_line(4,
                 "iter = " << iter << ", f(x) = " << fx
                           << ", |grad| = " << gn
                           << ", #PCG_iters = " << n_pcg_iters
                           << " -----> stopped because line search could not make any progress!");
        converged_ = false;
        break;
      }

      // set new CG tolerance
      //      double eta = std::min(0.1,gn);
      double eta = std::min(0.1,std::sqrt(gn)); // ToDo: is this one suggested in Nocedal better?
      //      pcg.setTolerance(eta);

      //      pcg.compute(H);
      //      if(pcg.info() != Eigen::Success)
      //      {
      //        std::cerr << "Warning: pcg.compute(H) did not succeed" << std::endl;
      //      }


      // is previous step a descent direction?
      //      double gdx = g.dot(dx);
      //      if(gdx < 0.0)
      //      {
      //        double s = -gdx/(dx.transpose()*H*dx);
      //        xt = s*dx;
      //        dx = pcg.solveWithGuess(-g, xt);
      //      }
      //      else
      //      {
      //        dx = pcg.solve(-g);
      //      }
      //
      bool cg_converged = false;

      // projeced CG
      {
        dx.setZero();
        r = H*dx+g;

        // project r  --> q
        v = ldlt.solve(_A*r);
        q = r - _A.transpose()*v;
        d = -q;

        double rtq = r.transpose()*q;
        std::cerr << "rtq = " << rtq << std::endl;

        for (int i = 0; i < max_pcg_iters; ++i)
        {
          Hd = H*d;

          double alpha = rtq/(d.transpose()*Hd);
          dx += alpha*d;
          r1 = r + alpha*Hd;

          // project r1  --> q1
          v = ldlt.solve(_A*r1);
          q1 = r1 - _A.transpose()*v;

          //          double beta = (r1.transpose()*q1).eval()/rtq;
          double beta = r1.transpose()*q1;
          beta /= rtq;
          d = -q1 + beta*d;

          q.swap(q1);
          r.swap(r1);

          // check convergence
          rtq = r.transpose()*q;

          if(rtq < eta || i == max_pcg_iters-1)
          {
            n_pcg_iters += i+1;
            constraint_violation = (_A*(x+dx) -_b).norm();
            std::cerr << "projected CG terminated with constraint violation " << constraint_violation << std::endl;
            cg_converged = !(i == max_pcg_iters-1);
            break;
          }
        }
      }

      // get maximal reasonable step
      double t_max  = std::min(1.0,
                               max_feasible_step_safety_factor_ * _problem->max_feasible_step(x.data(), dx.data()));

      // backtracking line search
      double  gdx = g.dot(dx);
      std::cerr << "gdx = " << gdx << std::endl;
      double t=t_max;
      xt = x + t*dx;
      double fxt = _problem->eval_f(xt.data());
      int iter_ls = 0;
      ls_succeeded = false;
      while(! (fxt <= fx + alpha_ls_*gdx*t) && iter_ls < max_iter_ls)
      {
        t *= beta_ls_;
        xt = x + t*dx;
        fxt = _problem->eval_f(xt.data());
        ++iter_ls;
      }

      // update x, fx and fx_old
      if( iter_ls < max_iter_ls)
      {
        x.swap(xt);
        fx_old = fx;
        fx = fxt;
        ls_succeeded = true;
      }

      DEB_line(4,
               "iter = " << iter << ", f(x) = " << fx << ", t = " << t
                         << " (tmax=" << t_max << "), " << "#ls = " << iter_ls
                         << ", |grad| = " << gn
                         << ", " << "PCG_tol = " << eta
                         << ", " << "PCG_iters = " << n_pcg_iters
                         << ", " << "PCG_converged = " << int(cg_converged) );
    }

    // store result
    _problem->store_result(x.data());

    // return success
    return converged_;
  }


  // solve with linear constraints
  // Warning: so far only feasible starting points with (_A*_problem->initial_x() == b) are supported!
  // It is also required that the constraints are linearly independent
  int solve_experimental(NProblemInterface* _problem, const SMatrixD& _A, const VectorD& _b)
  {
    DEB_time_func_def;
    converged_ = false;

    {
      std::cerr << "compute QR decomposition of A^T...\n";
      SMatrixD At = _A.transpose();
      At.makeCompressed();
      StopWatch swq; swq.start();
      Eigen::SparseQR<SMatrixD,Eigen::COLAMDOrdering<int> > sqr(At);
      std::cerr << "done! time = " << swq.stop()/1000.0 << "s\n";
    }

    // number of unknowns
    Eigen::Index n = _problem->n_unknowns();
    Eigen::Index m = _A.rows();
    DEB_line(2, "optimize via TruncatedNewtonPCG with " << (int)n << " unknowns and " << (int)m << " linear constraints");

    // Newton parameters
    const int    max_iters     = max_iters_;
    const double newton_tol = eps_;
    //    const double relative_improvement_thres = 1e-3;
    // CG parameters
//    const int max_pcg_iters = 2*(n-m);
    const int max_pcg_iters = 500;

    // allow minimal steps of eps_ls_
    const int max_iter_ls   = int ( std::log(eps_ls_)/std::log(beta_ls_));

    // initialize vectors of unknowns
    VectorD x(n);
    _problem->initial_x(x.data());

    // storage of update vector dx and gradient
    VectorD dx(n), g(n), g2(n), xt(n);
    dx.setZero();

    // hessian matrix
    SMatrixD H(n,n), H2(n,n), A2(m,n);

    // preconditioner
    VectorD D(n), D2(n);

    // storage for CG
    VectorD r(n), q(n), v(m), d(n), Hd(n), r1(n), q1(n);

    // cholesky for constraint projection
    Eigen::SimplicialLDLT<SMatrixD> ldlt;

    double constraint_violation(0.0);

    // statistics
    int n_pcg_iters(0);

    bool ls_succeeded = true;

    // get function value at current point
    double fx = _problem->eval_f(x.data());
    double fx_old = fx;

    for(int iter=0; iter<max_iters; ++iter)
    {
      // get gradient and Hessian
      _problem->eval_gradient(x.data(), g.data());
      _problem->eval_hessian(x.data(), H);

      if(iter == 0)
      {
        // setup preconditioner D2 = (|diag(H)|)^(-1/2)
        // and D = (|diag(H)|)^(1/2)
        for(int j=0; j<n ;++j)
        {
          D[j]  = std::sqrt(std::abs(H.coeffRef(j,j)));
          D2[j] = 1.0/D[j];
        }

        // check
        if(0)
        {
          H2 = D2.asDiagonal() * H * D2.asDiagonal();
          std::cerr << "************ Hessian diagonal: " << H.diagonal() << std::endl;
          std::cerr << "************ Preconditioned Hessian diagonal: " << H2.diagonal() << std::endl;
        }

        // prepare constraint projection
        DEB_line(2, "prepare constraint projection ..." );
        constraint_violation = (_A*x-_b).norm();
        DEB_line(2, "initial constraint violation = " << constraint_violation);
        A2 = _A*D2.asDiagonal();
        SMatrixD AAt = A2*A2.transpose();
        ldlt.compute(AAt);
        if(ldlt.info() != Eigen::Success)
        {
          for(unsigned int j=0; j<10; ++j)
          {
            double reg = 1e-8 * AAt.diagonal().sum() / double(m);
            DEB_line(2,
                     "Warning: LDLT factorization failed --> regularize AAt (could lead to reduced accuracy), reg=" << reg);
            for (Eigen::Index j = 0; j < AAt.rows(); ++j)
              AAt.coeffRef(j, j) += reg;

            ldlt.compute(AAt);
            if(ldlt.info() == Eigen::Success)
              break;
          }

        }
        DEB_line(2, "done!" );

      }

      // stopping criterion: gradient norm of eliminated problem (or equivalently norm of projected gradient)
      // project g  --> g2
      g2 = D2.asDiagonal()*g;

      v = ldlt.solve(A2*g2);
      g2 -= A2.transpose()*v;
      g2 = D.asDiagonal()*g2;

      double gn = g2.norm();
      if(gn < newton_tol)
      {
        DEB_line(4,
                 "iter = " << iter << ", f(x) = " << fx
                           << ", |projected grad| = " << gn
                           << ", #PCG_iters = " << n_pcg_iters
                           << " -----> converged to desired accuracy!");
        converged_ = true;
        break;
      }

      // stop if line search did not make any progress
      if(!ls_succeeded)
      {
        DEB_line(4,
                 "iter = " << iter << ", f(x) = " << fx
                           << ", |grad| = " << gn
                           << ", #PCG_iters = " << n_pcg_iters
                           << " -----> stopped because line search could not make any progress!");
        converged_ = false;
        break;
      }

      // set new CG tolerance
      //      double eta = std::min(0.1,gn);
      double eta = std::min(0.1,std::sqrt(gn)); // ToDo: is this one suggested in Nocedal better?
      //      pcg.setTolerance(eta);

      //      pcg.compute(H);
      //      if(pcg.info() != Eigen::Success)
      //      {
      //        std::cerr << "Warning: pcg.compute(H) did not succeed" << std::endl;
      //      }


      // is previous step a descent direction?
      //      double gdx = g.dot(dx);
      //      if(gdx < 0.0)
      //      {
      //        double s = -gdx/(dx.transpose()*H*dx);
      //        xt = s*dx;
      //        dx = pcg.solveWithGuess(-g, xt);
      //      }
      //      else
      //      {
      //        dx = pcg.solve(-g);
      //      }
      //
      bool cg_converged = false;

      // projeced CG
      {
        // precondition QP
        H2 = D2.asDiagonal() * H * D2.asDiagonal();
        g2 = D2.asDiagonal()*g;

        dx.setZero();
        r = H2*dx+g2;

        // norm of pre-conditioned gradient of QP
        double gn2 = g2.norm();

        double eta = std::min(0.1,std::sqrt(gn2));
        eta = 1e-12;

        // project r  --> q
        v = ldlt.solve(A2*r);
        q = r - A2.transpose()*v;
        d = -q;

        double rtq = r.transpose()*q;

        for (int i = 0; i < max_pcg_iters; ++i)
        {
          Hd = H2*d;

          double alpha = rtq/(d.transpose()*Hd);
          dx += alpha*d;
          r1 = r + alpha*Hd;

          // project r1  --> q1
          v = ldlt.solve(A2*r1);
          q1 = r1 - A2.transpose()*v;

          //          double beta = (r1.transpose()*q1).eval()/rtq;
          double beta = r1.transpose()*q1;
          beta /= rtq;
          d = -q1 + beta*d;

          q.swap(q1);
          r.swap(r1);

          // check convergence
          rtq = r.transpose()*q;

          std::cerr << "CG iter: rtq = " << rtq << ", eps = " << eta*gn2 << std::endl;

          if(rtq < eta*gn2 || i == max_pcg_iters-1)
          {
            n_pcg_iters += i+1;
            // transfrom dx to original coordinates
            dx = (D2.asDiagonal()*dx).eval();
            constraint_violation = (_A*(x+dx) -_b).norm();
            std::cerr << "projected CG terminated with constraint violation " << constraint_violation << " and r^Tq = " << rtq << std::endl;
            cg_converged = !(i == max_pcg_iters-1);
            break;
          }
        }
      }

      // get maximal reasonable step
      double t_max  = std::min(1.0,
                               max_feasible_step_safety_factor_ * _problem->max_feasible_step(x.data(), dx.data()));

      // backtracking line search
      double  gdx = g.dot(dx);
      //      std::cerr << "gdx = " << gdx << std::endl;
      double t=t_max;
      xt = x + t*dx;
      double fxt = _problem->eval_f(xt.data());
      int iter_ls = 0;
      ls_succeeded = false;
      while(! (fxt <= fx + alpha_ls_*gdx*t) && iter_ls < max_iter_ls)
      {
        t *= beta_ls_;
        xt = x + t*dx;
        fxt = _problem->eval_f(xt.data());
        ++iter_ls;
      }

      // update x, fx and fx_old
      if( iter_ls < max_iter_ls)
      {
        x.swap(xt);
        fx_old = fx;
        fx = fxt;
        ls_succeeded = true;
      }

      DEB_line(4,
               "iter = " << iter << ", f(x) = " << fx << ", t = " << t
                         << " (tmax=" << t_max << "), " << "#ls = " << iter_ls
                         << ", |proj grad| = " << gn
                         << ", " << "PCG_tol = " << eta
                         << ", " << "PCG_iters = " << n_pcg_iters
                         << ", " << "PCG_converged = " << int(cg_converged) );
    }

    // store result
    _problem->store_result(x.data());

    // return success
    return converged_;
  }

  // solve with linear constraints
  // Warning: so far only feasible starting points with (_A*_problem->initial_x() == b) are supported!
  // It is also required that the constraints are linearly independent
  int solve_reduced_system( NProblemInterface* _problem, const SMatrixD& _A, const VectorD& _b )
  {
    DEB_time_func_def;
    converged_ = false;

    // number of unknowns
    Eigen::Index n = _problem->n_unknowns();
    Eigen::Index m = _A.rows();

    std::cerr << "compute QR decomposition of A^T...\n";
    SMatrixD At = _A.transpose();
    At.makeCompressed();
    StopWatch swq; swq.start();
    Eigen::SparseQR<SMatrixD,Eigen::COLAMDOrdering<int> > sqr(At);
    int mq = sqr.rank();
    // #variables of reduced problem
    int nz = n - mq;
    std::cerr << "done! time = " << swq.stop()/1000.0 << "s\n";
    // get null-space basis Z
    SMatrixD Q;
    Q = sqr.matrixQ();
    SMatrixD Z(n, nz);
    Z = Q.middleCols(mq,nz); //ToDo: Is the performance of this operation ok?
    std::cerr << "||A*Z|| = " << (_A*Z).norm() << std::endl;

    DEB_line(2, "optimize via ReducedSystemNewtonPCG with " << (int)n << " unknowns and " << (int)m << " linear constraints (" << (int)mq << " linearly independent)");
    DEB_line(2, "----> reduced (unconstrained) problem has " << (int)nz << " unknowns");

    // Newton parameters
    const int    max_iters  = max_iters_;
    const double newton_tol = eps_; // norm of gradient of reduced (unconstrained) problem
    // CG parameters
//    const int    max_pcg_iters = n-m;
    const int    max_pcg_iters = n-m;
    const double pcg_tol = 1e-4;
    // allow only steps larger than 1e-8 in line search
    const int max_iter_ls   = int ( std::log(eps_)/std::log(beta_ls_));

    // initialize vectors of unknowns
    VectorD x(n);
    _problem->initial_x(x.data());

    double constraint_violation = (_A*x-_b).norm();
    DEB_line(2, "initial constraint violation = " << constraint_violation);
    if(constraint_violation > 1e-6)
      std::cerr << "Warning: this method expects a feasible starting point!" << std::endl;

    // storage of update vector dx and gradient
    VectorD g(n), dx(n), xn(n);

    // storage for reduced problem and CG
    VectorD xz(nz), cz(nz), rz(nz), gz(nz), dz(nz), ZtHzdz(nz), rz2(nz), gz2(nz),cz0(nz);
    xz.setZero();

    // hessian matrix
    SMatrixD H(n,n);

    // preconditioner
    VectorD W(n), Wi(n);

    // statistics
    int n_pcg_iters(0);

    double eta = pcg_tol;

    bool ls_succeeded = true;

    bool cg_converged = false;

    // get function value at current point
    double fx = _problem->eval_f(x.data());

    for(int iter=0; iter<max_iters; ++iter)
    {
      // get gradient and Hessian
      _problem->eval_gradient(x.data(), g.data());
      _problem->eval_hessian(x.data(), H);

      // diagonal preconditioner
      for (int j = 0; j < n; ++j)
      {
        W[j] = std::abs(H.coeffRef(j, j));
        Wi[j] = 1.0 / W[j];
//
//        // identity preconditioner
//        W[j] = 1.0;
//        Wi[j] = 1.0;
      }

      //--------------------------

      // solve Hz xz = -gz via CG
      cg_converged = false;
      xz.setZero();
      // project gradient
      cz = Z.transpose() * g;
      rz = cz; // Note: rz = Z.transpose()*H*Z*xz if xz \neq 0
      // apply preconditioner
      gz = Z.transpose() * (Wi.asDiagonal() * (Z * rz));
      dz = -gz;

      // check convergence
      double gzn = cz.norm();
      if(gzn < newton_tol)
      {
        constraint_violation = (_A*x-_b).norm();

        DEB_line(4,"converged" <<  ", f(x) = " << fx
                               << ", ||gz|| = " << gzn
                               << ", constraint_violation = " << constraint_violation
                               << ", " << "PCG_iters = " << n_pcg_iters );
        break;
      }

      // choose CG tolerance
      if(adaptive_tolerance_)
        eta = std::min(0.1, std::sqrt(gzn))*gzn;
      else
        eta = pcg_tol*gzn;

      cz0 = cz; // debug

      for (int iter_cg = 0; iter_cg < max_pcg_iters; ++iter_cg)
      {
        double rztgz = rz.transpose() * gz;
        ZtHzdz = Z.transpose()*(H * (Z * dz));
        double alpha = rztgz / (dz.transpose() * ZtHzdz);

        xz += alpha * dz;

        rz2 = rz + alpha * ZtHzdz;

        gz2 = Z.transpose() * (Wi.asDiagonal() * (Z * rz2));

        double rz2tgz2 = rz2.transpose() * gz2;
        double beta = rz2tgz2 / rztgz;

        dz = (-gz2 + beta * dz).eval();

        // update current gz and rz
        gz.swap(gz2);
        rz.swap(rz2);

        ++n_pcg_iters;

        // is this correct?
        if(rz.norm() < eta) break;

        if(0)
        {
          // debug output
          VectorD rr = Z.transpose()*H*Z*xz + cz0;
          std::cerr << "CG iter=" << iter_cg
                    << ", ||xz||=" << xz.norm()
                    << ", ||rz||=" << rz.norm()
                    << ", ||rr||=" << rr.norm()
                    << ", ||gz||=" << gz.norm()
                    << std::endl;
        }
      }

      // compute dx of full problem
      //      dx = (x0 + Z * xz) - x;
      dx = Z * xz;

      // get maximal reasonable step
      double t_max = std::min(1.0,
                              max_feasible_step_safety_factor_ * _problem->max_feasible_step(x.data(), dx.data()));

      // backtracking line search
      double gdx = g.dot(dx);
      double t = t_max;
      xn = x + t*dx;
      double fxn = _problem->eval_f(xn.data());
      int iter_ls = 0;
      ls_succeeded = false;
      while (!(fxn <= fx + alpha_ls_ * gdx * t) && iter_ls < max_iter_ls) {
        t *= beta_ls_;
        xn = x + t * dx;
        fxn = _problem->eval_f(xn.data());
        ++iter_ls;
      }

      // update x, fx, xz
      if (iter_ls < max_iter_ls)
      {
        x.swap(xn);
        fx = fxn;
        ls_succeeded = true;
        // update reduced variables
        //        xz = Z.transpose()*(x-x0);
      }
      else
      {
        std::cerr << "Warning: line search failed ---> skip update" << std::endl;
        break;
      }


      DEB_line(4,
               "iter = " << iter << ", f(x) = " << fx << ", t = " << t
                         << " (tmax=" << t_max << "), " << "#ls = " << iter_ls
                         << ", |reduced grad| = " << gzn
                         << ", " << "PCG_tol = " << eta
                         << ", " << "PCG_iters = " << n_pcg_iters
                         << ", " << "PCG_converged = " << int(cg_converged));


    }

    // store result
    _problem->store_result(x.data());

    // return success
    return converged_;
  }


  // solve with linear constraints
  // Warning: so far only feasible starting points with (_A*_problem->initial_x() == b) are supported!
  // It is also required that the constraints are linearly independent
  int solve_reduced_system_EigenCG( NProblemInterface* _problem, const SMatrixD& _A, const VectorD& _b )
  {
    DEB_time_func_def;
    converged_ = false;

    // number of unknowns
    Eigen::Index n = _problem->n_unknowns();
    Eigen::Index m = _A.rows();

    std::cerr << "compute QR decomposition of A^T...\n";
    StopWatch swq; swq.start();
    SMatrixD At = _A.transpose();
    At.makeCompressed();
    Eigen::SparseQR<SMatrixD,Eigen::COLAMDOrdering<int> > sqr(At);
    int mq = sqr.rank();
    // #variables of reduced problem
    int nz = n - mq;
    std::cerr << "done! time = " << swq.stop()/1000.0 << "s\n";
    // get null-space basis Z
    SMatrixD Q;
    Q = sqr.matrixQ();
    SMatrixD Z(n, nz);
    Z = Q.middleCols(mq,nz); //ToDo: Is the performance of this operation ok?
//    std::cerr << "||A*Z|| = " << (_A*Z).norm() << std::endl;
//    std::cerr << "||Z^T Z||^2 = " << std::pow((Z.transpose()*Z).norm(), 2) << std::endl;
//    std::cerr << "||Z Z^T||^2 = " << std::pow((Z*Z.transpose()).norm(), 2) << std::endl;

    DEB_line(2, "optimize via ReducedSystemNewtonPCGSimple with " << (int)n << " unknowns and " << (int)m << " linear constraints (" << (int)mq << " linearly independent)");
    DEB_line(2, "----> reduced (unconstrained) problem has " << (int)nz << " unknowns");

    // Newton parameters
    const int    max_iters  = max_iters_;
    const double newton_tol = eps_; // norm of gradient of reduced (unconstrained) problem
    // CG parameters
    const int    max_pcg_iters = n-m;
    const double pcg_tol = 1e-4;
    // allow only steps larger than eps_ls_ in line search
    const int max_iter_ls   = int ( std::log(eps_ls_)/std::log(beta_ls_));

    // initialize vectors of unknowns
    VectorD x(n);
    _problem->initial_x(x.data());

    double constraint_violation = (_A*x-_b).norm();
    DEB_line(2, "initial constraint violation = " << constraint_violation);
    if(constraint_violation > 1e-6)
      std::cerr << "Warning: this method expects a feasible starting point!" << std::endl;

    // storage of update vector dx and gradient
    VectorD g(n), dx(n), xn(n);

    // storage of update vector dx and gradient
    VectorD  dxz(nz), gz(nz);

    // hessian matrix
    SMatrixD H(n,n), Hz(nz,nz);

    // statistics
    int n_pcg_iters(0);

    // choose precision of CG solve
    double eta = pcg_tol;

    bool ls_succeeded = true;

    bool cg_converged = false;

    // get function value at current point
    double fx = _problem->eval_f(x.data());

    for(int iter=0; iter<max_iters; ++iter)
    {
      // get gradient and Hessian
      _problem->eval_gradient(x.data(), g.data());
      _problem->eval_hessian(x.data(), H);

      // get gradient and Hessian of reduced problem
      gz = Z.transpose() * g;
      Hz = Z.transpose()*H*Z;

      // check convergence
      double gzn = gz.norm();
      if(gzn < newton_tol)
      {
        constraint_violation = (_A*x-_b).norm();

        DEB_line(4,"converged" <<  ", f(x) = " << fx
                                        << ", ||gz|| = " << gzn
                                        << ", constraint_violation = " << constraint_violation
                                        << ", " << "PCG_iters = " << n_pcg_iters );
        break;
      }

      // choose CG tolerance
      if(adaptive_tolerance_)
        eta = std::min(0.1, std::sqrt(gzn));
      else
        eta = pcg_tol;

      Eigen::ConjugateGradient<SMatrixD, Eigen::Lower|Eigen::Upper> pcg;
//      Eigen::ConjugateGradient<SMatrixD, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double> > pcg;
      pcg.setMaxIterations(max_pcg_iters);
      pcg.setTolerance(eta);
      pcg.compute(Hz);
      dxz=pcg.solve(-gz);
      n_pcg_iters += pcg.iterations();
      cg_converged = (pcg.info() == Eigen::Success);

      if(0)
      { // debug
        std::cerr << "pcg residual=" << (Z.transpose() * H * Z * dxz + gz).norm() / gz.norm()
                  << std::endl;
      }

      // compute dx of full problem
      dx = Z * dxz;

      // get maximal reasonable step
      double t_max = std::min(1.0,
                              max_feasible_step_safety_factor_ * _problem->max_feasible_step(x.data(), dx.data()));

      // backtracking line search
      double gdx = g.dot(dx);
      double t = t_max;
      xn = x + t*dx;
      double fxn = _problem->eval_f(xn.data());
      int iter_ls = 0;
      ls_succeeded = false;
      while (!(fxn <= fx + alpha_ls_ * gdx * t) && iter_ls < max_iter_ls) {
        t *= beta_ls_;
        xn = x + t * dx;
        fxn = _problem->eval_f(xn.data());
        ++iter_ls;
      }

      // update x, fx, xz
      if (iter_ls < max_iter_ls)
      {
        x.swap(xn);
        fx = fxn;
        ls_succeeded = true;
        // update reduced variables
        //        xz = Z.transpose()*(x-x0);
      }
      else
      {
        std::cerr << "Warning: line search failed ---> terminate" << std::endl;
        break;
      }


      DEB_line(4,
               "iter = " << iter << ", f(x) = " << fx << ", t = " << t
                         << " (tmax=" << t_max << "), " << "#ls = " << iter_ls
                         << ", |reduced grad| = " << gzn
                         << ", " << "PCG_tol = " << eta
                         << ", " << "PCG_iters = " << n_pcg_iters
                         << ", " << "PCG_converged = " << int(cg_converged));


    }

    // store result
    _problem->store_result(x.data());

    // return success
    return converged_;
  }

  // solve with linear constraints
  // It is required that the constraints are linearly independent
  int solve_projected_normal_equation( NProblemInterface* _problem, const SMatrixD& _A, const VectorD& _b )
  {
//    DEB_time_func_def;
    DEB_enter_func;

    converged_ = false;
    feasible_solution_found_ = false;

    // number of unknowns
    Eigen::Index n = _problem->n_unknowns();
    Eigen::Index m = _A.rows();

    DEB_line_if(!silent_, 2, "optimize via TruncatedNewtonProjectedNormalEquationsPCG with " << (int)n << " unknowns and " << (int)m << " linear constraints");

    // Newton parameters
    const int    max_iters  = max_iters_;
    const double newton_tol = eps_; // norm of gradient of reduced (unconstrained) problem
    // CG parameters
//    const int    max_pcg_iters = n-m;
    const int    max_pcg_iters = max_pcg_iters_;
    const double pcg_tol = 1e-4;
    // allow only steps larger than eps_ls_ in line search
    const int max_iter_ls   = int ( std::log(eps_ls_)/std::log(beta_ls_));

    // initialize vectors of unknowns
    VectorD x(n);
    _problem->initial_x(x.data());


    double constraint_violation = (_A*x-_b).lpNorm<Eigen::Infinity>();
    if(constraint_violation < eps_constraints_violation_)
      feasible_solution_found_ = true;
    DEB_line_if(!silent_, 2, "initial inf-norm constraint violation = " << constraint_violation << " numerically feasible = " << int(feasible_solution_found_));

    // storage of update vector dx and gradient
    VectorD g(n), dx(n), xn(n), gz(n), dz(n);
    dx.setZero(); // required for warmstart

    // storage CG
    VectorD  r(n), v(n), q(n), p(n), r2(n), q2(n), Hp(n);

    // hessian matrix
    SMatrixD H(n,n);

    // Diagonal Preconditioner
    VectorD W(n), Wi(n);

    // cholesky decomposition of projection
//    Eigen::SimplicialLDLT<SMatrixD, Eigen::Lower, Eigen::MetisOrdering<int> > ldlt;
//    Eigen::SimplicialLDLT<SMatrixD, Eigen::Lower, Eigen::NaturalOrdering<int> > ldlt;
//    Eigen::SimplicialLDLT<SMatrixD> ldlt;
    Eigen::CholmodSupernodalLLT<SMatrixD> ldlt;
//    Eigen::UmfPackLU<SMatrixD> ldlt;
//    Eigen::CholmodSimplicialLLT<SMatrixD> ldlt;

    // prepare solver
    bool ldlt_initialized = false;

    // statistics
    int n_pcg_iters(0);

    // choose precision of CG solve
    double eta = pcg_tol;

    bool ls_succeeded = true;

    bool cg_converged = false;

    // get function value at current point
    double fx = _problem->eval_f(x.data());
    DEB_line_if(!silent_, 2, "initial objective value = " << fx );
    if(!std::isfinite(fx))
    {
      DEB_line_if(!silent_, 2, "ERROR: intial objective value is not finite ---> abort");
      return false;
    }

    // data for adaptive hessian update
    int count_hessian_skip = 0;
    bool hessian_updated = false;
    double alpha_line_search = 0.0;
    double rel_objective_decrease = 0.0;

    size_t iter = 0;
    for(; iter<max_iters; ++iter)
    {
      // get gradient and Hessian
      _problem->eval_gradient(x.data(), g.data());

      if(iter == 0 || // first iteration
         count_hessian_skip >= hessian_max_skips_ || // maximum number of Hessians have been skipped
         alpha_line_search < hessian_min_acceptable_alpha_ || // line search truncated too much
         rel_objective_decrease < hessian_min_acceptable_rel_objective_decrease_ // objective did not decrease sufficiently
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

      // pre-factor projection matrix and update preconditioner
      if(iter == 0 || (always_update_preconditioner_ && hessian_updated) )
      {
        // setup preconditioner W = (|diag(H)|)
        // and Wi = (|diag(H)|)^-1

        double minW = DBL_MAX;
        double maxW = 0.0;

        for(int j=0; j<n ;++j)
        {
          W[j]  = std::abs(H.coeffRef(j,j));
          Wi[j] = 1.0/W[j];

          minW = std::min(minW,W[j]);
          maxW = std::max(maxW,W[j]);

//          // no preconditioner
//          W[j]  = 1.0;
//          Wi[j] = 1.0;
        }

        // make sure that the preconditioner is numerically safe
        if(maxW == 0.0 || maxW > max_preconditioner_range_*minW)
        {
          double minW_new = maxW/max_preconditioner_range_;
          DEB_line_if(!silent_, 2,  "iter = " << iter << ", correct preconditioner range [" << minW << ", " << maxW << "] to [" << minW_new << ", " << maxW << "]")

          for(int j=0; j<n ;++j)
            if(W[j] < minW_new)
            {
              W[j]  = minW_new;
              Wi[j] = 1.0/W[j];
            }
        }

        // prepare constraint projection
        if(iter == 0)
          DEB_line_if(!silent_, 2, "prepare constraint projection ..." );

        SMatrixD AWiAt = _A*Wi.asDiagonal()*_A.transpose();

        // update factorization (if valid matrix)
        if(_A.rows() > 0 && _A.cols() > 0)
        {
          if (!ldlt_initialized)
          {
            ldlt.analyzePattern(AWiAt);
            ldlt_initialized = true;
          }
          ldlt.factorize(AWiAt);

          // ldlt.compute(AWiAt); // old update


          if (ldlt.info() != Eigen::Success)
          {
            for (unsigned int j = 0; j < 10; ++j)
            {
              double reg = 1e-8 * AWiAt.diagonal().sum() / double(m);
              DEB_line_if(!silent_, 2,
                          "Warning: LDLT factorization failed --> regularize ADAt (could lead to reduced accuracy), reg="
                                  << reg);
              for (Eigen::Index j = 0; j < AWiAt.rows(); ++j)
                AWiAt.coeffRef(j, j) += reg;  // operation is safe since all diagonal entries are nonzero!!!

              ldlt.compute(AWiAt);
              if (ldlt.info() == Eigen::Success)
                break;
            }
          }
        }
        if(iter == 0)
          DEB_line_if(!silent_, 2, "done!" );
      }

      // perform feasiblity step?
      bool feasibility_step_productive = false;
      if( constraint_violation >= eps_constraints_violation_)
      {
//        // optimize constraint violation in Hessian-diagonal norm
//        v  = ldlt.solve(_b-_A*x);
//        dz = Wi.asDiagonal()*_A.transpose()*v;

       // optimize full KKT-system but approximate Hessian only through diagonal
        if(_A.rows() > 0 && _A.cols() > 0)
          v  = ldlt.solve(_b-_A*x + _A*Wi.asDiagonal()*g);
        else
          v.setZero();

        dz = Wi.asDiagonal()*(_A.transpose()*v - g);

        // debug
//        std::cerr << "debug: constraint violation after full projection = " << (_A*(x+dz)-_b).norm() << std::endl;
        // truncate feasiblity step
        double t_max = std::min(1.0,
                                max_infeasibility_step_safety_factor_ * _problem->max_feasible_step(x.data(), dz.data()));
  //      std::cerr << "debug: t_max = " << t_max << std::endl;

        // new x
        xn = x + t_max*dz;

        double constraint_violation_old = constraint_violation;
        constraint_violation = (_A*xn-_b).lpNorm<Eigen::Infinity>();

        // is step productive? (can be numerically non-productive)
        if(constraint_violation < constraint_violation_old)
        {
          DEB_line_if(!silent_, 2, "iter = " << iter << ", feasibility step reducing constraint violation "
                                             << constraint_violation_old << " ---> " << constraint_violation);

          // update objective
          fx = _problem->eval_f(xn.data());
          if(std::isfinite(fx)) // valid objective function?
          {
            // update
            xn.swap(x);
            // update gradient
            _problem->eval_gradient(x.data(), g.data());
            feasibility_step_productive = true;
          }
          else
          {
            // revert
            fx = _problem->eval_f(x.data());
            constraint_violation = constraint_violation_old;

            DEB_line_if(!silent_, 2, "iter = " << iter << ", feasibility step resulted in non-finite objective value --> revert to previous x");
          }
        }
        else
        {
          // revert step (do not swap xn and x)
          constraint_violation = constraint_violation_old;

          DEB_line_if(!silent_, 2, "iter = " << iter << ", feasibility step was not productive --> revert to previous x");
        }

        // feasible solution found?
        if(constraint_violation < eps_constraints_violation_)
          feasible_solution_found_ = true;
      }
      // ---------------------------
      // Projected-Preconditioned-CG

      // choose starting point
      bool warmstart = false;
      double gdx2 = g.dot(dx);
      if(gdx2 < 0.0 && allow_warmstart_)
        warmstart = true;
      else
        dx.setZero();

      r = g;

      // project r  --> q
      if(_A.rows() > 0 && _A.cols() > 0)
        v  = ldlt.solve(_A*Wi.asDiagonal()*r);
      else
        v.setZero();
      gz = r - _A.transpose()*v;
      q  = Wi.asDiagonal()*gz;

      p = -q;

      double rtq = r.dot(q);

      // check convergence
      double gzn = gz.norm();
      if(gzn < newton_tol && constraint_violation < eps_constraints_violation_)
      {
        constraint_violation = (_A*x-_b).lpNorm<Eigen::Infinity>();
        feasible_solution_found_ = true;
        converged_ = true;
        reduced_gradient_norm_ = gzn;
        DEB_line_if(!silent_, 4,"converged" <<  ", f(x) = " << fx
                               << ", ||gz|| = " << gzn
                               << ", max_constraint_violation = " << constraint_violation
                               << ", " << "PCG_iters = " << n_pcg_iters );
        break;
      }

      // choose CG tolerance
      if(adaptive_tolerance_)
        eta = adaptive_tolerance_modifier_*std::min(0.1, std::sqrt(gzn))*gzn;
      else
        eta = pcg_tol*gzn;

      for(int pcg_iter=0; pcg_iter < max_pcg_iters; ++pcg_iter)
      {
        // compute norm of projected residual (in original norm)
        double rpn = (W.asDiagonal()*q).norm();

        // debug
        if(0)
        {
          VectorD rr = H * dx + g;
          if(_A.rows() > 0 && _A.cols() > 0)
            v = ldlt.solve(_A*Wi.asDiagonal()*rr);
          else
            v.setZero();
          VectorD rr2 = (rr - _A.transpose()*v);

          if(!silent_)
          std::cerr << "++++++++++++++++++++++++++" << std::endl
                    << "||P(Hdx+g)|| = " << rr2.norm() << std::endl
                    << "||qn||       = " << rpn << std::endl
                    << "sqrt(rtq)    = " << std::sqrt(rtq) << std::endl
                    << std::endl;
        }

        // stop if accuracy sufficient
        if(rpn < eta)
        {
          cg_converged = true;
          break;
        }

        // cache re-used quantities
        Hp = H*p;
        double pHp = p.dot(Hp);

        // handle directions of negative curvature (important for Hessians of non-convex problems)
        if( pHp <= 0.0)
        {
          if (pcg_iter == 0)
          {
            // use p = proj(-grad f)
            // TODO: determine better scaling of step
            dx = p;
            DEB_line_if(!silent_, 4,  "iter = " << iter << ", PCG found direction of negative curvature at pcg_iter = 0" );
          }
          else
          {
            // use current dx
            DEB_line_if(!silent_, 4, "iter = " << iter << ", PCG found direction of negative curvature at pcg_iter = " << pcg_iter );
          }
          // quit PCG iteration
          break;
        }

          // optimal step length along p (minimizer of quadratic)
        double alpha = rtq/pHp;

        // update dx
        dx += alpha*p;
        r2 = r + alpha*Hp;

        // project r2  --> q2
        if(_A.rows() > 0 && _A.cols() > 0)
          v = ldlt.solve(_A*Wi.asDiagonal()*r2);
        else
          v.setZero();
        q2 = Wi.asDiagonal()*(r2 - _A.transpose()*v);


        double rtq2 = r2.dot(q2);
        double beta = rtq2/rtq;
        p = -q2 + beta*p;

        // swap vectors
        q.swap(q2);
        // r.swap(r2);
        // constraint refinement
        r = r2 - _A.transpose()*v;
        // update rtq
        rtq = rtq2;

        // count number of pcg iters
        ++n_pcg_iters;
      }

      // get maximal reasonable step
      double t_max = std::min(1.0,
                              max_feasible_step_safety_factor_ * _problem->max_feasible_step(x.data(), dx.data()));

      // backtracking line search
      double gdx = g.dot(dx);
      double t = t_max;
      xn = x + t*dx;
      double fxn = _problem->eval_f(xn.data());
      int iter_ls = 0;
      ls_succeeded = false;
      while (!(fxn <= fx + alpha_ls_ * gdx * t) && iter_ls < max_iter_ls) {
        t *= beta_ls_;
        xn = x + t * dx;
        fxn = _problem->eval_f(xn.data());
        ++iter_ls;
      }

      // store line search truncation
      alpha_line_search = t;
      // store relative reduction
      rel_objective_decrease = std::abs((fx-fxn)/fx);

      // update x, fx, xz
      if (iter_ls < max_iter_ls)
      {
        if(allow_warmstart_)
          dx=xn-x;
        x.swap(xn);
        fx = fxn;
        ls_succeeded = true;
        // update constraint violation
        constraint_violation = (_A*x-_b).lpNorm<Eigen::Infinity>();
        feasible_solution_found_ = (constraint_violation < eps_constraints_violation_);

        // update reduced variables
        //        xz = Z.transpose()*(x-x0);
      }
      else
      {
        if(feasible_solution_found_ || !feasibility_step_productive)
        {
          DEB_line_if(!silent_, 4, "Warning: line search failed ---> terminate");
          DEB_line_if(!silent_, 4, "|g| = " << g.norm());
          DEB_line_if(!silent_, 4, "|dx| = " << dx.norm());
          DEB_line_if(!silent_, 4, "|W| = " << W.norm());
          DEB_line_if(!silent_, 4, "|Wi| = " << Wi.norm());
          DEB_line_if(!silent_, 4, "gdx = " << gdx);

          reduced_gradient_norm_ = gzn;

          break;
        }
        else
          DEB_line_if(!silent_, 4, "Warning: line search failed ---> continue since feasibility_step_productive");
      }

      DEB_line_if(!silent_, 4,
               "iter = " << iter << ", f(x) = " << fx << ", t = " << t
                         << " (tmax=" << t_max << "), " << "#ls = " << iter_ls
                         << ", constraint_violation = " << constraint_violation
                         << ", |reduced grad| = " << gzn
                         << ", " << "PCG_tol = " << eta
                         << ", " << "PCG_iters = " << n_pcg_iters
                         << ", " << "PCG_converged = " << int(cg_converged)
                         << ", " << "g*dx = " << gdx
                         << ", " << "warmstart = " << int(warmstart)
                         << ", " << "hessian_update = " << int(hessian_updated)
                         << ", " << "rel_obj_decrease = " << rel_objective_decrease
                         );

      if(std::abs(gdx) < eps_gdx_)
      {
        constraint_violation = (_A*x-_b).lpNorm<Eigen::Infinity>();
        reduced_gradient_norm_ = gzn;
        if(constraint_violation < eps_constraints_violation_)
        {
          converged_ = true;
          feasible_solution_found_ = true;
          DEB_line_if(!silent_, 4, "converged" << ", f(x) = " << fx
                                               << ", ||gz|| = " << gzn
                                               << ", " << "g*dx = " << gdx
                                               << ", max_constraint_violation = " << constraint_violation
                                               << ", " << "PCG_iters = " << n_pcg_iters);
          break;
        }
      }

    }
    n_iterations_used_ = iter; // TODO: or increase? may make more sense for incremental solves.

    // store result
    _problem->store_result(x.data());

    // compute dual variables
    if(compute_dual_variables_)
    {
      if(!ldlt_initialized)
      {
        std::cerr << "ERROR: dual variables cannot be computed since ldlt is not initialized!!!" << std::endl;
        return converged_;
      }

      // Variant I ---> compute dual variables with preconditioner (does not require additional factorization)
      if(_A.rows() > 0 && _A.cols() > 0)
        nue_ = ldlt.solve(-_A*Wi.asDiagonal()*(g+H*dx));
      else
        nue_.setZero();

//      // Variant II ---> no preconditioner to comptue dual variables (requires additional factorization)
//      ldlt.compute(_A*_A.transpose());
//      nue_ = ldlt.solve(-_A*(g+H*dx));

//      std::cerr << "Truncated Newton dual variables residual ||A^T nue + g + H dx|| = " << (_A.transpose()*nue_+g+H*dx).norm() << std::endl;
//      std::cerr << "||dx|| = " << dx.norm() << std::endl;
    }

    // return success
    return converged_;
  }


  int solve(NProblemInterface* _problem, std::vector<LinearConstraint>& _constraints)
  {
    // convert constraints
    SMatrixD A;
    VectorD b;
    LinearConstraintConverter::nsolver_to_eigen(_constraints, A, b, _problem->n_unknowns());

    return solve(_problem, A, b);
  }

  int solve_experimental(NProblemInterface* _problem, std::vector<LinearConstraint>& _constraints)
  {
    // convert constraints
    SMatrixD A;
    VectorD b;
    LinearConstraintConverter::nsolver_to_eigen(_constraints, A, b, _problem->n_unknowns());

    return solve_experimental(_problem, A, b);
  }

  int solve_reduced_system(NProblemInterface* _problem, std::vector<LinearConstraint>& _constraints)
  {
    // convert constraints
    SMatrixD A;
    VectorD b;
    LinearConstraintConverter::nsolver_to_eigen(_constraints, A, b, _problem->n_unknowns());

    return solve_reduced_system(_problem, A, b);
  }

  int solve_reduced_system_EigenCG(NProblemInterface* _problem, std::vector<LinearConstraint>& _constraints)
  {
    // convert constraints
    SMatrixD A;
    VectorD b;
    LinearConstraintConverter::nsolver_to_eigen(_constraints, A, b, _problem->n_unknowns());

    return solve_reduced_system_EigenCG(_problem, A, b);
  }

  int solve_projected_normal_equation(NProblemInterface* _problem, std::vector<LinearConstraint>& _constraints)
  {
    // convert constraints
    SMatrixD A;
    VectorD b;
    LinearConstraintConverter::nsolver_to_eigen(_constraints, A, b, _problem->n_unknowns());

    return solve_projected_normal_equation(_problem, A, b);
  }

  int solve_projected_normal_equation(NProblemInterface* _problem, std::vector<NConstraintInterface*>& _constraints)
  {
    // convert constraints
    SMatrixD A;
    VectorD b;
    LinearConstraintConverter::nsolver_to_eigen(_constraints, A, b, _problem->n_unknowns());

    return solve_projected_normal_equation(_problem, A, b);
  }

  bool converged() { return converged_; }
  bool feasible_solution_found() { return feasible_solution_found_; }
  void set_silent(const bool _silent) { silent_ = _silent;}

  void set_eps_constraints_violation(double eps){ eps_constraints_violation_ = eps; }

  double reduced_gradient_norm() { return reduced_gradient_norm_;};

  bool& always_update_preconditioner() { return always_update_preconditioner_;};

  double& adaptive_tolerance_modifier() { return adaptive_tolerance_modifier_;};

  double& tolerance_gdx() { return eps_gdx_;}

  int& max_iters()     { return max_iters_;}
  int& max_pcg_iters() { return max_pcg_iters_;}

  bool& allow_warmstart() { return allow_warmstart_;};

  double& max_feasible_step_safety_factor() { return max_feasible_step_safety_factor_;}

  // parameters to adaptively skip hessian matrix updates
  int&    hessian_max_skips() {return hessian_max_skips_;}
  double& hessian_min_acceptable_alpha() {return hessian_min_acceptable_alpha_;}
  double& hessian_min_acceptable_rel_objective_decrease() {return hessian_min_acceptable_rel_objective_decrease_;}

  double& max_preconditioner_range() { return max_preconditioner_range_;}

  double& alpha_ls() {return alpha_ls_;}
  double& beta_ls()  {return beta_ls_;}
  size_t n_iterations_used() const {return n_iterations_used_;}

  bool& compute_dual_variables() {return compute_dual_variables_;}

  VectorD& nue() {return nue_;}

private:

  double eps_;
  double eps_ls_;
  double eps_gdx_;
  int    max_iters_;
  int    max_pcg_iters_;
  double alpha_ls_;
  double beta_ls_;

  double eps_constraints_violation_ = 1e-6;

  double max_preconditioner_range_ = 1e8;

  double max_feasible_step_safety_factor_;
  double max_infeasibility_step_safety_factor_ = 0.6;

  // adaptively choose tolerance of CG optimization?
  bool adaptive_tolerance_;
  bool always_update_preconditioner_;
  bool allow_warmstart_;

  // parameters for adaptive hessian update
  int    hessian_max_skips_ = 10; // maximum number of skipped hessian updates
  double hessian_min_acceptable_alpha_ = 0.1; // minimal acceptable line-search step length to skip hessian update
  double hessian_min_acceptable_rel_objective_decrease_ = 0.1; // minimal acceptable relative objectitve descrease to skip hessian update

  double adaptive_tolerance_modifier_;

  bool converged_;
  bool feasible_solution_found_;
  bool silent_;

  double reduced_gradient_norm_;

  bool compute_dual_variables_ = false;

  // dual variables
  VectorD nue_;

  size_t n_iterations_used_ = 0;
};


//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // COMISO_TRUNCATEDNEWTONPCG_HH defined
//=============================================================================

