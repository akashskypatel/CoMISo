#pragma once
/*===========================================================================*\
 *                                                                           *
 *                            TruncatedNewtonPCG                             *
 *      Copyright (C) 2024 by Computer Graphics Group, University of Bern    *
 *                           http://cgg.unibe.ch                             *
 *                                                                           *
 *      Author: David Bommes                                                 *
 *                                                                           *
\*===========================================================================*/




//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>

//== INCLUDES =================================================================

#include <iomanip>
#include <CoMISo/Config/CoMISoDefines.hh>
#include <CoMISo/NSolver/NProblemInterface.hh>
#include <CoMISo/NSolver/NConstraintInterface.hh>
#include <CoMISo/NSolver/LinearConstraint.hh>

#include <Eigen/Dense>
#include <Eigen/Sparse>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================



class COMISODLLEXPORT TruncatedNewtonPCG
{
public:

  typedef Eigen::VectorXd             VectorD;
  typedef Eigen::SparseMatrix<double> SMatrixD;
  typedef Eigen::Triplet<double>      Triplet;

  struct OptimizerStatus
  {
    // check whether optimizer converged to a locally optimal point
    bool converged_to_local_optimum() const
    {
      return (     feasible
                && !negative_curvature_step
//                && line_search_t == 1.0
                && (newton_decrement_within_tolerance || projected_gradient_norm_within_tolerance) );
    }

    // check whether optimizer converged to a locally infeasible point
    bool converged_to_infeasible_point() const
    {
      return (!feasible && !feasibility_step_productive && (newton_decrement_within_tolerance || line_search_t == 0.0));
    }

    // function value at last iterate
    double fx = DBL_MAX;

    // was last iterate feasible w.r.t. eps_constraints_violation_, i.e. ||residual||_max <= eps_constraints_violation_
    bool feasible = false;
    // constraint violation ||residual||_max of last iterate
    double constraint_violation_inf_norm = DBL_MAX;
    // did feasibility step of last iterate improve feasbility?
    bool feasibility_step_productive = false;

    // was Newton decrement of last iterate within specified tolerance eps_gdx_?
    bool   newton_decrement_within_tolerance = false;
    // Newton decrement of last iterate
    double newton_decrement = DBL_MAX;

    // was projected gradient norm of last iterate within specified tolerance eps_?
    bool projected_gradient_norm_within_tolerance = false;
    // projected gradient norm of last iterate
    double projected_gradient_norm = DBL_MAX;

    // did last iterate perform a step into negative curvature direction?
    bool negative_curvature_step = false;

    // line search parameter of last step in [0,1]
    double line_search_t = 0.0;
    double line_search_t_max_feasible = 0.0;
    int    line_search_iterations = 0;

    // line search parameter of last step in [0,1]
    double line_search_t_inf = 0.0;
    double line_search_t_inf_max_feasible = 0.0;
    int    line_search_inf_iterations = 0;

    // number of refinement iters in current Newton iter
    int refinement_iters = 0;
    // total number of iterative refinement iterations
    int refinement_iters_total = 0;

    // ConjugateGradient data
    bool cg_converged = false;
    int  cg_iterations = 0;
    int  cg_iterations_total = 0;

    // Hessian updated in last iteration?
    bool hessian_updated = false;

    // total number of performed Newton iterations
    int n_newton_iters = 0;
    // total number of negative curvature iterations
    int n_negative_curvature_iters = 0;
  };


  /// Default constructor  (for default parameters see initialization below)
  TruncatedNewtonPCG()
  {}

  /// deprecated: old constructor provided for backwards compatibility
  TruncatedNewtonPCG(const double _eps, const double _eps_line_search = 1e-8,
                     const int _max_iters = 500, const double _alpha_ls = 0.2,
                     const double _beta_ls = 0.6)
          : eps_(_eps), eps_ls_(_eps_line_search), eps_gdx_(0.1), max_iters_(_max_iters), max_pcg_iters_(500),
            alpha_ls_(_alpha_ls), beta_ls_(_beta_ls), max_feasible_step_safety_factor_(0.5),
            adaptive_tolerance_(true), always_update_preconditioner_(true), allow_warmstart_(false), adaptive_tolerance_modifier_(1.0), silent_(false)
  {
  }

  // optimize unconstrained problem
  int solve(NProblemInterface* _problem);

  // optimize with linear constraints
  // Note: It is required that the constraints are linearly independent
  int solve( NProblemInterface* _problem, const SMatrixD& _A, const VectorD& _b );
  int solve( NProblemInterface* _problem, std::vector<LinearConstraint>& _constraints);
  int solve( NProblemInterface* _problem, std::vector<NConstraintInterface*>& _constraints);
  // deprecated naming: identical to above but kept for backward compatibility
  int solve_projected_normal_equation( NProblemInterface* _problem, const SMatrixD& _A, const VectorD& _b );
  int solve_projected_normal_equation(NProblemInterface* _problem, std::vector<LinearConstraint>& _constraints);
  int solve_projected_normal_equation(NProblemInterface* _problem, std::vector<NConstraintInterface*>& _constraints);

//  // solve with linear constraints
//  // Warning: so far only feasible starting points with (_A*_problem->initial_x() == b) are supported!
//  // It is also required that the constraints are linearly independent
//  int solve_experimental(NProblemInterface* _problem, const SMatrixD& _A, const VectorD& _b)
//  int solve_reduced_system( NProblemInterface* _problem, const SMatrixD& _A, const VectorD& _b )
//  int solve_reduced_system_EigenCG( NProblemInterface* _problem, const SMatrixD& _A, const VectorD& _b )
//  int solve_experimental(NProblemInterface* _problem, std::vector<LinearConstraint>& _constraints)
//  int solve_reduced_system(NProblemInterface* _problem, std::vector<LinearConstraint>& _constraints)
//  int solve_reduced_system_EigenCG(NProblemInterface* _problem, std::vector<LinearConstraint>& _constraints)


  // obtain deatiled information of optimization status
  OptimizerStatus& status() {return status_;}

  bool converged() { return status_.converged_to_local_optimum(); }
  bool feasible_solution_found() { return status_.feasible; }
  void set_silent(const bool _silent) { silent_ = _silent;}

  void set_eps_constraints_violation(double eps){ eps_constraints_violation_ = eps; }

  void set_eps_constraints_violation_desirable(double eps){ eps_constraints_violation_desirable_ = eps; }

  double reduced_gradient_norm() { return status_.projected_gradient_norm;};

  bool& always_update_preconditioner() { return always_update_preconditioner_;};

  double& adaptive_tolerance_modifier() { return adaptive_tolerance_modifier_;};

  double& tolerance_newton_decrement() { return eps_gdx_;}
  // deprecated
  double& tolerance_gdx() { return eps_gdx_;}

  double& tolerance_reduced_gradient_norm() { return eps_;}

  int& max_iters()     { return max_iters_;}
  int& max_pcg_iters() { return max_pcg_iters_;}

  bool& allow_warmstart() { return allow_warmstart_;};

  double& max_feasible_step_safety_factor() { return max_feasible_step_safety_factor_;}

  // parameters to adaptively skip hessian matrix updates
  int&    hessian_max_skips() {return hessian_max_skips_;}
  double& hessian_min_acceptable_alpha() {return hessian_min_acceptable_alpha_;}
  double& hessian_min_acceptable_rel_objective_decrease() {return hessian_min_acceptable_rel_objective_decrease_;}

  [[deprecated("Not used anymore, replaced by min_preconditioner_value and max_preconditioner_value")]]
  double& max_preconditioner_range() { return max_preconditioner_range_;}

  double& min_preconditioner_value() { return min_preconditioner_value_;}
  double& max_preconditioner_value() { return max_preconditioner_value_;}

  double& alpha_ls() {return alpha_ls_;}
  double& beta_ls()  {return beta_ls_;}
  size_t n_iterations_used() const {return n_iterations_used_;}

  bool line_search_feasibility_step() {return line_search_feasibility_step_;}

  bool& compute_dual_variables() {return compute_dual_variables_;}

  VectorD& nue() {return nue_;}

  double backtracking_line_search(NProblemInterface* _problem,
                           const VectorD& _x, const double _fx,
                           const VectorD& _dx,
                           const double _gdx, const double _t_max, const int _max_iter_ls,
                           VectorD& _x_new, double& _fx_new, int& _iter_ls) const;

  double line_search_negative_curvature( NProblemInterface* _problem,
                                         const VectorD& _x, const double _fx,
                                         const VectorD& _dx,
                                         const double _t_max, const int _max_iter_ls,
                                         VectorD& _x_new, double& _fx_new, int& _iter_ls) const;


  double backtracking_line_search_infeasible_merit_l1(NProblemInterface* _problem, const SMatrixD& _H,
                                               const SMatrixD& _A, const VectorD& _b,
                                               const VectorD& _x, const double& _fx,
                                               const VectorD& _g, VectorD& _dx,
                                               VectorD& _x_new, double& _fx_new,
                                               double& _mu_merit,
                                               const double _t_start, const int _max_ls_iters, int& _n_iters) const;

  double max_abs_cos_angle(const SMatrixD& _A, const VectorD& _A_inv_row_norm, const VectorD& _v) const;

  void print_iteration_data(const OptimizerStatus& _status) const;
  void print_summary(const OptimizerStatus& _status) const;


private:

  double eps_           = 1e-3;
  double eps_ls_        = 1e-8;
  double eps_gdx_       = 1e-3;
  int    max_iters_     = 500;
  int    max_pcg_iters_ = 500;
  double pcg_tolerance_ = 1e-4;

  double alpha_ls_      = 0.1;
  double beta_ls_       = 0.8;

  // parameters of L1 Merit function line search
  bool   line_search_feasibility_step_ = true;
  double rho_merit_     = 0.5;

  // max inf-norm constraint violation allowed for feasibility
  double eps_constraints_violation_ = 1e-6;
  // max inf-norm constraint violation above which feasibility steps are still performed
  double eps_constraints_violation_desirable_ = 1e-9;

  // deprecated parameter
  double max_preconditioner_range_ = 1e9;

  // set limits for preconditioner values
//  double min_preconditioner_value_ = 1e-12;
//  double max_preconditioner_value_ = 1e12;
  double min_preconditioner_value_ = 1e-6;
  double max_preconditioner_value_ = 1e6;

  double max_feasible_step_safety_factor_ = 0.5;
  double max_infeasibility_step_safety_factor_ = 0.6;

  // adaptively choose tolerance of CG optimization?
  bool adaptive_tolerance_ = true;
  bool always_update_preconditioner_ = true;
  bool allow_warmstart_ = false;

  // parameters for adaptive hessian update
  int    hessian_max_skips_ = 10; // maximum number of skipped hessian updates
  double hessian_min_acceptable_alpha_ = 0.1; // minimal acceptable line-search step length to skip hessian update
  double hessian_min_acceptable_rel_objective_decrease_ = 0.1; // minimal acceptable relative objective descrease to skip hessian update

  double adaptive_tolerance_modifier_ = 1.0;

  bool compute_dual_variables_ = false;

  // iterative refinement parameters
  bool enable_iterative_refinement_ = true;
  int  max_iterative_refinement_iters_ = 5;
  double iterative_refinement_cos_angle_threshold_ = 1e-12;
  // perform additional projection for Newton search direction (should not be necessary when iterative refinement is used)
  bool project_dx_before_update_ = false;

  // dual variables
  VectorD nue_;

  size_t n_iterations_used_ = 0;

  // Optimizer Status for last iterate
  OptimizerStatus status_;

  bool silent_ = false;

  // deprecated ---> remove and only use status_
  bool feasible_solution_found_ = false;
};


//=============================================================================
} // namespace COMISO
//=============================================================================

