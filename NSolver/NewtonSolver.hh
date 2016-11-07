//=============================================================================
//
//  CLASS NewtonSolver
//
//=============================================================================


#ifndef COMISO_NEWTONSOLVER_HH
#define COMISO_NEWTONSOLVER_HH

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>

//== INCLUDES =================================================================

#include <CoMISo/Config/CoMISoDefines.hh>
#include <CoMISo/Utils/StopWatch.hh>
#include "NProblemInterface.hh"
#include "NProblemGmmInterface.hh"

#if COMISO_SUITESPARSE_AVAILABLE
  #include <Eigen/UmfPackSupport>
#endif

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================

	      

/** \class NewtonSolver NewtonSolver.hh <ACG/.../NewtonSolver.hh>

    Brief Description.
  
    A more elaborate description follows.
*/
class COMISODLLEXPORT NewtonSolver
{
public:

  typedef Eigen::VectorXd             VectorD;
  typedef Eigen::SparseMatrix<double> SMatrixD;
  typedef Eigen::Triplet<double>      Triplet;

  /// Default constructor
  NewtonSolver(const double _eps = 1e-6, const int _max_iters = 200, const double _alpha_ls=0.2, const double _beta_ls = 0.6)
    : eps_(_eps), max_iters_(_max_iters), alpha_ls_(_alpha_ls), beta_ls_(_beta_ls), constant_hessian_structure_(false) {}

  /// Destructor
  ~NewtonSolver() {}

  // solve without linear constraints
  int solve(NProblemInterface* _problem)
  {
    SMatrixD A(0,_problem->n_unknowns());
    VectorD b(VectorD::Index(0));
    return solve(_problem, A, b);
  }

  // solve with linear constraints
  // Warning: so far only feasible starting points with (_A*_problem->initial_x() == b) are supported!
  // Extension to infeasible starting points is planned
  int solve(NProblemInterface* _problem, const SMatrixD& _A, const VectorD& _b)
  {
    // time solution procedure
    COMISO::StopWatch sw; sw.start();

    // number of unknowns
    int n = _problem->n_unknowns();
    // number of constraints
    int m = _b.size();

    std::cerr << "optmize via Newton with " << n << " unknowns and " << m << " linear constraints" << std::endl;

    // initialize vectors of unknowns
    VectorD x(n);
    _problem->initial_x(x.data());

    // storage of update vector dx and rhs of KKT system
    VectorD dx(n+m), rhs(n+m), g(n);
    rhs.setZero();

    // resize temp vector for line search (and set to x1 to approx Hessian correctly if problem is non-quadratic!)
    x_ls_ = x;

    // start with no regularization
    double regularize(0.0);
    int iter=0;
    while( iter < max_iters_)
    {
      // get Newton search direction by solving LSE
      factorize(_problem, _A, _b, x, regularize);

      // get rhs
      _problem->eval_gradient(x.data(), g.data());
      rhs.head(n) = -g;
      rhs.tail(m) = _b - _A*x;

      // solve KKT system
      dx = lu_solver_.solve(rhs);

      // get maximal reasonable step
      double t_max  = std::min(1.0, 0.5*_problem->max_feasible_step(x.data(), dx.data()));

      // perform line-search
      double newton_decrement(0.0);
      double fx(0.0);
      double t = backtracking_line_search(_problem, x, g, dx, newton_decrement, fx, t_max);

      x += dx.head(n)*t;

      std::cerr << "iter: " << iter
                << ", f(x) = " << fx
                << ", t = " << t << " (tmax=" << t_max << ")"
                << ", eps = [Newton decrement] = " << newton_decrement
                << ", constraint violation  = " << rhs.tail(m).norm()
                << std::endl;

      // converged?
      if(newton_decrement < eps_ || t == 0.0)
        break;

      ++iter;
    }

    // store result
    _problem->store_result(x.data());

    double solution_time = sw.stop();
    std::cerr << "Newton Method finished in " << solution_time/1000.0 << "s" << std::endl;

    // return success
    return 1;
  }

protected:

  void factorize(NProblemInterface* _problem, const SMatrixD& _A, const VectorD& _b, const VectorD& _x, double& _regularize)
  {
    const int n  = _problem->n_unknowns();
    const int m  = _A.rows();
    const int nf = n+m;

    // get hessian of quadratic problem
    SMatrixD H(n,n);
    _problem->eval_hessian(_x.data(), H);

    // set up KKT matrix
    // create sparse matrix
    std::vector< Triplet > trips;
    trips.reserve(H.nonZeros() + 2*_A.nonZeros());

    // add elements of H
    for (int k=0; k<H.outerSize(); ++k)
      for (SMatrixD::InnerIterator it(H,k); it; ++it)
        trips.push_back(Triplet(it.row(),it.col(),it.value()));

    // add elements of _A
    for (int k=0; k<_A.outerSize(); ++k)
      for (SMatrixD::InnerIterator it(_A,k); it; ++it)
      {
        // insert _A block below
        trips.push_back(Triplet(it.row()+n,it.col(),it.value()));

        // insert _A^T block right
        trips.push_back(Triplet(it.col(),it.row()+n,it.value()));
      }

    if(_regularize != 0.0)
      for( int i=0; i<m; ++i)
        trips.push_back(Triplet(n+i,n+i,_regularize));

    // create KKT matrix
    SMatrixD KKT(nf,nf);
    KKT.setFromTriplets(trips.begin(), trips.end());

    // compute LU factorization
    lu_solver_.compute(KKT);

    if(lu_solver_.info() != Eigen::Success)
    {
//      DEB_line(2, "Eigen::SparseLU reported problem: " << lu_solver_.lastErrorMessage());
//      DEB_line(2, "-> re-try with regularized constraints...");

      std::cerr << "Eigen::SparseLU reported problem while factoring KKT system: " << lu_solver_.lastErrorMessage() << std::endl;
//#if COMISO_SUITESPARSE_AVAILABLE
//      std::cerr << "Eigen::SparseLU reported problem while factoring KKT system. " << std::endl;
//#else
//      std::cerr << "Eigen::SparseLU reported problem while factoring KKT system: " << lu_solver_.lastErrorMessage() << std::endl;
//#endif

      // add more regularization
      if(_regularize == 0.0)
        _regularize = 1e-8;
      else
        _regularize *= 2.0;

      for( int i=0; i<m; ++i)
        trips.push_back(Triplet(n+i,n+i,_regularize));

      // create KKT matrix
      KKT.setFromTriplets(trips.begin(), trips.end());

      // compute LU factorization
      lu_solver_.compute(KKT);

//      IGM_THROW_if(lu_solver_.info() != Eigen::Success, QGP_BOUNDED_DISTORTION_FAILURE);
    }
  }

  double backtracking_line_search(NProblemInterface* _problem, VectorD& _x, VectorD& _g, VectorD& _dx, double& _newton_decrement, double& _fx, const double _t_start = 1.0)
  {
    int n = _x.size();

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

    std::cerr << "Warning: line search could not find a valid step within 100 iterations..." << std::endl;
    _fx = fx;
    return 0.0;
  }


  // deprecated function!
  // solve
  int solve(NProblemGmmInterface* _problem);

  // deprecated function!
  // solve specifying parameters
  int solve(NProblemGmmInterface* _problem, int _max_iter, double _eps)
  {
    max_iters_ = _max_iter;
    eps_ = _eps;
    return solve(_problem);
  }

  // deprecated function!
  bool& constant_hessian_structure() { return constant_hessian_structure_; }

protected:
  double* P(std::vector<double>& _v)
  {
    if( !_v.empty())
      return ((double*)&_v[0]);
    else
      return 0;
  }

private:

  double eps_;
  int    max_iters_;
//  double accelerate_;
  double alpha_ls_;
  double beta_ls_;

  VectorD x_ls_;

  // Sparse LU decomposition
  Eigen::SparseLU<SMatrixD> lu_solver_;
// ToDo: Check why UmfPack is not working!
//#if COMISO_SUITESPARSE_AVAILABLE
//  Eigen::UmfPackLU<SMatrixD> lu_solver_;
//#else
//  Eigen::SparseLU<SMatrixD> lu_solver_;
//#endif

  // deprecated
  bool   constant_hessian_structure_;
};


//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // ACG_NEWTONSOLVER_HH defined
//=============================================================================

