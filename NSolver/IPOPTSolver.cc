//=============================================================================
//
//  CLASS IPOPTSolver - IMPLEMENTATION
//
//=============================================================================

//== INCLUDES =================================================================

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_IPOPT_AVAILABLE
//=============================================================================


#include "IPOPTSolver.hh"

#include <Base/OutcomeUtils.hh>
#include <Base/Debug/DebOut.hh>


#include <IpTNLP.hpp>
#include <IpIpoptApplication.hpp>
#include <IpSolveStatistics.hpp>



//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION PROBLEM INSTANCE=========================================================


class NProblemIPOPT : public Ipopt::TNLP
{
public:

  // Ipopt Types
  typedef Ipopt::Number                    Number;
  typedef Ipopt::Index                     Index;
  typedef Ipopt::SolverReturn              SolverReturn;
  typedef Ipopt::IpoptData                 IpoptData;
  typedef Ipopt::IpoptCalculatedQuantities IpoptCalculatedQuantities;

  // sparse matrix and vector types
  typedef NConstraintInterface::SVectorNC SVectorNC;
  typedef NConstraintInterface::SMatrixNC SMatrixNC;
  typedef NProblemInterface::SMatrixNP    SMatrixNP;

  /** default constructor */
  NProblemIPOPT(NProblemInterface* _problem, const std::vector<NConstraintInterface*>& _constraints)
   : problem_(_problem), store_solution_(false) { split_constraints(_constraints); analyze_special_properties(_problem, _constraints);}

  /** default destructor */
  virtual ~NProblemIPOPT() {};

  /**@name Overloaded from TNLP */
  //@{
  /** Method to return some info about the nlp */
  virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style);

  /** Method to return the bounds for my problem */
  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u);

  /** Method to return the starting point for the algorithm */
  virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda);

  /** Method to return the objective value */
  virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

  /** Method to return the gradient of the objective */
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

  /** Method to return the constraint residuals */
  virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

  /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
  virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
                          Number* values);

  /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   */
  virtual bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess, Index* iRow,
                      Index* jCol, Number* values);

  //@}

  /** @name Solution Methods */
  //@{
  /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
  virtual void finalize_solution(SolverReturn status,
                                 Index n, const Number* x, const Number* z_L, const Number* z_U,
                                 Index m, const Number* g, const Number* lambda,
                                 Number obj_value,
                                 const IpoptData* ip_data,
                                 IpoptCalculatedQuantities* ip_cq);
  //@}

  // special properties of problem
  bool hessian_constant() const;
  bool jac_c_constant() const;
  bool jac_d_constant() const;

  bool&                 store_solution()  {return store_solution_; }
  std::vector<double>&  solution()        {return x_;}

private:
  /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   *
   */
  //@{
  //  MyNLP();
  NProblemIPOPT(const NProblemIPOPT&);
  NProblemIPOPT& operator=(const NProblemIPOPT&);
  //@}

  // split user-provided constraints into general-constraints and bound-constraints
  void split_constraints(const std::vector<NConstraintInterface*>& _constraints);

  // determine if hessian_constant, jac_c_constant or jac_d_constant
  void analyze_special_properties(const NProblemInterface* _problem, const std::vector<NConstraintInterface*>& _constraints);


protected:
  double* P(std::vector<double>& _v)
  {
    if( !_v.empty())
      return ((double*)&_v[0]);
    else
      return 0;
  }

private:

  // pointer to problem instance
  NProblemInterface* problem_;
  // reference to constraints vector
  std::vector<NConstraintInterface*> constraints_;
  std::vector<BoundConstraint*>      bound_constraints_;

  bool hessian_constant_;
  bool jac_c_constant_;
  bool jac_d_constant_;

  bool store_solution_;
  std::vector<double> x_;
};


//== CLASS DEFINITION PROBLEM INSTANCE=========================================================


class NProblemGmmIPOPT : public Ipopt::TNLP
{
public:

  // Ipopt Types
  typedef Ipopt::Number                    Number;
  typedef Ipopt::Index                     Index;
  typedef Ipopt::SolverReturn              SolverReturn;
  typedef Ipopt::IpoptData                 IpoptData;
  typedef Ipopt::IpoptCalculatedQuantities IpoptCalculatedQuantities;

  // sparse matrix and vector types
  typedef NConstraintInterface::SVectorNC SVectorNC;
  typedef NConstraintInterface::SMatrixNC SMatrixNC;
  typedef gmm::wsvector<double>           SVectorNP;
  typedef NProblemGmmInterface::SMatrixNP SMatrixNP;

  typedef gmm::array1D_reference<       double* > VectorPT;
  typedef gmm::array1D_reference< const double* > VectorPTC;

  typedef gmm::array1D_reference<       Index* > VectorPTi;
  typedef gmm::array1D_reference< const Index* > VectorPTCi;

  typedef gmm::linalg_traits<SVectorNP>::const_iterator SVectorNP_citer;
  typedef gmm::linalg_traits<SVectorNP>::iterator       SVectorNP_iter;

  /** default constructor */
  NProblemGmmIPOPT(NProblemGmmInterface* _problem, std::vector<NConstraintInterface*>& _constraints)
   : problem_(_problem), constraints_(_constraints), nnz_jac_g_(0), nnz_h_lag_(0) 
   {}

  /** default destructor */
  virtual ~NProblemGmmIPOPT() {};

  /**@name Overloaded from TNLP */
  //@{
  /** Method to return some info about the nlp */
  virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style);

  /** Method to return the bounds for my problem */
  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u);

  /** Method to return the starting point for the algorithm */
  virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda);

  /** Method to return the objective value */
  virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

  /** Method to return the gradient of the objective */
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

  /** Method to return the constraint residuals */
  virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

  /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
  virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
                          Number* values);

  /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   */
  virtual bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess, Index* iRow,
                      Index* jCol, Number* values);

  //@}

  /** @name Solution Methods */
  //@{
  /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
  virtual void finalize_solution(SolverReturn status,
                                 Index n, const Number* x, const Number* z_L, const Number* z_U,
                                 Index m, const Number* g, const Number* lambda,
                                 Number obj_value,
                                 const IpoptData* ip_data,
                                 IpoptCalculatedQuantities* ip_cq);
  //@}

private:
  /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   *
   */
  //@{
  //  MyNLP();
  NProblemGmmIPOPT(const NProblemGmmIPOPT&);
  NProblemGmmIPOPT& operator=(const NProblemGmmIPOPT&);
  //@}


private:

  // pointer to problem instance
  NProblemGmmInterface* problem_;
  // reference to constraints vector
  std::vector<NConstraintInterface*>& constraints_;

  int nnz_jac_g_;
  int nnz_h_lag_;

  // constant structure of jacobian_of_constraints and hessian_of_lagrangian
  std::vector<Index> jac_g_iRow_;
  std::vector<Index> jac_g_jCol_;
  std::vector<Index> h_lag_iRow_;
  std::vector<Index> h_lag_jCol_;

  // Sparse Matrix of problem (don't initialize every time!!!)
  SMatrixNP HP_;
};


//== IMPLEMENTATION ========================================================== 


// smart pointer to IpoptApplication to set options etc.
class IPOPTSolver::Impl 
{// Create an instance of the IpoptApplication
public:
  Impl() : app_(IpoptApplicationFactory()) {}

public:
  Ipopt::SmartPtr<Ipopt::IpoptApplication> app_;
};

// Constructor
IPOPTSolver::IPOPTSolver()
  : impl_(new Impl)
{

  // Switch to HSL if available in Comiso
#if COMISO_HSL_AVAILABLE
  impl_->app_->Options()->SetStringValue("linear_solver", "ma57");
#else
  impl_->app_->Options()->SetStringValue("linear_solver", "mumps");
#endif

#ifdef WIN32
  // Restrict memory to be able to run larger problems on windows
  // with the default mumps solver
  // TODO: find out what this does and whether it makes sense to do it
  impl_->app_->Options()->SetIntegerValue("mumps_mem_percent", 5);
#endif

  // set default parameters
  impl_->app_->Options()->SetIntegerValue("max_iter", 100);
  //  app->Options()->SetStringValue("derivative_test", "second-order");
  //  app->Options()->SetIntegerValue("print_level", 0);
  //  app->Options()->SetStringValue("expect_infeasible_problem", "yes");
}

IPOPTSolver::~IPOPTSolver()
{ delete impl_; }
//-----------------------------------------------------------------------------


void IPOPTSolver::solve(NProblemInterface* _problem, 
  const std::vector<NConstraintInterface*>& _constraints)
{
  DEB_enter_func;
  //----------------------------------------------------------------------------
  // 1. Create an instance of IPOPT NLP
  //----------------------------------------------------------------------------
  Ipopt::SmartPtr<Ipopt::TNLP> np = new NProblemIPOPT(_problem, _constraints);
  NProblemIPOPT* np2 = dynamic_cast<NProblemIPOPT*> (Ipopt::GetRawPtr(np));

  //----------------------------------------------------------------------------
  // 2. exploit special characteristics of problem
  //----------------------------------------------------------------------------

  DEB_out(2,"exploit detected special properties: ");
  if(np2->hessian_constant())
  {
    DEB_out(2,"*constant hessian* ");
    impl_->app_->Options()->SetStringValue("hessian_constant", "yes");
  }

  if(np2->jac_c_constant())
  {
    DEB_out(2, "*constant jacobian of equality constraints* ");
    impl_->app_->Options()->SetStringValue("jac_c_constant", "yes");
  }

  if(np2->jac_d_constant())
  {
    DEB_out(2, "*constant jacobian of in-equality constraints*");
    impl_->app_->Options()->SetStringValue("jac_d_constant", "yes");
  }
  DEB_out(2,"\n");

  //----------------------------------------------------------------------------
  // 3. solve problem
  //----------------------------------------------------------------------------

  // Initialize the IpoptApplication and process the options
  Ipopt::ApplicationReturnStatus status = impl_->app_->Initialize();
  if (status != Ipopt::Solve_Succeeded) 
    THROW_OUTCOME(IPOPT_INITIALIZATION_FAILED);

  status = impl_->app_->OptimizeTNLP( np);

  //----------------------------------------------------------------------------
  // 4. output statistics
  //----------------------------------------------------------------------------
  if (!(status == Ipopt::Solve_Succeeded || status == Ipopt::Solved_To_Acceptable_Level))
  {
    // TODO: we could trnslate these return codes, but will not do it for now
    //  enum ApplicationReturnStatus
    //    {
    //      Solve_Succeeded=0,
    //      Solved_To_Acceptable_Level=1,
    //      Infeasible_Problem_Detected=2,
    //      Search_Direction_Becomes_Too_Small=3,
    //      Diverging_Iterates=4,
    //      User_Requested_Stop=5,
    //      Feasible_Point_Found=6,
    //
    //      Maximum_Iterations_Exceeded=-1,
    //      Restoration_Failed=-2,
    //      Error_In_Step_Computation=-3,
    //      Maximum_CpuTime_Exceeded=-4,
    //      Not_Enough_Degrees_Of_Freedom=-10,
    //      Invalid_Problem_Definition=-11,
    //      Invalid_Option=-12,
    //      Invalid_Number_Detected=-13,
    //
    //      Unrecoverable_Exception=-100,
    //      NonIpopt_Exception_Thrown=-101,
    //      Insufficient_Memory=-102,
    //      Internal_Error=-199
    //    };
    //------------------------------------------------------

    THROW_OUTCOME(IPOPT_OPTIMIZATION_FAILED);
  }
  
  // Retrieve some statistics about the solve
  Ipopt::Index iter_count = impl_->app_->Statistics()->IterationCount();
  DEB_out(1,"\n\n*** IPOPT: The problem solved in " << iter_count << "iterations!\n");

  Ipopt::Number final_obj = impl_->app_->Statistics()->FinalObjective();
  DEB_out(1,"\n\n*** IPOPT: The final value of the objective function is "
    << final_obj << "\n");
}


//-----------------------------------------------------------------------------


void IPOPTSolver::solve(
      NProblemInterface*                        _problem,
      const std::vector<NConstraintInterface*>& _constraints,
      const std::vector<NConstraintInterface*>& _lazy_constraints,
      const double                              _almost_infeasible,
      const int                                 _max_passes        )
{
  DEB_enter_func;
  //----------------------------------------------------------------------------
  // 0. Initialize IPOPT Applicaiton
  //----------------------------------------------------------------------------
  
  StopWatch sw; sw.start();

  // Initialize the IpoptApplication and process the options
  Ipopt::ApplicationReturnStatus status;
  status = impl_->app_->Initialize();
  if (status != Ipopt::Solve_Succeeded)
    THROW_OUTCOME(IPOPT_INITIALIZATION_FAILED);

  bool feasible_point_found = false;
  int  cur_pass = 0;
  double acceptable_tolerance = 0.01; // hack: read out from ipopt!!!
  // copy default constraints
  std::vector<NConstraintInterface*> constraints = _constraints;
  std::vector<bool> lazy_added(_lazy_constraints.size(),false);

  // cache statistics of all iterations
  std::vector<int> n_inf;
  std::vector<int> n_almost_inf;

  while(!feasible_point_found && cur_pass <(_max_passes-1))
  {
    ++cur_pass;
    //----------------------------------------------------------------------------
    // 1. Create an instance of current IPOPT NLP
    //----------------------------------------------------------------------------
    Ipopt::SmartPtr<Ipopt::TNLP> np = new NProblemIPOPT(_problem, constraints);
    NProblemIPOPT* np2 = dynamic_cast<NProblemIPOPT*> (Ipopt::GetRawPtr(np));
    // enable caching of solution
    np2->store_solution() = true;

    //----------------------------------------------------------------------------
    // 2. exploit special characteristics of problem
    //----------------------------------------------------------------------------

    DEB_out(2, "detected special properties which will be exploit: ");
    if(np2->hessian_constant())
    {
      DEB_out(2, "*constant hessian* ");
      impl_->app_->Options()->SetStringValue("hessian_constant", "yes");
    }

    if(np2->jac_c_constant())
    {
      DEB_out(2, "*constant jacobian of equality constraints* ");
      impl_->app_->Options()->SetStringValue("jac_c_constant", "yes");
    }

    if(np2->jac_d_constant())
    {
      DEB_out(2, "*constant jacobian of in-equality constraints*");
      impl_->app_->Options()->SetStringValue("jac_d_constant", "yes");
    }
    DEB_out(2, "\n");

    //----------------------------------------------------------------------------
    // 3. solve problem
    //----------------------------------------------------------------------------
    status = impl_->app_->OptimizeTNLP( np);

    // check lazy constraints
    n_inf.push_back(0);
    n_almost_inf.push_back(0);
    feasible_point_found = true;
    for(unsigned int i=0; i<_lazy_constraints.size(); ++i)
      if(!lazy_added[i])
      {
        NConstraintInterface* lc = _lazy_constraints[i];

        double v = lc->eval_constraint(&(np2->solution()[0]));

        bool inf        = false;
        bool almost_inf = false;

        if(lc->constraint_type() == NConstraintInterface::NC_EQUAL)
        {
          v = std::abs(v);
          if(v>acceptable_tolerance)
            inf = true;
          else
            if(v>_almost_infeasible)
              almost_inf = true;
        }
        else
          if(lc->constraint_type() == NConstraintInterface::NC_GREATER_EQUAL)
          {
            if(v<-acceptable_tolerance)
              inf = true;
            else
              if(v<_almost_infeasible)
                almost_inf = true;
          }
          else
            if(lc->constraint_type() == NConstraintInterface::NC_LESS_EQUAL)
            {
              if(v>acceptable_tolerance)
                inf = true;
              else
                if(v>-_almost_infeasible)
                  almost_inf = true;
            }

        // infeasible?
        if(inf)
        {
          constraints.push_back(lc);
          lazy_added[i] = true;
          feasible_point_found = false;
          ++n_inf.back();
        }

        // almost violated or violated? -> add to constraints
        if(almost_inf)
        {
          constraints.push_back(lc);
          lazy_added[i] = true;
          ++n_almost_inf.back();
        }
      }
  }

  // no termination after max number of passes?
  if(!feasible_point_found)
  {
    ++cur_pass;

    DEB_warning(2, "*************** could not find feasible point after "
      << _max_passes-1 << " -> solving with all lazy constraints...\n");
    for(unsigned int i=0; i<_lazy_constraints.size(); ++i)
      if(!lazy_added[i])
        constraints.push_back(_lazy_constraints[i]);

    //----------------------------------------------------------------------------
    // 1. Create an instance of current IPOPT NLP
    //----------------------------------------------------------------------------
    Ipopt::SmartPtr<Ipopt::TNLP> np = new NProblemIPOPT(_problem, constraints);
    NProblemIPOPT* np2 = dynamic_cast<NProblemIPOPT*> (Ipopt::GetRawPtr(np));
    // enable caching of solution
    np2->store_solution() = true;

    //----------------------------------------------------------------------------
    // 2. exploit special characteristics of problem
    //----------------------------------------------------------------------------

    DEB_out(2, "exploit detected special properties: ");
    if(np2->hessian_constant())
    {
      DEB_out(2, "*constant hessian* ");
      impl_->app_->Options()->SetStringValue("hessian_constant", "yes");
    }

    if(np2->jac_c_constant())
    {
      DEB_out(2, "*constant jacobian of equality constraints* ");
      impl_->app_->Options()->SetStringValue("jac_c_constant", "yes");
    }

    if(np2->jac_d_constant())
    {
      DEB_out(2, "*constant jacobian of in-equality constraints*");
      impl_->app_->Options()->SetStringValue("jac_d_constant", "yes");
    }
    std::cerr << std::endl;

    //----------------------------------------------------------------------------
    // 3. solve problem
    //----------------------------------------------------------------------------
    status = impl_->app_->OptimizeTNLP( np);
  }

  const double overall_time = sw.stop()/1000.0;

  //----------------------------------------------------------------------------
  // 4. output statistics
  //----------------------------------------------------------------------------
  if (!(status == Ipopt::Solve_Succeeded || status == Ipopt::Solved_To_Acceptable_Level))
    THROW_OUTCOME(IPOPT_OPTIMIZATION_FAILED);

  // Retrieve some statistics about the solve
  Ipopt::Index iter_count = impl_->app_->Statistics()->IterationCount();
  DEB_out(1, "\n\n*** IPOPT: The problem solved in " 
    << iter_count << " iterations!\n");

  Ipopt::Number final_obj = impl_->app_->Statistics()->FinalObjective();
  DEB_out(1, "\n\n*** IPOPT: The final value of the objective function is "
    << final_obj << "\n");

  DEB_out(2, "############# IPOPT with lazy constraints statistics ###############\n");
  DEB_out(2, "overall time: " << overall_time << "s\n");
  DEB_out(2, "#passes     : " << cur_pass << "( of " << _max_passes << ")\n");
  for(unsigned int i=0; i<n_inf.size(); ++i)
    DEB_out(3, "pass " << i << " induced " << n_inf[i] 
      << " infeasible and " << n_almost_inf[i] << " almost infeasible\n")
}


//-----------------------------------------------------------------------------


void IPOPTSolver::solve(NProblemInterface*    _problem)
{
  std::vector<NConstraintInterface*> constraints;
  solve(_problem, constraints);
}


//-----------------------------------------------------------------------------


void IPOPTSolver::solve(NProblemGmmInterface* _problem, std::vector<NConstraintInterface*>& _constraints)
{
  DEB_enter_func;
  DEB_warning(1,"******NProblemGmmInterface is deprecated!!! -> use NProblemInterface *******");

  //----------------------------------------------------------------------------
  // 1. Create an instance of IPOPT NLP
  //----------------------------------------------------------------------------
  Ipopt::SmartPtr<Ipopt::TNLP> np = new NProblemGmmIPOPT(_problem, _constraints);

  //----------------------------------------------------------------------------
  // 2. solve problem
  //----------------------------------------------------------------------------

  // Initialize the IpoptApplication and process the options
  Ipopt::ApplicationReturnStatus status = impl_->app_->Initialize();
  if (status != Ipopt::Solve_Succeeded)
     THROW_OUTCOME(IPOPT_INITIALIZATION_FAILED);

  //----------------------------------------------------------------------------
  // 3. solve problem
  //----------------------------------------------------------------------------
  status = impl_->app_->OptimizeTNLP(np);

  //----------------------------------------------------------------------------
  // 4. output statistics
  //----------------------------------------------------------------------------
  if (!(status == Ipopt::Solve_Succeeded || status == Ipopt::Solved_To_Acceptable_Level))
    THROW_OUTCOME(IPOPT_OPTIMIZATION_FAILED);

  // Retrieve some statistics about the solve
  Ipopt::Index iter_count = impl_->app_->Statistics()->IterationCount();
  DEB_out(1,"\n\n*** IPOPT: The problem solved in " << iter_count << " iterations!\n");

  Ipopt::Number final_obj = impl_->app_->Statistics()->FinalObjective();
  DEB_out(1, "\n\n*** IPOPT: The final value of the objective function is " << final_obj << "\n");
}


//== IMPLEMENTATION PROBLEM INSTANCE==========================================================


void
NProblemIPOPT::
split_constraints(const std::vector<NConstraintInterface*>& _constraints)
{
  // split user-provided constraints into general-constraints and bound-constraints
  constraints_      .clear();       constraints_.reserve(_constraints.size());
  bound_constraints_.clear(); bound_constraints_.reserve(_constraints.size());

  for(unsigned int i=0; i<_constraints.size(); ++i)
  {
    BoundConstraint* bnd_ptr = dynamic_cast<BoundConstraint*>(_constraints[i]);

    if(bnd_ptr)
      bound_constraints_.push_back(bnd_ptr);
    else
      constraints_.push_back(_constraints[i]);
  }
}


//-----------------------------------------------------------------------------


void
NProblemIPOPT::
analyze_special_properties(const NProblemInterface* _problem, const std::vector<NConstraintInterface*>& _constraints)
{
  hessian_constant_ = true;
  jac_c_constant_   = true;
  jac_d_constant_   = true;

  if(!_problem->constant_hessian())
    hessian_constant_ = false;

  for(unsigned int i=0; i<_constraints.size(); ++i)
  {
    if(!_constraints[i]->constant_hessian())
      hessian_constant_ = false;

    if(!_constraints[i]->constant_gradient())
    {
      if(_constraints[i]->constraint_type() == NConstraintInterface::NC_EQUAL)
        jac_c_constant_ = false;
      else
        jac_d_constant_ = false;
    }

    // nothing else to check?
    if(!hessian_constant_ && !jac_c_constant_ && !jac_d_constant_)
      break;
  }

  //hessian of Lagrangian is only constant, if all hessians of the constraints are zero (due to lambda multipliers)
  if(!jac_c_constant_ || !jac_d_constant_)
    hessian_constant_ = false;
}


//-----------------------------------------------------------------------------


bool NProblemIPOPT::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                         Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // number of variables
  n = problem_->n_unknowns();

  // number of constraints
  m = constraints_.size();

  // get non-zeros of hessian of lagrangian and jacobi of constraints
  nnz_jac_g = 0;
  nnz_h_lag = 0;

  // get nonzero structure
  std::vector<double> x(n);
  problem_->initial_x(P(x));

  // nonzeros in the jacobian of C_ and the hessian of the lagrangian
  SMatrixNP HP;
  SVectorNC g;
  SMatrixNC H;
  problem_->eval_hessian(P(x), HP);

  // get nonzero structure of hessian of problem
  for(int i=0; i<HP.outerSize(); ++i)
    for (SMatrixNP::InnerIterator it(HP,i); it; ++it)
      if(it.row() >= it.col())
        ++nnz_h_lag;

  // get nonzero structure of constraints
  for( int i=0; i<m; ++i)
  {
    constraints_[i]->eval_gradient(P(x),g);

    nnz_jac_g += g.nonZeros();

    // count lower triangular elements
    constraints_[i]->eval_hessian (P(x),H);

    SMatrixNC::iterator m_it = H.begin();
    for(; m_it != H.end(); ++m_it)
      if( m_it.row() >= m_it.col())
        ++nnz_h_lag;
  }

  // We use the standard fortran index style for row/col entries
  index_style = C_STYLE;

  return true;
}


//-----------------------------------------------------------------------------


bool NProblemIPOPT::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{
  DEB_enter_func;
  // check dimensions
  DEB_warning_if(( n != (Index)problem_->n_unknowns() ), 1,
    "IPOPT #unknowns != n " << n << problem_->n_unknowns() );
  DEB_warning_if(( m != (Index)constraints_.size() ), 1, 
    "Warning: IPOPT #constraints != m " << m << constraints_.size() );


  // first clear all variable bounds
  for( int i=0; i<n; ++i)
  {
    // x_l[i] = Ipopt::nlp_lower_bound_inf;
    // x_u[i] = Ipopt::nlp_upper_bound_inf;

    x_l[i] = -1.0e19;
    x_u[i] =  1.0e19;
  }

  // iterate over bound constraints and set them
  for(unsigned int i=0; i<bound_constraints_.size(); ++i)
  {
    if((Index)(bound_constraints_[i]->idx()) < n)
    {
      switch(bound_constraints_[i]->constraint_type())
      {
      case NConstraintInterface::NC_LESS_EQUAL:
      {
        x_u[bound_constraints_[i]->idx()] = bound_constraints_[i]->bound();
      }break;

      case NConstraintInterface::NC_GREATER_EQUAL:
      {
        x_l[bound_constraints_[i]->idx()] = bound_constraints_[i]->bound();
      }break;

      case NConstraintInterface::NC_EQUAL:
      {
        x_l[bound_constraints_[i]->idx()] = bound_constraints_[i]->bound();
        x_u[bound_constraints_[i]->idx()] = bound_constraints_[i]->bound();
      }break;
      }
    }
    else
      DEB_warning(2, "invalid bound constraint in IPOPTSolver!!!")
  }

  // set bounds for constraints
  for( int i=0; i<m; ++i)
  {
    // enum ConstraintType {NC_EQUAL, NC_LESS_EQUAL, NC_GREATER_EQUAL};
    switch(constraints_[i]->constraint_type())
    {
      case NConstraintInterface::NC_EQUAL         : g_u[i] = 0.0   ; g_l[i] =  0.0   ; break;
      case NConstraintInterface::NC_LESS_EQUAL    : g_u[i] = 0.0   ; g_l[i] = -1.0e19; break;
      case NConstraintInterface::NC_GREATER_EQUAL : g_u[i] = 1.0e19; g_l[i] =  0.0   ; break;
      default                                     : g_u[i] = 1.0e19; g_l[i] = -1.0e19; break;
    }
  }

  return true;
}


//-----------------------------------------------------------------------------


bool NProblemIPOPT::get_starting_point(Index n, bool init_x, Number* x,
                               bool init_z, Number* z_L, Number* z_U,
                               Index m, bool init_lambda,
                               Number* lambda)
{
  // get initial value of problem instance
  problem_->initial_x(x);

  return true;
}


//-----------------------------------------------------------------------------


bool NProblemIPOPT::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  // return the value of the objective function
  obj_value = problem_->eval_f(x);
  return true;
}


//-----------------------------------------------------------------------------


bool NProblemIPOPT::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  problem_->eval_gradient(x, grad_f);

  return true;
}


//-----------------------------------------------------------------------------


bool NProblemIPOPT::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  // evaluate all constraint functions
  for( int i=0; i<m; ++i)
    g[i] = constraints_[i]->eval_constraint(x);

  return true;
}


//-----------------------------------------------------------------------------


bool NProblemIPOPT::eval_jac_g(Index n, const Number* x, bool new_x,
                       Index m, Index nele_jac, Index* iRow, Index *jCol,
                       Number* values)
{
  DEB_enter_func;
  if (values == NULL)
  {
    // get x for evaluation (arbitrary position should be ok)
    std::vector<double> x_rnd(problem_->n_unknowns(), 0.0);

    int gi = 0;
    SVectorNC g;
    for( int i=0; i<m; ++i)
    {
      constraints_[i]->eval_gradient(&(x_rnd[0]), g);
      SVectorNC::InnerIterator v_it(g);
      for( ; v_it; ++v_it)
      {
        iRow[gi] = i;
        jCol[gi] = v_it.index();
        ++gi;
      }
    }
  }
  else
  {
    // return the values of the jacobian of the constraints

    // return the structure of the jacobian of the constraints
    // global index
    int gi = 0;
    SVectorNC g;

    for( int i=0; i<m; ++i)
    {
      constraints_[i]->eval_gradient(x, g);

      SVectorNC::InnerIterator v_it(g);

      for( ; v_it; ++v_it)
      {
        values[gi] = v_it.value();
        ++gi;
      }
    }

    DEB_warning_if((gi != nele_jac), 1, 
      "number of non-zeros in Jacobian of C is incorrect: "
                << gi << " vs " << nele_jac)
  }

  return true;
}


//-----------------------------------------------------------------------------


bool NProblemIPOPT::eval_h(Index n, const Number* x, bool new_x,
                   Number obj_factor, Index m, const Number* lambda,
                   bool new_lambda, Index nele_hess, Index* iRow,
                   Index* jCol, Number* values)
{
  DEB_enter_func;
  if (values == NULL)
  {
    // return structure

    // get x for evaluation (arbitrary position should be ok)
    std::vector<double> x_rnd(problem_->n_unknowns(), 0.0);

     // global index
     int gi = 0;
     // get hessian of problem
     SMatrixNP HP;
     problem_->eval_hessian(&(x_rnd[0]), HP);

     for(int i=0; i<HP.outerSize(); ++i)
       for (SMatrixNP::InnerIterator it(HP,i); it; ++it)
       {
         // store lower triangular part only
         if(it.row() >= it.col())
         {
           //         it.value();
           iRow[gi] = it.row();
           jCol[gi] = it.col();
           ++gi;
         }
       }

    // Hessians of Constraints
    for(unsigned int j=0; j<constraints_.size(); ++j)
    {
      SMatrixNC H;
      constraints_[j]->eval_hessian(&(x_rnd[0]), H);

      SMatrixNC::iterator m_it  = H.begin();
      SMatrixNC::iterator m_end = H.end();

      for(; m_it != m_end; ++m_it)
      {
        // store lower triangular part only
        if( m_it.row() >= m_it.col())
        {
          iRow[gi] = m_it.row();
          jCol[gi] = m_it.col();
          ++gi;
        }
      }
    }

    // error check
    DEB_warning_if(( gi != nele_hess), 1,
      "number of non-zeros in Hessian of Lagrangian is incorrect while indexing: "
                << gi << " vs " << nele_hess )
  }
  else
  {
    // return values.

    // global index
    int gi = 0;
    // get hessian of problem
    SMatrixNP HP;
    problem_->eval_hessian(x, HP);

    for(int i=0; i<HP.outerSize(); ++i)
      for (SMatrixNP::InnerIterator it(HP,i); it; ++it)
      {
        // store lower triangular part only
        if(it.row() >= it.col())
        {
          values[gi] = obj_factor*it.value();
          ++gi;
        }
      }

    // Hessians of Constraints
    for(unsigned int j=0; j<constraints_.size(); ++j)
    {
      SMatrixNC H;
      constraints_[j]->eval_hessian(x, H);

      SMatrixNC::iterator m_it  = H.begin();
      SMatrixNC::iterator m_end = H.end();

      for(; m_it != m_end; ++m_it)
      {
        // store lower triangular part only
        if( m_it.row() >= m_it.col())
        {
          values[gi] = lambda[j]*(*m_it);
          ++gi;
        }
      }
    }

    // error check
    DEB_warning_if(( gi != nele_hess), 1,
      "number of non-zeros in Hessian of Lagrangian is incorrect2: "
                << gi << " vs " << nele_hess )
  }
  return true;
}


//-----------------------------------------------------------------------------


//inline double _QNT(double x) 
//{
//	return double(float(x));
//}

//double _QNT(const double x)
//    {
//    // clear the 12 least significant mantissa bits to reduce noise
//    //const double fact = pow(2., 41); 
//
//	const double fact = pow(2., 37); 
//    int i;
//    double m = frexp(x, &i);
//    m *= fact;
//    int sgn_x = m < 0 ? -1 : 1;
//    m = sgn_x * floor(fabs(m));
//    m /= fact;
//    double xq = ldexp(m, i);
//    return xq;
//    }

double _QNT(const double x) { return x; }

void NProblemIPOPT::finalize_solution(SolverReturn status,
                              Index n, const Number* x, const Number* z_L, const Number* z_U,
                              Index m, const Number* g, const Number* lambda,
                              Number obj_value,
                              const IpoptData* ip_data,
                              IpoptCalculatedQuantities* ip_cq)
{
  DEB_enter_func;
	DEB_out(1, "Quantanizing the IPOPT solution\n");
	std::vector<double> x_qnt(n);
    for( Index i=0; i<n; ++i)
      x_qnt[i] = _QNT(x[i]);



  // problem knows what to do
  problem_->store_result(&x_qnt[0]);

  if(store_solution_)
  {
    x_.resize(n);
    for( Index i=0; i<n; ++i)
      x_[i] = x_qnt[i];
  }
}


//-----------------------------------------------------------------------------


bool NProblemIPOPT::hessian_constant() const
{
  return hessian_constant_;
}


//-----------------------------------------------------------------------------


bool NProblemIPOPT::jac_c_constant() const
{
  return jac_c_constant_;
}


//-----------------------------------------------------------------------------


bool NProblemIPOPT::jac_d_constant() const
{
  return jac_d_constant_;
}


//== IMPLEMENTATION PROBLEM INSTANCE==========================================================


bool NProblemGmmIPOPT::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                         Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // number of variables
  n = problem_->n_unknowns();

  // number of constraints
  m = constraints_.size();

  // get nonzero structure
  std::vector<double> x(n);
  problem_->initial_x(&(x[0]));
  // ToDo: perturb x

  // nonzeros in the jacobian of C_ and the hessian of the lagrangian
  SMatrixNP HP;
  SVectorNC g;
  SMatrixNC H;
  problem_->eval_hessian(&(x[0]), HP);
  nnz_jac_g = 0;
  nnz_h_lag = 0;

  // clear old data
  jac_g_iRow_.clear();
  jac_g_jCol_.clear();
  h_lag_iRow_.clear();
  h_lag_jCol_.clear();

  // get non-zero structure of initial hessian
  // iterate over rows
  for( int i=0; i<n; ++i)
  {
    SVectorNP& ri = HP.row(i);

    SVectorNP_citer v_it  = gmm::vect_const_begin(ri);
    SVectorNP_citer v_end = gmm::vect_const_end  (ri);

    for(; v_it != v_end; ++v_it)
    {
      // store lower triangular part only
      if( i >= (int)v_it.index())
      {
        h_lag_iRow_.push_back(i);
        h_lag_jCol_.push_back(v_it.index());
        ++nnz_h_lag;
      }
    }
  }


  // get nonzero structure of constraints
  for( int i=0; i<m; ++i)
  {
    constraints_[i]->eval_gradient(&(x[0]),g);
    constraints_[i]->eval_hessian (&(x[0]),H);

    // iterate over sparse vector
    SVectorNC::InnerIterator v_it(g);
    for(; v_it; ++v_it)
    {
      jac_g_iRow_.push_back(i);
      jac_g_jCol_.push_back(v_it.index());
      ++nnz_jac_g;
    }

    // iterate over superSparseMatrix
    SMatrixNC::iterator m_it  = H.begin();
    SMatrixNC::iterator m_end = H.end();
    for(; m_it != m_end; ++m_it)
      if( m_it.row() >= m_it.col())
      {
        h_lag_iRow_.push_back(m_it.row());
        h_lag_jCol_.push_back(m_it.col());
        ++nnz_h_lag;
      }
  }

  // store for error checking...
  nnz_jac_g_ = nnz_jac_g;
  nnz_h_lag_ = nnz_h_lag;

  // We use the standard fortran index style for row/col entries
  index_style = C_STYLE;

  return true;
}


//-----------------------------------------------------------------------------


bool NProblemGmmIPOPT::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{
  // first clear all variable bounds
  for( int i=0; i<n; ++i)
  {
    // x_l[i] = Ipopt::nlp_lower_bound_inf;
    // x_u[i] = Ipopt::nlp_upper_bound_inf;

    x_l[i] = -1.0e19;
    x_u[i] =  1.0e19;
  }

  // set bounds for constraints
  for( int i=0; i<m; ++i)
  {
    // enum ConstraintType {NC_EQUAL, NC_LESS_EQUAL, NC_GREATER_EQUAL};
    switch(constraints_[i]->constraint_type())
    {
      case NConstraintInterface::NC_EQUAL         : g_u[i] = 0.0   ; g_l[i] =  0.0   ; break;
      case NConstraintInterface::NC_LESS_EQUAL    : g_u[i] = 0.0   ; g_l[i] = -1.0e19; break;
      case NConstraintInterface::NC_GREATER_EQUAL : g_u[i] = 1.0e19; g_l[i] =  0.0   ; break;
      default                                     : g_u[i] = 1.0e19; g_l[i] = -1.0e19; break;
    }
  }

  return true;
}


//-----------------------------------------------------------------------------


bool NProblemGmmIPOPT::get_starting_point(Index n, bool init_x, Number* x,
                               bool init_z, Number* z_L, Number* z_U,
                               Index m, bool init_lambda,
                               Number* lambda)
{
  // get initial value of problem instance
  problem_->initial_x(x);

  return true;
}


//-----------------------------------------------------------------------------


bool NProblemGmmIPOPT::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  // return the value of the objective function
  obj_value = problem_->eval_f(x);
  return true;
}


//-----------------------------------------------------------------------------


bool NProblemGmmIPOPT::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  problem_->eval_gradient(x, grad_f);

  return true;
}


//-----------------------------------------------------------------------------


bool NProblemGmmIPOPT::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  // evaluate all constraint functions
  for( int i=0; i<m; ++i)
    g[i] = constraints_[i]->eval_constraint(x);

  return true;
}


//-----------------------------------------------------------------------------


bool NProblemGmmIPOPT::eval_jac_g(Index n, const Number* x, bool new_x,
                       Index m, Index nele_jac, Index* iRow, Index *jCol,
                       Number* values)
{
  DEB_enter_func;
  if (values == NULL)
  {
    // return the (cached) structure of the jacobian of the constraints
    gmm::copy(jac_g_iRow_, VectorPTi(iRow, jac_g_iRow_.size()));
    gmm::copy(jac_g_jCol_, VectorPTi(jCol, jac_g_jCol_.size()));
  }
  else
  {
    // return the values of the jacobian of the constraints

    // return the structure of the jacobian of the constraints
    // global index
    int gi = 0;
    SVectorNC g;

    for( int i=0; i<m; ++i)
    {
      constraints_[i]->eval_gradient(x, g);

      SVectorNC::InnerIterator v_it(g);

      for( ; v_it; ++v_it)
      {
        if(gi < nele_jac)
          values[gi] = v_it.value();
        ++gi;
      }
    }

    DEB_warning_if(( gi != nele_jac), 1,
      "number of non-zeros in Jacobian of C is incorrect: "
       << gi << " vs " << nele_jac)
  }

  return true;
}


//-----------------------------------------------------------------------------


bool NProblemGmmIPOPT::eval_h(Index n, const Number* x, bool new_x,
                   Number obj_factor, Index m, const Number* lambda,
                   bool new_lambda, Index nele_hess, Index* iRow,
                   Index* jCol, Number* values)
{
  DEB_enter_func;
  if (values == NULL)
  {
    // return the (cached) structure of the hessian
    gmm::copy(h_lag_iRow_, VectorPTi(iRow, h_lag_iRow_.size()));
    gmm::copy(h_lag_jCol_, VectorPTi(jCol, h_lag_jCol_.size()));
  }
  else
  {
    // return values.

    // global index
    int gi = 0;

    // get hessian of problem
    problem_->eval_hessian(x, HP_);

    for( int i=0; i<n; ++i)
    {
      SVectorNP& ri = HP_.row(i);

      SVectorNP_citer v_it  = gmm::vect_const_begin(ri);
      SVectorNP_citer v_end = gmm::vect_const_end  (ri);

      for(; v_it != v_end; ++v_it)
      {
        // store lower triangular part only
        if( i >= (int)v_it.index())
        {
          if( gi < nele_hess)
            values[gi] = obj_factor*(*v_it);
          ++gi;
        }
      }
    }

    // Hessians of Constraints
    for(unsigned int j=0; j<constraints_.size(); ++j)
    {
      SMatrixNC H;
      constraints_[j]->eval_hessian(x, H);

      SMatrixNC::iterator m_it  = H.begin();
      SMatrixNC::iterator m_end = H.end();

      for(; m_it != m_end; ++m_it)
      {
        // store lower triangular part only
        if( m_it.row() >= m_it.col())
        {
          if( gi < nele_hess)
            values[gi] = lambda[j]*(*m_it);
          ++gi;
        }
      }
    }

    // error check
    DEB_warning_if(( gi != nele_hess), 1, 
      "number of non-zeros in Hessian of Lagrangian is incorrect: "
        << gi << " vs " << nele_hess);
  }
  return true;
}


//-----------------------------------------------------------------------------


void NProblemGmmIPOPT::finalize_solution(SolverReturn status,
                              Index n, const Number* x, const Number* z_L, const Number* z_U,
                              Index m, const Number* g, const Number* lambda,
                              Number obj_value,
                              const IpoptData* ip_data,
                              IpoptCalculatedQuantities* ip_cq)
{
  // problem knows what to do
  problem_->store_result(x);
}



//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // COMISO_IPOPT_AVAILABLE
//=============================================================================
