//=============================================================================
//
//  CLASS GUROBISolver - IMPLEMENTATION
//
//=============================================================================

//== INCLUDES =================================================================

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_GUROBI_AVAILABLE
//=============================================================================
 #include <gurobi_c++.h>


#include "GUROBISolver.hh"
#include <CoMISo/Utils/StopWatch.hh>

#include <Base/Utils/OutcomeUtils.hh>
#include <stdexcept>

#include <Base/Debug/DebOut.hh>

DEB_module("NSLV")


//== NAMESPACES ===============================================================

namespace COMISO {

//== IMPLEMENTATION ========================================================== 



GUROBISolver::
GUROBISolver()
{
}

//-----------------------------------------------------------------------------
//#define TRACE_GUROBI_INPUT(DESCR, EXPR) \
//	std::cout << DESCR << ": " << EXPR << std::endl

#define TRACE_GUROBI_INPUT(DESCR, EXPR) DEB_out(7, DESCR << ":" << EXPR << "\n")

#define _QNT(X) X
//inline double _QNT(double x) 
//{
//	return double(float(x));
//}

//#define TRACE_GUROBI_INPUT(DESCR, EXPR)

// ********** SOLVE **************** //
void GUROBISolver::solve(
  NProblemInterface*                  _problem,
  std::vector<NConstraintInterface*>& _constraints,
  std::vector<PairIndexVtype>&        _discrete_constraints,
  const double                        /*_time_limit*/)
{
  DEB_enter_func;
  try
  {
    //----------------------------------------------
    // 0. set up environment
    //----------------------------------------------

    GRBEnv   env   = GRBEnv();
    GRBModel model = GRBModel(env);

    //model.getEnv().set(GRB_DoubleParam_TimeLimit, _time_limit);

    // stop when a solution is found
    model.getEnv().set(GRB_IntParam_SolutionLimit, 2);

    //----------------------------------------------
    // 1. allocate variables
    //----------------------------------------------

    // determine variable types: 0->real, 1->integer, 2->bool
    std::vector<char> vtypes(_problem->n_unknowns(),0);
    for(unsigned int i=0; i<_discrete_constraints.size(); ++i)
      switch(_discrete_constraints[i].second)
      {
        case Integer: vtypes[_discrete_constraints[i].first] = 1; break;
        case Binary : vtypes[_discrete_constraints[i].first] = 2; break;
        default     : break;
      }

    // GUROBI variables
    std::vector<GRBVar> vars;
    // first all
    for( int i=0; i<_problem->n_unknowns(); ++i)
      switch(vtypes[i])
      {
        case 0 : vars.push_back( model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS) ); break;
        case 1 : vars.push_back( model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER   ) ); break;
        case 2 : vars.push_back( model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_BINARY    ) ); break;
      }


    // Integrate new variables
    model.update();

    //----------------------------------------------
    // 2. setup constraints
    //----------------------------------------------

    // get zero vector
    std::vector<double> x(_problem->n_unknowns(), 0.0);

    for(unsigned int i=0; i<_constraints.size();  ++i)
    {
      TRACE_GUROBI_INPUT("Constraint index: ", i);
      DEB_warning_if(!_constraints[i]->is_linear(),1,
        "GUROBISolver received a problem with non-linear constraints!!!");

      GRBLinExpr lin_expr;
      NConstraintInterface::SVectorNC gc;
      _constraints[i]->eval_gradient(P(x), gc);

      TRACE_GUROBI_INPUT("Constraint type: ", 
        (int)_constraints[i]->constraint_type());

      NConstraintInterface::SVectorNC::InnerIterator v_it(gc);
      for(; v_it; ++v_it)
      {
        //        lin_expr += v_it.value()*vars[v_it.index()];
        lin_expr += vars[v_it.index()]*_QNT(v_it.value());
        TRACE_GUROBI_INPUT("Constraint linear coeff", _QNT(v_it.value())); 
        TRACE_GUROBI_INPUT("Constraint linear var idx", v_it.index()); 
      }

      double b = _QNT(_constraints[i]->eval_constraint(P(x)));
      TRACE_GUROBI_INPUT("Constraint B", b); 

      switch(_constraints[i]->constraint_type())
      {
      case NConstraintInterface::NC_EQUAL         : model.addConstr(lin_expr + b == 0); break;
      case NConstraintInterface::NC_LESS_EQUAL    : model.addConstr(lin_expr + b <= 0); break;
      case NConstraintInterface::NC_GREATER_EQUAL : model.addConstr(lin_expr + b >= 0); break;
      }
    }
    model.update();

    //----------------------------------------------
    // 3. setup energy
    //----------------------------------------------

    DEB_warning_if(!_problem->constant_hessian(), 1, 
      "GUROBISolver received a problem with non-constant hessian!!!");

    GRBQuadExpr objective;

    // add quadratic part
    NProblemInterface::SMatrixNP H;
    _problem->eval_hessian(P(x), H);
    for( int i=0; i<H.outerSize(); ++i)
    {
      for (NProblemInterface::SMatrixNP::InnerIterator it(H,i); it; ++it)
      {
        objective += 0.5*_QNT(it.value())*vars[it.row()]*vars[it.col()];
        TRACE_GUROBI_INPUT("Objective quadratic term coeff", _QNT(it.value())); 
        TRACE_GUROBI_INPUT("Objective quadratic term var i idx", it.row()); 
        TRACE_GUROBI_INPUT("Objective quadratic term var j idx", it.col()); 
      }
    }

    // add linear part
    std::vector<double> g(_problem->n_unknowns());
    _problem->eval_gradient(P(x), P(g));
    for(unsigned int i=0; i<g.size(); ++i)
    {
      objective += _QNT(g[i])*vars[i];
      TRACE_GUROBI_INPUT("Objective linear term coeff", _QNT(g[i])); 
      TRACE_GUROBI_INPUT("Objective linear term var idx", i); 
    }
    // add constant part
    const double c = _QNT(_problem->eval_f(P(x)));
    objective += c;
    TRACE_GUROBI_INPUT("Objective constant term coeff", c); 

    model.set(GRB_IntAttr_ModelSense, 1);
    model.setObjective(objective);
    model.update();


    //----------------------------------------------
    // 4. solve problem
    //----------------------------------------------


    if (solution_input_path_.empty())
    {
      //if (!problem_env_output_path_.empty())
      //{
      //  std::cout << "Writing problem's environment into file \"" << problem_env_output_path_ << "\"." << std::endl;
      //  model.getEnv().writeParams(problem_env_output_path_);
      //}
      //if (!problem_output_path_.empty())
      //{
      //  std::cout << "Writing problem into file \"" << problem_output_path_ << "\"." << std::endl;
      //  GurobiHelper::outputModelToMpsGz(model, problem_output_path_);
      //}

      model.optimize();
    }
    else
    {
      DEB_out(2, "Reading solution from file \"" << solution_input_path_ << "\"\n")
    }

    //----------------------------------------------
    // 5. store result
    //----------------------------------------------

    if (solution_input_path_.empty())
    {
      // store computed result
      for(unsigned int i=0; i<vars.size(); ++i)
        x[i] = vars[i].get(GRB_DoubleAttr_X);
    }
    //else
    //{
    //    std::cout << "Loading stored solution from \"" << solution_input_path_ << "\"." << std::endl;
    //    // store loaded result
    //    const size_t oldSize = x.size();
    //    x.clear();
    //    GurobiHelper::readSolutionVectorFromSOL(x, solution_input_path_);
    //    if (oldSize != x.size()) {
    //        std::cerr << "oldSize != x.size() <=> " << oldSize << " != " << x.size() << std::endl;
    //        throw std::runtime_error("Loaded solution vector doesn't have expected dimension.");
    //    }
    //}

    _problem->store_result(P(x));

    // ObjVal is only available if the optimize was called.
    DEB_out_if(solution_input_path_.empty(), 2,
        "GUROBI Objective: " << model.get(GRB_DoubleAttr_ObjVal) << "\n");
  }
  catch (GRBException& e)
  {
    DEB_warning(1, "Error code = " << e.getErrorCode() << "[" << e.getMessage() << "]\n");
    
    // NOTE: we could propagate e.getMessage() either using std::exception, or a specialized Reform exception type
    switch ( e.getErrorCode() )
    {
    // TODO: Find where these are defined and use the correct name for 10010 
    case 10010: THROW_OUTCOME(GUROBI_LICENCE_MODEL_TOO_LARGE);
    default: THROW_OUTCOME(UNSPECIFIED_GUROBI_EXCEPTION);
    }
  }
}


//-----------------------------------------------------------------------------

#if 0
bool
GUROBISolver::
solve(NProblemInterface*                        _problem,
      const std::vector<NConstraintInterface*>& _constraints,
      const std::vector<NConstraintInterface*>& _lazy_constraints,
      std::vector<PairIndexVtype>&              _discrete_constraints,   // discrete constraints
      const double                              _almost_infeasible,
      const int                                 _max_passes,
      const double                              _time_limit,
      const bool                                _silent)
{
  DEB_enter_func;
//  // hack! solve with all constraints
//  std::vector<NConstraintInterface*> all_constraints;
//  std::copy(_constraints.begin(),_constraints.end(),std::back_inserter(all_constraints));
//  std::copy(_lazy_constraints.begin(),_lazy_constraints.end(),std::back_inserter(all_constraints));
//
//  return solve(_problem, all_constraints, _discrete_constraints, _time_limit);
  StopWatch sw; sw.start();

  bool feasible_point_found = false;
  int  cur_pass = 0;
  double acceptable_tolerance = 0.01; // hack: read out from ipopt!!!
  std::vector<bool> lazy_added(_lazy_constraints.size(),false);

  // cache statistics of all iterations
  std::vector<int> n_inf;
  std::vector<int> n_almost_inf;

  try
  {
    //----------------------------------------------
    // 0. set up environment
    //----------------------------------------------

    GRBEnv   env   = GRBEnv();
    GRBModel model = GRBModel(env);

    model.getEnv().set(GRB_DoubleParam_TimeLimit, _time_limit);


    //----------------------------------------------
    // 1. allocate variables
    //----------------------------------------------

    // determine variable types: 0->real, 1->integer, 2->bool
    std::vector<char> vtypes(_problem->n_unknowns(),0);
    for(unsigned int i=0; i<_discrete_constraints.size(); ++i)
      switch(_discrete_constraints[i].second)
      {
        case Integer: vtypes[_discrete_constraints[i].first] = 1; break;
        case Binary : vtypes[_discrete_constraints[i].first] = 2; break;
        default     : break;
      }

    // GUROBI variables
    std::vector<GRBVar> vars;
    // first all
    for( int i=0; i<_problem->n_unknowns(); ++i)
      switch(vtypes[i])
      {
        case 0 : vars.push_back( model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS) ); break;
        case 1 : vars.push_back( model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER   ) ); break;
        case 2 : vars.push_back( model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_BINARY    ) ); break;
      }


    // Integrate new variables
    model.update();

    //----------------------------------------------
    // 2. setup constraints
    //----------------------------------------------

    // get zero vector
    std::vector<double> x(_problem->n_unknowns(), 0.0);

    for(unsigned int i=0; i<_constraints.size();  ++i)
    {
      add_constraint_to_model(_constraints[i], model, vars, P(x));
    }
    model.update();

    //----------------------------------------------
    // 3. setup energy
    //----------------------------------------------

    DEB_warning_if(!_problem->constant_hessian(), 1,
      "GUROBISolver received a problem with non-constant hessian!!!\n");

    GRBQuadExpr objective;

    // add quadratic part
    NProblemInterface::SMatrixNP H;
    _problem->eval_hessian(P(x), H);
    for( int i=0; i<H.outerSize(); ++i)
      for (NProblemInterface::SMatrixNP::InnerIterator it(H,i); it; ++it)
        objective += 0.5*it.value()*vars[it.row()]*vars[it.col()];


    // add linear part
    std::vector<double> g(_problem->n_unknowns());
    _problem->eval_gradient(P(x), P(g));
    for(unsigned int i=0; i<g.size(); ++i)
      objective += g[i]*vars[i];

    // add constant part
    objective += _problem->eval_f(P(x));

    model.set(GRB_IntAttr_ModelSense, 1);
    model.setObjective(objective);


    //----------------------------------------------
    // 4. iteratively solve problem
    //----------------------------------------------
    bool solution_found = false;

    while(!feasible_point_found && cur_pass <(_max_passes-1))
    {
      ++cur_pass;
      //----------------------------------------------------------------------------
      // 1. optimize current Model and get result
      //----------------------------------------------------------------------------

      model.update();
      model.optimize();
      // store computed result
      if(model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
        for(unsigned int i=0; i<vars.size(); ++i)
          x[i] = vars[i].get(GRB_DoubleAttr_X);

      //----------------------------------------------------------------------------
      // 2. Check lazy constraints
      //----------------------------------------------------------------------------
      // check lazy constraints
      n_inf.push_back(0);
      n_almost_inf.push_back(0);
      feasible_point_found = true;
      for(unsigned int i=0; i<_lazy_constraints.size(); ++i)
        if(!lazy_added[i])
        {
          NConstraintInterface* lc = _lazy_constraints[i];

          double v = lc->eval_constraint(P(x));

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
                add_constraint_to_model( lc, model, vars, P(x));
                lazy_added[i] = true;
                feasible_point_found = false;
                ++n_inf.back();
              }
              else // almost violated or violated? -> add to constraints
                if(almost_inf)
                {
                  add_constraint_to_model( lc, model, vars, P(x));
                  lazy_added[i] = true;
                  ++n_almost_inf.back();
                }
        }
    }

    // no termination after max number of passes?
    if(!feasible_point_found)
    {
      ++cur_pass;

      DEB_out(1, "*************** could not find feasible point after " 
        << _max_passes-1 << " -> solving with all lazy constraints...\n");
      for(unsigned int i=0; i<_lazy_constraints.size(); ++i)
        if(!lazy_added[i])
        {
          add_constraint_to_model( _lazy_constraints[i], model, vars, P(x));
        }

      model.update();
      model.optimize();

      // store computed result
      if(model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
        for(unsigned int i=0; i<vars.size(); ++i)
          x[i] = vars[i].get(GRB_DoubleAttr_X);
    }

    const double overall_time = sw.stop()/1000.0;

    //----------------------------------------------------------------------------
    // 4. output statistics
    //----------------------------------------------------------------------------
    // Make this code DEB_only
    DEB_out(2,"############# GUROBI with lazy constraints statistics ###############\n");
    DEB_out(2, "overall time: " << overall_time << "s\n");
    DEB_out(2, "#passes     : " << cur_pass << "( of " << _max_passes << ")\n");
    for(unsigned int i=0; i<n_inf.size(); ++i)
      DEB_out(2, "pass " << i << " induced " << n_inf[i]
        << " infeasible and " << n_almost_inf[i] << " almost infeasible\n");
    DEB_out(2, "GUROBI Objective: " << model.get(GRB_DoubleAttr_ObjVal) << "\n");

    //----------------------------------------------
    // 5. store result
    //----------------------------------------------
    if(model.get(GRB_IntAttr_Status) != GRB_OPTIMAL)
      return false;
    else
    {
      _problem->store_result(P(x));
      return true;
    }
  }
  catch(GRBException& e)
  {
    DEB_warning(1, 
      "Error code = " << e.getErrorCode() << "[" << e.getMessage() << "\n");
    return false;
  }
  catch(...)
    {
    DEB_warning)1,"Exception during optimization\n");
    return false;
  }

  return false;

//    //----------------------------------------------
//    // 4. iteratively solve problem
//    //----------------------------------------------
//    IloBool solution_found = IloFalse;
//
//    while(!feasible_point_found && cur_pass <(_max_passes-1))
//    {
//      ++cur_pass;
//      //----------------------------------------------------------------------------
//      // 1. Create an instance of current Model and optimize
//      //----------------------------------------------------------------------------
//
//      if(!_silent)
//        std::cerr << "cplex -> setup IloCPlex...\n";
//      IloCplex cplex(model);
//      cplex.setParam(IloCplex::TiLim, _time_limit);
//      // silent mode?
//      if(_silent)
//        cplex.setOut(env_.getNullStream());
//
//      if(!_silent)
//        std::cerr << "cplex -> optimize...\n";
//      solution_found = cplex.solve();
//
//      if(solution_found != IloFalse)
//      {
//        for(unsigned int i=0; i<vars.size(); ++i)
//          x[i] = cplex.getValue(vars[i]);
//
//        _problem->store_result(P(x));
//      }
//      else continue;
//
//      //----------------------------------------------------------------------------
//      // 2. Check lazy constraints
//      //----------------------------------------------------------------------------
//      // check lazy constraints
//      n_inf.push_back(0);
//      n_almost_inf.push_back(0);
//      feasible_point_found = true;
//      for(unsigned int i=0; i<_lazy_constraints.size(); ++i)
//        if(!lazy_added[i])
//        {
//          NConstraintInterface* lc = _lazy_constraints[i];
//
//          double v = lc->eval_constraint(P(x));
//
//          bool inf        = false;
//          bool almost_inf = false;
//
//          if(lc->constraint_type() == NConstraintInterface::NC_EQUAL)
//          {
//            v = std::abs(v);
//            if(v>acceptable_tolerance)
//              inf = true;
//            else
//              if(v>_almost_infeasible)
//                almost_inf = true;
//          }
//          else
//            if(lc->constraint_type() == NConstraintInterface::NC_GREATER_EQUAL)
//            {
//              if(v<-acceptable_tolerance)
//                inf = true;
//              else
//                if(v<_almost_infeasible)
//                  almost_inf = true;
//            }
//            else
//              if(lc->constraint_type() == NConstraintInterface::NC_LESS_EQUAL)
//              {
//                if(v>acceptable_tolerance)
//                  inf = true;
//                else
//                  if(v>-_almost_infeasible)
//                    almost_inf = true;
//              }
//
//          // infeasible?
//          if(inf)
//          {
//            add_constraint_to_model( lc, vars, model, env_, P(x));
//            lazy_added[i] = true;
//            feasible_point_found = false;
//            ++n_inf.back();
//          }
//          else // almost violated or violated? -> add to constraints
//            if(almost_inf)
//            {
//              add_constraint_to_model( lc, vars, model, env_, P(x));
//              lazy_added[i] = true;
//              ++n_almost_inf.back();
//            }
//        }
//    }
//
//    // no termination after max number of passes?
//    if(!feasible_point_found)
//    {
//      ++cur_pass;
//
//      std::cerr << "*************** could not find feasible point after " << _max_passes-1 << " -> solving with all lazy constraints..." << std::endl;
//      for(unsigned int i=0; i<_lazy_constraints.size(); ++i)
//        if(!lazy_added[i])
//        {
//          add_constraint_to_model( _lazy_constraints[i], vars, model, env_, P(x));
//        }
//
//      //----------------------------------------------------------------------------
//      // 1. Create an instance of current Model and optimize
//      //----------------------------------------------------------------------------
//
//      IloCplex cplex(model);
//      cplex.setParam(IloCplex::TiLim, _time_limit);
//      // silent mode?
//      if(_silent)
//        cplex.setOut(env_.getNullStream());
//      solution_found = cplex.solve();
//
//      if(solution_found != IloFalse)
//      {
//        for(unsigned int i=0; i<vars.size(); ++i)
//          x[i] = cplex.getValue(vars[i]);
//
//        _problem->store_result(P(x));
//      }
//    }
//
//    const double overall_time = sw.stop()/1000.0;
//
//    //----------------------------------------------------------------------------
//    // 4. output statistics
//    //----------------------------------------------------------------------------
////    if (solution_found != IloFalse)
////    {
////      // Retrieve some statistics about the solve
////      Ipopt::Index iter_count = app_->Statistics()->IterationCount();
////      printf("\n\n*** IPOPT: The problem solved in %d iterations!\n", iter_count);
////
////      Ipopt::Number final_obj = app_->Statistics()->FinalObjective();
////      printf("\n\n*** IPOPT: The final value of the objective function is %e.\n", final_obj);
////    }
//
//    std::cerr <<"############# CPLEX with lazy constraints statistics ###############" << std::endl;
//    std::cerr << "overall time: " << overall_time << "s" << std::endl;
//    std::cerr << "#passes     : " << cur_pass << "( of " << _max_passes << ")" << std::endl;
//    for(unsigned int i=0; i<n_inf.size(); ++i)
//      std::cerr << "pass " << i << " induced " << n_inf[i] << " infeasible and " << n_almost_inf[i] << " almost infeasible" << std::endl;
//
//    return (solution_found != IloFalse);
}

#endif//0

//-----------------------------------------------------------------------------


void
GUROBISolver::
add_constraint_to_model(COMISO::NConstraintInterface* _constraint, GRBModel& _model, std::vector<GRBVar>& _vars, double * _x)
{
  DEB_enter_func;
  DEB_warning_if(!_constraint->is_linear(), 1,
     "GUROBISolver received a problem with non-linear constraints!!!\n");

  GRBLinExpr lin_expr;
  NConstraintInterface::SVectorNC gc;
  _constraint->eval_gradient(_x, gc);

  NConstraintInterface::SVectorNC::InnerIterator v_it(gc);
  for(; v_it; ++v_it)
    lin_expr += _vars[v_it.index()]*v_it.value();

  double b = _constraint->eval_constraint(_x);

  switch(_constraint->constraint_type())
  {
    case NConstraintInterface::NC_EQUAL         : _model.addConstr(lin_expr + b == 0); break;
    case NConstraintInterface::NC_LESS_EQUAL    : _model.addConstr(lin_expr + b <= 0); break;
    case NConstraintInterface::NC_GREATER_EQUAL : _model.addConstr(lin_expr + b >= 0); break;
  }
}

//-----------------------------------------------------------------------------


void
GUROBISolver::
set_problem_output_path( const std::string &_problem_output_path)
{
  problem_output_path_ = _problem_output_path;
}


//-----------------------------------------------------------------------------


void
GUROBISolver::
set_problem_env_output_path( const std::string &_problem_env_output_path)
{
  problem_env_output_path_ = _problem_env_output_path;
}


//-----------------------------------------------------------------------------


void
GUROBISolver::
set_solution_input_path(const std::string &_solution_input_path)
{
  solution_input_path_ = _solution_input_path;
}


//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // COMISO_GUROBI_AVAILABLE
//=============================================================================

