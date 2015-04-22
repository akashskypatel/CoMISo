//=============================================================================
//
//  CLASS DOCloudSolver - IMPLEMENTATION
//
//=============================================================================

// TODO: this uses CBC for the MPS file export; consider implementing our own 
// MPS/LP export to remove the dependency.

//== INCLUDES =================================================================

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_DOCLOUD_AVAILABLE

//=============================================================================
#include "DOCloudSolver.hh"

#include <Base/Debug/DebTime.hh>
#include <Base/Utils/OutcomeUtils.hh>

// For Branch and bound
#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcStrategy.hpp"
#include "CbcHeuristicLocal.hpp"
#include "OsiClpSolverInterface.hpp"

#include <stdexcept>
#include <stdio.h>

DEB_module("DOCloudSolver")

#define CBC_INFINITY COIN_DBL_MAX

//== NAMESPACES ===============================================================

namespace COMISO {

//== IMPLEMENTATION ==========================================================
#define TRACE_CBT(DESCR, EXPR) DEB_line(7, DESCR << ": " << EXPR)
#define P(X) ((X).data())

namespace {

void solve_impl(
  NProblemInterface*                        _problem,
  const std::vector<NConstraintInterface*>& _constraints,
  const std::vector<PairIndexVtype>&        _discrete_constraints,
  const double                              /*_time_limit*/
)
{
  DEB_enter_func;
  DEB_warning_if(!_problem->constant_hessian(), 1,
    "DOCloudSolver received a problem with non-constant hessian!");
  DEB_warning_if(!_problem->constant_gradient(), 1,
    "DOCloudSolver received a problem with non-constant gradient!");

  const int n_rows = _constraints.size(); // Constraints #
  const int n_cols = _problem->n_unknowns(); // Unknowns #

  // expand the variable types from discrete mtrx array
  std::vector<VariableType> var_type(n_cols, Real);
  for (size_t i = 0; i < _discrete_constraints.size(); ++i)
    var_type[_discrete_constraints[i].first] = _discrete_constraints[i].second;

  //----------------------------------------------
  // set up constraints
  //----------------------------------------------
  std::vector<double> x(n_cols, 0.0); // solution

  CoinPackedMatrix mtrx(false, 0, 0);// make a row-ordered, empty mtrx
  mtrx.setDimensions(0, n_cols);

  std::vector<double> row_lb(n_rows, -CBC_INFINITY);
  std::vector<double> row_ub(n_rows,  CBC_INFINITY);

  for (size_t i = 0; i < _constraints.size();  ++i)
  {
    DEB_error_if(!_constraints[i]->is_linear(), "constraint " << i <<
                 " is non-linear");

    NConstraintInterface::SVectorNC gc;
    _constraints[i]->eval_gradient(P(x), gc);

    // the setup below appears inefficient, the interface is rather restrictive
    CoinPackedVector lin_expr;
    lin_expr.reserve(gc.size());
    for (NConstraintInterface::SVectorNC::InnerIterator v_it(gc); v_it; ++v_it)
      lin_expr.insert(v_it.index(), v_it.value());
    mtrx.appendRow(lin_expr);

    const double b = -_constraints[i]->eval_constraint(P(x));
    switch (_constraints[i]->constraint_type())
    {
    case NConstraintInterface::NC_EQUAL         :
      row_lb[i] = row_ub[i] = b;
      break;
    case NConstraintInterface::NC_LESS_EQUAL    :
      row_ub[i] = b;
      break;
    case NConstraintInterface::NC_GREATER_EQUAL :
      row_lb[i] = b;
      break;
    }

    TRACE_CBT("Constraint " << i << " is of type " <<
              _constraints[i]->constraint_type() << "; RHS val", b);
  }

  //----------------------------------------------
  // setup energy
  //----------------------------------------------

  std::vector<double> objective(n_cols);
  _problem->eval_gradient(P(x), P(objective));
  TRACE_CBT("Objective linear term", objective);

  const double c = _problem->eval_f(P(x));
  // ICGB: Where does objective constant term go in CBC?
  // MCM: I could not find this either: It is possible that the entire model
  // can be translated to accomodate the constant (if != 0). Ask DB!
  DEB_error_if(c > FLT_EPSILON, "Ignoring a non-zero constant objective term: "
               << c);
  TRACE_CBT("Objective constant term", c);

  // CBC Problem initialize
  OsiClpSolverInterface si;
  si.setHintParam(OsiDoReducePrint, true, OsiHintTry);

  const std::vector<double> col_lb(n_cols, -CBC_INFINITY);
  const std::vector<double> col_ub(n_cols,  CBC_INFINITY);
  si.loadProblem(mtrx, P(col_lb), P(col_ub), P(objective), P(row_lb), P(row_ub));

  // set variables
  for (int i = 0; i < n_cols; ++i)
  {
    switch (var_type[i])
    {
    case Real:
      si.setContinuous(i);
      break;
    case Integer:
      si.setInteger(i);
      break;
    case Binary :
      si.setInteger(i);
      si.setColBounds(i, 0.0, 1.0);
      break;
    }
  }

  // TODO: make this accessible through the DEB system instead
  volatile static bool dump_problem = true; // change on run-time if necessary
  if (dump_problem)
  {
    static int n_mps_dumps = 0;
    char filename[64];
    sprintf_s(filename, "DOCloud_problem_dump_%04i", n_mps_dumps); 
    si.writeMps(filename); //output problem as .MPS
  }

  THROW_OUTCOME(TODO); // unimplemented
  //THROW_OUTCOME_if(solution == nullptr, UNSPECIFIED_CBC_EXCEPTION);
  //_problem->store_result(solution);
}

}// namespace 

#undef TRACE_CBT
#undef P

void DOCloudSolver::solve(
  NProblemInterface*                        _problem,
  const std::vector<NConstraintInterface*>& _constraints,
  const std::vector<PairIndexVtype>&        _discrete_constraints,
  const double                              _time_limit
)
{
  DEB_enter_func;
  try
  {
    solve_impl(_problem, _constraints, _discrete_constraints, _time_limit);
  }
  catch (CoinError& ce)
  {
    DEB_warning(1, "CoinError code = " << ce.message() << "]\n");
    THROW_OUTCOME(TODO);
  }
}


//=============================================================================
} // namespace COMISO
//=============================================================================

#endif // COMISO_DOCLOUD_AVAILABLE
//=============================================================================

