//=============================================================================
//
//  CLASS LazyConstraintSolver
//
//=============================================================================


#ifndef COMISO_LAZYCONSTRAINTSOLVER_HH
#define COMISO_LAZYCONSTRAINTSOLVER_HH


//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
//#if COMISO_OSQP_AVAILABLE // TODO

//== INCLUDES =================================================================

#include <CoMISo/Config/CoMISoDefines.hh>
#include <vector>
#include <string>
#include "NProblemInterface.hh"
#include "NConstraintInterface.hh"

#include <Base/Debug/DebOut.hh>

//== FORWARDDECLARATIONS ======================================================


//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================

namespace detail
{
  enum feasibility { feasible, infeasible, almost_infeasible };

  feasibility get_feasibility(NConstraintInterface* c, double* x, double _acceptable_tolerance, double _almost_infeasible_threshold)
  {
    auto v = c->eval_constraint(x);

    if (c->constraint_type() == NConstraintInterface::NC_EQUAL)
    {
      if (std::abs(v) < _acceptable_tolerance)
        return feasible;
      else
        return infeasible;
    }
    else if (c->constraint_type() == NConstraintInterface::NC_LESS_EQUAL)
    {
      if (v >= -_acceptable_tolerance)
        return infeasible;
      else if (v >= -_almost_infeasible_threshold)
        return almost_infeasible;
      else
        return feasible;
    }
    else if (c->constraint_type() == NConstraintInterface::NC_GREATER_EQUAL)
    {
      if (v <= _acceptable_tolerance)
        return infeasible;
      else if (v <= _almost_infeasible_threshold)
        return almost_infeasible;
      else return feasible;
    }
    else
    {
      DEB_error("Unknown constraint type");
      return infeasible;
    }
  }
}

/// SolveFunction should be callable with two arguments: NProblemInterface* and const std::vector<NConstraintInterface*>
/// Result function should return double* to solution
template <typename SolveFunction, typename ResultFunction>
auto solve_with_lazy_constraints(SolveFunction& _solve,
                                 ResultFunction& _get_result,
                                 NProblemInterface* _problem,
                                 const std::vector<NConstraintInterface*>& _initial_constraints,
                                 const std::vector<NConstraintInterface*>& _lazy_constraints,
                                 double _acceptable_tolerance = 1e-8,
                                 double _almost_infeasible_threshold = 0.5,
                                 int _max_passes = 5,
                                 bool _final_step_with_all_constraints = false) -> decltype(_solve(_problem, _initial_constraints))
{
  DEB_enter_func;

  auto res = _solve(_problem, _initial_constraints);

  std::vector<NConstraintInterface*> constraints = _initial_constraints;
  std::vector<bool> added(_lazy_constraints.size(), false);

  std::vector<int> n_infeasible;
  std::vector<int> n_almost_infeasible;

  for (int pass = 0; pass < _max_passes; ++pass)
  {
    n_infeasible.push_back(0);
    n_almost_infeasible.push_back(0);

    double* solution_x = _get_result();

    for (size_t i = 0; i < _lazy_constraints.size(); ++i)
    {
      if (added[i])
        continue; // we already added this constraint

      auto f = detail::get_feasibility(_lazy_constraints[i], solution_x, _acceptable_tolerance, _almost_infeasible_threshold);
      if (f == detail::infeasible)
        ++n_infeasible.back();
      else if (f == detail::almost_infeasible)
        ++n_almost_infeasible.back();

      if (f != detail::feasible)
      {
        constraints.push_back(_lazy_constraints[i]);
        added[i] = true;
      }
    }

    if (n_infeasible.back() == 0)
      break; // if nothing is infeasible we are done

    res = _solve(_problem, constraints);
//    _almost_infeasible_threshold *= 2;
  }


  // Retrieve some statistics about the solve
  DEB_line(4, "############# lazy constraints statistics ###############");
  DEB_line(4, _lazy_constraints.size() << " lazy constraints in input.");
  DEB_line(4, "#passes     : " << n_infeasible.size() << "( of " << _max_passes << ")");
  for(int i=0; i<n_infeasible.size(); ++i)
    DEB_line(5, "pass " << i+1 << " induced " << n_infeasible[i]
      << " infeasible and " << n_almost_infeasible[i] << " almost infeasible");

  return res;
}


//=============================================================================
} // namespace COMISO

//=============================================================================
//#endif // COMISO_OSQP_AVAILABLE
//=============================================================================
#endif // COMISO_OSQPSOLVER_HH defined
//=============================================================================

