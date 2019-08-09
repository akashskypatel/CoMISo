//=============================================================================
//
//  CLASS IPOPTSolverLean
//
//=============================================================================


#ifndef COMISO_IPOPTLEANSOLVER_HH
#define COMISO_IPOPTLEANSOLVER_HH

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_IPOPT_AVAILABLE

//== INCLUDES =================================================================
#include <CoMISo/Config/CoMISoDefines.hh>
#include <vector>
#include <cstddef>


//== NAMESPACES ===============================================================
namespace COMISO {

//== FORWARDDECLARATIONS ======================================================
class NProblemInterface;
class NConstraintInterface;

//== CLASS DEFINITION =========================================================

/** \class IPOPTSolverLean
    Solver for Interior Point optimization problems.

    Solves an interior point problem, given an NProblemInterface
    instance and optionally a set of constraints as well as "lazy
    constraints" via NConstraintInterface.

    Lazy constraints are not active while the initial solution to the
    problem is computed. After the first solution is found, the lazy
    constraints are checked and added to the set of active constraints
    if they are violated. This process is then repeated until all
    constraints are satisfied OR a maximum number of solution attempts
    has been reached. In that case the optimization is started once
    more, with all lazy constraints active.
*/
class COMISODLLEXPORT IPOPTSolverLean
{
public:

  IPOPTSolverLean();
  ~IPOPTSolverLean();

  // *********** OPTIONS **************//

  /*!
  Set the maximum number of iterations
  */
  void set_max_iterations(const int _max_iterations);
  int max_iterations() const;

  /*!  Set the threshold on the lazy inequality constraint to decide
  if we are near the constraint boundary.
  */
  void set_almost_infeasible_threshold(const double _alm_infsb_thrsh);
  double almost_infeasible_threshold() const;

  /*!
  Set the max number of incremental lazy constraint iterations before switching
  to the fully constrained problem.
  \note The default value is 5.
  */
  void set_incremental_lazy_constraint_max_iteration_number
    (const int _incr_lazy_cnstr_max_iter_nmbr);
  int incremental_lazy_constraint_max_iteration_number() const;

  /*
  Turn on/off solving the fully constraint problem after exhausting the
  incremental lazy constraint iterations.

  \note The default value of this is true.
  */
  void set_enable_all_lazy_contraints(const bool _enbl_all_lzy_cnstr);
  bool enable_all_lazy_contraints() const;

  // ********** SOLVE **************** //

  //! Solve a problem instance with an optional set of constraints.
  //! \throws Outcome
  void solve
  (NProblemInterface* _problem,
   const std::vector<NConstraintInterface*>& _constraints = {});

  //! Same as above with additional lazy constraints that are only
  //! added iteratively to the problem if not satisfied.
  //! \throws Outcome
  void solve
  (NProblemInterface* _problem,
   const std::vector<NConstraintInterface*>& _constraints,
   const std::vector<NConstraintInterface*>& _lazy_constraints);

  void solve(NProblemInterface* _problem);

  //! Get the computed solution energy
  double energy();

private:
  class Impl;
  Impl* impl_;

  // inhibit copy
  IPOPTSolverLean(const IPOPTSolverLean&);
  IPOPTSolverLean& operator=(const IPOPTSolverLean&);
};

} // namespace COMISO

//=============================================================================
#endif // COMISO_IPOPT_AVAILABLE
//=============================================================================
#endif // COMISO_IPOPTLEANSOLVER_HH defined
//=============================================================================
