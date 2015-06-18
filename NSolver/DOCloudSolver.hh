//=============================================================================
//
//  CLASS DOCloudSolver
//
//=============================================================================
#ifndef COMISO_DOCloudSolver_HH
#define COMISO_DOCloudSolver_HH

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_DOCLOUD_AVAILABLE

//== INCLUDES =================================================================

#include <CoMISo/Config/CoMISoDefines.hh>
#include <vector>
#include <string>
#include "NProblemInterface.hh"
#include "NConstraintInterface.hh"
#include "VariableType.hh"
#include "GurobiHelper.hh"

//== FORWARDDECLARATIONS ======================================================
class GRBModel;
class GRBVar;

//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================

/**
    Solver interface for the IBM Decision Optimization Cloud.

    \todo A more elaborate description.
*/
class COMISODLLEXPORT DOCloudSolver
{
public:
  /// Default constructor
  DOCloudSolver() {}

  /// Destructor
  ~DOCloudSolver() {}

  static void set_api_key(const char* _api_key);

  // ********** SOLVE **************** //
  //! \throws Outcome
  void solve(
    NProblemInterface* _problem, // problem instance
    const std::vector<NConstraintInterface*>& _constraints, // linear constraints
    const std::vector<PairIndexVtype>& _discrete_constraints, // discrete constraints
    const double _time_limit = 60); // time limit in seconds

  //! \throws Outcome
  void solve(
    NProblemInterface* _problem, // problem instance
    const std::vector<NConstraintInterface*>& _constraints, // linear constraints
    const double _time_limit = 60) // time limit in seconds
  {
    std::vector<PairIndexVtype> dc;
    return solve(_problem, _constraints, dc, _time_limit);
  }
};


//=============================================================================
} // namespace COMISO

//=============================================================================
#endif // COMISO_DOCLOUD_AVAILABLE
//=============================================================================
#endif // COMISO_DOCloudSolver_HH
//=============================================================================

