//=============================================================================
//
//  CLASS IPOPTSolver
//
//=============================================================================


#ifndef COMISO_IPOPTSOLVER_HH
#define COMISO_IPOPTSOLVER_HH


//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_IPOPT_AVAILABLE

//== INCLUDES =================================================================

#include <CoMISo/Config/CoMISoDefines.hh>
#include <vector>
#include <cstddef>
#include <gmm/gmm.h>
#include "NProblemGmmInterface.hh"
#include "NProblemInterface.hh"
#include "NConstraintInterface.hh"
#include "BoundConstraint.hh"

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================


/** \class NewtonSolver NewtonSolver.hh <ACG/.../NewtonSolver.hh>

    Brief Description.
  
    A more elaborate description follows.
*/
class COMISODLLEXPORT IPOPTSolver
{
public:
   
  /// Default constructor
  IPOPTSolver();
 
  /// Destructor
  ~IPOPTSolver();

  // ********** SOLVE **************** //
  //! \throws Outcome
  void solve(NProblemInterface*    _problem, const std::vector<NConstraintInterface*>& _constraints);

  // same as above with additional lazy constraints that are only added iteratively to the problem if not satisfied
  //! \throws Outcome
  void solve(NProblemInterface*                        _problem,
            const std::vector<NConstraintInterface*>& _constraints,
            const std::vector<NConstraintInterface*>& _lazy_constraints,
            const double                              _almost_infeasible = 0.5,
            const int                                 _max_passes        = 5   );

  // for convenience, if no constraints are given
  //! \throws Outcome
  void solve(NProblemInterface*    _problem);

  // deprecated interface for backwards compatibility
  //! \throws Outcome
  void solve(NProblemGmmInterface* _problem, 
    std::vector<NConstraintInterface*>& _constraints);

  // ********* CONFIGURATION ********************* //
  // access the ipopt-application (for setting parameters etc.)
  // examples: app().Options()->SetIntegerValue("max_iter", 100);
  //           app().Options()->SetStringValue("derivative_test", "second-order");
  //Ipopt::IpoptApplication& app() {return (*app_); }


protected:
  double* P(std::vector<double>& _v)
  {
    if( !_v.empty())
      return ((double*)&_v[0]);
    else
      return 0;
  }

private:
  class Impl;
  Impl* impl_;
};


//=============================================================================
} // namespace COMISO

//=============================================================================
#endif // COMISO_IPOPT_AVAILABLE
//=============================================================================
#endif // ACG_IPOPTSOLVER_HH defined
//=============================================================================

