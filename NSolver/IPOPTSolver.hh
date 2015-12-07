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


//== NAMESPACES ===============================================================

namespace COMISO {

//== FORWARDDECLARATIONS ======================================================
class NProblemGmmInterface; // deprecated
class NProblemInterface;
class NConstraintInterface;

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
  void solve(NProblemInterface* _problem, 
    const std::vector<NConstraintInterface*>& _constraints);

  //! Same as above with additional lazy constraints that are only added iteratively to the problem if not satisfied
  //! \throws Outcome
  void solve(NProblemInterface* _problem,
    const std::vector<NConstraintInterface*>& _constraints,
    const std::vector<NConstraintInterface*>& _lazy_constraints,
    const double _almost_infeasible = 0.5,
    const int _max_passes = 5   );

  // for convenience, if no constraints are given
  //! \throws Outcome
  void solve(NProblemInterface* _problem);

  // deprecated interface for backwards compatibility
  //! \throws Outcome
  void solve(NProblemGmmInterface* _problem, 
    std::vector<NConstraintInterface*>& _constraints);

  //! Get the computed solution energy 
  double energy();

private:
  class Impl;
  Impl* impl_;

  // inhibit copy
  IPOPTSolver(const IPOPTSolver&);
  IPOPTSolver& operator=(const IPOPTSolver&);
};


//=============================================================================
} // namespace COMISO

//=============================================================================
#endif // COMISO_IPOPT_AVAILABLE
//=============================================================================
#endif // ACG_IPOPTSOLVER_HH defined
//=============================================================================

