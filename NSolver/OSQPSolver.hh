//=============================================================================
//
//  CLASS OSQPSolver
//
//=============================================================================


#ifndef COMISO_OSQPSOLVER_HH
#define COMISO_OSQPSOLVER_HH


//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
//#if COMISO_OSQP_AVAILABLE // TODO

//== INCLUDES =================================================================

#include <CoMISo/Config/CoMISoDefines.hh>
#include <vector>
#include <string>
#include "NProblemInterface.hh"
#include "NConstraintInterface.hh"

//== FORWARDDECLARATIONS ======================================================


//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================

/** \class OSQP Solver OSQPSolver.hh

    Solver for quadratic problem with linear equality and linear inequality constraints
    based on OSQP.
*/
class COMISODLLEXPORT OSQPSolver
{
public:
  // ********** SOLVE **************** //
  bool solve(NProblemInterface*                        _problem,                // problem instance
             const std::vector<NConstraintInterface*>& _constraints            // linear constraints
             );

// TODO:
//  // same as above with additional lazy constraints that are only added iteratively to the problem if not satisfied
//  bool solve(NProblemInterface*                        _problem,
//                    const std::vector<NConstraintInterface*>& _constraints,
//                    const std::vector<NConstraintInterface*>& _lazy_constraints)


protected:

private:

  void regularize_hessian(NProblemInterface::SMatrixNP& _H);

  NProblemInterface::SMatrixNP get_hessian(NProblemInterface* _problem);
  Eigen::VectorXd get_linear_energy_coefficients(NProblemInterface* _problem);

  void get_constraints(int _n_cols, const std::vector<NConstraintInterface*>& _constraints, COMISO::NProblemInterface::SMatrixNP& _C, Eigen::VectorXd& _lower_bounds, Eigen::VectorXd& _upper_bounds);
};



//=============================================================================
} // namespace COMISO

//=============================================================================
//#endif // COMISO_OSQP_AVAILABLE
//=============================================================================
#endif // COMISO_OSQPSOLVER_HH defined
//=============================================================================

