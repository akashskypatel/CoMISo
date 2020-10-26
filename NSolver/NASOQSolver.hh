//=============================================================================
//
//  CLASS NASOQSolver
//
//=============================================================================


#ifndef COMISO_NASOQSOLVER_HH
#define COMISO_NASOQSOLVER_HH


//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
//#if COMISO_NASOQ_AVAILABLE // TODO

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

/** \class NASOQ Solver NASOQSolver.hh

    Solver for quadratic problem with linear equality and linear inequality constraints
    based on NASOQ.
*/
class COMISODLLEXPORT NASOQSolver
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
  Eigen::VectorXd get_q(NProblemInterface* _problem);

  void get_equality_constraints(int _n_cols, const std::vector<NConstraintInterface*>& _constraints, COMISO::NProblemInterface::SMatrixNP& _A, Eigen::VectorXd& _rhs);
  void get_inequality_constraints(int _n_cols, const std::vector<NConstraintInterface*>& _constraints, COMISO::NProblemInterface::SMatrixNP& _C, Eigen::VectorXd& _rhs);
};



//=============================================================================
} // namespace COMISO

//=============================================================================
//#endif // COMISO_NASOQ_AVAILABLE
//=============================================================================
#endif // COMISO_NASOQSOLVER_HH defined
//=============================================================================

