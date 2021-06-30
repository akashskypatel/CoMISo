// Copyright 2021 Autodesk, Inc. All rights reserved.

#include "EigenLSQConstrainedSolverT.cc"

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#if (COMISO_EIGEN3_AVAILABLE)
//== INCLUDES =================================================================

namespace COMISO
{

template class EigenLSQConstrainedSolverT<1>;
template class EigenLSQConstrainedSolverT<2>;
template class EigenLSQConstrainedSolverT<3>;

} // namespace COMISO

//=============================================================================
#endif // COMISO_EIGEN3_AVAILABLE
//=============================================================================
