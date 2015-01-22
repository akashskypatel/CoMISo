// (C) Copyright 2014 by Autodesk, Inc.
//
// The information contained herein is confidential, proprietary
// to Autodesk,  Inc.,  and considered a trade secret as defined
// in section 499C of the penal code of the State of California.
// Use of  this information  by  anyone  other  than  authorized
// employees of Autodesk, Inc.  is granted  only under a written
// non-disclosure agreement,  expressly  prescribing  the  scope
// and manner of such use.

#ifndef GMMTYPES_HH_INCLUDED
#define GMMTYPES_HH_INCLUDED

#include <gmm/gmm_matrix.h>

namespace COMISO_GMM
{  

// Matrix Types
typedef gmm::col_matrix< gmm::wsvector<double> > WSColMatrix;
typedef gmm::row_matrix< gmm::wsvector<double> > WSRowMatrix;
typedef gmm::col_matrix< gmm::rsvector<double> > RSColMatrix;
typedef gmm::csc_matrix<double> CSCMatrix;

}//namespace COMISO_GMM

#endif//GMMTYPES_HH_INCLUDED
