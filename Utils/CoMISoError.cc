// (C) Copyright 2015 by Autodesk, Inc.
//
// The information contained herein is confidential, proprietary
// to Autodesk,  Inc.,  and considered a trade secret as defined
// in section 499C of the penal code of the State of California.
// Use of  this information  by  anyone  other  than  authorized
// employees of Autodesk, Inc.  is granted  only under a written
// non-disclosure agreement,  expressly  prescribing  the  scope
// and manner of such use.

#include "CoMISoError.hh"

namespace COMISO { 

static const char* ERROR_MESSAGE[] =
{
  #define DEFINE_ERROR(CODE, MSG) MSG,
  #include "CoMISoErrorInc.hh"
  #undef DEFINE_ERROR
};

const char* Error::message() const { return ERROR_MESSAGE[idx_]; }

}//namespace COMISO
