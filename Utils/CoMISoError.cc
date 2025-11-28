// (C) Copyright 2015 by Autodesk, Inc.

#include "CoMISoError.hh"

namespace COMISO { 

static const char* ERROR_MESSAGE[] =
{
  #define DEFINE_ERROR(CODE, MSG) MSG,
  #include "CoMISoErrorInc.hh"
  #undef DEFINE_ERROR
};

Error::Error(const int _idx) : Error(_idx, ERROR_MESSAGE[_idx]) {}
const char* Error::message() const { return ERROR_MESSAGE[idx_]; }

}//namespace COMISO
