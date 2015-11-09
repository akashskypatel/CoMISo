// (C) Copyright 2015 by Autodesk, Inc.
//
// The information contained herein is confidential, proprietary
// to Autodesk,  Inc.,  and considered a trade secret as defined
// in section 499C of the penal code of the State of California.
// Use of  this information  by  anyone  other  than  authorized
// employees of Autodesk, Inc.  is granted  only under a written
// non-disclosure agreement,  expressly  prescribing  the  scope
// and manner of such use.

#ifndef COMISO_ERROR_HH_INCLUDED
#define COMISO_ERROR_HH_INCLUDED

#include <Base/Utils/BaseError.hh>
#include <CoMISo/Config/CoMISoDefines.hh>

namespace COMISO {

class COMISODLLEXPORT Error : public Base::Error
{
public:
  enum Index
  {
  #define DEFINE_ERROR(CODE, MSG) CODE,
  #include <CoMISo/Utils/CoMISoErrorInc.hh>
  #undef DEFINE_ERROR
  };

public:
  //! Constructor.
  Error(const Index _idx) : Base::Error((int)_idx) {}

  //! Return the error message 
  virtual const char* message() const;

protected:
  Error(const int _idx) : Base::Error(_idx) {}
};

}//namespace COMISO

#define COMISO_THROW(INDEX) { THROW_ERROR_MODULE(COMISO, INDEX); }
#define COMISO_THROW_if(COND, INDEX) { if (COND) COMISO_THROW(INDEX); }

#define COMISO_THROW_TODO(MSG) { THROW_ERROR_TODO_MODULE(COMISO, MSG); }
#define COMISO_THROW_TODO_if(COND, MSG) { if (COND) COMISO_THROW_TODO(MSG); }

#endif//COMISO_ERROR_HH_INCLUDED
