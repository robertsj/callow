//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DBC.hh
 * \author robertsj
 * \date   Sep 13, 2012
 * \brief  Design by Contract macros
 */
//---------------------------------------------------------------------------//

#ifndef callow_DBC_HH_
#define callow_DBC_HH_

#include "callow_config.hh"
#include "CallowException.hh"

#ifdef CALLOW_ENABLE_DBC

#define Assert(c)     if (!(c)) throw callow::CallowException( __LINE__, __FILE__,#c)
#define IsValid(obj)  Assert((obj) != NULL && (obj)->is_valid())
#define Require(c)    Assert(c)
#define Ensure(c)     Assert(c)

#else

#define Assert(c)   ((void) 0)
#define IsValid(c)  ((void) 0)
#define Require(c)  ((void) 0)
#define Ensure(c)   ((void) 0)

#endif

#define Insist(c,m)   if (!(c)) {std::cerr << m << std::endl; throw callow::CallowException( __LINE__, __FILE__,#c);}

#endif /* callow_DBC_HH_ */
