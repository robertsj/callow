//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Initialization.hh
 * \author robertsj
 * \date   Sep 18, 2012
 * \brief  Initialization class definition.
 */
//---------------------------------------------------------------------------//

#ifndef INITIALIZATION_HH_
#define INITIALIZATION_HH_

#include "callow_config.hh"

/// Initialize external packages, if enabled
inline void callow_initialize(int argc, char** argv)
{
#ifdef CALLOW_ENABLE_PETSC
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
#endif
#ifdef CALLOW_ENABLE_SLEPC
  SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
#endif
}

/// Finalize external packages, if enabled
inline void callow_finalize()
{
#ifdef CALLOW_ENABLE_PETSC
  PetscFinalize();
#endif
#ifdef CALLOW_ENABLE_SLEPC
  SlepcFinalize();
#endif
}

#endif /* INITIALIZATION_HH_ */
