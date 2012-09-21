//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PetscSolver.i.hh
 * \author robertsj
 * \date   Sep 20, 2012
 * \brief  PetscSolver.i class definition.
 */
//---------------------------------------------------------------------------//

#ifndef PETSCSOLVER_I_HH_
#define PETSCSOLVER_I_HH_

#include "matrix/Matrix.hh"
#include <cmath>
#include <cstdio>

namespace callow
{

//---------------------------------------------------------------------------//
// SOLVE
//---------------------------------------------------------------------------//

inline void PetscSolver::
solve_impl(const Vector<PetscScalar> &b, Vector<PetscScalar> &x)
{
  PetscErrorCode ierr;
  ierr = KSPSolve(d_petsc_solver, const_cast<Vector<PetscScalar>* >(&b)->petsc_vector(), x.petsc_vector());
  Insist(!ierr, "Error in KSPSolve.");
}

} // end namespace callow

#endif /* PETSCSOLVER_I_HH_ */
