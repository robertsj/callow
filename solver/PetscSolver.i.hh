//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PetscSolver.i.hh
 * \author robertsj
 * \date   Sep 20, 2012
 * \brief  PetscSolver.i class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
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
// CONSTRUCTOR & DESTRUCTOR
//---------------------------------------------------------------------------//

PetscSolver::PetscSolver(const double  atol,
                         const double  rtol,
                         const int     maxit)
  : LinearSolver<PetscScalar>(atol, rtol, maxit, "solver_petsc")
{
  // Create the KSP object.
  PetscErrorCode ierr = KSPCreate(PETSC_COMM_SELF, &d_petsc_solver);
  Insist(!ierr, "Error creating KSP object.");

   // Set the operator.
   KSPSetOperators(d_solver, d_operator, d_operator, SAME_NONZERO_PATTERN);
}

PetscSolver::~PetscSolver()
{
  KSPDestroy(&d_petsc_solver);
}

//---------------------------------------------------------------------------//
// SOLVE
//---------------------------------------------------------------------------//

inline void PetscSolver::
solve_impl(const Vector<PetscScalar> &b, Vector<PetscScalar> &x)
{
  PetscErrorCode ierr;
  ierr = KSPSolve(d_petsc_solver, b.petsc_vector(), x.petsc_vector());
  Insist(!ierr, "Error in KSPSolve.");
}

} // end namespace callow

#endif /* PETSCSOLVER_I_HH_ */
