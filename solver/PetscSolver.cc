//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PetscSolver.cc
 * \brief  PetscSolver 
 * \author Jeremy Roberts
 * \date   Sep 20, 2012
 */
//---------------------------------------------------------------------------//

#include "PetscSolver.hh"

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

}

PetscSolver::~PetscSolver()
{
  KSPDestroy(&d_petsc_solver);
}


} // end namespace callow

//---------------------------------------------------------------------------//
//              end of file PetscSolver.cc
//---------------------------------------------------------------------------//
