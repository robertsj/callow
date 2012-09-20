//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PetscSolver.hh
 * \author robertsj
 * \date   Sep 20, 2012
 * \brief  PetscSolver class definition.
 */
//---------------------------------------------------------------------------//

#ifndef PETSCSOLVER_HH_
#define PETSCSOLVER_HH_

#include "LinearSolver.hh"

namespace callow
{

/**
 *  \class PetscSolver
 *  \brief Uses PETSc to solve a system
 *
 *
 */
class PetscSolver: public LinearSolver<PetscScalar>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef LinearSolver<PetscScalar> Base;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  PetscSolver(const double atol, const double rtol, const int maxit);

  virtual ~PetscSolver();

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  // expose base class members
  using LinearSolver<PetscScalar>::status;
  using LinearSolver<PetscScalar>::d_absolute_tolerance;
  using LinearSolver<PetscScalar>::d_relative_tolerance;
  using LinearSolver<PetscScalar>::d_maximum_iterations;
  using LinearSolver<PetscScalar>::d_L1_residual;
  using LinearSolver<PetscScalar>::d_L2_residual;
  using LinearSolver<PetscScalar>::d_LI_residual;
  using LinearSolver<PetscScalar>::d_number_iterations;
  using LinearSolver<PetscScalar>::d_A;
  using LinearSolver<PetscScalar>::d_PL;
  using LinearSolver<PetscScalar>::d_PR;

  // petsc solver type
  KSP d_petsc_solver;

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL LINEAR SOLVERS MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  /**
   *  \param b  right hand side
   *  \param x  unknown vector
   */
  void solve_impl(const Vector<PetscScalar> &b, Vector<PetscScalar> &x);

};

} // end namespace callow

// Inline member definitions
#include "PetscSolver.i.hh"

#endif /* PETSCSOLVER_HH_ */
