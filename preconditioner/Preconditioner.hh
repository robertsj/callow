//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Preconditioner.hh
 * \brief  Preconditioner 
 * \author Jeremy Roberts
 * \date   Sep 18, 2012
 */
//---------------------------------------------------------------------------//

#ifndef PRECONDITIONER_HH_
#define PRECONDITIONER_HH_

#include "callow_config.hh"
#include "matrix/MatrixBase.hh"
#include "vector/Vector.hh"
#include <string>

namespace callow
{

/*!
 *  \class Preconditioner
 *  \brief Defines a preconditioner for linear solves
 *
 *  Consider the linear system
 *  \f[
 *      \mathbf{A}x = b \, .
 *  \f]
 *  When using iterative methods to solve this system,
 *  one can often make the system easier to solve.  Suppose
 *  we define an operator \f$ \mathbf{P} \f$ such that
 *  \f$ \mathbf{P}^{-1} \mathbf{A} \approx \mathbf{I} \f$.
 *  If the action of \f$ \mathbf{P}^{-1} \f$ is relatively
 *  easy to compute (as compared to inverting \mathbf{A}),
 *  \f$ \mathbf{P} \f$ is a good preconditioner.
 *
 *  We apply a precondition on the left to obtain the
 *  modified system
 *  \f[
 *      \mathbf{P^{-1} A}x = \mathbf{P}^{-1} b \, ,
 *  \f]
 *  or on the right to get
 *  \f[
 *      \mathbf{AP^{-1}} \overbrace{\mathbf{P}x}^{y} =  b \, ,
 *  \f]
 *  following which we solve
 *  \f[
 *      x = \mathbf{P}^{-1} y \, .
 *  \f]
 *
 *  Within callow, the Jacobi and ILU(0) preconditioners are
 *  available along with user-defined shell preconditioners.
 *  If built with PETSc, all preconditioners are available (to PETSc)
 *  as shells.  Otherwise, the user can set PETSc preconditioners
 *  with PetscSolver parameters.
 */
template <class T>
class Preconditioner
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef SP<Preconditioner<T> >              SP_preconditioner;
  typedef typename MatrixBase<T>::SP_matrix   SP_matrix;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  Preconditioner(std::string name = "preconditioner"){};

  virtual ~Preconditioner(){};

#ifdef CALLOW_ENABLE_PETSC
  /**
   *  set the PETSc preconditioner and do other setup
   *
   *  this should only be called by PetscSolver
   */
  void set_petsc_pc(PC pc);
  /// return petsc preconditioner
  PC petsc_pc() {return d_petsc_pc;}
#endif

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL PRECONDITIONERS MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  /// Solve Px = b
  virtual void apply(Vector<T> &b, Vector<T> &x) = 0;

protected:

  /// pc name
  std::string d_name;

#ifdef CALLOW_ENABLE_PETSC
  /// PETSc preconditioner
  PC d_petsc_pc;
#endif

};

#ifdef CALLOW_ENABLE_PETSC
// this is the function petsc actual calls; internally, it redirects
// to our own operation.  all callow preconditioners are views by
// petsc as shells.
PetscErrorCode pc_apply_wrapper(PC pc, Vec b, Vec x);
#endif

} // end namespace callow

#include "Preconditioner.i.hh"

#endif // PRECONDITIONER_HH_ 

//---------------------------------------------------------------------------//
//              end of file Preconditioner.hh
//---------------------------------------------------------------------------//
