//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PCShell.i.hh
 * \author robertsj
 * \date   Sep 20, 2012
 * \brief  PCShell.i class definition.
 */
//---------------------------------------------------------------------------//

#ifndef PCSHELL_I_HH_
#define PCSHELL_I_HH_

namespace callow
{

template <class T>
PCShell<T>::PCShell(std::string name)
  : Preconditioner<T>(name)
{

#ifdef CALLOW_ENABLE_PETSC
  PetscErrorCode ierr;
  // set the PC context
  ierr = PCShellSetContext(d_petsc_pc, this);
  Insist(!ierr, "Error setting shell preconditioner context.");
  // set the PC operator
  ierr = PCShellSetApply(d_petsc_pc, pc_apply_wrapper);
  Insist(!ierr, "Error setting shell preconditioner operator.");
  // set the PC name for good measure
  ierr = PCShellSetName(d_petsc_pc, d_name.c_str());
  Insist(!ierr, "Error within-group preconditioner name.");
#endif

}

#ifdef CALLOW_ENABLE_PETSC

inline PetscErrorCode pc_apply_wrapper(PC pc, Vec b, Vec x)
{
  // get the context and cast
  PetscErrorCode ierr;
  void *context;
  ierr = PCShellGetContext(pc, &context); CHKERRQ(ierr);
  PCShell<PetscScalar> *foo = (PCShell<PetscScalar> *) context;
  // wrap the petsc vectors
  Vector<PetscScalar> B(b);
  Vector<PetscScalar> X(x);
  // call the actual apply operator.
  foo->apply(B, X);
  return ierr;
}

#endif

} // end namespace callow

#endif /* PCSHELL_I_HH_ */
