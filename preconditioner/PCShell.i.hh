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

template <class T>
PCShell::PCShell
{

#ifdef CALLOW_ENABLE_PETSC

  // Set the PC context
  ierr = PCShellSetContext(d_PC, this);
  Insist(!ierr, "Error setting shell preconditioner context.");

  // Set the PC operator
  ierr = PCShellSetApply(d_PC, apply_petsc);
  Insist(!ierr, "Error setting shell preconditioner operator.");

  // Set the PC name for good measure
  ierr = PCShellSetName(wg_pc, "DSA Preconditioner");
  Insist(!ierr, "Error within-group preconditioner name.");

#endif
}

#ifdef CALLOW_ENABLE_PETSC

inline PetscErrorCode apply_petsc(PC pc, Vec b, Vec x);
{
  // Get the context and cast as InnerGMRES pointer.
  PetscErrorCode ierr;
  void *ctx;
  ierr = PCShellGetContext(pc, &ctx); CHKERRQ(ierr);
  callow::PCShell *foo = (callow::PCShell *) ctx;
  // Call the actual apply operator.
  return foo->apply(b, x);
}

#endif

#endif /* PCSHELL_I_HH_ */
