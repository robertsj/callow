//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PCShell.hh
 * \author robertsj
 * \date   Sep 20, 2012
 * \brief  PCShell class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef PCSHELL_HH_
#define PCSHELL_HH_

#include "Preconditioner.hh"
#include "matrix/Matrix.hh"

namespace callow
{

/*!
 *  \class PCShell
 *  \brief Applies a shell preconditioner
 *
 *  A shell preconditioner allows the user to define their own
 *  preconditioning processes that are potentially matrix free.
 *
 */
template <class T>
class PCShell: public Preconditioner<T>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef Preconditioner<T>                   Base;
  typedef typename Base::SP_preconditioner    SP_preconditioner;
  typedef typename Matrix<T>::SP_matrix       SP_matrix;
  typedef typename Vector<T>::SP_vector       SP_vector;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /// Construct a shell preconditioner
  PCShell();

  /// Virtual destructor
  virtual ~PCShell(){};

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL PRECONDITIONERS MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  /// Solve Px = b
  virtual void apply(Vector<T> &b, Vector<T> &x) = 0;

protected:

#ifdef CALLOW_ENABLE_PETSC
  /// PETSc preconditioner
  PC d_petsc_pc;
#endif

};

#ifdef CALLOW_ENABLE_PETSC
/// Wrapper apply for PETSc
PetscErrorCode apply_petsc(PC pc, Vec x, Vec y);
#endif

} // end namespace callow

#include "PCShell.i.hh"



#endif /* PCSHELL_HH_ */
