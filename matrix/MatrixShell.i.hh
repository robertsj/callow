//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MatrixShell.i.hh
 * \author robertsj
 * \date   Sep 20, 2012
 * \brief  MatrixShell.i class definition.
 */
//---------------------------------------------------------------------------//

#ifndef callow_MATRIXSHELL_I_HH_
#define callow_MATRIXSHELL_I_HH_

namespace callow
{

template <class T>
MatrixShell<T>::MatrixShell(void* context)
  : d_context(context)
{
  /* ... */
}

template <class T>
MatrixShell<T>::MatrixShell(void* context, const int m, const int n)
  : d_context(context)
{
  set_size(m, n);
}

template <class T>
MatrixShell<T>::~MatrixShell()
{
  if (d_is_ready)
  {
#ifdef CALLOW_ENABLE_PETSC
    // destroy the petsc matrix.  note, since we constructed it
    // using our own arrays, those still need to be deleted.
    MatDestroy(&d_petsc_matrix);
#endif
  }
}

} // end namespace callow


#endif /* callow_MATRIXSHELL_I_HH_ */
