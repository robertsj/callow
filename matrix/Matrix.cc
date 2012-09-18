//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Matrix.cc
 * \author robertsj
 * \date   Sep 13, 2012
 * \brief  Matrix class definition.
 */
//---------------------------------------------------------------------------//

#include "Matrix.hh"
#include <iostream>

namespace callow
{

// Instantiations
#ifdef CALLOW_ENABLE_PETSC
template class Matrix<PetscScalar>;
#else
template class Matrix<float>;
template class Matrix<double>;
#endif

} // end namespace callow

