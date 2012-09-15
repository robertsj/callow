//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   matrix_fixture.hh
 * \brief  matrix_fixture 
 * \author Jeremy Roberts
 * \date   Sep 13, 2012
 */
//---------------------------------------------------------------------------//
#ifndef MATRIX_FIXTURE_HH_
#define MATRIX_FIXTURE_HH_

#include "matrix/Matrix.hh"

namespace callow
{

template<class T>
typename Matrix<T>::SP_matrix test_matrix_1(int n = 5)
{

  typename Matrix<T>::SP_matrix A;
  A = new Matrix<T>(n, n);
  A->preallocate(3);

  double l = -0.20;
  double d = 0.50;
  double u = -0.20;
//  double l = -1;
//  double d = 2;
//  double u = -1;

  for (int row = 0; row < n; row++)
  {
    if (row == 0)
    {
      int c[] =
      { 0, 1 };
      double v[] =
      { d, u };
      A->add_row(row, c, v, 2);
    }
    else if (row == n - 1)
    {
      int c[] =
        { n - 2, n - 1 };
      double v[] =
        { l, d };
      A->add_row(row, c, v, 2);
    }
    else
    {
      int c[] =
        { row - 1, row, row + 1 };
      double v[] =
        { l, d, u };
      A->add_row(row, c, v);
    }
  }
  A->assemble();
  return A;
}

} // end namespace detran

#endif // MATRIX_FIXTURE_HH_ 
//---------------------------------------------------------------------------//
//              end of file matrix_fixture.hh
//---------------------------------------------------------------------------//
