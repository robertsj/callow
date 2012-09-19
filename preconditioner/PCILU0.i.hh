//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PCILU0.i.hh
 * \brief  PCILU0.i 
 * \author Jeremy Roberts
 * \date   Sep 18, 2012
 */
//---------------------------------------------------------------------------//

#ifndef PCILU0_I_HH_
#define PCILU0_I_HH_

namespace callow
{

template <class T>
PCILU0<T>::PCILU0(SP_matrix A)
{
  // preconditions
  Require(A);
  Require(A->number_rows() == A->number_columns());

  // copy A
  d_P(new Matrix<T>(*A));

  // working array
  int* iw = new int(A->number_columns());
  for (int k = 0; k < A->number_columns(); ++k) iw[k] = 0;

  // Get the CSR structure
  int *ia = d_P->row_pointers();
  int *ja = d_P->column_indices();
  int *uptr = d_P->diagonal_indices();
  T   *a  = d_P->value();

  // main loop
  for (int k = 0; k < A->number_columns(); ++k)
  {
    int j1 = ia[k];
    int j2 = ia[k+1];
    for (int j = j1; j < j2; ++j)
      iw[ja[j]] = j;
    j = j1;
    while (j < j2)
    {
      int jrow = ja[j];
      // if not diagonal
      if (jrow >= k)
        break;
      else
      {
        // compute the multiplier for the jrow
        a[j] *= a[uptr[jrow]];
        // perform linear combination
        for (int jj = uptr[jrow]+1; jj < ia[jrow+1]; ++jj)
        {
          int jw = iw[ja[jj]];
          if (jw != 0) a[jw] -= t1 * a[jj];
        }
        ++j;
      }
    }
    // store pointer to diagonal
    uptr[k] = j;
    if ( (jrow != k) or (a[j] == 0.0) )
    {
      THROW("ZERO PIVOT IN ILU(0)");
    }
    a[j] = 1.0 / a[j];
    // reset iw to zero
    for (int i = j1; i < j2; ++i) iw[ja[i]] = 0;
  }

  delete [] w;

}

template <class T>
void PCILU0<T>::apply(Vector<T> &x)
{

}



} // end namespace detran

#endif // PCILU0_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file PCILU0.i.hh
//---------------------------------------------------------------------------//
