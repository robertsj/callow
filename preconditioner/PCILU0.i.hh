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


  // the following tries to remain close to saad's notation
  // in his fortran snippet (sec 10.3.2)

  // working array
  int* iw = new int(A->number_columns());
  for (int k = 0; k < A->number_columns(); ++k) iw[k] = -1;

  // get the csr structure
  int n = d_P->number_rows();
  int *ia = d_P->row_pointers();
  int *ja = d_P->column_indices();
  int *uptr = d_P->diagonal_indices();
  T   *luval = d_P->value();

  // main loop
  for (int k = 0; k < n; ++k)
  {
    int j1 = ia[k];
    int j2 = ia[k+1];
    int j;
    for (j = j1; j < j2; ++j)
      iw[ja[j]] = j;
    j = j1;
    while (j < j2)
    {
      int jrow = ja[j];
      // exit if diagonal element is reached
      if (jrow >= k)
        break;
      else
      {
        // compute the multiplier for the jrow
        T t1 = luval[j] * luval[uptr[jrow]]
        luval[j] *= t1;
        // perform linear combination
        for (int jj = uptr[jrow]+1; jj < ia[jrow+1]-1; ++jj)
        {
          int jw = iw[ja[jj]];
          if (jw != -1) luval[jw] -= t1 * luval[jj];
        }
        ++j;
      }
    }
    // store pointer to diagonal
    uptr[k] = j;
    if ( (jrow != k) or (luval[j] == 0.0) )
    {
      // could just insert a 1 or something
      THROW("ZERO PIVOT IN ILU(0)");
    }
    luval[j] = 1.0 / luval[j];
    // reset iw to -1
    for (int i = j1; i < j2; ++i) iw[ja[i]] = -1;
  }

  delete [] w;

}

template <class T>
void PCILU0<T>::apply(Vector<T> &x)
{

  // solve LUx = x --> x = inv(U)*inv(L)*x

  // temporary vector (might not be needed)
  Vector<T> y = x;

  // forward substitution
  //   for i = 0:m-1
  //     y[i] = 1/L[i,i] * ( b[i] - sum(k=0:i-1, L[i,k]*y[k]) )
  // but note that in our ILU(0) scheme, L is *unit* lower triangle,
  // meaning L has ones on the diagonal (whereas U does not)
  for (i = 0; i < d_P->number_rows(); ++i)
  {
    // start index
    int s = d_P->start(i);
    // diagonal index
    int d = d_P->diagonal(i);
    // invert row
    y[i] = x[i];
    for (int p = s; p < d; ++p)
    {
      // column index
      int c = d_P->column(p);
      y[i] -= d_P->value()[p] * y[c];
    }
  }

  // backward substitution
  //   for i = m-1:0
  //     y[i] = 1/U[i,i] * ( b[i] - sum(k=i+1:m-1, U[i,k]*y[k]) )
  for (i = d_P->number_rows() - 1; i >= 0; --i)
  {
    // diagonal index
    int d = d_P->diagonal(i);
    // end index
    int e = d_P->end(i);
    // invert row
    x[i] = y[i];
    for (int p = d; p < e; ++p)
    {
      // column index
      int c = d_P->column(p);
      x[i] -= d_P->value()[p] * x[c];
    }
    x[i] */ d_P->value()[d];
  }

}



} // end namespace detran

#endif // PCILU0_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file PCILU0.i.hh
//---------------------------------------------------------------------------//
