//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Matrix.i.hh
 * \author robertsj
 * \date   Sep 13, 2012
 * \brief  Matrix.i class definition.
 */
//---------------------------------------------------------------------------//

#ifndef MATRIX_I_HH_
#define MATRIX_I_HH_

#include "utils/DBC.hh"
#include <algorithm>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>

namespace callow
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR & DESTRUCTOR
//---------------------------------------------------------------------------//

template <class T>
Matrix<T>::Matrix()
  : d_allocated(false)
{
  /* ... */
}

template <class T>
Matrix<T>::Matrix(const int m, const int n)
  : MatrixBase<T>(m, n)
  , d_allocated(false)
{
  /* ... */
}

template <class T>
Matrix<T>::Matrix(const int m, const int n, const int nnz)
  : MatrixBase<T>(m, n)
  , d_allocated(false)
{
  preallocate(nnz);
}

template <class T>
Matrix<T>::~Matrix()
{
  if (d_is_ready)
  {
    delete [] d_value;
    delete [] d_column_indices;
    delete [] d_row_pointers;
    delete [] d_diagonal;
  }
}

//---------------------------------------------------------------------------//
// PREALLOCATION
//---------------------------------------------------------------------------//

template <class T>
inline void Matrix<T>::preallocate(const int nnz)
{
  // preconditions
  Require(d_sizes_set);
  Require(!d_allocated);
  Require(nnz > 0);
  // set number of nonzeros, preallocate triplets, and initialize counter
  d_nnz = nnz;
  d_aij = new triplet[nnz];
  for (int i = 0; i < nnz; i++)
  {
    d_aij[i].i = 0;
    d_aij[i].j = 0;
    d_aij[i].v = 0.0;
  }
  d_counter = 0;
  d_allocated = true;
}

//---------------------------------------------------------------------------//
// ASSEMBLING
//---------------------------------------------------------------------------//

/*
 *  We construct in COO format, ie with (i,j,v) triplets.  This makes
 *  adding entries much easier.  Now, we need to order everything
 *  so that the resulting CSR storage has for all i-->j, with j in
 *  order.  The diagonal pointer is also set.  If (i,i,v) doesn't exist,
 *  we insert it, since that's needed for Gauss-Seidel, preconditioning,
 *  and other things.
 */
template <class T>
inline void Matrix<T>::assemble()
{
  // preconditions
  Require(!d_is_ready);
  Require(d_allocated);

  // sort the triplets by row
  std::cout << " sorting rows..." << std::endl;
  sort_triplets(d_aij, d_counter);

  // sort within row by column
  std::cout << " sorting cols..." << std::endl;
  int row_s = 0;
  int row_e = 0;
  for (int i = 0; i < d_m; i++)
  {
    bool diag = false;
    row_e = row_s;
    Assert(row_e < d_nnz);
    while (d_aij[row_e].i == i)
    {
      // if we find a diagonal entry or it shouldn't exist, set it to true
      if (d_aij[row_e].j == i or i >= d_n) diag = true;
      // don't save the zero
      if (d_aij[row_e].v == 0.0) d_counter--;
      if (++row_e == d_nnz) break;
    }
    // add entry for diagonal, even if zero.
    if (!diag) d_counter++;
    // sort this row (pointer, count, flag for column)
    sort_triplets(&d_aij[row_s], row_e-row_s, false);
    row_s = row_e;
  }
  std::cout << " creating csr..." << std::endl;

  // allocate
  d_nnz            = d_counter;
  d_row_pointers   = new int[d_m + 1];
  d_column_indices = new int[d_nnz];
  d_value          = new T[d_nnz];
  d_diagonal       = new int[d_m];

  // fill the csr storage
  int aij_idx = 0;
  int val_idx = 0;
  d_row_pointers[0] = 0;
  for (int i = 0; i < d_m; i++)
  {
    bool diag = false;
    int row_count = 0;
    while (d_aij[aij_idx].i == i)
    {
      Assert(val_idx < d_nnz);
      // lower diagonal
      if (d_aij[aij_idx].j < i)
      {
        d_value[val_idx] = d_aij[aij_idx].v;
        d_column_indices[val_idx] = d_aij[aij_idx].j;
      }
      // diagonal, if present
      else if (d_aij[aij_idx].j == i)
      {
        d_value[val_idx] = d_aij[aij_idx].v;
        d_column_indices[val_idx] = d_aij[aij_idx].j;
        d_diagonal[i] = val_idx;
        diag = true;
      }
      // upper diagonal
      else
      {
        // add diagonal if it wasn't done and it should exist on this row
        if (!diag and i < d_n)
        {
          d_value[val_idx] = 0.0;
          d_column_indices[val_idx] = i;
          d_diagonal[i] = val_idx;
          diag = true;
          ++val_idx;
        }
        Assert(val_idx < d_nnz);
        d_value[val_idx] = d_aij[aij_idx].v;
        d_column_indices[val_idx] = d_aij[aij_idx].j;
      }
      ++val_idx;
      ++aij_idx;
      ++row_count;
      if (aij_idx >= d_counter) break;
    }
    d_row_pointers[i + 1] = row_count + d_row_pointers[i];
  }
  // delete the coo storage and specify the matrix is set to use
  delete [] d_aij;
  d_is_ready = true;
}

template <class T>
inline void Matrix<T>::sort_triplets(triplet* aij, const int n, bool flag)
{
  const int MAX_LEVELS = 300;
  triplet piv;
  int beg[MAX_LEVELS];
  int end[MAX_LEVELS];
  int i = 0, L, R, swap;
  beg[0] = 0;
  end[0] = n;
  while (i >= 0)
  {
    L = beg[i];
    R = end[i] - 1;
    if (L < R)
    {
      piv = aij[L];
      while (L < R)
      {
        if (flag) // by row
        {
          while (aij[R].i >= piv.i && L < R) R--;
          if (L < R) aij[L++] = aij[R];
          while (aij[L].i <= piv.i && L < R) L++;
          if (L < R) aij[R--] = aij[L];
        }
        else // by column
        {
          while (aij[R].j >= piv.j && L < R) R--;
          if (L < R) aij[L++] = aij[R];
          while (aij[L].j <= piv.j && L < R) L++;
          if (L < R) aij[R--] = aij[L];
        }
      }
      aij[L    ] = piv;
      beg[i + 1] = L + 1;
      end[i + 1] = end[i];
      end[i++  ] = L;
      if (end[i] - beg[i] > end[i - 1] - beg[i - 1])
      {
        swap       = beg[i];
        beg[i    ] = beg[i - 1];
        beg[i - 1] = swap;
        swap       = end[i];
        end[i    ] = end[i - 1];
        end[i - 1] = swap;
      }
    }
    else
    {
      i--;
    }
  }

}

//---------------------------------------------------------------------------//
// INDEXING
//---------------------------------------------------------------------------//

template <class T>
inline int Matrix<T>::start(const int i) const
{
  Require(d_is_ready);
  Require(i >= 0 and i < d_m);
  return d_row_pointers[i];
}
template <class T>
inline int Matrix<T>::diagonal(const int i) const
{
  Require(d_is_ready);
  Require(i >= 0 and i < d_m);
  return d_diagonal[i];
}

template <class T>
inline int Matrix<T>::end(const int i) const
{
  Require(d_is_ready);
  Require(i >= 0 and i < d_m);
  return d_row_pointers[i + 1];
}

template <class T>
inline int Matrix<T>::column(const int p) const
{
  Require(d_is_ready);
  Require(p >= 0 and p < d_nnz);
  return d_column_indices[p];
}

template <class T>
inline T Matrix<T>::operator[](const int p) const
{
  Require(d_is_ready);
  Require(p >= 0 and p < d_nnz);
  return d_value[p];
}

//---------------------------------------------------------------------------//
// ACCESS
//---------------------------------------------------------------------------//

template <class T>
inline T Matrix<T>::operator()(const int i, const int j) const
{
  Require(d_is_ready);
  Require(i >= 0 and i < d_m);
  Require(j >= 0 and j < d_n);
  // loop through elements in this row
  for (int k = d_row_pointers[i]; k < d_row_pointers[i + 1]; k++)
  {
    // if j is a column index, return the value
    if (d_column_indices[k] == j) return d_value[k];
  }
  // otherwise, this is a zeros
  return 0;
}

//---------------------------------------------------------------------------//
// MULTIPLY
//---------------------------------------------------------------------------//

template <class T>
inline void Matrix<T>::multiply(const Vector<T> &x, Vector<T> &y)
{
  Require(d_is_ready);
  Require(x.size() == d_n);
  Require(y.size() == d_m);
  // clear the output vector
  y.scale(0);
  // row pointer
  int p = 0;
  // for all rows
  #pragma omp for private(p)
  for (int i = 0; i < d_m; i++)
  {
    T temp = y[i];
    // for all columns
    for (p = d_row_pointers[i]; p < d_row_pointers[i + 1]; p++)
    {
      int j = d_column_indices[p];
      temp += x[j] * d_value[p];
    }
    y[i] = temp;
  }
}

// \todo good threading option needed
template <class T>
inline void Matrix<T>::multiply_transpose(const Vector<T> &x, Vector<T> &y)
{
  Require(d_is_ready);
  Require(x.size() == d_m);
  Require(y.size() == d_n);
  // clear the output vector
  y.scale(0);
  // for all rows (now columns)
  for (int i = 0; i < d_m; i++)
  {
    // for all columns (now rows)
    for (int p = d_row_pointers[i]; p < d_row_pointers[i + 1]; p++)
    {
      int j = d_column_indices[p];
      y[j] += x[i] * d_value[p];
    }
  }
}

//---------------------------------------------------------------------------//
// INSERTING VALUES
//---------------------------------------------------------------------------//

template <class T>
inline bool Matrix<T>::insert(int i, int j, T v)
{
  Require(!d_is_ready);
  Require(d_allocated);
  Require(i >= 0 and i < d_m);
  Require(j >= 0 and j < d_n);
  // return if storage unavailable
  if (d_counter >= d_nnz) return false;
  // otherwise, add the entry
  d_aij[d_counter].i = i;
  d_aij[d_counter].j = j;
  d_aij[d_counter].v = v;
  ++d_counter;
  return true;
}

template <class T>
inline bool Matrix<T>::insert(int i, int *j, T *v, int n)
{
  Require(!d_is_ready);
  Require(d_allocated);
  Require(i >= 0 and i < d_m);
  // return if storage unavailable
  if (d_counter + n >= d_nnz) return false;
  // otherwise, add the entries
  for (int jj = 0; jj < n; ++jj)
  {
    Require(j[jj] >= 0 and j[jj] < d_n);
    d_aij[d_counter].i = i;
    d_aij[d_counter].j = j[jj];
    d_aij[d_counter].v = v[jj];
    ++d_counter;
  }
  return true;
}

template <class T>
inline bool Matrix<T>::insert(int *i, int j, T *v, int n)
{
  Require(!d_is_ready);
  Require(d_allocated);
  Require(j >= 0 and j < d_n);
  // return if storage unavailable
  if (d_counter + n >= d_nnz) return false;
  // otherwise, add the entries
  for (int ii = 0; ii < n; ++ii)
  {
    Require(i[ii] >= 0 and i[ii] < d_m);
    d_aij[d_counter].i = i[ii];
    d_aij[d_counter].j = j;
    d_aij[d_counter].v = v[ii];
    ++d_counter;
  }
  return true;
}

template <class T>
inline bool Matrix<T>::insert(int *i, int *j, T *v, int n)
{
  Require(!d_is_ready);
  Require(d_allocated);
  // return if storage unavailable
  if (d_counter + n >= d_nnz) return false;
  // otherwise, add the entries
  for (int k = 0; k < n; ++k)
  {
    Require(i[k] >= 0 and i[k] < d_m);
    Require(j[k] >= 0 and j[k] < d_n);
    d_aij[d_counter].i = i[k];
    d_aij[d_counter].j = j[k];
    d_aij[d_counter].v = v[k];
    ++d_counter;
  }
  return true;
}

//---------------------------------------------------------------------------//
// IO
//---------------------------------------------------------------------------//

template <class T>
inline void Matrix<T>::display() const
{
  Require(d_is_ready);
  printf(" CSR matrix \n");
  printf(" ---------------------------\n");
  printf("      number rows = %5i \n",   d_m);
  printf("   number columns = %5i \n",   d_n);
  printf("      stored size = %5i \n\n", d_nnz);
  printf("\n");
  if (d_m > 20 or d_n > 20)
  {
    printf("  *** matrix not printed for m or n > 20 *** ");
    return;
  }
  for (int i = 0; i < d_m; i++)
  {
    printf(" row  %3i | ", i);
    for (int p = d_row_pointers[i]; p < d_row_pointers[i + 1]; p++)
    {
      int j = d_column_indices[p];
      T   v = d_value[p];
      printf(" %3i (%13.6e)", j, v);
    }
    printf("\n");
  }
  printf("\n");
}

} // end namespace callow


#endif /* MATRIX_I_HH_ */
