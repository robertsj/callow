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

}

template <class T>
Matrix<T>::Matrix(const int m, const int n, const int nnz_per_row)
  : MatrixBase<T>(m, n)
  , d_allocated(false)
{
  preallocate(nnz_per_row);
}

template <class T>
Matrix<T>::Matrix(const int m, const int n, int* nnz_per_row)
  : MatrixBase<T>(m, n)
  , d_allocated(false)
{
  preallocate(nnz_per_row);
}

template <class T>
Matrix<T>::~Matrix()
{
  if (d_allocated)
  {
    delete [] d_value;
    delete [] d_column_indices;
    delete [] d_row_pointers;
    delete [] d_nnz_per_row;
    delete [] d_count_per_row;
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
  Require(p >= 0 and p < d_total_nnz);
  return d_column_indices[p];
}

template <class T>
inline T Matrix<T>::operator[](const int p) const
{
  Require(d_is_ready);
  Require(p >= 0 and p < d_total_nnz);
  return d_value[p];
}

//---------------------------------------------------------------------------//
// PREALLOCATION
//---------------------------------------------------------------------------//

template <class T>
inline void Matrix<T>::preallocate(const int nnz_per_row)
{
  Require(d_sizes_set);
  Require(nnz_per_row > 0);
  Require(nnz_per_row <= d_m);
  // total nonzeros
  d_total_nnz = nnz_per_row * d_m;
  d_nnz_per_row = new int[d_m];
  for (int i = 0; i < d_m; i++) d_nnz_per_row[i] = nnz_per_row;
  preallocate();
}


template <class T>
inline void Matrix<T>::preallocate(int* nnz_per_row)
{
  Require(d_sizes_set);
  // total nonzeros
  d_total_nnz = 0;
  for (int i = 0; i < d_m; i++)
  {
    int nnz = nnz_per_row[i];
    Assert(nnz >= 0 and nnz < d_m);
    d_total_nnz += nnz;
  }
  Assert(d_total_nnz > 0);
  d_nnz_per_row = nnz_per_row;
  preallocate();
}

template <class T>
inline void Matrix<T>::preallocate()
{
  d_value          = new T[d_total_nnz];
  d_column_indices = new int[d_total_nnz];
  d_row_pointers   = new int[d_m + 1];
  d_count_per_row  = new int[d_m];
  // initialize value and columns to zero
  for (int i = 0; i < d_total_nnz; i++) d_value[i] = 0;
  for (int i = 0; i < d_total_nnz; i++) d_column_indices[i] = 0;
  // set row pointers
  d_row_pointers[0] = 0;
  for (int i = 1; i <= d_m; i++)
  {
    d_row_pointers[i] = d_row_pointers[i - 1] + d_nnz_per_row[i - 1];
    d_count_per_row[i - 1] = 0;
  }
  d_allocated = true;
}

//---------------------------------------------------------------------------//
// ORDERING AND SUCH
//---------------------------------------------------------------------------//

/*
 *  We want this matrix to be available for use in Jacobia and
 *  Gauss-Seidel iteration.  To facilitate this, we need to
 *  know something about the L, D, and U structure.  This
 *  routine goes through, sorts the rows by column, truncates
 *  zeros, and sets a pointer to the diagonal index.  This can be used
 *  to iterate on L, D, or U.
 *
 *  This DIES if the central diagonal does not exist!!  Note
 *
 */
template <class T>
inline void Matrix<T>::assemble()
{
  Require(d_allocated);
  d_diagonal = new int[d_m];
  std::string msg = "The diagonal has to be included, even if zero: ";
  std::stringstream buffer;
  bool have_diag = true;
  for (int i = 0; i < d_m; i++)
  {
    d_diagonal[i] = 0;
    int j_last = 0;
    if (d_m == d_n) have_diag = false;
    for (int p = d_row_pointers[i]; p < d_row_pointers[i + 1]; ++p)
    {
      int j = d_column_indices[p];
      Insist( (j >= j_last) or j == 0,
             "Column indices must be monotonic increasing within a row.");
      if (j == i and !have_diag)
      {
        d_diagonal[i] = p;
        have_diag = true;
      }
      j_last = j;
    }
    if (!have_diag and d_m == d_n) buffer << " Row " << i << std::endl;
    Insist(have_diag, msg + buffer.str());
  }
  d_is_ready = true;
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

  // for all rows
  for (int i = 0; i < d_m; i++)
  {
    // for all columns
    for (int p = d_row_pointers[i]; p < d_row_pointers[i + 1]; p++)
    {
      int j = d_column_indices[p];
      y[i] += x[j] * d_value[p];
    }
  }

}

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
inline void Matrix<T>::add_single(int i, int j, T v)
{
  Require(!d_is_ready);
  Require(d_allocated);
  Require(i >= 0 and i < d_m);
  Require(j >= 0 and j < d_n);
  Require(d_count_per_row[i] < d_nnz_per_row[i]);
  int offset = d_count_per_row[i];
  d_column_indices[d_row_pointers[i] + offset] = j;
  d_value[d_row_pointers[i] + offset] = v;
  ++d_count_per_row[i];
}

template <class T>
inline void Matrix<T>::add_row(int i, int *j, T *v, int n)
{
  Require(!d_is_ready);
  Require(d_allocated);
  Require(i >= 0 and i < d_m);
  Require(d_count_per_row[i] == 0);
  // allow adding only partial rows
  int j_size = d_nnz_per_row[i];
  if (n != 0) j_size = n;
  for (int k = 0; k < j_size; k++)
  {
    Assert(j[k] < d_n);
    d_column_indices[d_row_pointers[i] + k] = j[k];
    d_value[d_row_pointers[i] + k] = v[k];
  }
  // but disallow further additions
  d_count_per_row[i] = d_nnz_per_row[i];
}

template <class T>
inline void Matrix<T>::add_csr(int *i, int *j, T *v)
{
  Require(!d_is_ready);
  Require(d_allocated);
  // do copies so that client is responsible for their
  // input arrays
  for (int k = 0; k <= d_m; k++) d_row_pointers[k] = i[k];
  for (int k = 0; k < d_total_nnz; k++)
  {
    Assert(j[k] < d_n);
    d_column_indices[k] = j[k];
    d_value[k] = v[k];
  }
}

//---------------------------------------------------------------------------//
// IO
//---------------------------------------------------------------------------//

template <class T>
inline void Matrix<T>::display() const
{
  Require(d_is_ready);
  std::cout << " CSR matrix " << std::endl;
  std::cout << " ---------------------------" << std::endl;
  std::cout << "      number rows = " << d_m << std::endl;
  std::cout << "   number columns = " << d_m << std::endl;
  std::cout << "   allocated size = " << d_total_nnz << std::endl;
  std::cout << std::endl;
  for (int i = 0; i < d_m; i++)
  {
    std::cout << " row " << i << " | ";
    for (int p = d_row_pointers[i]; p < d_row_pointers[i + 1]; p++)
    {
      int j = d_column_indices[p];
      T   v = d_value[p];
      std::cout << " " << j << "(" << v << ")";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

//---------------------------------------------------------------------------//
// QUERY
//---------------------------------------------------------------------------//

template <class T>
inline int Matrix<T>::nnz(const int i) const
{
  Require(d_allocated);
  Require(i >= 0 and i < d_m);
  return d_nnz_per_row[i];
}

template <class T>
inline int Matrix<T>::nnz() const
{
  Require(d_allocated);
  return d_total_nnz;
}

} // end namespace callow


#endif /* MATRIX_I_HH_ */
