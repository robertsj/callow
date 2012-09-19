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
Matrix<T>::Matrix(const int m, const int n, const int nnzrow)
  : MatrixBase<T>(m, n)
  , d_allocated(false)
{
  preallocate(nnzrow);
}

template <class T>
Matrix<T>::Matrix(Matrix<T> &A)
  : MatrixBase<T>(A.number_rows(), A.number_columns())
{
  d_nnz = A.number_nonzeros();
  d_row_pointers = new int[d_m + 1];
  d_column_indices = new int[d_nnz];
  d_value = new T[d_nnz];
  d_diagonal = new int[d_m];
  d_row_pointers[0] = 0.0;
  for (int i = 0; i < d_m; ++i)
  {
    d_row_pointers[i+1] = A.row_pointer()[i+1];
    d_diagonal[i] = A.diagonal_indices()[i];
  }
  for (int i = 0; i < d_nnz; ++i)
  {
    d_column_indices[i] = A.column_indices()[i+1];
    d_value[i] = A.value()[i];
  }
  d_allocated = true;
  d_is_ready = true;
}

template <class T>
Matrix<T>::~Matrix()
{
  if (d_is_ready)
  {
#ifdef CALLOW_ENABLE_PETSC
    // destroy the petsc matrix.  note, since we constructed it
    // using our own arrays, those still need to be deleted.
    MatDestroy(&d_petsc_matrix);
#endif
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
inline void Matrix<T>::preallocate(const int nnzrow)
{
  // preconditions
  Require(d_sizes_set);
  Require(!d_allocated);
  Require(nnzrow > 0);
  // set number of nonzeros, preallocate triplets, and initialize counter
  d_aij.resize(d_m, std::vector<triplet_T>(nnzrow));
  d_counter.resize(d_m, 0);
  d_allocated = true;
}

template <class T>
inline void Matrix<T>::preallocate(int* nnzrows)
{
  // preconditions
  Require(d_sizes_set);
  Require(!d_allocated);
  // set number of nonzeros, preallocate triplets, and initialize counter
  d_aij.resize(d_m);
  for (int i = 0; i < d_m; i++)
  {
    Require(nnzrows[i] >= 0);
    d_aij[i].resize(nnzrows[i]);
  }
  d_counter.resize(d_m, 0);
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

  typedef typename std::vector<triplet_T>::iterator it_T;

  // sort the coo structure
  d_nnz = 0;
  for (int i = 0; i < d_m; i++)
  {
//    std::cout << " before sort" << std::endl;
//    for (int j = 0; j < d_aij[i].size(); ++j)
//      std::cout << d_aij[i][j] << std::endl;

    // do the sort
    std::sort(d_aij[i].begin(), d_aij[i].end(), compare_triplet<T>);

//    std::cout << " after sort" << std::endl;
//    for (int j = 0; j < d_aij[i].size(); ++j)
//      std::cout << d_aij[i][j] << std::endl;

    // remove empty entries (if not entered, i,j==-1)
    while ( (d_aij[i].end()-1)->j == -1 ) d_aij[i].pop_back();

//    std::cout << " after pop" << std::endl;
//    for (int j = 0; j < d_aij[i].size(); ++j)
//      std::cout << d_aij[i][j] << std::endl;

    // find and/or insert the diagonal
    int j = 0;
    for (int k = 0; k < d_aij[i].size(); ++k)
      if (d_aij[i][k].j <= d_aij[i][k].i) ++j;
    if (j) --j;
    if (d_aij[i][j].i != d_aij[i][j].j and i < d_n)
      d_aij[i].insert(d_aij[i].begin()+j, triplet_T(i, i, 0.0));


//    std::cout << " after diag" << std::endl;
//    for (int j = 0; j < d_aij[i].size(); ++j)
//      std::cout << d_aij[i][j] << std::endl;

    // update the total nonzero count
    d_nnz += d_aij[i].size();
  }

  // allocate
  d_row_pointers   = new int[d_m + 1];
  d_column_indices = new int[d_nnz];
  d_value          = new T[d_nnz];
  d_diagonal       = new int[d_m];

  // fill the csr storage
  int aij_idx = 0;
  int val_idx = 0;
  d_row_pointers[0] = 0;

  // pointer to value index
  int p = 0;
  // for all rows
  for (int i = 0; i < d_m; i++)
  {
    // for all columns in the row
    for (int j = 0; j < d_aij[i].size(); ++j, ++p)
    {
      d_column_indices[p] = d_aij[i][j].j;
      d_value[p] = d_aij[i][j].v;
      // store diagonal index
      if (d_column_indices[p] == i) d_diagonal[i] = p;
    }
    d_row_pointers[i + 1] = p;
    // delete this row of coo storage
    d_aij[i].clear();
  }
  // delete the coo storage and specify the matrix is set to use
  d_aij.clear();

#ifdef CALLOW_ENABLE_PETSC
  PetscErrorCode ierr;
  ierr = MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, d_m, d_n,
                                   d_row_pointers, d_column_indices,
                                   d_value, &d_petsc_matrix);
  Assert(!ierr);
#endif
  d_is_ready = true;
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
  if (i == j)
    return d_value[d_diagonal[i]];
  else if (i < j)
  {
    for (int k = d_row_pointers[i]; k < d_diagonal[i]; k++)
      if (d_column_indices[k] == j) return d_value[k];
  }
  else
  {
    for (int k = d_diagonal[i] + 1; k < d_row_pointers[i + 1]; k++)
      if (d_column_indices[k] == j) return d_value[k];
  }
  // otherwise, this is a zeros
  return 0.0;
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
#ifdef CALLOW_ENABLE_PETSC_OPS
  MatMult(d_petsc_matrix, const_cast<Vector<T>* >(&x)->petsc_vector(), y.petsc_vector());
#else
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
#endif
}

// \todo good threading option needed
template <class T>
inline void Matrix<T>::multiply_transpose(const Vector<T> &x, Vector<T> &y)
{
  Require(d_is_ready);
  Require(x.size() == d_m);
  Require(y.size() == d_n);
#ifdef CALLOW_ENABLE_PETSC_OPS
  MatMultTranspose(d_petsc_matrix, const_cast<Vector<T>* >(&x)->petsc_vector(), y.petsc_vector());
#else
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
#endif
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
  if (d_counter[i] >= d_aij[i].size()) return false;
  // otherwise, add the entry
  d_aij[i][d_counter[i]].i = i;
  d_aij[i][d_counter[i]].j = j;
  d_aij[i][d_counter[i]].v = v;
  ++d_counter[i];
  return true;
}

template <class T>
inline bool Matrix<T>::insert(int i, int *j, T *v, int n)
{
  Require(!d_is_ready);
  Require(d_allocated);
  Require(i >= 0 and i < d_m);
  // return if storage unavailable
  //std::cout << " i = " << i << " count = " << d_counter[i] << " n = " << n << " size = " << d_aij[i].size() << std::endl;
  if (d_counter[i] + n >  d_aij[i].size()) return false;
  // otherwise, add the entries
  for (int jj = 0; jj < n; ++jj)
  {
    Require(j[jj] >= 0 and j[jj] < d_n);
    d_aij[i][d_counter[i]].i = i;
    d_aij[i][d_counter[i]].j = j[jj];
    d_aij[i][d_counter[i]].v = v[jj];
    ++d_counter[i];
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
  for (int ii = 0; ii < n; ++ii)
  {
    Require(i[ii] >= 0 and i[ii] < d_m);
    if (d_counter[i[ii]] + 1 > d_aij[i[ii]].size())
      return false;
  }
  // otherwise, add the entries
  for (int ii = 0; ii < n; ++ii)
  {
    d_aij[i[ii]][d_counter[i[ii]]].i = i[ii];
    d_aij[i[ii]][d_counter[i[ii]]].j = j;
    d_aij[i[ii]][d_counter[i[ii]]].v = v[ii];
    ++d_counter[i[ii]];
  }
  return true;
}

template <class T>
inline bool Matrix<T>::insert(int *i, int *j, T *v, int n)
{
  Require(!d_is_ready);
  Require(d_allocated);
  // return if storage unavailable
  // \todo this assumes one entry per row---fix
  for (int k = 0; k < n; ++k)
  {
    Require(i[k] >= 0 and i[k] < d_m);
    Require(j[k] >= 0 and j[k] < d_n);
    if (d_counter[i[k]] + 1 > d_aij[i[k]].size())
      return false;
  }
  // otherwise, add the entries
  for (int k = 0; k < n; ++k)
  {
    d_aij[i[k]][d_counter[i[k]]].i = i[k];
    d_aij[i[k]][d_counter[i[k]]].j = j[k];
    d_aij[i[k]][d_counter[i[k]]].v = v[k];
    ++d_counter[i[k]];
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
