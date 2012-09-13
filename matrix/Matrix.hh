//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Matrix.hh
 * \author robertsj
 * \date   Sep 13, 2012
 * \brief  Matrix class definition.
 */
//---------------------------------------------------------------------------//

#ifndef callow_MATRIX_HH_
#define callow_MATRIX_HH_

#include "MatrixBase.hh"

namespace callow
{

/*!
 *  \class Matrix
 *  \brief CRS matrix
 *
 *  Example:
 *
 *   | 7  0  0  2 |
 *   | 0  2  0  4 |
 *   | 1  0  0  0 |
 *   | 3  8  0  6 |
 *
 * value           = [7 2 2 4 1 3 8 6]
 * column indices  = [0 3 1 3 0 0 1 3]
 * row pointers    = [0 2 4 5 8]
 *
 * We're keeping this simple.  The user can define the number
 * of nonzeros for each row using one value per row or a
 * specific value for each row.  The size CANNOT change after.
 *
 * Inserting values can be done one at a time, by row, or by
 * inserting three arrays representing the csr format
 *
 * WARNING: For now, this really simple implementation assumes
 *          that entries are entered in increasing column order.
 *
 *          Also, it assumes that the diagonal entry exists.
 */
template <class T>
class Matrix: public MatrixBase<T>
{

public:

  //---------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //---------------------------------------------------------------------------//

  Matrix();
  Matrix(const int m, const int n);
  Matrix(const int m, const int n, const int nnz_per_row);
  Matrix(const int m, const int n, int* nnz_per_row);
  virtual ~Matrix();

  //---------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //---------------------------------------------------------------------------//

  // create memory; sizes must be set first
  void preallocate(const int nnz_per_row);
  void preallocate(int* nnz_per_row);

  // adding
  void add_single(int  i, int  j, T  v);
  void add_row(int i, int *j, T *v, int n = 0);
  void add_csr(int *i, int *j, T *v);

  // inquiry
  int number_rows() const { return d_m; }
  int number_columns() const { return d_n; }
  int number_nonzeros() const { return d_total_nnz; }

  int start(const int i) const;
  int diagonal(const int i) const;
  int end(const int i) const ;
  int column(const int p) const;
  T   operator[](const int p) const;

  // access to elements. returns 0 if ij not nonzero.
  T operator()(const int i, const int j) const;

  // get underlying storage and indexing. careful.
  T* value() {return d_value;}
  int* column_indices() {return d_column_indices;}
  int* row_pointer() {return d_row_pointers;}

  // is memory allocated?
  bool allocated() const { return d_allocated; }

  //---------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MATRIX OBJECTS MUST IMPLEMENT THESE
  //---------------------------------------------------------------------------//

  // postprocess storage
  void assemble();
  // action y <-- A * x
  void multiply(const Vector<T> &x,  Vector<T> &y);
  // action y <-- A' * x
  void multiply_transpose(const Vector<T> &x, Vector<T> &y);
  // pretty print to screen
  void display() const;

private:

  //---------------------------------------------------------------------------//
  // DATA
  //---------------------------------------------------------------------------//

  /// expose base members
  using MatrixBase<T>::d_m;
  using MatrixBase<T>::d_n;
  using MatrixBase<T>::d_sizes_set;
  using MatrixBase<T>::d_is_ready;

  /// matrix elements
  T* d_value;
  /// column indices
  int* d_column_indices;
  /// row pointers
  int* d_row_pointers;

  /// total number of nonzeros
  int d_total_nnz;
  /// number per row
  int* d_nnz_per_row;
  /// count per row. this makes sure we don't add more than we allocated.
  int* d_count_per_row;
  int* d_diagonal;
  /// are we allocated?
  bool d_allocated;

  //---------------------------------------------------------------------------//
  // IMPLEMENTATION
  //---------------------------------------------------------------------------//

  /// internal preallocation
  void preallocate();
  /// nnz per row i
  int nnz(const int i) const;
  /// total nnz
  int nnz() const;


};

} // end namespace callow

// Inline members
#include "Matrix.i.hh"

#endif /* callow_MATRIX_HH_ */
