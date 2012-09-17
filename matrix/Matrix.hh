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
#include <vector>

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
 * We're keeping this simple to use.  The user must provide
 * the number of nonzeros for the matrix.  The actual number
 * can be less than this, but no more (unless we implement
 * a size increase mechanism later)
 *
 * During the construction process,
 * a COO (row, column, value) format is used.
 * This allows the user to add entries
 * one at a time, a row at a time, a column at a time, or
 * several rcv triples at a time. At the construction process,
 * the  storage is streamlined into CSR format, ensuring that
 * entries are stored by row and column and that a diagonal
 * entry exists.  That latter is required for things like
 * the \ref Jacobi or \ref GaussSeidel solvers, along with
 * certain preconditioner types.
 *
 */
template <class T>
class Matrix: public MatrixBase<T>
{

public:

  //---------------------------------------------------------------------------//
  // TYPEDEFS
  //---------------------------------------------------------------------------//

  typedef SP<Matrix<T> >  SP_matrix;
  typedef struct{ int i; int j; T v;} triplet;

  //---------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //---------------------------------------------------------------------------//

  // construction with deferred sizing and allocation
  Matrix();
  // construction with sizing but deferred allocation
  Matrix(const int m, const int n);
  // construction with sizing and allocation
  Matrix(const int m, const int n, const int nnz);
  // destructor
  virtual ~Matrix();

  //---------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //---------------------------------------------------------------------------//

  /// create memory; sizes must be set first
  void preallocate(const int nnz);
  /// add one value (return false if can't add)
  bool insert(int  i, int  j, T  v);
  /// add n values to a row  (return false if can't add)
  bool insert(int  i, int *j, T *v, int n);
  /// add n values to a column  (return false if can't add)
  bool insert(int *i, int  j, T *v, int n);
  /// add n triplets  (return false if can't add)
  bool insert(int *i, int *j, T *v, int n);
  /// number of rows
  int number_rows() const { return d_m; }
  /// number of columns
  int number_columns() const { return d_n; }
  /// number of nonzeros
  int number_nonzeros() const { return d_nnz; }
  /// starting index for a row
  int start(const int i) const;
  /// diagonal index for a row
  int diagonal(const int i) const;
  /// ending index for a row
  int end(const int i) const ;
  /// column index from cardinal index
  int column(const int p) const;
  /// value at a cardinal index
  T   operator[](const int p) const;
  /// value at ij and returns 0 if not present
  T operator()(const int i, const int j) const;
  // get underlying storage and indexing. careful!
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
  /// pointer to diagonal index of each row
  int* d_diagonal;
  /// number of nonzeros
  int d_nnz;
  /// are we allocated?
  bool d_allocated;
  // temporaries
  triplet* d_aij;
  // counts entries added
  int d_counter;

  //---------------------------------------------------------------------------//
  // IMPLEMENTATION
  //---------------------------------------------------------------------------//

  /// internal preallocation
  void preallocate();
  /// sort portion of triplets
  void sort_triplets(triplet* aij, const int n, bool flag=true);

};

} // end namespace callow

// Inline members
#include "Matrix.i.hh"

#endif /* callow_MATRIX_HH_ */
