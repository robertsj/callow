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
#include <ostream>
#include <vector>
#include <string>

namespace callow
{

/// Structure for storing COO matrix
template <class T>
struct triplet
{
  int i;
  int j;
  T   v;
  // i, j initialized to -1 to ensure they get sorted out
  triplet(int ii=-1, int jj=-1, T vv=0.0)
    : i(ii), j(jj), v(vv)
  {}
};

template <class T>
inline bool compare_triplet(const triplet<T> &x, const triplet<T> &y)
{
  // leave -1's on the right
  if (x.j == -1) return false;
  if (y.j == -1) return true;
  // otherwise, order from low to high j
  return x.j < y.j;
}

template <class T>
inline std::ostream& operator<< (std::ostream &out, triplet<T> s)
{
  out << "(i = " << s.i << ", j = " << s.j << ", v = " << s.v << ")";
  return out;
}

/*!
 *  \class Matrix
 *  \brief CRS matrix
 *
 * This is the base matrix class used within callow.  It implements
 * a compressed row storage (CRS) matrix.  An example of what this
 * means is as follows.
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
 * We're keeping this simple to use.  The use must specify
 * the number of nonzero entries per row, either via one
 * value applied to all rows or an array of values for
 * each row.  The reason the row storage is needed rather
 * than say the total number of nonzero entries is to make
 * construction easier.  During construction, the matrix
 * is stored (temporarily) in coordinate (COO) format, i.e.
 * a list of (i, j, value) triplets.  If this were stored
 * in one monolithic array, the construction of the CSR
 * structure would require we sort all the triplets by
 * row and then by column (since we want explicit access
 * to L, D, and U).  Initial testing proved that sorting
 * is just too time consuming, even with an n*log(n) method
 * like quicksort.
 *
 * Consequently, we do keep (i, j, value) triplets, but they
 * are stored by row, and to store by row, we need an initial
 * guess of how many entries there are.  For now, the size
 * can not be increased, but that would not be too difficult.
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
  typedef triplet<T>      triplet_T;

  //---------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //---------------------------------------------------------------------------//

  // construction with deferred sizing and allocation
  Matrix();
  // construction with sizing but deferred allocation
  Matrix(const int m, const int n);
  // construction with sizing and allocation
  Matrix(const int m, const int n, const int nnz);
  // copy constructor
  Matrix(Matrix<T> &A);
  // destructor
  virtual ~Matrix();

  //---------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //---------------------------------------------------------------------------//

  /// allocate using constant row size
  void preallocate(const int nnz_row);
  /// allocate using variable row size
  void preallocate(int *nnz_rows);

  /// add one value (return false if can't add)
  bool insert(int  i, int  j, T  v);
  /// add n values to a row  (return false if can't add)
  bool insert(int  i, int *j, T *v, int n);
  /// add n values to a column  (return false if can't add)
  bool insert(int *i, int  j, T *v, int n);
  /// add n triplets  (return false if can't add)
  bool insert(int *i, int *j, T *v, int n);


  /// starting index for a row
  int start(const int i) const;
  /// diagonal index for a row
  int diagonal(const int i) const;
  /// ending index for a row
  int end(const int i) const ;
  /// column index from cardinal index
  int column(const int p) const;

  /// value at a cardinal index
  T operator[](const int p) const;
  /// value at ij and returns 0 if not present
  T operator()(const int i, const int j) const;

  // get underlying storage and indexing. careful!
  T*   values()    {return d_values;}
  int* columns()   {return d_columns;}
  int* rows()      {return d_rows;}
  int* diagonals() {return d_diagonals;}

  /// number of nonzeros
  int number_nonzeros() const { return d_nnz; }
  /// is memory allocated?
  bool allocated() const {return d_allocated;}
  /// print (i, j, v) to ascii file with 1-based indexing for matlab
  void print_matlab(std::string filename = "matrix.out") const;

  //---------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MATRICES MUST IMPLEMENT
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

#ifdef CALLOW_ENABLE_PETSC
  using MatrixBase<T>::d_petsc_matrix;
#endif
  /// matrix elements
  T* d_values;
  /// column indices
  int* d_columns;
  /// row pointers
  int* d_rows;
  /// pointer to diagonal index of each row
  int* d_diagonals;
  /// number of nonzeros
  int d_nnz;
  /// are we allocated?
  bool d_allocated;
  // temporaries [number of rows][nonzeros per row]
  std::vector<std::vector<triplet_T> > d_aij;
  // counts entries added per row
  std::vector<int> d_counter;

  //---------------------------------------------------------------------------//
  // IMPLEMENTATION
  //---------------------------------------------------------------------------//

  /// internal preallocation
  void preallocate();

};



} // end namespace callow

// Inline members
#include "Matrix.i.hh"

#endif /* callow_MATRIX_HH_ */
