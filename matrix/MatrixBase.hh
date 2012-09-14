//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MatrixBase.hh
 * \author robertsj
 * \date   Sep 13, 2012
 * \brief  MatrixBase class definition.
 */
//---------------------------------------------------------------------------//

#ifndef callow_MATRIXBASE_HH_
#define callow_MATRIXBASE_HH_

#include "utils/SP.hh"
#include "vector/Vector.hh"

namespace callow
{

template <class T>
class MatrixBase
{

public:

  //---------------------------------------------------------------------------//
  // TYPEDEFS
  //---------------------------------------------------------------------------//

  typedef SP<MatrixBase<T> >    SP_matrix;

  //---------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //---------------------------------------------------------------------------//

  MatrixBase()
    : d_m(0)
    , d_n(0)
    , d_sizes_set(false)
    , d_is_ready(false)
  {
    /* ... */
  }

  MatrixBase(const int m, const int n)
    : d_m(m)
    , d_n(n)
    , d_is_ready(false)
  {
    Require(m > 0 and n > 0);
    d_sizes_set = true;
  }

  virtual ~MatrixBase(){}

  //---------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //---------------------------------------------------------------------------//

  void set_size(const int m, const int n)
  {
    Require(m > 0 and n > 0);
    d_m = m;
    d_n = n;
    d_sizes_set = true;
  }

  int number_rows() const { return d_m; }
  int number_columns() const { return d_n; }
  bool is_ready() const { return d_is_ready; }

  //---------------------------------------------------------------------------//
  // ABSTRACT INTERFACE
  //---------------------------------------------------------------------------//

  // postproces storage
  virtual void assemble() = 0;
  // action y <-- A * x
  virtual void multiply(const Vector<T> &x,  Vector<T> &y) = 0;
  // action y <-- A' * x
  virtual void multiply_transpose(const Vector<T> &x, Vector<T> &y) = 0;

  virtual void display() const = 0;



protected:

  //---------------------------------------------------------------------------//
  // DATA
  //---------------------------------------------------------------------------//

  /// number of rows
  int d_m;
  /// number of columns
  int d_n;
  /// are m and n set?
  bool d_sizes_set;
  /// am i good to go?
  bool d_is_ready;
};


} // end namespace callow

#endif /* callow_MATRIXBASE_HH_ */
