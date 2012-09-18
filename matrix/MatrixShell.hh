//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MatrixShell.hh
 * \author robertsj
 * \date   Sep 18, 2012
 * \brief  MatrixShell class definition.
 */
//---------------------------------------------------------------------------//

#ifndef callow_MATRIXSHELL_HH_
#define callow_MATRIXSHELL_HH_

namespace callow
{

/*!
 *  \class MatrixShell
 *  \brief Defines a matrix free operator
 *
 */
template <class T>
class MatrixShell: public MatrixBase<T>
{

public:

  //---------------------------------------------------------------------------//
  // TYPEDEFS
  //---------------------------------------------------------------------------//

  typedef SP<MatrixShell<T> >    SP_matrix;

  //---------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //---------------------------------------------------------------------------//

  MatrixShell()
    : d_m(0)
    , d_n(0)
    , d_sizes_set(false)
    , d_is_ready(false)
  {
    /* ... */
  }

  MatrixShell(const int m, const int n)
    : d_m(m)
    , d_n(n)
    , d_is_ready(false)
  {
    Require(m > 0 and n > 0);
    d_sizes_set = true;
  }

  virtual ~MatrixShell(){}

  //---------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MATRICES MUST IMPLEMENT
  //---------------------------------------------------------------------------//

  // default shell assemble does nothing
  void assemble()
  {
    /* ... */
  }
  // default shell display gives just the sizes
  void display() const
  {
    std::cout << "MatrixShell:" << std::endl
              << "  # rows = " << d_m << std::endl
              << "  # cols = " << d_n << std::endl
              << std::endl;
  }
  // the client must implement the action y <-- A * x
  virtual void multiply(const Vector<T> &x,  Vector<T> &y) = 0;
  // the client must implement the action y <-- A' * x
  virtual void multiply_transpose(const Vector<T> &x, Vector<T> &y) = 0;

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


#endif /* callow_MATRIXSHELL_HH_ */
