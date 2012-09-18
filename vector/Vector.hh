//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Vector.hh
 * \author robertsj
 * \date   Sep 13, 2012
 * \brief  Vector class definition.
 */
//---------------------------------------------------------------------------//

#ifndef callow_VECTOR_HH_
#define callow_VECTOR_HH_

#include "callow_config.hh"
#include "utils/SP.hh"

namespace callow
{

/*!
 *  \class Vector
 *  \brief Dense vector object
 *
 *
 */
template <class T>
class Vector
{

public:

  //---------------------------------------------------------------------------//
  // ENUMERATIONS
  //---------------------------------------------------------------------------//

  enum vec_norm_types
  {
    L1, L2, LINF, END_VEC_NORM_TYPES
  };

  //---------------------------------------------------------------------------//
  // TYPEDEFS
  //---------------------------------------------------------------------------//

  typedef SP<Vector<T> >    SP_vector;

  //---------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //---------------------------------------------------------------------------//

  Vector();
  Vector(const int n, T v = 0);
  Vector(const Vector &x);
  Vector(Vector &x);
  virtual ~Vector();
  void resize(const int n);

  //---------------------------------------------------------------------------//
  // ACCESS
  //---------------------------------------------------------------------------//

  const T& operator[](const int i) const;
  T& operator[](const int i);
#ifdef CALLOW_ENABLE_PETSC
  Vec petsc_vector() {return d_petsc_vector;}
#endif

  //---------------------------------------------------------------------------//
  // VECTOR OPERATIONS
  //---------------------------------------------------------------------------//

  // scalar
  T dot(const Vector<T>& x);
  T norm(const int type = L2);
  T norm_residual(const Vector<T>& x, const int type = L2);

  // element-wise operations
  void set(const T v);
  void scale(const T v);
  void add(const Vector<T>& x);
  void subtract(const Vector<T>& x);
  void multiply(const Vector<T>& x);
  void divide(const Vector<T>& x);
  void copy(const Vector<T>& x);
  void add_a_times_x(const T a, const Vector<T>& x);

  //---------------------------------------------------------------------------//
  // QUERY
  //---------------------------------------------------------------------------//

  int size() const { return d_size; }
  void display() const;

protected:

  //---------------------------------------------------------------------------//
  // DATA
  //---------------------------------------------------------------------------//

  int d_size;
  T* d_value;
#ifdef CALLOW_ENABLE_PETSC
  Vec d_petsc_vector;
#endif

};

} // end namespace callow

#include "Vector.i.hh"

#endif /* callow_VECTOR_HH_ */
