//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Vector.i.hh
 * \author robertsj
 * \date   Sep 13, 2012
 * \brief  Vector inline member definitions
 */
//---------------------------------------------------------------------------//

#ifndef callow_VECTOR_I_HH_
#define callow_VECTOR_I_HH_

#include "utils/DBC.hh"
#include <cmath>
#include <iostream>

namespace callow
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR & DESTRUCTOR
//---------------------------------------------------------------------------//

template <class T>
Vector<T>::Vector()
  : d_size(0)
{
  /* ... */
}

template <class T>
Vector<T>::Vector(const int n, T v)
  : d_size(n),
    d_value(new T[n])
{
  // Preconditions
  Require(d_value);
  // Initialize value
  set(v);
}

template <class T>
Vector<T>::~Vector()
{
  if (d_size) delete [] d_value;
}

//---------------------------------------------------------------------------//
// ACCESS
//---------------------------------------------------------------------------//

template <class T>
inline const T& Vector<T>::operator[](const int i) const
{
  Require(i >= 0);
  Require(i < d_size);
  return d_value[i];
}

template <class T>
inline T& Vector<T>::operator[](const int i)
{
  Require(i >= 0);
  Require(i < d_size);
  return d_value[i];
}

//---------------------------------------------------------------------------//
// VECTOR OPERATIONS
//---------------------------------------------------------------------------//

template <class T>
inline T Vector<T>::dot(const Vector<T>& x)
{
  Require(d_size == x.size());
  T val = 0;
  for (int i = 0; i < d_size; i++)
    val += d_value[i] * x[i];
}

template <class T>
inline T Vector<T>::norm(const int type)
{
  T val = 0.0;
  if (type == L1)
  {
    for (int i = 0; i < d_size; i++)
      val += std::abs(d_value[i]);
  }
  else if (type == L2)
  {
    for (int i = 0; i < d_size; i++)
      val += d_value[i] * d_value[i];
    val = std::sqrt(val);
  }
  else if (type == LINF)
  {
    for (int i = 0; i < d_size; i++)
      val = std::max(val, std::abs(d_value[i]));
  }
  return val;
}

template <class T>
inline T Vector<T>::norm_residual(const Vector<T>& x, const int type)
{
  Require(d_size == x.size());
  T val = 0.0;
  if (type == L1)
  {
    for (int i = 0; i < d_size; i++)
      val += std::abs(d_value[i] - x[i]);
  }
  else if (type == L2)
  {
    for (int i = 0; i < d_size; i++)
      val += (d_value[i] - x[i])*(d_value[i] - x[i]);
    val = std::sqrt(val);
  }
  else if (type == LINF)
  {
    for (int i = 0; i < d_size; i++)
      val = std::max(val, std::abs(d_value[i] - x[i]));
  }
  return val;
}

template <class T>
inline void Vector<T>::set(const T v)
{
  for (int i = 0; i < d_size; i++)
    d_value[i] = v;
}

template <class T>
inline void Vector<T>::scale(const T v)
{
  for (int i = 0; i < d_size; i++)
    d_value[i] *= v;
}

template <class T>
inline void Vector<T>::add(const Vector &x)
{
  Require(x.size() == d_size);
  for (int i = 0; i < d_size; i++)
    d_value[i] += x[i];
}

template <class T>
inline void Vector<T>::subtract(const Vector &x)
{
  Require(x.size() == d_size);
  for (int i = 0; i < d_size; i++)
    d_value[i] -= x[i];
}

template <class T>
inline void Vector<T>::multiply(const Vector &x)
{
  Require(x.size() == d_size);
  for (int i = 0; i < d_size; i++)
    d_value[i] *= x[i];
}

template <class T>
inline void Vector<T>::divide(const Vector &x)
{
  Require(x.size() == d_size);
  for (int i = 0; i < d_size; i++)
    d_value[i] /= x[i];
}

//---------------------------------------------------------------------------//
// IO
//---------------------------------------------------------------------------//

template <class T>
void Vector<T>::display() const
{
  std::cout << " Vector " << std::endl;
  std::cout << " ---------------------------" << std::endl;
  std::cout << "      number rows = " << d_size << std::endl;
  std::cout << std::endl;
  for (int i = 0; i < d_size; i++)
  {
    std::cout << " row " << i << " | " << d_value[i] << std::endl;
  }
  std::cout << std::endl;

}

} // end namespace callow

#endif /* callow_VECTOR_I_HH_ */
