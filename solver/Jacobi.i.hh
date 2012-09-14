//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Jacobi.i.hh
 * \brief  Jacobi inline member definitions
 * \author Jeremy Roberts
 * \date   Sep 13, 2012
 */
//---------------------------------------------------------------------------//

#ifndef JACOBI_I_HH_
#define JACOBI_I_HH_

#include "matrix/Matrix.hh"
#include <cmath>
#include <cstdio>

namespace callow
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR & DESTRUCTOR
//---------------------------------------------------------------------------//

template <class T>
Jacobi<T>::Jacobi(const double  atol,
                  const double  rtol,
                  const int     maxit)
  : LinearSolver<T>(atol, rtol, maxit)
{
  /* ... */
}

//---------------------------------------------------------------------------//
// SOLVE
//---------------------------------------------------------------------------//


template <class T>
inline int Jacobi<T>::solve(const Vector<T> &b, Vector<T> &x)
{
  Require(x.size() == b.size());
  Require(x.size() == d_A->number_rows());
  Insist(dynamic_cast< Matrix<T>* >(d_A.bp()),
    "Need an explicit matrix for use with Jacobi iteration");
  typename Matrix<T>::SP_matrix A = d_A;

  typedef Vector<T> Vec;

  // temporary storage and pointers for swapping
  Vec temp(x.size(), 0.0);
  Vec* x0 = &x;
  Vec* x1 = &temp;
  Vec* swap;

  // iteration count
  int &iteration = d_number_iterations;

  // compute initial residual Ax - b and its norm
  A->multiply((*x0), (*x1));
  x1->scale(d_omega);
  double r0 = x1->norm_residual(b, Vec::L2);
  d_L2_residual[0] = r0;
  if (r0 < d_absolute_tolerance) return Base::SUCCESS;
  if (d_monitor)
    printf("iteration: %5i    residual: %12.8e \n", 0, r0);


  // perform iterations
  int status = Base::MAXIT;
  for (int iteration = 1; iteration < d_maximum_iterations; ++iteration)
  {

    //---------------------------------------------------//
    // compute X1 <-- inv(D)*(L+U)*X0 + inv(D)*b
    //---------------------------------------------------//

    T* a = A->value();
    for (int i = 0; i < A->number_rows(); i++)
    {
      T v = 0;
      int p = A->start(i);
      int d = A->diagonal(i);
      // L * X0
      for (; p < d; ++p)
        v += a[p] * (*x0)[A->column(p)];
      ++p; // skip diagonal
      // U * X0
      for (; p < A->end(i); ++p)
        v += a[p] * (*x0)[A->column(p)];
      (*x1)[i] = (b[i] - v) / a[d];
    }
    a = 0;

    //---------------------------------------------------//
    // compute residual norm
    //---------------------------------------------------//

    double r = x1->norm_residual(*x0, Vec::L2);
    d_L2_residual[iteration] = r;

    //---------------------------------------------------//
    // swap pointers
    //---------------------------------------------------//
    swap = x0;
    x0   = x1;
    x1   = swap;

    //---------------------------------------------------//
    // check convergence
    //---------------------------------------------------//

    if (d_monitor)
      printf("iteration: %5i    residual: %12.8e \n", iteration, r);

    if (r < std::max(d_relative_tolerance * r0, d_absolute_tolerance))
    {
      if (d_monitor)
        printf("*** Jacobi converged in %5i iterations with a residual of %12.8e \n", iteration, r );
      status = Base::SUCCESS;
      break;
    }
    else if (r - d_L2_residual[iteration - 1] > 0.0 and iteration > 1)
    {
      if (d_monitor)
        printf("*** Jacobi diverged at iteration %5i with a residual of %12.8e \n", iteration, r );
      status = Base::DIVERGE;
      break;
    }
  }
  if (status == Base::MAXIT and d_monitor)
    printf("*** Jacobi did not converge within the maximum number of iterations \n");

  // copy into the solution vector if needed
  if (x0 != &x) x = *x0;

  return status;
}

} // end namespace callow
#endif // JACOBI_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file Jacobi.i.hh
//---------------------------------------------------------------------------//
