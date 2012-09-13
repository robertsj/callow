//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Richardson.i.hh
 * \author robertsj
 * \date   Sep 13, 2012
 * \brief  Richardson inline member definitions
 */
//---------------------------------------------------------------------------//

#ifndef RICHARDSON_I_HH_
#define RICHARDSON_I_HH_

#include <cmath>
#include <cstdio>

namespace callow
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR & DESTRUCTOR
//---------------------------------------------------------------------------//

template <class T>
Richardson<T>::Richardson(const double  atol,
                          const double  rtol,
                          const int     maxit,
                          const double  omega)
  : LinearSolver<T>(atol, rtol, maxit)
  , d_omega(omega)
{

}

//---------------------------------------------------------------------------//
// SOLVE
//---------------------------------------------------------------------------//


template <class T>
inline int Richardson<T>::solve(const Vector<T> &b, Vector<T> &x)
{
  Require(x.size() == b.size());
  Require(x.size() == d_A->number_rows());

  typedef Vector<T> Vec;

  // scale rhs by relaxation factor
  Vec B(b.size(), 0.0);
  B.add(b);
  B.scale(d_omega);

  // temporary storage and pointers for swapping
  Vec temp(x.size(), 0.0);
  Vec* x0 = &x;
  Vec* x1 = &temp;
  Vec* swap;

  // iteration count
  int &iteration = d_number_iterations;

  // compute initial residual w(Ax - b) and its norm
  d_A->multiply((*x0), (*x1));
  x1->scale(d_omega);
  double r0 = x1->norm_residual(B, Vec::L2);
  d_L2_residual[0] = r0;
  if (r0 < d_absolute_tolerance) return Base::SUCCESS;

  // perform iterations
  int status = Base::MAXIT;
  for (int iteration = 1; iteration < d_maximum_iterations; ++iteration)
  {

    //---------------------------------------------------//
    // compute X1 <-- (I - w*A) * X0 + w*b
    //---------------------------------------------------//

    // X1 <-- A * X0
    d_A->multiply((*x0), (*x1));
    // X1 <-- w * X1 = w * A * X0
    x1->scale(d_omega);
    // X1 <-- X1 - X0 =  (A - I) * X0
    x1->subtract(*x0);
    // X1 <-- X1 - b = (A - I) * X0 - b
    x1->subtract(B);
    // X1 <-- -X1 = (I - A) * X0 + b
    x1->scale(-1);

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
        printf("*** Richardson converged in %5i iterations with a residual of %12.8e \n", iteration, r );
      status = Base::SUCCESS;
      break;
    }
    else if (r - d_L2_residual[iteration - 1] > 0.0)
    {
      if (d_monitor)
        printf("*** Richardson diverged \n");
      status = Base::DIVERGE;
      break;
    }
  }
  if (status == Base::MAXIT and d_monitor)
    printf("*** Richardson did not converge within the maximum number of iterations \n");

  return status;
}

} // end namespace callow



#endif /* RICHARDSON_I_HH_ */
