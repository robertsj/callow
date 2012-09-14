//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GMRES.i.hh
 * \author robertsj
 * \date   Sep 14, 2012
 * \brief  GMRES inline member definitions
 */
//---------------------------------------------------------------------------//

#ifndef GMRES_I_HH_
#define GMRES_I_HH_

#include "matrix/Matrix.hh"
#include <cmath>
#include <cstdio>

namespace callow
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR & DESTRUCTOR
//---------------------------------------------------------------------------//

template <class T>
GMRES<T>::GMRES(const double  atol,
                const double  rtol,
                const int     maxit,
                const int     restart)
  : LinearSolver<T>(atol, rtol, maxit, "GMRES")
  , d_restart(restart)
{
  Insist(d_restart > 2, "Need a restart of > 2");
  d_H = new T[restart * restart];
}

//---------------------------------------------------------------------------//
// SOLVE
//---------------------------------------------------------------------------//


template <class T>
inline void GMRES<T>::solve_impl(const Vector<T> &b, Vector<T> &x)
{

  typedef Vector<T> Vec;

  // temporary storage and pointers for swapping
  Vec r(x.size(), 0.0);

  // krylov subspace
  std::vector<Vec>  v(d_restart, Vec(x.size(), 0.0);

  // perform iterations
  iteration = 0;
  for (int outer = 1; outer <= d_maximum_iterations; ++iteration)
  {

    // r <-- inv(P) * (b - A * x0)
    d_A->multiply((*x0), r);
    r.subtract(b);
    r.scale(-1.0);

    // v(1) <-- r / ||r||_2
    T norm_r = r.norm(Vec::L2);
    v[1].copy(r);
    v[1].scale(1.0/norm_r);

    // perform inners
    for (int i = 0; i < d_restart; i++)
    {



    }


  }

  // copy into the solution vector if needed
  if (x0 != &x) x.copy(*x0);

  return;
}

} // end namespace callow




#endif /* GMRES_I_HH_ */
