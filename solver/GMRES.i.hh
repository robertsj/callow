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
#include <iostream>
using std::cout;
using std::endl;

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
  , d_reorthog(1)
{
  Insist(d_restart > 2, "Need a restart of > 2");
  d_H = new T*[(restart + 1)];
  for (int i = 0; i <= d_restart; i++)
  {
    d_H[i] = new T[restart];
    for (int j = 0; j < d_restart; j++)
      d_H[i][j] = 0.0;
  }
  d_c.resize(restart + 1);
  d_s.resize(restart + 1);
}

template <class T>
GMRES<T>::~GMRES()
{
  for (int i = 0; i <= d_restart; i++)
  {
    delete [] d_H[i];
  }
  delete [] d_H;
}
//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

template <class T>
inline void GMRES<T>::solve_impl(const Vector<T> &b, Vector<T> &x)
{
 cout << " solving gmres " << endl;
  typedef Vector<T> Vec_T;

  // krylov basis
  std::vector<Vec_T>  v(d_restart, Vec_T(x.size(), 0.0));

  // residual
  Vec_T r(x.size(), 0.0);

  // vector such that x = V*y
  Vec_T y(d_restart, 0.0);

  // vector such that g(1:k) = R*y --> x = V*inv(R)*g and |g(k+1)| is the residual
  Vec_T g(d_restart + 1, 0.0);

  // initialize c and s
  d_c.set(0.0);
  d_s.set(0.0);

  // outers
  int iteration = 0;
  bool done = false;
  // each iteration is an action of A
  while (!done and iteration < d_maximum_iterations)
  {
    // intialize H
    for (int i = 0; i < d_restart; i++)
    {
      v[i].set(0.0);
    }
    g.set(0.0);
    //y.set(0.0);

    // compute residual
    d_A->multiply(x, r);
    r.subtract(b);
    r.scale(-1.0);
    T rho = r.norm(Vec_T::L2);

    // check initial outer residual
    if (iteration == 0)
    {
      if (monitor_init(rho)) return;
      ++iteration;
    }
    else
    {
      if (d_monitor_output) cout << "restarting..." << endl;
    }

//    else
//    {
//      if (monitor(iteration, rho)) return;
//    }
    //iteration++;

    // initial krylov vector
    v[0].copy(r);
    v[0].scale(1.0 / rho);
    g[0] = rho;

    // inner iterations (of size restart)
    int k = 0;
    for (; k < d_restart - 1; ++k, ++iteration)
    {

      if (iteration >= d_maximum_iterations)
      {
        done = true;
        break;
      }

      // perform modified gram-schmidt
      d_A->multiply(v[k], v[k+1]);
      T norm_Av = v[k+1].norm();
      for (int j = 0; j <= k; ++j)
      {
        d_H[j][k] = v[k+1].dot(v[j]);
        v[k+1].add_a_times_x(-d_H[j][k], v[j]);
      }
      d_H[k+1][k] = v[k+1].norm(Vec_T::L2);
      T norm_Av_2 = d_H[k+1][k];

      // optional reorthogonalization
      if ( (d_reorthog == 1 and norm_Av + 0.001 * norm_Av_2 == norm_Av) or (d_reorthog == 2) )
      {
        cout << " reorthog ... " << endl;
        for (int j = 0; j < k; ++j)
        {
          T hr = v[j].dot(v[k+1]);
          d_H[j][k] += hr;
          v[k+1].add_a_times_x(-hr, v[j]);
        }
        d_H[k+1][k] = v[k+1].norm();
      }

      // careful of happy break down
      if (d_H[k+1][k] != 0.0)
      {
        v[k+1].scale(1.0/d_H[k+1][k]);
      }
      else
      {
        std::printf("happy breakdown for k = %5i (iteration = %5i) \n",
                    k, iteration);
      }

      // going to triangle form incrementally
      if (k > 0) apply_givens(k);
      double nu = std::sqrt(d_H[k][k]*d_H[k][k] + d_H[k+1][k]*d_H[k+1][k]);
      d_c[k] =  d_H[k  ][k] / nu;
      d_s[k] = -d_H[k+1][k] / nu;
      d_H[k  ][k] = d_c[k] * d_H[k][k] - d_s[k]*d_H[k+1][k];
      d_H[k+1][k] = 0.0;
      T g_0 = d_c[k]*g[k] - d_s[k]*g[k+1];
      T g_1 = d_s[k]*g[k] + d_c[k]*g[k+1];
      g[k  ] = g_0;
      g[k+1] = g_1;

      // residual
      rho = std::abs(g_1);
      if (monitor(iteration, rho))
      {
        done = true;
        break;
      }

    } // end inners

    // compute y
    compute_y(y, g, k);

    // update x = x0 + v[0]*y[0] + ...
    for (int i = 0; i < k; ++i)
    {
      x.add_a_times_x(y[i], v[i]);
    }


  } // end outers

  d_A->multiply(x, r);
  r.subtract(b);
  r.scale(-1.0);
  T resid = r.norm(Vec_T::L2);
  printf(" final residual norm is actually: %12.8e \n", resid);

  return;

}

template <class T>
inline void GMRES<T>::apply_givens(const int k)
{
  Require(k < d_restart);
  for (int i = 0; i < k; ++i)
  {
    T g_0 = d_c[i]*d_H[i][k] - d_s[i]*d_H[i+1][k];
    T g_1 = d_s[i]*d_H[i][k] + d_c[i]*d_H[i+1][k];
    d_H[i][k]   = g_0;
    d_H[i+1][k] = g_1;
  }
}

template <class T>
inline void GMRES<T>::compute_y(Vector<T> &y, const Vector<T> &g, const int k)
{
  // H is [m+1][m], though we may have only need for k*k
  // solves H[0:k][0:k]*y[0:k] = g[0:k]
  //  but only for H_ij for j>=i, i.e. upper triangle
  Require(k < d_restart);
  Require(y.size() >= k);
  Require(g.size() >= k);

  for (int i = k - 1; i >= 0; --i)
  {
    y[i] = g[i];
    for (int j = i + 1; j < k; ++j)
    {
      y[i] -= d_H[i][j] * y[j];
    }
    Assert(d_H[i][i] != 0.0);
    y[i] /= d_H[i][i];
  }

}

template <class T>
inline void GMRES<T>::initialize_H()
{
  for (int i = 0; i <= d_restart; i++)
    for (int j = 0; j < d_restart; j++)
      d_H[i][j] = 0.0;
}

} // end namespace callow




#endif /* GMRES_I_HH_ */
