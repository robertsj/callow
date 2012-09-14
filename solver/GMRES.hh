//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GMRES.hh
 * \author robertsj
 * \date   Sep 14, 2012
 * \brief  GMRES class definition.
 */
//---------------------------------------------------------------------------//

#ifndef GMRES_HH_
#define GMRES_HH_

#include "LinearSolver.hh"

namespace callow
{

/**
 *  \class GMRES
 *  \brief Uses preconditioned GMRES(m) iteration to solve a system
 *
 *  GMRES seeks to find the best solution \f$ x \f$
 *  to the linear system
 *  \f[
 *     \mathbf{A}x = b
 *  \f]
 *  such that \f$ x \in \mathcal{K}_n \f$, where the
 *  Krylov subspace is defined
 *  \f[
 *     \mathcal{K}_n \equiv
 *       [b, \mathbf{A}b, \mathbf{A}^2 b, \ldots, \mathbf{A}^{n-1} b] \, .
 *  \f]
 *  Specifically, GMRES finds \f$ x_n \f$ that satisfies
 *  \f[
 *     \min_{x_n \in \mathcal{K}_n} ||b-\mathbf{A}x_n||_2 \, ,
 *  \f]
 *  i.e. it finds the least squares fit within the current Krylov
 *  subspace at every iteration, stopping when that residual
 *  is small enough.
 *
 */
template<class T>
class GMRES: public LinearSolver<T>
{

public:

  typedef LinearSolver<T> Base;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  GMRES(const double atol, const double rtol, const int maxit,
        const int restart);

  virtual ~GMRES()
  {
    delete [] d_H;
  }

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  // expose base class members
  using LinearSolver<T>::status;
  using LinearSolver<T>::d_absolute_tolerance;
  using LinearSolver<T>::d_relative_tolerance;
  using LinearSolver<T>::d_maximum_iterations;
  using LinearSolver<T>::d_L1_residual;
  using LinearSolver<T>::d_L2_residual;
  using LinearSolver<T>::d_LI_residual;
  using LinearSolver<T>::d_number_iterations;
  using LinearSolver<T>::d_A;
  using LinearSolver<T>::d_P;

  /// maximum size of krylov subspace
  int d_restart;

  /// upper hessenberg (treated as dense)
  T d_H;

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL LINEAR SOLVERS MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  /**
   *  \param b  right hand side
   *  \param x  unknown vector
   */
  void solve_impl(const Vector<T> &b, Vector<T> &x);

};

} // end namespace callow

// Inline member definitions
#include "GMRES.i.hh"

#endif /* GMRES_HH_ */
