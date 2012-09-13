//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LinearSolver.hh
 * \author robertsj
 * \date   Sep 13, 2012
 * \brief  LinearSolver class definition.
 */
//---------------------------------------------------------------------------//

#ifndef LINEARSOLVER_HH_
#define LINEARSOLVER_HH_

#include "matrix/MatrixBase.hh"
#include <vector>

namespace callow
{

/**
 *  \class LinearSolver
 *  \brief Base class for iterative linear solvers
 *
 *  We solve problems of the form
 *  \f[
 *      \mathbf{A}x = b
 *  \f]
 *  via iterative methods.  A system is "solved" when the
 *  norm of the residual is small enough or some maximum
 *  iteration count is reached.  The residual is defined
 *  \f[
 *      \mathbf{A}x - b
 *  \f]
 *  and its norm is
 *  \f[
 *    r =  || \mathbf{A}x - b ||
 *  \f]
 *  By default, the L2 norm is used, though L1 and Linf
 *  are also recorded can can be used.  This represents
 *  the absolute norm.  Sometimes it makes more sense to
 *  check the residual with respect to the initial norm,
 *  in which case the relative norm is
 *  \f[
 *      r_n / r_0 = || \mathbf{A}x^{n} - b || / || \mathbf{A}x^{0} - b ||
 *  \f]
 *  The iteration terminates when
 *  \f[
 *      r_n < \mathrm{max} ( \tau_{\mathrm{rel}} r_0, \tau_{\mathrm{abs}} )
 *  \f]
 *
 *  Currently planned are
 *    - Jacobi
 *    - Gauss-Seidel
 *    - SOR
 *    - GMRES(m)
 *
 *  possibly with preconditioning
 */
template<class T>
class LinearSolver
{

public:

  //-------------------------------------------------------------------------//
  // ENUMERATIONS
  //-------------------------------------------------------------------------//

  enum status
  {
    SUCCESS, MAXIT, DIVERGE
  };

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  LinearSolver(const double atol, const double rtol, const int maxit)
    : d_absolute_tolerance(atol)
    , d_relative_tolerance(rtol)
    , d_maximum_iterations(maxit)
    , d_L1_residual(maxit + 1, 0)
    , d_L2_residual(maxit + 1, 0)
    , d_LI_residual(maxit + 1, 0)
    , d_number_iterations(0)
    , d_monitor(false)
  {
    Require(d_absolute_tolerance > 0.0);
    Require(d_relative_tolerance > 0.0);
    Require(d_maximum_iterations > 0.0);
  }

  virtual ~LinearSolver(){}

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /**
   *  \param A  linear operator
   *  \param P  optional preconditioning process
   */
  void set_operators(MatrixBase<T> *A, MatrixBase<T> *P = 0)
  {
    Require(A);
    d_A = A;
    Ensure(d_A->number_rows() == d_A->number_columns());
  }

  void set_tolerances(const double atol, const double rtol, const int maxit)
  {
    d_absolute_tolerance = atol;
    d_relative_tolerance = rtol;
    d_maximum_iterations = maxit;
  }

  void set_monitor(const bool v)
  {
    d_monitor = v;
  }

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL LINEAR SOLVERS MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  /**
   *  \param b  right hand side
   *  \param x  unknown vector
   */
  virtual int solve(const Vector<T> &b, Vector<T> &x) = 0;

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  double d_absolute_tolerance;
  double d_relative_tolerance;
  int    d_maximum_iterations;
  std::vector<double> d_L1_residual;
  std::vector<double> d_L2_residual;
  std::vector<double> d_LI_residual;
  int    d_number_iterations;
  MatrixBase<T>* d_A;
  MatrixBase<T>* d_P;
  bool d_monitor;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

};

} // end namespace callow

#endif /* LINEARSOLVER_HH_ */
