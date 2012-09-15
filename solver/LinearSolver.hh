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
#include "utils/SP.hh"
#include <cstdio>
#include <string>
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
 *  possibly with preconditioning.
 *
 *  Note, some linear solvers require that the matrix provides L, U, and
 *  D operations.  Here, we simply require those solvers to have a \ref
 *  Matrix operator or subclasses so that the elements can be accessed
 *  directly.
 *
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
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef SP<LinearSolver<T> >                SP_solver;
  typedef typename MatrixBase<T>::SP_matrix   SP_matrix;
  typedef typename Vector<T>::SP_vector       SP_vector;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  LinearSolver(const double atol,
               const double rtol,
               const int maxit,
               std::string name = "solver")
    : d_absolute_tolerance(atol)
    , d_relative_tolerance(rtol)
    , d_maximum_iterations(maxit)
    , d_L1_residual(maxit + 1, 0)
    , d_L2_residual(maxit + 1, 0)
    , d_LI_residual(maxit + 1, 0)
    , d_number_iterations(0)
    , d_monitor_output(false)
    , d_monitor_diverge(true)
    , d_name(name)
  {
    Require(d_absolute_tolerance > 0.0);
    Require(d_relative_tolerance > 0.0);
    Require(d_maximum_iterations >= 0);
  }

  virtual ~LinearSolver(){}

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /**
   *  \param A  linear operator
   *  \param P  optional preconditioning process
   */
  void set_operators(SP_matrix A, SP_matrix P = SP_matrix(0))
  {
    Require(A);
    d_A = A;
    Ensure(d_A->number_rows() == d_A->number_columns());
  }

  /**
   *  \param atol   absolute tolerance (||r_n|| < atol)
   *  \param rtol   relative tolerance (||r_n|| < rtol * ||r_0||)
   *  \param maxit  maximum iterations (n < maxit)
   */
  void set_tolerances(const double atol, const double rtol, const int maxit)
  {
    d_absolute_tolerance = atol;
    d_relative_tolerance = rtol;
    d_maximum_iterations = maxit;
    Require(d_absolute_tolerance > 0.0);
    Require(d_relative_tolerance > 0.0);
    Require(d_maximum_iterations >= 0);
  }

  /**
   *  Print residual norms and other diagonostic information.
   *
   *  \param v  monitor via stdout
   */
  void set_monitor_output(const bool v)
  {
    d_monitor_output = v;
  }

  /**
   *  Turn on monitoring of diverging iterations.
   */
  void set_monitor_diverge(const bool v)
  {
    d_monitor_diverge = v;
  }

  /**
   *  \param b  right hand side
   *  \param x  unknown vector
   */
  int solve(const Vector<T> &b, Vector<T> &x)
  {
    Require(x.size() == b.size());
    Require(x.size() == d_A->number_rows());

    d_status = MAXIT;
    solve_impl(b, x);
    if (d_status ==  MAXIT)
    {
      printf("*** %s did not converge within the maximum number of iterations\n",
             d_name.c_str());
    }
    return d_status;
  }

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  std::string d_name;
  double d_absolute_tolerance;
  double d_relative_tolerance;
  int    d_maximum_iterations;
  std::vector<double> d_L1_residual;
  std::vector<double> d_L2_residual;
  std::vector<double> d_LI_residual;
  int    d_number_iterations;
  SP_matrix d_A;
  SP_matrix d_P;
  bool d_monitor_output;
  bool d_monitor_diverge;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  // print out iteration and residual for initial
  bool monitor_init(T r)
  {
    d_L2_residual[0] = r;
    if (d_monitor_output) printf("iteration: %5i    residual: %12.8e \n", 0, r);
    if (r < d_absolute_tolerance)
    {
      printf("*** %s converged in %5i iterations with a residual of %12.8e \n",
             d_name.c_str(), 0, r );
      d_status = SUCCESS;
      return true;
    }
    return false;
  }

  // print out iteration and residual
  bool monitor(int it, T r)
  {
    d_number_iterations = it;
    d_L2_residual[it] = r;
    if (d_monitor_output) printf("iteration: %5i    residual: %12.8e \n", it, r);
    Assert(it > 0);
    if (r < std::max(d_relative_tolerance * d_L2_residual[it - 1],
                                d_absolute_tolerance))
    {
      printf("*** %s converged in %5i iterations with a residual of %12.8e \n",
             d_name.c_str(), it, r );
      d_status = SUCCESS;
      return true;
    }
    else if (d_monitor_diverge and it >  1 and r - d_L2_residual[it - 1] > 0.0)
    {
      printf("*** %s diverged \n", d_name.c_str());
      d_status = DIVERGE;
      return true;
    }
    return false;
  }

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL LINEAR SOLVERS MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  virtual void solve_impl(const Vector<T> &b, Vector<T> &x) = 0;

private:

  int d_status;

};

} // end namespace callow

#endif /* LINEARSOLVER_HH_ */
