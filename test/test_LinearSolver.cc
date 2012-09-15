//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_LinearSolver.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Test of LinearSolver class and its subclasses
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST              \
        FUNC(test_Richardson)  \
        FUNC(test_Jacobi)      \
        FUNC(test_GaussSeidel) \
        FUNC(test_GMRES)

#include "TestDriver.hh"
// solvers
#include "solver/Richardson.hh"
#include "solver/Jacobi.hh"
#include "solver/GaussSeidel.hh"
#include "solver/GMRES.hh"
//
#include "matrix/Matrix.hh"
#include "test/matrix_fixture.hh"
#include <iostream>
using std::cout;
using std::endl;

using namespace callow;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//---------------------------------------------------------------------------//
// TEST DEFINITIONS
//---------------------------------------------------------------------------//

// system size
int n = 100000;

// create two vectors
Vector<double> X(n, 0.0);
Vector<double> B(n, 1.0);

// convergence
double abstol = 1e-10;
double reltol = 1e-10;
int    maxit  = 1000;

int test_Richardson(int argc, char *argv[])
{
  Richardson<double> solver(abstol, reltol, maxit, 1);
  solver.set_operators(test_matrix_1<double>(n));
  solver.set_monitor_output(true);
  solver.set_monitor_diverge(false);
  int status = solver.solve(B, X);
//  TEST(status == 0);
//  TEST(soft_equiv(X[0],  5.0, 1e-9));
//  TEST(soft_equiv(X[1],  7.5, 1e-9));
//  TEST(soft_equiv(X[24], 9.9999995530, 1e-9));
  return 0;
}

int test_Jacobi(int argc, char *argv[])
{
  Jacobi<double> solver(abstol, reltol, maxit);
  solver.set_operators(test_matrix_1<double>(n));
  solver.set_monitor_output(true);
  solver.set_monitor_diverge(false);
  int status = solver.solve(B, X);
//  TEST(status == 0);
//  TEST(soft_equiv(X[0],  5.0, 1e-9));
//  TEST(soft_equiv(X[1],  7.5, 1e-9));
//  TEST(soft_equiv(X[24], 9.9999995530, 1e-9));
  return 0;
}

int test_GaussSeidel(int argc, char *argv[])
{
  GaussSeidel<double> solver(abstol, reltol, maxit);
  solver.set_operators(test_matrix_1<double>(n));
  solver.set_monitor_output(true);
  solver.set_monitor_diverge(false);

  int status = solver.solve(B, X);
  //X.display();
//  TEST(status == 0);
//  TEST(soft_equiv(X[0],  5.0, 1e-9));
//  TEST(soft_equiv(X[1],  7.5, 1e-9));
//  TEST(soft_equiv(X[24], 9.9999995530, 1e-9));
  return 0;
}

int test_GMRES(int argc, char *argv[])
{
  GMRES<double> solver(abstol, reltol, 200, 30);
  solver.set_operators(test_matrix_1<double>(n));
  solver.set_monitor_output(true);
  solver.set_monitor_diverge(true);
  int status = solver.solve(B, X);
//  GMRES<double> solver2(abstol, reltol, 10, 10 );
//  solver2.set_operators(test_matrix_1<double>(n));
//  solver2.set_monitor_output(true);
//  solver2.set_monitor_diverge(false);
//  status = solver2.solve(B, X);

 // X.display();
//  TEST(status == 0);
//  TEST(soft_equiv(X[0],  5.0, 1e-9));
//  TEST(soft_equiv(X[1],  7.5, 1e-9));
//  TEST(soft_equiv(X[24], 9.999999552965, 1e-6));
  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_Matrix.cc
//---------------------------------------------------------------------------//
