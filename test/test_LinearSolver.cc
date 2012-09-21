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
        FUNC(test_GMRES)       \
        FUNC(test_PetscSolver)

#include "TestDriver.hh"
#include "utils/Initialization.hh"
// solvers
#include "solver/Richardson.hh"
#include "solver/Jacobi.hh"
#include "solver/GaussSeidel.hh"
#include "solver/GMRES.hh"
#ifdef CALLOW_ENABLE_PETSC
#include "solver/PetscSolver.hh"
#endif
// pc
#include "preconditioner/PCJacobi.hh"
#include "preconditioner/PCILU0.hh"
//
#include "test/matrix_fixture.hh"
#include <iostream>

using std::cout;
using std::endl;

using namespace callow;

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

//---------------------------------------------------------------------------//
// TEST DEFINITIONS
//---------------------------------------------------------------------------//

// system size
int n = 20;

// convergence
double abstol = 1e-13;
double reltol = 1e-13;
int    maxit  = 1000;

// reference
double X_ref[] = {5.459135698786395e+00, 8.236037378063020e+00,
                  9.648531187020394e+00, 1.036694341855932e+01,
                  1.073221653274062e+01, 1.091771229837346e+01,
                  1.101148972685055e+01, 1.105810668452673e+01,
                  1.107978760425382e+01, 1.108701173915058e+01,
                  1.108356356535519e+01, 1.106847349927400e+01,
                  1.103582874555248e+01, 1.097247463295948e+01,
                  1.085272174937750e+01, 1.062793308617560e+01,
                  1.020677234863001e+01, 9.418093128951856e+00,
                  7.941390927380783e+00, 5.176556370952343e+00};

int test_Richardson(int argc, char *argv[])
{
  Vector<double> X(n, 0.0);
  Vector<double> B(n, 1.0);
  Richardson<double> solver(abstol, reltol, maxit, 1.0);
  solver.set_operators(test_matrix_1<double>(n));
  solver.set_monitor_output(true);
  solver.set_monitor_diverge(false);
  int status = solver.solve(B, X);
  TEST(status == 0);
  for (int i = 0; i < 20; ++i)
  {
    TEST(soft_equiv(X[i],  X_ref[i], 1e-9));
  }
  return 0;
}

int test_Jacobi(int argc, char *argv[])
{
  Vector<double> X(n, 0.0);
  Vector<double> B(n, 1.0);
  Jacobi<double> solver(abstol, reltol, maxit);
  solver.set_operators(test_matrix_1<double>(n));
  solver.set_monitor_output(true);
  solver.set_monitor_diverge(false);
  int status = solver.solve(B, X);
  TEST(status == 0);
  for (int i = 0; i < 20; ++i)
  {
    TEST(soft_equiv(X[i],  X_ref[i], 1e-9));
  }
  return 0;
}

int test_GaussSeidel(int argc, char *argv[])
{
  Vector<double> X(n, 0.0);
  Vector<double> B(n, 1.0);
  GaussSeidel<double> solver(abstol, reltol, maxit);
  solver.set_operators(test_matrix_1<double>(n));
  solver.set_monitor_output(true);
  solver.set_monitor_diverge(false);
  int status = solver.solve(B, X);
  TEST(status == 0);
  for (int i = 0; i < 20; ++i)
  {
    TEST(soft_equiv(X[i],  X_ref[i], 1e-9));
  }
  return 0;
}

int test_GMRES(int argc, char *argv[])
{
  Vector<double> X(n, 0.0);
  Vector<double> B(n, 1.0);
  GMRES<double> solver(abstol, reltol, maxit, 10);
  solver.set_operators(test_matrix_1<double>(n));
  solver.set_monitor_output(true);
  solver.set_monitor_diverge(false);
  int status = solver.solve(B, X);
  TEST(status == 0);
  for (int i = 0; i < 20; ++i)
  {
    TEST(soft_equiv(X[i],  X_ref[i], 1e-9));
  }
  return 0;
}

int test_PetscSolver(int argc, char *argv[])
{
#ifdef CALLOW_ENABLE_PETSC
  Vector<double> X(n, 0.0);
  Vector<double> B(n, 1.0);
  PetscSolver solver(abstol, reltol, maxit);
  solver.set_operators(test_matrix_1<double>(n));
  solver.set_monitor_output(true);
  solver.set_monitor_diverge(false);
  int status = solver.solve(B, X);
  X.display();
  //TEST(status == 0);
  for (int i = 0; i < 20; ++i)
  {
    TEST(soft_equiv(X[i],  X_ref[i], 1e-9));
  }
#endif
  return 0;
}



//---------------------------------------------------------------------------//
//              end of test_Matrix.cc
//---------------------------------------------------------------------------//
