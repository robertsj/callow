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
#include "utils/Initialization.hh"
// solvers
#include "solver/Richardson.hh"
#include "solver/Jacobi.hh"
#include "solver/GaussSeidel.hh"
#include "solver/GMRES.hh"
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
int n = 50;

// convergence
double abstol = 1e-10;
double reltol = 1e-10;
int    maxit  = 1000;

int test_Richardson(int argc, char *argv[])
{
  Vector<double> X(n, 0.0);
  Vector<double> B(n, 1.0);
  Richardson<double> solver(abstol, reltol, maxit, 1);
  solver.set_operators(test_matrix_1<double>(n));
  solver.set_monitor_output(true);
  solver.set_monitor_diverge(false);
  int status = solver.solve(B, X);
  TEST(status == 0);
  TEST(soft_equiv(X[0],  5.0, 1e-9));
  TEST(soft_equiv(X[1],  7.5, 1e-9));
  TEST(soft_equiv(X[24], 9.9999995530, 1e-9));
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
  TEST(soft_equiv(X[0],  5.0, 1e-9));
  TEST(soft_equiv(X[1],  7.5, 1e-9));
  TEST(soft_equiv(X[24], 9.9999995530, 1e-9));
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
  TEST(soft_equiv(X[0],  5.0, 1e-9));
  TEST(soft_equiv(X[1],  7.5, 1e-9));
  TEST(soft_equiv(X[24], 9.9999995530, 1e-9));
  return 0;
}

int test_GMRES(int argc, char *argv[])
{
  typename Matrix<double>::SP_matrix A = test_matrix_2<double>(4);
  typename PCJacobi<double>::SP_preconditioner PJ(new PCJacobi<double>(A));
  typename PCILU0<double>::SP_preconditioner   PI(new PCILU0<double>(A));

  GMRES<double> solver(0, 1e-10, 10000, 30);
  solver.set_operators(A);
  //solver.set_left_pc(PI);
  solver.set_monitor_output(true);
  //solver.set_monitor_diverge(true);

  Vector<double> X(A->number_rows(), 1.0);
  Vector<double> B(A->number_rows(), 1.0);
  int status = solver.solve(B, X);
  A->print_matlab("A.out");
  B.print_matlab("B.out");
  X.print_matlab("X.out");


  typename Matrix<double>::SP_matrix Q;
  Q = new Matrix<double>(3,3);
  Q->preallocate(3);
  Q->insert(0, 0, 1);
  Q->insert(0, 1, 2);
  Q->insert(1, 0, 3);
  Q->insert(1, 1, 4);
  Q->insert(1, 2, 5);
  Q->insert(2, 1, 6);
  Q->insert(2, 2, 7);
  Q->assemble();
  PCILU0<double> PQ(Q);
  Q->print_matlab("q.out");


//  GMRES<double> solver2(abstol, reltol, 10, 10 );
//  solver2.set_operators(test_matrix_1<double>(n));
//  solver2.set_monitor_output(true);
//  solver2.set_monitor_diverge(false);
//  status = solver2.solve(B, X);
//  TEST(status == 0);
  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_Matrix.cc
//---------------------------------------------------------------------------//
