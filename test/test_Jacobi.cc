//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Jacobi.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Test of Jacobi class.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST            \
        FUNC(test_Jacobi)

#include "TestDriver.hh"
#include "solver/Jacobi.hh"
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

int test_Jacobi(int argc, char *argv[])
{
  typedef Matrix<double>        Mat;
  typedef Vector<double>        Vec;
  typedef Jacobi<double>        Solver;


  typename Matrix<double>::SP_matrix A;
  int n = 5;
  A = test_matrix_1<double>(n);

  // Create two vectors
  Vector<double> X(n, 0.0);
  Vector<double> B(n, 1.0);

//  19.2181034482759
//  30.3965517241379
//  33.2500000000000
//  28.2758620689655
//  16.6379310344828

  // Create solver
  Solver solver(0.000001, 0.000001, 10);

  solver.set_operators(A);
  solver.set_monitor(true);
  int status = solver.solve(B, X);
  cout << " Status = " << status << endl;

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_Matrix.cc
//---------------------------------------------------------------------------//
