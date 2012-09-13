//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Richardson.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Test of Richardson class.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST            \
        FUNC(test_Richardson)

#include "TestDriver.hh"
#include "solver/Richardson.hh"
#include "matrix/Matrix.hh"
#include <iostream>
using std::cout;
using std::endl;

// Setup
using namespace callow;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//---------------------------------------------------------------------------//
// TEST DEFINITIONS
//---------------------------------------------------------------------------//

int test_Richardson(int argc, char *argv[])
{
  typedef Matrix<double>        Mat;
  typedef Vector<double>        Vec;
  typedef Richardson<double>    Solver;

  int n = 5;

  // Create test matrix
  Mat* A(new Mat(n, n));
  A->preallocate(3);

  for (int row = 0; row < n; row++)
  {
    if (row == 0)
    {
      int c[]    = { 0, 1 };
      double v[] = {0.4, -0.22};
      A->add_row(row, c, v, 2);
    }
    else if (row == n - 1)
    {
      int c[]    = { n - 2, n - 1};
      double v[] = {    -0.2,    0.4};
      A->add_row(row, c, v, 2);
    }
    else
    {
      int c[]    = { row - 1, row, row + 1};
      double v[] = {      -0.2, 0.4,      -0.22};
      A->add_row(row, c, v);
    }
  }
  A->assemble();

  // Create two vectors
  Vector<double> X(n, 0.0);
  Vector<double> B(n, 1.0);

//  19.2181034482759
//  30.3965517241379
//  33.2500000000000
//  28.2758620689655
//  16.6379310344828

  // Create solver
  Solver solver(0.000001, 0.000001, 1000, 2.4);
  solver.set_operators(A);
  solver.set_monitor(true);
  int status = solver.solve(B, X);
  cout << " Status = " << status << endl;

  A->display();
  X.display();
  B.display();

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_Matrix.cc
//---------------------------------------------------------------------------//
