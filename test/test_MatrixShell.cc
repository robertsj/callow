//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_MatrixShell.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Test of Matrix class.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST         \
        FUNC(test_MatrixShell)

#include "TestDriver.hh"
#include "matrix/Matrix.hh"
#include "utils/Initialization.hh"
#include "test/matrixshell_fixture.hh"
#include <iostream>
using std::cout;
using std::endl;

// Setup
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

// Test of basic public interface
int test_MatrixShell(int argc, char *argv[])
{
  typedef TestMatrixShell<double> Mat_T;
  typedef Vector<double>          Vec_T;

  // n * n
  {
    // Create test matrix
    int n = 5;
    Mat_T A(n);
    A.display();
    TEST(A.number_rows()    == n);
    TEST(A.number_columns() == n);



    // Create two vectors
    Vec_T X(n, 0.0);
    Vec_T Y(n, 0.0);

    for (int i = 0; i < X.size(); i++)
    {
      double value = i * i;
      X[i] = value;
    }

    // Multiply
    A.multiply(X, Y);
    double ref[] =
    {-0.210000000000000, -1.340000000000000, -4.090000000000000, -8.660000000000000, -9.800000000000001};

    for (int i = 0; i < Y.size(); i++)
    {
      TEST(soft_equiv(Y[i], ref[i]));
    }

    // Transpose
    A.multiply_transpose(Y, X);
    double ref2[] =
    {0.373000000000000, 1.532100000000000, 4.058400000000001, 7.148900000000000, 6.718600000000000};

    for (int i = 0; i < X.size(); i++)
    {
      TEST(soft_equiv(X[i], ref2[i]));
    }

  } // end n * n

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_MatrixShell.cc
//---------------------------------------------------------------------------//
