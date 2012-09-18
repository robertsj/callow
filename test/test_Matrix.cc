//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Matrix.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Test of Matrix class.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST         \
        FUNC(test_Matrix)

#include "TestDriver.hh"
#include "matrix/Matrix.hh"
#include "utils/Initialization.hh"
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
int test_Matrix(int argc, char *argv[])
{
  typedef Matrix<double> Mat;
  typedef Vector<double> Vec;

  // n * n
  {
    // Create test matrix
    int n = 5;
    Mat A(n, n);
    TEST(A.number_rows()    == n);
    TEST(A.number_columns() == n);
    A.preallocate(3);

    for (int row = n - 1; row >= 0; row--)
    {
      if (row == 0)
      {
        int c[]    = { 0, 1  };
        double v[] = {-2, 1.1};
        TEST(A.insert(row, c, v, 2));
      }
      else if (row == A.number_rows() - 1)
      {
        int c[] = {A.number_rows() - 2, A.number_rows() - 1};
        double v[] = {1, -2};
        TEST(A.insert(row, c, v, 2));
      }
      else
      {
        int c[] = {row - 1, row, row + 1};
        double v[] = {1, -2, 1.1};
        TEST(A.insert(row, c, v, 3));
      }
    }
    A.assemble();

    // Create two vectors
    Vector<double> X(n, 0.0);
    Vector<double> Y(n, 0.0);

    for (int i = 0; i < X.size(); i++)
    {
      double value = i * i;
      X[i] = value;
    }

    // Multiply
    A.multiply(X, Y);
    double ref[] =
    { 1.1, 2.4, 2.9, 3.6, -23.0 };
    for (int i = 0; i < Y.size(); i++)
    {
      TEST(soft_equiv(Y[i], ref[i]));
    }

    // Transpose
    A.multiply_transpose(Y, X);
    double ref2[] =
    { 0.2, -0.69, 0.44, -27.01, 49.96 };
    for (int i = 0; i < X.size(); i++)
    {
      TEST(soft_equiv(X[i], ref2[i]));
    }

    // Indexing
    double *a = A.value();
    for (int i = 0; i < A.number_rows(); i++)
    {
      std::cout << i << " " << A.diagonal(i) << std::endl;
//      for (int p = A.start(i); p < A.diagonal(i); p++)
//      {
//        std::cout << " " << A.column(p) << "(" << a[p] << ")";
//      }
//      std::cout << endl;
    }

    A.display();
    X.display();
    Y.display();

  } // end n * n

  // m * n
  {
    int m = 3;
    int n = 2;

    // Create test matrix
    Mat A(m, n);
    TEST(A.number_rows()    == m);
    TEST(A.number_columns() == n);
    A.preallocate(2);
    /*!
     *  1 2
     *  3 4
     *  5 6
     */
    A.insert(0, 0, 1.0);
    A.insert(0, 1, 2.0);
    A.insert(1, 0, 3.0);
    A.insert(1, 1, 4.0);
    A.insert(2, 0, 5.0);
    A.insert(2, 1, 6.0);
    A.assemble();

    // Create two vectors
    Vector<double> X(n, 1.0);
    Vector<double> Y(m, 0.0);

    // Multiply
    A.multiply(X, Y);
    double ref[] =
    { 3.0, 7.0, 11.0 };
    for (int i = 0; i < Y.size(); i++)
    {
      TEST(soft_equiv(Y[i], ref[i]));
    }

    // Transpose
    A.multiply_transpose(Y, X);
    return 0;
    double ref2[] =
    { 79.0, 100.0 };
    for (int i = 0; i < X.size(); i++)
    {
      TEST(soft_equiv(X[i], ref2[i]));
    }


    A.display();
    X.display();
    Y.display();

  } // end m * n

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_Matrix.cc
//---------------------------------------------------------------------------//
