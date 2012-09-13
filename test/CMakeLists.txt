
# Vector
ADD_EXECUTABLE(test_Vector              test_Vector.cc)
TARGET_LINK_LIBRARIES(test_Vector       callow )
ADD_TEST(test_Vector                    test_Vector 0)

# Matrix
ADD_EXECUTABLE(test_Matrix              test_Matrix.cc)
TARGET_LINK_LIBRARIES(test_Matrix       callow )
ADD_TEST(test_Matrix                    test_Matrix 0)

## Solvers

# RICHARDSON
ADD_EXECUTABLE(test_Richardson          test_Richardson.cc)
TARGET_LINK_LIBRARIES(test_Richardson   callow )
ADD_TEST(test_Richardson                test_Richardson 0)

# JACOBI

# GAUSS-SEIDEL