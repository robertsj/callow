# Set source
SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(SOLVER_SRC
  ${SRC_DIR}/Richardson.cc
  ${SRC_DIR}/Jacobi.cc
  ${SRC_DIR}/GaussSeidel.cc
  ${SRC_DIR}/GMRES.cc
  ${SRC_DIR}/PetscSolver.cc
  PARENT_SCOPE
)

# Install
INSTALL_FILES(/include/callow/solver .hh)


