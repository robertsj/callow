# Set source
SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(SOLVER_SRC
  ${SRC_DIR}/Richardson.cc
  ${SRC_DIR}/Jacobi.cc
  PARENT_SCOPE
)

# Install
INSTALL_FILES(/include/callow/solver .hh)


