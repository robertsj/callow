# Set source
SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(MATRIX_SRC
  ${SRC_DIR}/Matrix.cc
  PARENT_SCOPE
)

# Install
INSTALL_FILES(/include/callow/matrix .hh)

