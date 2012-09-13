# Set source
SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(UTILS_SRC
  ${SRC_DIR}/CallowException.cc
  PARENT_SCOPE
)

# Install
INSTALL_FILES(/include/callow/utils .hh)
