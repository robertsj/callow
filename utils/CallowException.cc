//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CallowException.cc
 * \author robertsj
 * \date   Sep 13, 2012
 * \brief  CallowException member definitions
 */
//---------------------------------------------------------------------------//

#include "CallowException.hh"
#include <sstream>

namespace callow
{

std::string itoa(int i)
{
  std::stringstream out;
  out << i;
  return out.str();
}

std::string dtoa(double d)
{
  std::stringstream out;
 out << d;
 return out.str();
}

std::string CallowException::d_prepend = "callow exception";

CallowException::CallowException()
{
    d_message = d_prepend;
}

CallowException::CallowException(int line, std::string file, std::string msg)
{
    d_message = d_prepend + "\n"
                          + "           on line: " + itoa(line) + "\n"
                          + "           in file: " + file + "\n"
                          + "           message: " + msg;
}

const char* CallowException::what() const throw()
{
    return d_message.c_str();
}

CallowException::~CallowException() throw()
{
  /* ... */
}

} // end namespace callow

