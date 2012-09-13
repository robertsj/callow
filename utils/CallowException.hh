//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CallowException.hh
 * \author robertsj
 * \date   Sep 13, 2012
 * \brief  CallowException class definition
 */
//---------------------------------------------------------------------------//

#ifndef callow_CALLOWEXCEPTION_HH_
#define callow_CALLOWEXCEPTION_HH_

#include <iostream>
#include <exception>
#include <string>

namespace callow
{

/*!
 *  A generic mechanism to manually manage exceptions
 */
class CallowException: public std::exception
{

public:

  /// Constructs a new CallowException with the default message.
  CallowException();

  /*!
   *  \brief Constructs a new GenException with a provided message
   *  \param line   line of code erring
   *  \param file   file in which error occurs
   *  \param msg    the message
   */
  CallowException(int line, std::string file, std::string msg);

  /*!
   *  \brief Returns the error message associated with this CallowException
   *  \return the message
   */
  virtual const char* what() const throw ();

  /*!
   *  \brief Destroys this GenException.
   */
  virtual ~CallowException() throw ();


protected:

  /// The message associated with this exception.
  std::string d_message;

  /// A string to prepend to all message of this class.
  static std::string d_prepend;

};

} // end namespace callow

/// Easy macro for throwing exceptions.
#define THROW(m) throw callow::CallowException(__LINE__,__FILE__,m);

#endif /* callow_CALLOWEXCEPTION_HH_ */
