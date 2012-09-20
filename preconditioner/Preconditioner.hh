//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Preconditioner.hh
 * \brief  Preconditioner 
 * \author Jeremy Roberts
 * \date   Sep 18, 2012
 */
//---------------------------------------------------------------------------//

#ifndef PRECONDITIONER_HH_
#define PRECONDITIONER_HH_

#include "callow_config.hh"
#include "matrix/MatrixBase.hh"
#include "vector/Vector.hh"
#include <string>

namespace callow
{

/*!
 *  \class Preconditioner
 *  \brief Defines a preconditioner for linear solves
 *
 *  Consider the linear system
 *  \f[
 *      \mathbf{A}x = b \, .
 *  \f]
 *  When using iterative methods to solve this system,
 *  one can often make the system easier to solve.  Suppose
 *  we define an operator \f$ \mathbf{P} \f$ such that
 *  \f$ \mathbf{P}^{-1} \mathbf{A} \approx \mathbf{I} \f$.
 *  If the action of \f$ \mathbf{P}^{-1} \f$ is relatively
 *  easy to compute (as compared to inverting \mathbf{A}),
 *  \f$ \mathbf{P} \f$ is a good preconditioner.
 *
 *  We apply a precondition on the left to obtain the
 *  modified system
 *  \f[
 *      \mathbf{P^{-1} A}x = \mathbf{P}^{-1} b \, ,
 *  \f]
 *  or on the right to get
 *  \f[
 *      \mathbf{AP^{-1}} \overbrace{\mathbf{P}x}^{y} =  b \, ,
 *  \f]
 *  following which we solve
 *  \f[
 *      x = \mathbf{P}^{-1} y \, .
 *  \f]
 *
 *  Within callow, the Jacobi and ILU(0) preconditioners are
 *  available.
 */
template <class T>
class Preconditioner
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef SP<Preconditioner<T> >              SP_preconditioner;
  typedef typename MatrixBase<T>::SP_matrix   SP_matrix;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  Preconditioner(std::string name = "preconditioner"){};

  virtual ~Preconditioner(){};

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL PRECONDITIONERS MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  /// Solve Px = b
  virtual void apply(Vector<T> &b, Vector<T> &x) = 0;

protected:

  /// pc name
  std::string d_name;

};

} // end namespace callow

#endif // PRECONDITIONER_HH_ 

//---------------------------------------------------------------------------//
//              end of file Preconditioner.hh
//---------------------------------------------------------------------------//
