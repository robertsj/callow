//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PCILU0.hh
 * \brief  PCILU0 
 * \author Jeremy Roberts
 * \date   Sep 18, 2012
 */
//---------------------------------------------------------------------------//

#ifndef PCILU0_HH_
#define PCILU0_HH_


#include "Preconditioner.hh"
#include "matrix/Matrix.hh"

namespace callow
{

/*!
 *  \class PCILU0
 *  \brief Implements the ILU(0) preconditioner
 *
 *  Following Saad, ILU(0) is defined
 *  \code
 *    for i = 2, n
 *      for k = 1, i - 1
 *        for (i, k) in nonzeros of lower A
 *          A(i,k) = A(i,k)/A(k,k)
 *          for j = k + 1 .. n
 *            for (i, j) in nonzeros of upper A
 *              A(i, j) = A(i, j) - A(i, k)*A(k, k)
 *            end
 *          end
 *        end
 *      end
 *    end
 *  \endcode
 */
template <class T>
class PCILU0: public Preconditioner<T>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef Preconditioner<T>                   Base;
  typedef typename Base::SP_preconditioner    SP_preconditioner;
  typedef typename Matrix<T>::SP_matrix       SP_matrix;
  typedef typename Vector<T>::SP_vector       SP_vector;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /// Construct a Jacobi preconditioner for the explicit matrix A
  PCILU0(SP_matrix A);

  /// Virtual destructor
  virtual ~PCILU0(){};

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL PRECONDITIONERS MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  /// Solve Px = y
  void apply(Vector<T> &x);

protected:

  /// ILU decomposition of A
  SP_matrix d_P;

};

} // end namespace callow

#include "PCILU0.i.hh"

#endif // PCILU0_HH_ 

//---------------------------------------------------------------------------//
//              end of file PCILU0.hh
//---------------------------------------------------------------------------//
