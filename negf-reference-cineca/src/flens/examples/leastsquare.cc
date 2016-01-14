// This example shows the usage of various least square solvers
// (at the moment the LAPACK functions GELS and GELSS).
// GELSS can also be used to calculate the pseudo-inverse of a matrix.
//
// For the differences of these methods see:
// http://en.wikipedia.org/wiki/Linear_least_squares
// http://www.netlib.org/lapack/lug/node27.html
//
// 2007, Georg Holzmann

#include "flens/flens.h"
#include <iostream>

using namespace std;
using namespace flens;

typedef GeMatrix<FullStorage<float, ColMajor> > DEMatrix;

int main(int argc, char *argv[])
{
  cout << "\n## OVERDETERMINED SYSTEM ##\n\n";

  DEMatrix A(3,2);
  DEMatrix B(3,2);

  A = 1, 2,
      3, 4,
      5, 6;
  B = 1, 4,
      2, 5,
      3, 6;

  cout << "A: " << A << "B: " << B << endl;

  // calc with gels
  flens::ls(flens::NoTrans, A, B);
  cout << "GELS result: " << B(_(1,2),_) << "Residual: " << B(3,_) << endl;

  // ATTENTION: we have to reinizialize the values, because
  //            ls() changes the content of A and B !
  A = 1, 2,
      3, 4,
      5, 6;
  B = 1, 4,
      2, 5,
      3, 6;

  // calc with gelss
  flens::lss(A, B);
  cout << "GELSS result: " << B(_(1,2),_) << "Residual: " << B(3,_) << endl;


  cout << "\n## UNDERDETERMINED SYSTEM ##\n\n";

  A.resize(2,3);
  B.resize(3,2);

  A = 1, 3, 5,
      2, 4, 6;
  B =  1,   4,
       2,   5,
      -99, -99; // these two values are only for the result storage

  cout << "A: " << A << "B: " << B(_(1,2),_) << endl;

  // calc with gels
  flens::ls(flens::NoTrans, A, B);
  cout << "GELS result: " << B(_(1,2),_) << "Residual: " << B(3,_) << endl;

  A = 1, 3, 5,
      2, 4, 6;
  B =  1,   4,
       2,   5,
      -99, -99;

  // calc with gelss
  flens::lss(A, B);
  cout << "GELSS result: " << B(_(1,2),_) << "Residual: " << B(3,_) << endl;

  return 0;
}
