#include <cassert>
#include <cmath>
#include <flens/flens.h>
#include <iostream>

#include <scifi/poisson1d.h>
#include <scifi/gaussseidel.h>
#include <scifi/restriction.h>
#include <scifi/prolongation.h>

using namespace flens;
using namespace std;

//
// Assignment 6:  edit this file
//


//-- Problem setup -------------------------------------------------------------

//----> right-hand side for -u'' = f
void
setRHS(DenseVector<Array<double> > &f)
{
    //-- BEGIN
    // your code here

    //    -> use  f(x) = pi^2 * sin(pi*x) 
    //-- END
}

//----> exact solution
void
exactSolution(DenseVector<Array<double> > &u)
{
    //-- BEGIN
    // your code here

    //    -> exact solution:  u(x) = sin(pi*x) 
    //-- END
}

//-- Multigrid -----------------------------------------------------------------

//----> data structures
const long long l = 4;
const int n = (1<<l) - 1;  // n = 2^l -1
DenseVector<Array<double> > u[l+1], f[l+1], r[l+1];
Poisson1D                   A[l+1];

//----> init data structures
void
initMultigrid()
{
    typedef DenseVector<Array<double> > Vec;

    //-- BEGIN
    // your code here

    // NOTE:
    //    -> initialize u[], f[], r[], A[]
    //       make sure they are of correct size!
    
    // example:  u[3] = Vec(_(a,b));
    //    -> initialize u[3] with vector of range [a,...,b]
    //-- END
}

//----> multigrid method
Restriction  R;
Prolongation P;

template <typename S>
void
multigrid(int l, int v1, int v2)
{
    //-- BEGIN
    // your code here
    
    // implement the multigrid method with 15 lines of code
    
    //-- END
    
    // NOTES:
    //    -> for l==2 solve exactly
    //    -> v1 number of pre-smoothing steps
    //    -> v2 number of post-smoothing steps
}

//-- Statistics ----------------------------------------------------------------

DenseVector<Array<double> > e(_(0,n+1));

void
printStat(int iteration)
{
    double rNorm = 0;
    double error = 0;
    
    //-- BEGIN
    // your code here
    
    //    -> compute l2-norm of residual
    //    -> compute l2-norm of error
    
    //-- END

    cout.width(3);
    cout << iteration << ") | ";
    
    cout.precision(12);
    cout.setf(std::ios::fixed);
    cout.width(17);
    cout << rNorm << " | ";
    
    cout.precision(12);
    cout.setf(std::ios::fixed);
    cout.width(16);
    cout << error << " | " << endl;
}


int
main()
{
    initMultigrid();
    setRHS(f[l]);
    exactSolution(e);
    
    //-- BEGIN
    // your code here
    
    //    -> init r[l] with residual for initial guess u[l]
    //    -> start multigrid with v1 = 0 and v2 = 1  (why?)
    //    -> loop multigrid with v1 = 1 and v2 = 1
    //    -> after each iteration call printStat(k)
    //       where k is the iteration counter
    
    //-- END
}
