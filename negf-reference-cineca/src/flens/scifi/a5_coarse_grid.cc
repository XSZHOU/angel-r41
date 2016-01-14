#include <cassert>
#include <cmath>
#include <iostream>

#include <flens/flens.h>
#include <scifi/poisson1d.h>
#include <scifi/gaussseidel.h>
#include <scifi/restriction.h>
#include <scifi/prolongation.h>

using namespace flens;
using namespace std;

//
// Assignment 5:  this file
//

int
main()
{
    //-- fine grid -------------------------------------------------------------
    int    N = 63;
    double h = 1./(N+1);
    
    Poisson1D                   A(N+1);
    DenseVector<Array<double> > u(_(0,N+1)), f(_(0,N+1)), r(_(0,N+1));
    GaussSeidel<Poisson1D>      S(A, f, 1.);

    //-- coarse grid -------------------------------------------------------------
    Restriction  R;
    Prolongation P;
    // int n = ?;
    
    //      -> define variables A_c, u_c, f_c, S_c of appropriate type and size
    
    //-- initial guess with two frequencies ------------------------------------
    
    //      -> initialize u with sin(27*pi*x) + sin(3*pi*x)

    //-- smooth on fine grid ---------------------------------------------------
    
    //      -> apply v1 smoothing steps (choose different v1)
    //      -> compute residual 

    //-- solve on coarse grid --------------------------------------------------

    //      -> restrict residual to coarse grid
    //      -> this restriction becomes new right-hand side f_c
    //      -> solve system of lin. eq.   A_c * u_c = f_c
    
    //-- apply correction ------------------------------------------------------

    //      -> update u with prolongation of u_c

    //-- smooth again on fine grid ---------------------------------------------

    //      -> perform another v2 smoothing steps (choose differen v2)
    //      -> compute residual

}

