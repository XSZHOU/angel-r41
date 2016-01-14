#include <cassert>
#include <cmath>
#include <iostream>

#include <flens/flens.h>
#include <scifi/poisson1d.h>
#include <scifi/gaussseidel.h>

using namespace flens;
using namespace std;

//
// Assignment 2:  edit file poisson1d.h
//


int
main()
{
    int N = 32;
    //int    N = 63;
    double h = 1./(N+1);
    
    
    Poisson1D                   A(N+1);
    DenseVector<Array<double> > u(_(0,N+1)), f(_(0,N+1)), r(_(0,N+1));
    GaussSeidel<Poisson1D>      S(A, f, 1.);
    
    // initial guess for u consisting of two frequencies
    for (int k=1; k<=N; ++k) {
        double x = k*h;
        u(k) = sin(27*M_PI*x) + sin(3*M_PI*x);
    }

    SnapShot initial("A2_gs_initial");
    initial(0) << u << endl;
   
    SnapShot snap("A2_gs_res");
    for (int k=1; k<=30; ++k) {
        r = f - A*u;
        snap(k) << r << endl;
     
        u = S*u;
    }
}
