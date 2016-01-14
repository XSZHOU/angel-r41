#include <cassert>
#include <cmath>
#include <iostream>

#include <flens/flens.h>
#include <scifi/poisson1d.h>
#include <scifi/gaussseidel.h>

using namespace flens;
using namespace std;

//
// Assignment 1:  edit file gaussseidel.h
//

int
main()
{
    Poisson1D                   A(8);
    DenseVector<Array<double> > x(_(0,8)), f(_(0,8));
    GaussSeidel<Poisson1D>      S(A, f, 1.3);
    
    
    // fuer verschieden x testen
    
    f = 1, 1, 1, 1, 1, 1, 1;                            
    
    SnapShot snap("A1_gs");
    for (int k=1; k<=50; ++k) {
        x = S*x;
        snap(k) << x;
    }
}
