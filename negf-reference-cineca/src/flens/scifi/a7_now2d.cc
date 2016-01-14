#include <cassert>
#include <cmath>
#include <flens/flens.h>
#include <iostream>

#include <examples/poisson2d.h>

using namespace flens;
using namespace std;

/*******************************************************************************

   Edit poisson2d.h

*******************************************************************************/

void
setRHS(GridVector2D &f)
{
    GridVector2D::Grid &F = f.grid;

    int n = F.lastRow()-1;
    double h = 1./(n+1);
    
    for (int i=1; i<=n; ++i) {
        for (int j=1; j<=n; ++j) {
            double x = i*h;
            double y = j*h;
            
            F(i,j) = 2*y*(1-y) + 2*x*(1-x);
        }
    }
}

int
main()
{
    int N = 15;
    Poisson2D A(N+1);
    Restriction R;
    Prolongation P;
    
    GridVector2D u(N+1), f(N+1), e(N+1), r(N+1);
    GaussSeidel<Poisson2D, GridVector2D> S(A, f, 1.7);


    int n = 7;
    Poisson2D A_c(n+1);
 
    GridVector2D u_c(n+1), f_c(n+1);
    GaussSeidel<Poisson2D, GridVector2D> S_c(A_c, f_c, 1.7);
 

    // init f
    setRHS(f);

    // smooth 5 times
    for (int k=1; k<=5; ++k) {
        u = S*u;
    }
    // compute residual
    r = f - A*u;
    cout << "r = " << norm(r) << "  ";

    // restrict
    f_c = R*r;
    
    // smooth on coarse grid
    for (int k=1; k<=10; ++k) {
        u_c = S_c*u_c;
    }
    
    // update u with prolonagtion of u_c
    u += P*u_c;
    r = f - A*u;
    cout << "r = " << norm(r) << "  ";

    // smooth on fine grid
    for (int k=1; k<=5; ++k) {
        u = S*u;
    }
    r = f - A*u;
    cout << "r = " << norm(r) << "  ";
}
