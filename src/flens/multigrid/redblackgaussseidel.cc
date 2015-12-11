#include <cassert>
#include <redblackgaussseidel.h>

namespace flens {

RedBlackGaussSeidel::RedBlackGaussSeidel(const DistributedVector &_x, double _omega)
    : x(_x), omega(_omega)
{    
}

void
mv(double alpha, const RedBlackGaussSeidel &A, const DistributedVector &x,
   double beta, DistributedVector &y)
{
    assert(alpha==1.);
    assert(beta==0.);
    assert(&x==&y);
    
    mv(A, y);
}

void
mv(const RedBlackGaussSeidel &A, DistributedVector &y)
{
    const DistributedVector::Grid &F = A.x.grid;
    DistributedVector::Grid       &U = y.grid;
    
    double c         = A.omega/4.;
    double hh        = 1./((y.N+1)*(y.N+1));


    // red nodes
    for (int i=y.firstRow; i<=y.m; ++i) {
        for (int j=y.firstCol; j<=y.n; ++j) {
            
            if ((i+j)%2==1) {
                continue;
            }
             
            U(i,j) = (1-A.omega)*U(i,j)
                   + c*(F(i,j)*hh + U(i-1,j) + U(i+1,j) + U(i,j-1) + U(i,j+1));        
        }
    }

    y.updateBoundary();

    // black nodes

    for (int i=y.firstRow; i<=y.m; ++i) {
        for (int j=y.firstCol; j<=y.n; ++j) {
            
            if ((i+j)%2==0) {
                continue;
            }
             
            U(i,j) = (1-A.omega)*U(i,j)
                   + c*(F(i,j)*hh + U(i-1,j) + U(i+1,j) + U(i,j-1) + U(i,j+1));        
        }
    }

    y.updateBoundary();    
}   
 
} // namespace flens
