#include <cmath>

#include <transitions.h>
#include <distributedvector.h>
#include <poisson2d.h>
#include <redblackgaussseidel.h>

using namespace flens;

double
norm(const DistributedVector &x)
{
    double h = 1./(x.N+1);
    
    // compute (frobenius) norm
    // -> local norm
    double r = 0;
    for (int i=x.firstRow; i<=x.m; ++i) {
        for (int j=x.firstCol; j<=x.n; ++j) {
            r += x.grid(i,j)*x.grid(i,j);
        }
    }

    // -> global norm
    double R;
    x.mpiCart.comm.Allreduce(&r, &R, 1, MPI::DOUBLE, MPI::SUM);
    return sqrt(h*R);
}

double
errorNorm(const DistributedVector &e, const DistributedVector &u)
{
    double h = 1./(e.N+1);
    
    double r = 0;
    for (int i=e.firstRow; i<=e.m; ++i) {
        for (int j=e.firstCol; j<=e.n; ++j) {
            r += (e.grid(i,j)-u.grid(i,j))*(e.grid(i,j)-u.grid(i,j));
        }
    }
    
    double R;
    u.mpiCart.comm.Allreduce(&r, &R, 1, MPI::DOUBLE, MPI::SUM);
    return sqrt(h*R);
}

void
setRHS(DistributedVector &f)
{
    double h = 1./(f.N+1);
    
    for (int i=f.firstRow; i<=f.m; ++i) {
        for (int j=f.firstCol; j<=f.n; ++j) {
            double x = (i+f.m0)*h;
            double y = (j+f.n0)*h;
            
            f.grid(i,j) =  2*y*(1-y) + 2*x*(1-x);
        }
    }
}

void
setExactSolution(DistributedVector &e)
{
    double h = 1./(e.N+1);
    
    for (int i=e.firstRow; i<=e.m; ++i) {
        for (int j=e.firstCol; j<=e.n; ++j) {
            double x = (i+e.m0)*h;
            double y = (j+e.n0)*h;
            
            e.grid(i,j) =  x*(1-x)*y*(1-y);
        }
    }
}

//------------------------------------------------------------------------------

const long l = 9;
const int  n = (1<<l) - 1;  // n = 2^l -1

Restriction R;
Prolongation P;

DistributedVector  u[l+1], f[l+1], r[l+1];
Poisson2D A[l+1];

void
initMultigrid(const MpiCart &mpiCart)
{
    typedef DistributedVector Vec;

    int N = n;
    for (int i=l; i>=1; --i, N/=2) {
        u[i] = Vec(mpiCart, N+1);
        f[i] = Vec(mpiCart, N+1);
        r[i] = Vec(mpiCart, N+1);
        A[i] = Poisson2D(N+1);
    }
    for (int i=l; i>=3; --i) {
        f[i-1] = R*f[i];
    }
}

namespace flens {

void
residual(const DistributedVector &f,
         const Poisson2D &A,
         const DistributedVector &u,
         DistributedVector &r)
{
    int rhh = A.rh*A.rh;
 
    const DistributedVector::Grid &F = f.grid;
    const DistributedVector::Grid &U = u.grid;
    DistributedVector::Grid       &R = r.grid;
    
    for (int i=f.firstRow; i<=f.m; ++i) {
        for (int j=f.firstCol; j<=f.n; ++j) {
            R(i,j) = F(i,j) - rhh*(4*U(i,j) -U(i+1,j) -U(i-1,j) -U(i,j+1) -U(i,j-1));
        }
    }
}

} // namespace flens

void
printStat(int k,
          const DistributedVector &r,
          const DistributedVector &u,
          const DistributedVector &e)
{
    double rNorm = norm(r);
    double eNorm = errorNorm(e,u);

    if ((r.mpiCart.row==0) && (r.mpiCart.col==0)) {
        std::cout.width(3);
        std::cout << k << ") | ";
    
        std::cout.precision(12);
        std::cout.setf(std::ios::fixed);
        std::cout.width(17);
        std::cout << rNorm << " | ";
    
        std::cout.precision(12);
        std::cout.setf(std::ios::fixed);
        std::cout.width(16);
        std::cout << eNorm << " | " << std::endl;
    }
}

void
multigrid(int l, int v1, int v2)
{
    RedBlackGaussSeidel S(f[l], 1.);
    if (l==2) {
        for (int k=1; k<=1000; ++k) {
            u[l] = S * u[l];
        }
    } else {
        for (int v=1; v<=v1; ++v) {
            u[l] = S * u[l];
        }
        r[l] = f[l] - A[l]*u[l];
        
        f[l-1] = R * r[l];

        u[l-1].grid = 0;
        multigrid(l-1,v1,v2);
        
        u[l] += P * u[l-1];
        
        for (int v=1; v<=v2; ++v) {
            u[l] = S * u[l];
        }
    }
}

void
fullMultigrid()
{
    RedBlackGaussSeidel S(f[2], 1.);
    for (int k=1; k<=1000; ++k) {
        u[2] = S * u[2];
    }
    
    for (int i=2; i<l; ++i) {
        u[i+1] += P * u[i];
        multigrid(i+1, 1, 1);
    }
}

//------------------------------------------------------------------------------

using namespace std;

int
main(int argc, char *argv[])
{
    MPI::Init(argc, argv);
    
    MpiCart mpiCart(1,1);

    initMultigrid(mpiCart);

    DistributedVector exact(mpiCart, n+1);
    setRHS(f[l]);
    setExactSolution(exact);
    
    r[l] = f[l] - A[l]*u[l];
    printStat(0, r[l], u[l], exact);
    
    for (int k=1; k<=20; ++k) {
        multigrid(l, 1, 1);
        printStat(k, r[l], u[l], exact);
    }
    
    MPI::Finalize();
    
    return 0;	
}
