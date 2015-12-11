#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iostream>
#include <mpi.h>

#include <flens/flens.h>

//------------------------------------------------------------------------------

using namespace flens;
using namespace std;


//------------------------------------------------------------------------------

enum CartDir { North=0, East=1, South=2, West=3 };

struct MpiCart
{
    MpiCart(int _numRows, int _numCols)
        : numRows(_numRows), numCols(_numCols)
    {
        int  dims[2];
        bool periods[2];

        // create mpi_cart:  2-dimensional, 2 x 2, non-periodic
        dims[0]    = numRows;
        dims[1]    = numCols;
        periods[0] = false;
        periods[1] = false;
        comm = MPI::COMM_WORLD.Create_cart(2, dims, periods, true);
        
        // find your position in the cart
        int rank = MPI::COMM_WORLD.Get_rank();
        int coords[2];
        comm.Get_coords(rank, 2, coords);
        row = coords[0];
        col = coords[1];
    
        // find your neighbors
        enum CartDir { North=0, East=1, South=2, West=3 };
        comm.Shift(0, 1, neighbors[North], neighbors[South]);
        comm.Shift(1, 1, neighbors[West],  neighbors[East]);
    }
    
    MpiCart(const MpiCart &rhs)
        : comm(rhs.comm),
          numRows(rhs.numRows), numCols(rhs.numCols),
          row(rhs.row), col(rhs.col)
    {
        neighbors[North] = rhs.neighbors[North];
        neighbors[East]  = rhs.neighbors[East];
        neighbors[South] = rhs.neighbors[South];
        neighbors[West]  = rhs.neighbors[West];
    }
    
    MPI::Cartcomm comm;
    int numRows, numCols;
    int row, col;
    int neighbors[4];
};

struct DistributedVector
{
    typedef GeMatrix<FullStorage<double, RowMajor> > Grid;

    DistributedVector(const MpiCart &_mpiCart, int rh)
        : mpiCart(_mpiCart), N(rh-1)
    {
        // global matrix has size (N+1) x (N+1)
        
        // compute m and n, such that local matrix has size (-1,m+1) x (-1,n+1)
        m = (N+1) / mpiCart.numRows - 1;
        n = (N+1) / mpiCart.numCols - 1;
        grid.resize(_(-1,m+1),_(-1,n+1));
        
        // m0, n0: row and col offset
        m0 = (m+1)*mpiCart.row;
        n0 = (n+1)*mpiCart.col;
        
        // firstRow, firstCol:
        // -> non-zero value of local grid is stored in (firstRow, firstCol)
        firstRow = (mpiCart.row==0) ? 1 : 0;
        firstCol = (mpiCart.col==0) ? 1 : 0;
                     
        // types to send rows and cols
        MpiRow = MPI::DOUBLE.Create_contiguous(n+1); 
        MpiRow.Commit();
            
        MpiCol = MPI::DOUBLE.Create_vector(m+1, 1, grid.leadingDimension()); 
        MpiCol.Commit();
    }
    
    void
    updateBoundary()
    {
        MPI::Request req[8];
        MPI::Status  stat[8];
    
        req[0] = mpiCart.comm.Isend(&grid(0,n), 1, MpiCol, mpiCart.neighbors[East],  0);  
        req[1] = mpiCart.comm.Isend(&grid(m,0), 1, MpiRow, mpiCart.neighbors[South], 0);    
        req[2] = mpiCart.comm.Isend(&grid(0,0), 1, MpiCol, mpiCart.neighbors[West],  0);  
        req[3] = mpiCart.comm.Isend(&grid(0,0), 1, MpiRow, mpiCart.neighbors[North], 0);    

        req[4] = mpiCart.comm.Irecv(&grid(  0,n+1), 1, MpiCol, mpiCart.neighbors[East],  0);
        req[5] = mpiCart.comm.Irecv(&grid(m+1,  0), 1, MpiRow, mpiCart.neighbors[South], 0);
        req[6] = mpiCart.comm.Irecv(&grid(  0, -1), 1, MpiCol, mpiCart.neighbors[West],  0);
        req[7] = mpiCart.comm.Irecv(&grid( -1,  0), 1, MpiRow, mpiCart.neighbors[North], 0);
        
        MPI::Request::Waitall(8, req, stat);
    }
    
    MpiCart         mpiCart;
    int             N;
    int             m, n;
    int             m0, n0;
    int             firstRow, firstCol;
    Grid            grid;
    MPI::Datatype   MpiRow, MpiCol;
};

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

void
gaussSeidelRedBlack(double omega, const DistributedVector &f, DistributedVector &u)
{
    const DistributedVector::Grid &F = f.grid;
    DistributedVector::Grid       &U = u.grid;
    
    double c         = omega/4.;
    double hh        = 1./((u.N+1)*(u.N+1));


    // red nodes
    for (int i=u.firstRow; i<=u.m; ++i) {
        for (int j=u.firstCol; j<=u.n; ++j) {
            
            if ((i+j)%2==1) {
                continue;
            }
             
            U(i,j) = (1-omega)*U(i,j)
                   + c*(F(i,j)*hh + U(i-1,j) + U(i+1,j) + U(i,j-1) + U(i,j+1));        
        }
    }

    u.updateBoundary();

    // black nodes

    for (int i=u.firstRow; i<=u.m; ++i) {
        for (int j=u.firstCol; j<=u.n; ++j) {
            
            if ((i+j)%2==0) {
                continue;
            }
             
            U(i,j) = (1-omega)*U(i,j)
                   + c*(F(i,j)*hh + U(i-1,j) + U(i+1,j) + U(i,j-1) + U(i,j+1));        
        }
    }

    u.updateBoundary();
}

void
residualPoisson(const DistributedVector &f,
                const DistributedVector &u,
                DistributedVector &r)
{
    int rhh = (f.N+1)*(f.N+1);
 
    const DistributedVector::Grid &F = f.grid;
    const DistributedVector::Grid &U = u.grid;
    DistributedVector::Grid       &R = r.grid;

    // compute R(i,j)
    
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
printStat(int k,
          const DistributedVector &r,
          const DistributedVector &u,
          const DistributedVector &e)
{
    double rNorm = norm(r);
    double eNorm = errorNorm(e,u);

    if ((r.mpiCart.row==0) && (r.mpiCart.col==0)) {
        cout.width(3);
        cout << k << ") | ";
    
        cout.precision(12);
        cout.setf(std::ios::fixed);
        cout.width(17);
        cout << rNorm << " | ";
    
        cout.precision(12);
        cout.setf(std::ios::fixed);
        cout.width(16);
        cout << eNorm << " | " << endl;
    }
}

int
main(int argc, char *argv[])
{
    MPI::Init(argc, argv);
    
    MpiCart mpiCart(2,2);

    DistributedVector f(mpiCart, 8),
                      e(mpiCart, 8),
                      u(mpiCart, 8),
                      r(mpiCart, 8);
    
    setRHS(f);
    setExactSolution(e);

    
    residualPoisson(f, u, r);    
    printStat(0, r, u, e);
    for (int k=1; k<=10; ++k) {
        gaussSeidelRedBlack(1., f, u);
        
        residualPoisson(f, u, r);
        printStat(k, r, u, e);
    }

    MPI::Finalize();
    return 0;
}
