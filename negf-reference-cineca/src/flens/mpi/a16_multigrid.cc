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

enum CartDir { North=0, East=1, South=2, West=3};

struct MpiCart
{
    MpiCart()
    {}
    
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
        comm.Shift(0, 1, neighbors[North], neighbors[South]);
        comm.Shift(1, 1, neighbors[West],  neighbors[East]);
    }
    
    MpiCart(const MpiCart &rhs)
        : comm(rhs.comm),
          numRows(rhs.numRows), numCols(rhs.numCols),
          row(rhs.row), col(rhs.col)
    {
        for (int k=0; k<8; ++k) {
            neighbors[k] = rhs.neighbors[k];
        }
    }

    int
    cartRank(int row, int col)
    {
        if ((row<0) || (col<0) || (row==numRows) || (col==numCols)) {
            return MPI_PROC_NULL;
        }
        int coords[2];
        coords[0] = row;
        coords[1] = col;
        return comm.Get_cart_rank(coords);
    }
    
    MPI::Cartcomm comm;
    int numRows, numCols;
    int row, col;
    int neighbors[4];
};

struct DistributedVector
{
    typedef GeMatrix<FullStorage<double, RowMajor> > Grid;

    DistributedVector()
    {
    }

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
        // MpiRow = MPI::DOUBLE.Create_contiguous(n+1); 
        MpiRow = MPI::DOUBLE.Create_contiguous(n+2); 
        MpiRow.Commit();
            
        MpiCol = MPI::DOUBLE.Create_vector(m+1, 1, grid.leadingDimension()); 
        MpiCol.Commit();
    }
    
    void
    updateBoundary()
    {
        MPI::Request req[8];
        MPI::Status  stat[8];
    
        req[ 0] = mpiCart.comm.Isend(&grid(0,n), 1, MpiCol, mpiCart.neighbors[East],  0);  
        req[ 1] = mpiCart.comm.Isend(&grid(m,0), 1, MpiRow, mpiCart.neighbors[South], 0);    
        req[ 2] = mpiCart.comm.Isend(&grid(0,0), 1, MpiCol, mpiCart.neighbors[West],  0);  
        req[ 3] = mpiCart.comm.Isend(&grid(0,0), 1, MpiRow, mpiCart.neighbors[North], 0);    

        req[ 4] = mpiCart.comm.Irecv(&grid(  0,n+1), 1, MpiCol, mpiCart.neighbors[East],  0);
        req[ 5] = mpiCart.comm.Irecv(&grid(m+1,  0), 1, MpiRow, mpiCart.neighbors[South], 0);
        req[ 6] = mpiCart.comm.Irecv(&grid(  0, -1), 1, MpiCol, mpiCart.neighbors[West],  0);
        req[ 7] = mpiCart.comm.Irecv(&grid( -1,  0), 1, MpiRow, mpiCart.neighbors[North], 0);
        
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
    
    for (int i=f.firstRow; i<=f.m; ++i) {
        for (int j=f.firstCol; j<=f.n; ++j) {
            R(i,j) = F(i,j) - rhh*(4*U(i,j) -U(i+1,j) -U(i-1,j) -U(i,j+1) -U(i,j-1));
        }
    }
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

// f_c = R*r
void
restriction(const DistributedVector &r, DistributedVector &f_c)
{
    const DistributedVector::Grid &R   = r.grid;
    DistributedVector::Grid       &F_C = f_c.grid;
    
    int I = 2*f_c.firstRow;
    for (int i=f_c.firstRow; i<=f_c.m; ++i, I+=2) {
        int J = 2*f_c.firstCol;
        for (int j=f_c.firstCol; j<=f_c.n; ++j, J+=2) {
            F_C(i,j) = (                 R(I-1,J)
                        + R(I  ,J-1) + 4*R(I  ,J) + R(I  ,J+1)
                                        +R(I+1,J))/8; 
        }
    }
}

// u = u + P*u_c
void
prolongation(const DistributedVector &u_c, DistributedVector &u)
{
      
    const DistributedVector::Grid &U_C   = u_c.grid;
    DistributedVector::Grid       &U     = u.grid;

    int I = 2*u_c.firstRow;
    for (int i=u_c.firstRow; i<=u_c.m; ++i, I+=2) {
        int J = 2*u_c.firstCol;
        for (int j=u_c.firstCol; j<=u_c.n; ++j, J+=2) {
            U(I-1,J-1)+=.25*U_C(i,j);  U(I-1,J)+=.5*U_C(i,j); U(I-1,J+1)+=.25*U_C(i,j);
            U(I  ,J-1)+=.50*U_C(i,j);  U(I  ,J)+=   U_C(i,j); U(I  ,J+1)+=.50*U_C(i,j);
            U(I+1,J-1)+=.25*U_C(i,j);  U(I+1,J)+=.5*U_C(i,j); U(I+1,J+1)+=.25*U_C(i,j);
        }
    }
    
    int i = u_c.m+1; 
    I = 2*i;
    int J = 2*u_c.firstCol;
    for (int j=u_c.firstCol; j<=u_c.n; ++j, J+=2) {
        U(I-1,J-1)+=.25*U_C(i,j);  U(I-1,J)+=.5*U_C(i,j); U(I-1,J+1)+=.25*U_C(i,j);
        U(I  ,J-1)+=.50*U_C(i,j);  U(I  ,J)+=   U_C(i,j); U(I  ,J+1)+=.50*U_C(i,j);
    }


    int j = u_c.n+1; 
    J = 2*j;
    I = 2*u_c.firstRow;
    for (int i=u_c.firstRow; i<=u_c.m; ++i, I+=2) {
        U(I-1,J-1)+=.25*U_C(i,j);  U(I-1,J)+=.5*U_C(i,j);
        U(I  ,J-1)+=.50*U_C(i,j);  U(I  ,J)+=   U_C(i,j);
        U(I+1,J-1)+=.25*U_C(i,j);  U(I+1,J)+=.5*U_C(i,j);
    }

    i = u_c.m+1; 
    j = u_c.n+1; 
    J = 2*j;
    I = 2*i;

    U(I-1,J-1)+=.25*U_C(i,j);  U(I-1,J)+=.5*U_C(i,j);
    U(I  ,J-1)+=.50*U_C(i,j);  U(I  ,J)+=   U_C(i,j);
    
    u.updateBoundary();
}

void
write(const DistributedVector &v, const string &name)
{
    ostringstream s;
    s << name << "_local_" << v.mpiCart.row << "_" << v.mpiCart.col << ".dat";
    ofstream out(s.str().c_str());
    out << v.N << " " << v.N << endl;
    out << v.m-v.firstRow+1 << " " << v.n-v.firstCol+1<< endl; 
    out << v.grid(_(v.firstRow, v.m),_(v.firstCol, v.n));
}

void
read(DistributedVector::Grid &grid, const string &name, int numProcX, int numProcY)
{
    int row = 1;
    for (int px=0; px<numProcX; ++px) {
        int col = 1;
	    for (int py=0; py<numProcY; ++py) {
	        ostringstream s;
	        s << name << "_local_" << px << "_" << py << ".dat";
	        ifstream in(s.str().c_str());
	        int N, M;
            in >> N >> M;
	        if ((px==0) && (py==0)) {
		        grid.resize(N+2,M+2,0,0);
	        }
            int r; int c;
            in >> r; in >> c;
	        for (int i=0; i<r; ++i) {
                for (int j=0; j<c; ++j) {
                    in >> grid(row+i,col+j);
		        }
	        }
            col += c;
	        if (py==numProcY-1) {
	            row += r;
            }
	    }
    }
}

//------------------------------------------------------------------------------

const long l = 9;
const int  n = (1<<l) - 1;  // n = 2^l -1

DistributedVector  u[l+1], f[l+1], r[l+1];

void
initMultigrid(const MpiCart &mpiCart)
{
    typedef DistributedVector Vec;

    int N = n;
    for (int i=l; i>=1; --i, N/=2) {
        u[i] = Vec(mpiCart, N+1);
        f[i] = Vec(mpiCart, N+1);
        r[i] = Vec(mpiCart, N+1);
    }
    for (int i=l; i>=3; --i) {
        restriction(f[i], f[i-1]);
    }
}

void
multigrid(int l, int v1, int v2)
{
    if (l==2) {
        for (int k=1; k<=1000; ++k) {
            gaussSeidelRedBlack(1., f[l], u[l]);
        }
    } else {
        for (int v=1; v<=v1; ++v) {
            gaussSeidelRedBlack(1., f[l], u[l]);
        }
        residualPoisson(f[l], u[l], r[l]);

        restriction(r[l], f[l-1]);


        u[l-1].grid = 0;
        multigrid(l-1,v1,v2);
        
        prolongation(u[l-1], u[l]);
        
        
        for (int v=1; v<=v2; ++v) {
            gaussSeidelRedBlack(1., f[l], u[l]);
        }
    }
}

void
fullMultigrid()
{
    for (int k=1; k<=1000; ++k) {
        gaussSeidelRedBlack(1., f[2], u[2]);
    }
    
    for (int i=2; i<l; ++i) {
        prolongation(u[i], u[i+1]);
        multigrid(i+1, 1, 1);
    }
}

//------------------------------------------------------------------------------


int
main(int argc, char *argv[])
{
    MPI::Init(argc, argv);
    
    MpiCart mpiCart(1,1);

    initMultigrid(mpiCart);
    
    DistributedVector e(mpiCart, n+1);
    setRHS(f[l]);
    setExactSolution(e);
    
    residualPoisson(f[l], u[l], r[l]);  
    printStat(0, r[l], u[l], e);
    
    fullMultigrid();
    printStat(1, r[l], u[l], e);
    
    for (int k=2; k<=20; ++k) {
        multigrid(l,1,1);
        printStat(k, r[l], u[l], e);
    }
    MPI::Finalize();
}
