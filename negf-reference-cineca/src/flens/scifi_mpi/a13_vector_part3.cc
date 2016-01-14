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
    return sqrt(R);
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

    // ...

    u.updateBoundary();

    // black nodes

    // ...

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

int
main(int argc, char *argv[])
{
    MPI::Init(argc, argv);
    
    MpiCart mpiCart(2,2);

    DistributedVector f(mpiCart, 8),
                      e(mpiCart, 8),
                      u(mpiCart, 8);
    
    setRHS(f);
    setExactSolution(e);
    
    gaussSeidelRedBlack(1., f, u);
    
    // print one local matrix (for instance that for rank==0)
    if ((mpiCart.row==0) && (mpiCart.col==0)) {
        cout << "f = " << f.grid << endl;
        cout << "e = " << e.grid << endl;
        cout << "u = " << u.grid << endl;
    }

    write(f, "f");
    write(e, "e");
    write(u, "u");
    
    MPI::COMM_WORLD.Barrier();    

    DistributedVector::Grid ff, ee, uu;
    if ((mpiCart.row==0) && (mpiCart.col==0)) {
        read(ff, "f", 2, 2);
        read(ee, "e", 2, 2);
        read(uu, "u", 2, 2);
        cout << ff << endl;
        cout << ee << endl;
        cout << uu << endl;
    }
    
    MPI::Finalize();
    return 0;
}
