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

int
main(int argc, char *argv[])
{
    MPI::Init(argc, argv);

    GeMatrix<FullStorage<double, RowMajor> > A(_(-1,4),_(-1,4));
    
    A = 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 5, 0,
        0, 2, 0, 0, 5, 0,
        0, 2, 0, 0, 5, 0,
        0, 2, 3, 3, 3, 0,
        0, 0, 0, 0, 0, 0;
    
    // types to send rows and cols
    MPI::Datatype   MpiRow, MpiCol;
    
    MpiRow = MPI::DOUBLE.Create_contiguous(4); 
    MpiRow.Commit();
            
    MpiCol = MPI::DOUBLE.Create_vector(4, 1, A.leadingDimension()); 
    MpiCol.Commit();
    
    // create mpi_cart:  2-dimensional, 2 x 2, non-periodic
    MpiCart mpiCart(2,2);

    // print one local matrix (for instance that for rank==0)
    if ((mpiCart.row==0) && (mpiCart.col==0)) {
        cout << A << endl;

        cout << "pos:    " << mpiCart.row << ", " << mpiCart.col << endl;
        cout << "North = " << mpiCart.neighbors[North] << endl;
        cout << "East  = " << mpiCart.neighbors[East] << endl;
        cout << "South = " << mpiCart.neighbors[South] << endl;
        cout << "West  = " << mpiCart.neighbors[West] << endl;
        
    }

    // update boundaries 
    MPI::Request req[8];
    MPI::Status  stat[8];
    
    req[0] = mpiCart.comm.Isend(&A(0,3), 1, MpiCol, mpiCart.neighbors[East],  0);  
    req[1] = mpiCart.comm.Isend(&A(3,0), 1, MpiRow, mpiCart.neighbors[South], 0);    
    req[2] = mpiCart.comm.Isend(&A(0,0), 1, MpiCol, mpiCart.neighbors[West],  0);  
    req[3] = mpiCart.comm.Isend(&A(0,0), 1, MpiRow, mpiCart.neighbors[North], 0);    

    req[4] = mpiCart.comm.Irecv(&A( 0, 4), 1, MpiCol, mpiCart.neighbors[East],  0);
    req[5] = mpiCart.comm.Irecv(&A( 4, 0), 1, MpiRow, mpiCart.neighbors[South], 0);
    req[6] = mpiCart.comm.Irecv(&A( 0,-1), 1, MpiCol, mpiCart.neighbors[West],  0);
    req[7] = mpiCart.comm.Irecv(&A(-1, 0), 1, MpiRow, mpiCart.neighbors[North], 0);
        
    MPI::Request::Waitall(8, req, stat);
    
    // compute (frobenius) norm
    // -> local norm
    double r = 0;
    for (int i=0; i<=3; ++i) {
        for (int j=0; j<=3; ++j) {
            r += A(i,j)*A(i,j);
        }
    }

    // -> global norm
    double R;
    mpiCart.comm.Allreduce(&r, &R, 1, MPI::DOUBLE, MPI::SUM);
    r = sqrt(R);
    
    if ((mpiCart.row==0) && (mpiCart.col==0)) {
        // should be 21.633308
        cout << "norm = " << r << endl;
    }
    
    // store local matrix
    ostringstream s;
	s << "A_local_" << mpiCart.row << "_" << mpiCart.col << ".dat";
	ofstream file(s.str().c_str());
	file << A << endl;
    
    MPI::Finalize();
    return 0;
}
