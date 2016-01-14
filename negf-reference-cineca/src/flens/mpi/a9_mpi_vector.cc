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
    
    
    int  dims[2];
    bool periods[2];

    dims[0]    = 2;
    dims[1]    = 2;
    periods[0] = false;
    periods[1] = false;
    MPI::Cartcomm comm = MPI::COMM_WORLD.Create_cart(2, dims, periods, true);
    
    // find your position in the cart
    int rank = MPI::COMM_WORLD.Get_rank();
    int coords[2];
    comm.Get_coords(rank, 2, coords);
    int row = coords[0];
    int col = coords[1];
    
    // find your neighbors
    int neighbors[4];
    enum CartDir { North=0, East=1, South=2, West=3 };
    comm.Shift(0, 1, neighbors[North], neighbors[South]);
    comm.Shift(1, 1, neighbors[West],  neighbors[East]);

    // print one local matrix (for instance that for rank==0)
    if (rank==0) {
        cout << A << endl;

        cout << "rank:   " << rank << endl;
        cout << "pos:    " << row << ", " << col << endl;
        cout << "North = " << neighbors[North] << endl;
        cout << "East  = " << neighbors[East] << endl;
        cout << "South = " << neighbors[South] << endl;
        cout << "West  = " << neighbors[West] << endl;
        
    }

    // update boundaries 
    MPI::Request req[8];
    MPI::Status  stat[8];
    
    req[0] = comm.Isend(&A(0,3), 1, MpiCol, neighbors[East],  0);  
    req[1] = comm.Isend(&A(3,0), 1, MpiRow, neighbors[South], 0);    
    req[2] = comm.Isend(&A(0,0), 1, MpiCol, neighbors[West],  0);  
    req[3] = comm.Isend(&A(0,0), 1, MpiRow, neighbors[North], 0);    

    req[4] = comm.Irecv(&A( 0, 4), 1, MpiCol, neighbors[East],  0);
    req[5] = comm.Irecv(&A( 4, 0), 1, MpiRow, neighbors[South], 0);
    req[6] = comm.Irecv(&A( 0,-1), 1, MpiCol, neighbors[West],  0);
    req[7] = comm.Irecv(&A(-1, 0), 1, MpiRow, neighbors[North], 0);
        
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
    comm.Allreduce(&r, &R, 1, MPI::DOUBLE, MPI::SUM);
    r = sqrt(R);
    
    if (rank==0) {
        // should be 21.633308
        cout << "norm = " << r << endl;
    }
    
    // store local matrix
    ostringstream s;
	s << "A_local_" << row << "_" << col << ".dat";
	ofstream file(s.str().c_str());
	file << A << endl;
    
    MPI::Finalize();
    return 0;
}
