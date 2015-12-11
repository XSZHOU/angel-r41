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

    /*
     * this functions create and return MPI::Datatypes
     *  
     *  MPI::Datatype
     *  MPI::DOUBLE.Create_contiguous(int len);    
     *
     *  MPI::Datatype
     *  MPI::DOUBLE.Create_vector(int len, int blocklen, int stride);
     *
     *  returns stride between rows of A:
     *      A.leadingDimension()
     */
    
    // MpiRow = ...
    MpiRow.Commit();
            
    // MpiCol = ...
    MpiCol.Commit();
    
    // create mpi_cart:  2-dimensional, 2 x 2, non-periodic
    int  dims[2];
    bool periods[2];

    // dims[0]  = ...
    // ...
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
    // find neighbors[North], neighbors[East], neighbors[South], neighbors[West]
    //     comm.Shift(dim, disp, source, destination);

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
    
    // what is this line doing???
    req[0] = comm.Isend(&A(0,3), 1, MpiCol, neighbors[East],  0);  
    //  ... also send to other directions    

    req[4] = comm.Irecv(&A( 0, 4), 1, MpiCol, neighbors[East],  0);
    // ... also receive from others
        
    // wait for all requests:
    //    MPI::Request::Waitall(int numReq, MPI::Request req[], MPI::Status stat[]);
    
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
