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
    MPI_Datatype   MpiRow, MpiCol;
    
    MPI_Type_contiguous(4, MPI_DOUBLE, &MpiRow); 
    MPI_Type_commit(&MpiRow);
            
    MPI_Type_vector(4, 1, A.leadingDimension(), MPI_DOUBLE, &MpiCol); 
    MPI_Type_commit(&MpiCol);
    
    // create mpi_cart:  2-dimensional, 2 x 2, non-periodic
    int dims[2], periods[2];

    dims[0]    = 2;
    dims[1]    = 2;
    periods[0] = 0;
    periods[1] = 0;
    
    MPI_Comm comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm);
    
    // find your position in the cart
    int rank = MPI::COMM_WORLD.Get_rank();
    int coords[2];
    MPI_Cart_coords(comm, rank, 2, coords);
    int row = coords[0];
    int col = coords[1];
    
    // find your neighbors
    int neighbors[4];
    enum CartDir { North=0, East=1, South=2, West=3 };
    MPI_Cart_shift(comm, 0, 1, &neighbors[North], &neighbors[South]);
    MPI_Cart_shift(comm, 1, 1, &neighbors[West],  &neighbors[East]);

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
    MPI_Request req[8];
    MPI_Status  stat[8];
    
    MPI_Isend(&A(0,3), 1, MpiCol, neighbors[East],  0, comm, &req[0]);  
    MPI_Isend(&A(3,0), 1, MpiRow, neighbors[South], 0, comm, &req[1]);    
    MPI_Isend(&A(0,0), 1, MpiCol, neighbors[West],  0, comm, &req[2]);  
    MPI_Isend(&A(0,0), 1, MpiRow, neighbors[North], 0, comm, &req[3]);    

    MPI_Irecv(&A( 0, 4), 1, MpiCol, neighbors[East],  0, comm, &req[4]);
    MPI_Irecv(&A( 4, 0), 1, MpiRow, neighbors[South], 0, comm, &req[5]);
    MPI_Irecv(&A( 0,-1), 1, MpiCol, neighbors[West],  0, comm, &req[6]);
    MPI_Irecv(&A(-1, 0), 1, MpiRow, neighbors[North], 0, comm, &req[7]);
        
    MPI_Waitall(8, req, stat);
    
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
    MPI_Allreduce(&r, &R, 1, MPI_DOUBLE, MPI_SUM, comm);
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
	
	// free types
	MPI_Type_free(&MpiRow);
	MPI_Type_free(&MpiCol);
    
    MPI::Finalize();
    return 0;
}
