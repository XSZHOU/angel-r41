#include <flens/flens.h>
#include <distributedvector.h>

namespace flens {

DistributedVector::DistributedVector()
{
}

DistributedVector::DistributedVector(const MpiCart &_mpiCart, int rh)
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
    MpiRow = MPI::DOUBLE.Create_contiguous(n+2); 
    MpiRow.Commit();
        
    MpiCol = MPI::DOUBLE.Create_vector(m+1, 1, grid.leadingDimension()); 
    MpiCol.Commit();
}
    
void
DistributedVector::updateBoundary()
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
    
void
DistributedVector::updateEdges()
{
    MPI::Request req[2];
    MPI::Status  stat[2];
    
    req[0] = mpiCart.comm.Isend(&grid(0,0), 1, MPI::DOUBLE, mpiCart.neighbors[NorthWest],  0);  
    req[1] = mpiCart.comm.Irecv(&grid(m+1,n+1), 1, MPI::DOUBLE, mpiCart.neighbors[SouthEast],  0);
    
    MPI::Request::Waitall(2, req, stat);
}

//------------------------------------------------------------------------------

void
copy(const DistributedVector &x, DistributedVector &y)
{
    y.operator=(x);
}

void
axpy(double alpha, const DistributedVector &x, DistributedVector &y)
{
    assert(alpha==1.);
    
    y.grid += x.grid;
}


} // namespace flens
