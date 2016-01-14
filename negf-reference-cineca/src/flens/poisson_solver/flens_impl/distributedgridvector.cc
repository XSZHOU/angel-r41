#include <poisson_solver/flens_impl/distributedgridvector.h>

namespace flens {

//== DistributedGridVector2D ===================================================

DistributedGridVector2D::DistributedGridVector2D()
    : rh(0)
{
}

DistributedGridVector2D::DistributedGridVector2D(const MpiCart &_mpiCart,
                                                 int _rh)
    : mpiCart(_mpiCart), rh(_rh), N(rh-1)
{
    // global matrix has size (0..N+1) x (0..N+1)
    // compute m and n, such that local matrix has size (-1,m+1) x (-1,n+1)
    m = (N+1) / mpiCart.numRows - 1;
    n = (N+1) / mpiCart.numCols - 1;
    grid.resize(_(-1,m+1),_(-1,n+1));

    // m0, n0: row and col offset
    i0 = (m+1)*mpiCart.row;
    j0 = (n+1)*mpiCart.col;

    // firstRow, firstCol:
    // -> non-zero value of local grid is stored in (firstRow, firstCol)
    firstRow = (mpiCart.row==0) ? 1 : 0;
    firstCol = (mpiCart.col==0) ? 1 : 0;

    // types to send rows and cols
    MpiRow = MPI::DOUBLE.Create_contiguous(n+1);
    MpiRow.Commit();

    MpiCol = MPI::DOUBLE.Create_vector(m+2, 1, grid.leadingDimension());
    MpiCol.Commit();

    MpiGrid = MPI::DOUBLE.Create_vector(m+2, n+2, grid.leadingDimension());
    MpiGrid.Commit();
}

DistributedGridVector2D &
DistributedGridVector2D::operator=(double value)
{
    localGrid() = value;
    return *this;
}

DistributedGridVector2D::ConstLocalGrid
DistributedGridVector2D::localGrid() const
{
    return grid(_(firstRow-1,m+1),_(firstCol-1,n+1), firstRow-1, firstCol-1);
}

DistributedGridVector2D::LocalGrid
DistributedGridVector2D::localGrid()
{
    return grid(_(firstRow-1,m+1),_(firstCol-1,n+1), firstRow-1, firstCol-1);
}

void
DistributedGridVector2D::setGhostNodes()
{
    MPI::Request req[8];
    MPI::Status  stat[8];

    const int *neighbors = mpiCart.neighbors;
    req[0] = mpiCart.comm.Isend(&grid(0,n), 1, MpiCol, neighbors[East],  0);
    req[1] = mpiCart.comm.Isend(&grid(0,0), 1, MpiRow, neighbors[South], 0);
    req[2] = mpiCart.comm.Isend(&grid(0,0), 1, MpiCol, neighbors[West],  0);
    req[3] = mpiCart.comm.Isend(&grid(m,0), 1, MpiRow, neighbors[North], 0);

    req[4] = mpiCart.comm.Irecv(&grid(  0,n+1), 1, MpiCol, neighbors[East],  0);
    req[5] = mpiCart.comm.Irecv(&grid( -1,  0), 1, MpiRow, neighbors[South], 0);
    req[6] = mpiCart.comm.Irecv(&grid(  0, -1), 1, MpiCol, neighbors[West],  0);
    req[7] = mpiCart.comm.Irecv(&grid(m+1,  0), 1, MpiRow, neighbors[North], 0);

    MPI::Request::Waitall(8, req, stat);
}

//------------------------------------------------------------------------------

std::ostream &
operator<<(std::ostream &out, const DistributedGridVector2D &v)
{
    DistributedGridVector2D::Grid global, local;
    MpiCart mpiCart = v.mpiCart;
    int N = v.N;

    int coords[2];
    int rank;
    if ((mpiCart.row==0) && (mpiCart.col==0)) {
        global.resize(_(0,N+1),_(0,N+1));
        global(_(0,v.m+1),_(0,v.n+1)) = v.grid(_(0,v.m+1),_(0,v.n+1));
        local.resize(_(-1, v.m+1),_(-1,v.n+1));

        for (int i=0; i<mpiCart.numRows; ++i) {
            for (int j=0; j<mpiCart.numCols; ++j) {
                if ((i==0) && (j==0)) {
                    continue;
                }
                coords[0] = i;
                coords[1] = j;
                rank = mpiCart.comm.Get_cart_rank(coords);
                mpiCart.comm.Recv(&local(0,0), 1, v.MpiGrid, rank,  0);
                int i0 = i*(v.m+1), i1 = (i+1)*(v.m+1);
                int j0 = j*(v.n+1), j1 = (j+1)*(v.n+1);
                global(_(i0,i1), _(j0,j1)) = local(_(0,v.m+1),_(0,v.n+1));
            }
        }
        out << global;
    } else {
        coords[0] = 0;
        coords[1] = 0;
        rank = mpiCart.comm.Get_cart_rank(coords);
        mpiCart.comm.Send(&v.grid(0, 0), 1, v.MpiGrid, rank,  0);
    }
    return out;
}

} // namespace flens
