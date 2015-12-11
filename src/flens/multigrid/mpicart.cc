#include <mpicart.h>

namespace flens {
    
MpiCart::MpiCart()
{
}

MpiCart::MpiCart(int _numRows, int _numCols)
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


    neighbors[NorthEast] = cartRank(row-1, col+1);
    neighbors[SouthEast] = cartRank(row+1, col+1);
    neighbors[SouthWest] = cartRank(row+1, col-1);
    neighbors[NorthWest] = cartRank(row-1, col-1);

}

MpiCart::MpiCart(const MpiCart &rhs)
    : comm(rhs.comm), numRows(rhs.numRows), numCols(rhs.numCols),
      row(rhs.row), col(rhs.col)
{
    for (int k=0; k<8; ++k) {
        neighbors[k] = rhs.neighbors[k];
    }
}

int
MpiCart::cartRank(int row, int col)
{
    if ((row<0) || (col<0) || (row==numRows) || (col==numCols)) {
        return MPI_PROC_NULL;
    }
    int coords[2];
    coords[0] = row;
    coords[1] = col;
    return comm.Get_cart_rank(coords);
}

} // namespace flens
