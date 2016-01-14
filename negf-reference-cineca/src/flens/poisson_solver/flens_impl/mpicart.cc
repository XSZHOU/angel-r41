#include <poisson_solver/flens_impl/mpicart.h>
#include <cstdlib>
#include <iostream>

namespace flens {

MpiCart::MpiCart()
{
}

MpiCart::MpiCart(int argc, char **argv)
    : numRows(1), numCols(1)
{
    MPI::Init(argc, argv);
    if (argc>=3) {
        numRows = std::atoi(argv[1]);
        numCols = std::atoi(argv[2]);
    } else {
        if (MPI::COMM_WORLD.Get_rank()==0) {
            std::cerr << "usage: " << argv[0] << " numRows numCols" << std::endl;
        }
        std::exit(0);
    }

    // create mpi_cart:  2-dimensional, non-periodic
    int  dims[2];
    bool periods[2];

    dims[0]    = numRows;
    dims[1]    = numCols;
    periods[0] = false;
    periods[1] = false;
    comm = MPI::COMM_WORLD.Create_cart(2, dims, periods, true);

    // find your position in the cart
    int rank = comm.Get_rank();
    int coords[2];
    comm.Get_coords(rank, 2, coords);
    row = coords[0];
    col = coords[1];

    // find your neighbors
    enum CartDir { North=0, East=1, South=2, West=3 };
    comm.Shift(0, 1, neighbors[South], neighbors[North]);
    comm.Shift(1, 1, neighbors[West],  neighbors[East]);
}

MpiCart::MpiCart(const MpiCart &rhs)
    : comm(rhs.comm),
      numRows(rhs.numRows), numCols(rhs.numCols),
      row(rhs.row), col(rhs.col)
{
    neighbors[North] = rhs.neighbors[North];
    neighbors[East]  = rhs.neighbors[East];
    neighbors[South] = rhs.neighbors[South];
    neighbors[West]  = rhs.neighbors[West];
}

} // namespace flens
