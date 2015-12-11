#ifndef MULTIGRID_MPICART_H
#define MULTIGRID_MPICART_H 1

#include <mpi.h>

namespace flens {
    
enum CartDir { North=0, East=1, South=2, West=3,
               NorthEast=4, SouthEast=5, SouthWest=6, NorthWest=7};
    
struct MpiCart
{
    MpiCart();

    MpiCart(int _numRows, int _numCols);

    MpiCart(const MpiCart &rhs);

    int
        cartRank(int row, int col);

    MPI::Cartcomm comm;
    int numRows, numCols;
    int row, col;
    int neighbors[8];
};
    
} // namespace flens

#endif // MULTIGRID_MPICART_H
