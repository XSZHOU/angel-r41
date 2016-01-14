#ifndef POISSON_SOLVER_FLENS_IMPL_MPICART_H
#define POISSON_SOLVER_FLENS_IMPL_MPICART_H 1

#include <mpi.h>

namespace flens {

enum CartDir { North=0, East=1, South=2, West=3 };

class MpiCart
{
    public:
        MpiCart();

        MpiCart(int argc, char **argv);

        MpiCart(const MpiCart &rhs);

        MPI::Cartcomm comm;
        int numRows, numCols;
        int row, col;
        int neighbors[4];
};

} // namespace flens

#endif // POISSON_SOLVER_FLENS_IMPL_MPICART_H
