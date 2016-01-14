#ifndef POISSON_FLENS_DISTRIBUTEDGRIDVECTOR_H
#define POISSON_FLENS_DISTRIBUTEDGRIDVECTOR_H 1

#include <flens/flens.h>
#include <poisson_solver/flens_impl/mpicart.h>

namespace flens {

//== DistributedGridVector2D ===================================================

class DistributedGridVector2D;

template <>
struct TypeInfo<DistributedGridVector2D>
{
    typedef DistributedGridVector2D Impl;
    typedef double                  ElementType;
};

//------------------------------------------------------------------------------

class DistributedGridVector2D
    : public Vector<DistributedGridVector2D>
{
    public:
        typedef GeMatrix<FullStorage<double, RowMajor> > Grid;
        typedef Grid::ConstView                          ConstLocalGrid;
        typedef Grid::View                               LocalGrid;

        DistributedGridVector2D();

        DistributedGridVector2D(const MpiCart &_mpiCart, int rh);

        DistributedGridVector2D &
        operator=(double value);

        template <typename RHS>
        DistributedGridVector2D &
        operator=(const Vector<RHS> &rhs);

        template <typename RHS>
        DistributedGridVector2D &
        operator+=(const Vector<RHS> &rhs);

        template <typename RHS>
        DistributedGridVector2D &
        operator-=(const Vector<RHS> &rhs);

        ConstLocalGrid
        localGrid() const;

        LocalGrid
        localGrid();

        void
        setGhostNodes();

        MpiCart         mpiCart;
        int             rh, N;
        int             m, n;
        int             i0, j0;
        int             firstRow, firstCol;
        Grid            grid;
        MPI::Datatype   MpiRow, MpiCol, MpiGrid;
};

//------------------------------------------------------------------------------

std::ostream &
operator<<(std::ostream &out, const DistributedGridVector2D &v);

} // namespace flens

#include <poisson_solver/flens_impl/distributedgridvector.tcc>

#endif // POISSON_SOLVER_FLENS_IMPL_DISTRIBUTEDGRIDVECTOR_H
