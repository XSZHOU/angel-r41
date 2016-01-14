#ifndef MULTIGRID_DISTRIBUTED_VECTOR_H
#define MULTIGRID_DISTRIBUTED_VECTOR_H 1

#include <flens/matvec.h>
#include <flens/fullstorage.h>
#include <flens/generalmatrix.h>
#include <mpicart.h>

namespace flens {

DefaultTypeInfo(DistributedVector,double);

struct DistributedVector
    : Vector<DistributedVector>
{
    typedef GeMatrix<FullStorage<double, RowMajor> > Grid;

    DistributedVector();

    DistributedVector(const MpiCart &_mpiCart, int rh);
        
    template <typename RHS>
        DistributedVector &
        operator=(const Vector<RHS> & rhs);
        
    template <typename RHS>
        DistributedVector &
        operator+=(const Vector<RHS> & rhs);
    
    void
    updateBoundary();
    
    void
    updateEdges();
    
    //---
    
    MpiCart         mpiCart;
    int             N;
    int             m, n;
    int             m0, n0;
    int             firstRow, firstCol;
    Grid            grid;
    MPI::Datatype   MpiRow, MpiCol;
};

} // namespace flens

#include <distributedvector.tcc>

#endif // MULTIGRID_DISTRIBUTED_VECTOR_H
