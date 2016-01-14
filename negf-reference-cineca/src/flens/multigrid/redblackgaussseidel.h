#ifndef MULTIGRID_REDBLACKGAUSSSEIDEL_H
#define MULTIGRID_REDBLACKGAUSSSEIDEL_H 1

#include <flens/symmetricmatrix.h>
#include <distributedvector.h>

namespace flens {
    
DefaultTypeInfo(RedBlackGaussSeidel, double);

struct RedBlackGaussSeidel
    : public SymmetricMatrix<RedBlackGaussSeidel>
{   
    RedBlackGaussSeidel(const DistributedVector &x, double omega);
    
    const DistributedVector &x;
    double omega;
};

void
mv(double alpha, const RedBlackGaussSeidel &A, const DistributedVector &x,
   double beta, DistributedVector &y);

void
mv(const RedBlackGaussSeidel &A, DistributedVector &y);
  
} // namespace flens

#endif // MULTIGRID_REDBLACKGAUSSSEIDEL_H