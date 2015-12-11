#ifndef MULTIGRID_POISSON2D_H
#define MULTIGRID_POISSON2D_H 1

#include <flens/matvec.h>
#include <flens/symmetricmatrix.h>

namespace flens {

DefaultTypeInfo(Poisson2D,double)

struct Poisson2D
    : public SymmetricMatrix<Poisson2D>
{
    Poisson2D();
    
    Poisson2D(int _rh);

    int rh;
};

} // namespace flens

#endif // MULTIGRID_POISSON2D_H