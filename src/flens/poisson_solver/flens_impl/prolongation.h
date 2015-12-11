#ifndef POISSON_SOLVER_FLENS_IMPL_PROLONGATION_H
#define POISSON_SOLVER_FLENS_IMPL_PROLONGATION_H 1

#include <flens/flens.h>

namespace flens {

//-- dirichlet poisson 1d ------------------------------------------------------

class Prolongation;

template <>
struct TypeInfo<Prolongation>
{
    typedef Prolongation Impl;
    typedef double       ElementType;
};

class Prolongation
    : public GeneralMatrix<Prolongation>
{
};

//-- neumann poisson 1d --------------------------------------------------------

class Prolongation_NBC;

template <>
struct TypeInfo<Prolongation_NBC>
{
    typedef Prolongation_NBC Impl;
    typedef double           ElementType;
};

class Prolongation_NBC
    : public GeneralMatrix<Prolongation_NBC>
{
};

} // namespace flens

#endif // POISSON_SOLVER_FLENS_IMPL_PROLONGATION_H
