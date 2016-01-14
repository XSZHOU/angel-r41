#ifndef POISSON_SOLVER_FLENS_IMPL_RESTRICTION_H
#define POISSON_SOLVER_FLENS_IMPL_RESTRICTION_H 1

#include <flens/flens.h>

namespace flens {

//-- dirichlet poisson 1d ------------------------------------------------------

class Restriction;

template <>
struct TypeInfo<Restriction>
{
    typedef Restriction Impl;
    typedef double      ElementType;
};

class Restriction
    : public GeneralMatrix<Restriction>
{
};

//-- neumann poisson 1d --------------------------------------------------------

class Restriction_NBC;

template <>
struct TypeInfo<Restriction_NBC>
{
    typedef Restriction_NBC Impl;
    typedef double          ElementType;
};

class Restriction_NBC
    : public GeneralMatrix<Restriction_NBC>
{
};

} // namespace flens

#endif // POISSON_SOLVER_FLENS_IMPL_RESTRICTION_H
