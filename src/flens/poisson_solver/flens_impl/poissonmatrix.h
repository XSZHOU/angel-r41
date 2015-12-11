#ifndef POISSON_SOLVER_FLENS_IMPL_POISSONMATRIX_H
#define POISSON_SOLVER_FLENS_IMPL_POISSONMATRIX_H 1

#include <flens/flens.h>

namespace flens {

//-- DirichletPoisson1D --------------------------------------------------------

class DirichletPoisson1D;

template <>
struct TypeInfo<DirichletPoisson1D>
{
    typedef DirichletPoisson1D Impl;
    typedef double             ElementType;
};

class DirichletPoisson1D
    : public SymmetricMatrix<DirichletPoisson1D>
{
    public:
        DirichletPoisson1D();

        DirichletPoisson1D(int _rh);

        int rh;
};

//-- NeumannPoisson1D --------------------------------------------------------

class NeumannPoisson1D;

template <>
struct TypeInfo<NeumannPoisson1D>
{
    typedef NeumannPoisson1D Impl;
    typedef double           ElementType;
};

class NeumannPoisson1D
    : public GeneralMatrix<NeumannPoisson1D>
{
    public:
        NeumannPoisson1D();

        NeumannPoisson1D(int _rh);

        int rh;
};

//-- DirichletPoisson2D --------------------------------------------------------

class DirichletPoisson2D;

template <>
struct TypeInfo<DirichletPoisson2D>
{
    typedef DirichletPoisson2D Impl;
    typedef double             ElementType;
};

class DirichletPoisson2D
    : public SymmetricMatrix<DirichletPoisson2D>
{
    public:
        DirichletPoisson2D();

        DirichletPoisson2D(int _rh);

        int rh;
};

//-- NeumannPoisson2D ----------------------------------------------------------

class NeumannPoisson2D;

template <>
struct TypeInfo<NeumannPoisson2D>
{
    typedef NeumannPoisson2D Impl;
    typedef double           ElementType;
};

class NeumannPoisson2D
    : public SymmetricMatrix<NeumannPoisson2D>
{
    public:
        NeumannPoisson2D();

        NeumannPoisson2D(int _rh);

        int rh;
};

} // namespace flens

#endif // POISSON_SOLVER_FLENS_IMPL_POISSONMATRIX_H
