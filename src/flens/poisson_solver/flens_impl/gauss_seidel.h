#ifndef POISSON_SOLVER_FLENS_IMPL_GAUSS_SEIDEL_H
#define POISSON_SOLVER_FLENS_IMPL_GAUSS_SEIDEL_H 1

#include <flens/flens.h>

namespace flens {

//-- GaussSeidel ---------------------------------------------------------------

template <typename M, typename V>
class GaussSeidel;

template <typename M, typename V>
struct TypeInfo<GaussSeidel<M,V> >
{
    typedef GaussSeidel<M,V> Impl;
    typedef double           ElementType;
};

template <typename M, typename V>
class GaussSeidel
    : public GeneralMatrix<GaussSeidel<M,V> >
{
    public:
        typedef M MatrixType;
        typedef V VectorType;

        GaussSeidel(const M &_A, const V &_f);

        const M &A;
        const V &f;
};

//-- GaussSeidelRedBlack -------------------------------------------------------

template <typename M, typename V>
class GaussSeidelRedBlack;

template <typename M, typename V>
struct TypeInfo<GaussSeidelRedBlack<M,V> >
{
    typedef GaussSeidelRedBlack<M,V> Impl;
    typedef double                   ElementType;
};

template <typename M, typename V>
class GaussSeidelRedBlack
    : public GeneralMatrix<GaussSeidelRedBlack<M,V> >
{
    public:
        typedef M MatrixType;
        typedef V VectorType;

        GaussSeidelRedBlack(const M &_A, const V &_f);

        const M &A;
        const V &f;
};

} // namespace flens

#include <poisson_solver/flens_impl/gauss_seidel.tcc>

#endif // POISSON_SOLVER_FLENS_IMPL_GAUSS_SEIDEL_H
