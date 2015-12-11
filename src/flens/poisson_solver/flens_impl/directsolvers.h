#ifndef POISSON_SOLVER_FLENS_IMPL_DIRECTSOLVERS_H
#define POISSON_SOLVER_FLENS_IMPL_DIRECTSOLVERS_H 1

#include <fftw3.h>

namespace flens {

//-- stationary iterative solver -----------------------------------------------

template <typename Method>
class StationaryIterativeSolver
{
    public:
        typedef typename Method::MatrixType MatrixType;
        typedef typename Method::VectorType VectorType;

        StationaryIterativeSolver(const MatrixType &A, const VectorType &f,
                                  VectorType &_u);

        void
        solve();

    private:
        Method            S;
        const MatrixType  &A;
        const VectorType  &f;
        VectorType        &u, r;
};

//-- fast poisson solver -------------------------------------------------------

template <typename MatrixType, typename VectorType>
class FastPoissonSolver
{
    public:
        FastPoissonSolver(const MatrixType &A, const VectorType &f,
                          VectorType &_u);

        ~FastPoissonSolver();

        void
        solve();

    private:
        const VectorType &f;
        VectorType       &u;
        fftw_plan        plan;
};

//-- cholesky factorization ----------------------------------------------------

template <typename MatrixType, typename VectorType>
class Cholesky
{
    public:
        Cholesky(const MatrixType &A, const VectorType &f,
                 VectorType &_u);

        void
        solve();

    private:
        const VectorType &f;
        VectorType       &u;
};

} // namespace flens

#include <poisson_solver/flens_impl/directsolvers.tcc>

#endif // POISSON_SOLVER_FLENS_IMPL_DIRECTSOLVERS_H
