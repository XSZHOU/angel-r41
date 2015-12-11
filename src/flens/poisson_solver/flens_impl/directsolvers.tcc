#include <iostream>
#include <limits>

namespace flens {

//-- stationary iterative solver -----------------------------------------------

template <typename Meth>
StationaryIterativeSolver<Meth>::StationaryIterativeSolver(const MatrixType &_A,
                                                           const VectorType &_f,
                                                           VectorType &_u)
    : S(_A, _f), A(_A), f(_f), u(_u), r(f)
{
}

template <typename Method>
void
StationaryIterativeSolver<Method>::solve()
{
    double eps = std::numeric_limits<double>::epsilon();

    r = f - A*u;
    for (int i=1; (i<=5000) && (normL2(r)>eps); ++i) {
        u = S*u;
        r = f - A*u;
    }
}

//-- fast poisson solver -------------------------------------------------------

template <typename MatrixType, typename VectorType>
FastPoissonSolver<MatrixType, VectorType>::~FastPoissonSolver()
{
    // fftw_destroy_plan(plan);
}

} // namespace flens
