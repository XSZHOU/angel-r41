#ifndef POISSON_SOLVER_KERNEL_IMPL_FASTPOISSONSOLVERS_H
#define POISSON_SOLVER_KERNEL_IMPL_FASTPOISSONSOLVERS_H 1

#include <fftw3.h>

namespace flens {

//-- dirichlet poisson 1d ------------------------------------------------------

template <typename U>
    fftw_plan
    dp1d_fastpoissonsolver_init(int rh, DenseVector<U> &u,
                                unsigned fftw_flags = FFTW_ESTIMATE);


template <typename F, typename U>
    void
    dp1d_fastpoissonsolver(int rh, fftw_plan &plan,
                           const DenseVector<F> &f, DenseVector<U> &u);


//-- dirichlet poisson 2d ------------------------------------------------------

template <typename U>
    fftw_plan
    dp2d_fastpoissonsolver_init(int rh, GeMatrix<U> &u,
                                unsigned fftw_flags = FFTW_ESTIMATE);


template <typename F, typename U>
    void
    dp2d_fastpoissonsolver(int rh, fftw_plan &plan,
                           const GeMatrix<F> &f, GeMatrix<U> &u);

} // namespace flens

#include <poisson_solver/kernel_impl/fastpoissonsolvers.tcc>

#endif // POISSON_SOLVER_KERNEL_IMPL_FASTPOISSONSOLVERS_H
