#ifndef POISSON_SOLVER_KERNEL_IMPL_DP_BLAS_H
#define POISSON_SOLVER_KERNEL_IMPL_DP_BLAS_H 1

namespace flens {

//-- dirichlet poisson 2d ------------------------------------------------------

template <typename X>
    void
    dp2d_copy(const GeMatrix<X> &x,  GeMatrix<X> &y);

template <typename X>
    void
    dp2d_scal(double alpha, GeMatrix<X> &x);

template <typename X, typename Y>
    double
    dp2d_dot(int rh, const GeMatrix<X> &x, const GeMatrix<Y> &y);

template <typename T, typename X, typename Y>
    void
    dp2d_mv(int rh, T alpha, const GeMatrix<X> &x, T beta, GeMatrix<Y> &y);

} // namespace flens

#include <poisson_solver/kernel_impl/dp_blas.tcc>

#endif // POISSON_SOLVER_KERNEL_IMPL_DP_BLAS_H
