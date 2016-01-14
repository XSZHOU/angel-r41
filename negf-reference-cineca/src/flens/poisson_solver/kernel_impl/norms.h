#ifndef POISSON_SOLVER_KERNEL_IMPL_NORMS_H
#define POISSON_SOLVER_KERNEL_IMPL_NORMS_H 1

namespace flens {

//-- dirichlet poisson 1d ------------------------------------------------------

template <typename U>
    double
    dp1d_normInf(const DenseVector<U> &u);

template <typename U>
    double
    dp1d_norm2sqr(const DenseVector<U> &u);

//-- dirichlet poisson 2d ------------------------------------------------------

template <typename U>
    double
    dp2d_normInf(const GeMatrix<U> &u);

template <typename U>
    double
    dp2d_norm2sqr(const GeMatrix<U> &u);

} // namespace flens

#include <poisson_solver/kernel_impl/norms.tcc>

#endif // POISSON_SOLVER_KERNEL_IMPL_NORMS_H
