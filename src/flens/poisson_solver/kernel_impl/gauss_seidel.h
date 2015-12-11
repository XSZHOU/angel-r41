#ifndef POISSON_SOLVER_KERNEL_IMPL_GAUSS_SEIDEL_H
#define POISSON_SOLVER_KERNEL_IMPL_GAUSS_SEIDEL_H 1

namespace flens {

//-- dirichlet poisson 1d ------------------------------------------------------

template <typename F, typename U>
    void
    dp1d_gauss_seidel(int rh, const DenseVector<F> &f, DenseVector<U> &u);

//-- dirichlet poisson 2d ------------------------------------------------------

template <typename F, typename U>
    void
    dp2d_gauss_seidel_red(int rh, const GeMatrix<F> &f, GeMatrix<U> &u);


template <typename F, typename U>
    void
    dp2d_gauss_seidel_black(int rh, const GeMatrix<F> &f, GeMatrix<U> &u);

//-- poisson 2d ----------------------------------------------------------------

template <typename F, typename U>
    void
    p2d_gauss_seidel_red(int rh, const GeMatrix<F> &f, GeMatrix<U> &u);

template <typename F, typename U>
    void
    p2d_gauss_seidel_black(int rh, const GeMatrix<F> &f, GeMatrix<U> &u);

} // namespace flens

#include <poisson_solver/kernel_impl/gauss_seidel.tcc>

#endif // POISSON_SOLVER_KERNEL_IMPL_GAUSS_SEIDEL_H
