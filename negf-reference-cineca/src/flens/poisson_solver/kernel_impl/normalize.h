#ifndef POISSON_SOLVER_KERNEL_IMPL_NORMALIZE_H
#define POISSON_SOLVER_KERNEL_IMPL_NORMALIZE_H 1

namespace flens {

//-- neumann poisson 1d --------------------------------------------------------

template <typename U>
    void
    np1d_normalize(DenseVector<U> &u);

//-- neumann poisson 2d --------------------------------------------------------

template <typename U>
    void
    np2d_bc(GeMatrix<U> &u);

template <typename U>
    void
    np2d_normalize(GeMatrix<U> &u);

} // namespace flens

#include <poisson_solver/kernel_impl/normalize.tcc>

#endif // POISSON_SOLVER_KERNEL_IMPL_NORMALIZE_H
