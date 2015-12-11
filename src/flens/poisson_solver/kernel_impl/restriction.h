#ifndef POISSON_SOLVER_KERNEL_IMPL_RESTRICTION_H
#define POISSON_SOLVER_KERNEL_IMPL_RESTRICTION_H 1

namespace flens {

//-- dirichlet poisson 1d ------------------------------------------------------

template <typename V, typename VC>
    void
    dp1d_restriction(const DenseVector<V> &v, DenseVector<VC> &vc);

//-- neumann poisson 1d --------------------------------------------------------

template <typename V, typename VC>
    void
    np1d_restriction(const DenseVector<V> &v, DenseVector<VC> &vc);

//-- dirichlet poisson 2d ------------------------------------------------------

template <typename V, typename VC>
    void
    dp2d_restriction_hw(const GeMatrix<V> &v, GeMatrix<VC> &vc);

//-- neumann poisson 2d --------------------------------------------------------

template <typename V, typename VC>
    void
    np2d_restriction_hw(const GeMatrix<V> &v, GeMatrix<VC> &vc);

} // namespace flens

#include <poisson_solver/kernel_impl/restriction.tcc>

#endif // POISSON_SOLVER_KERNEL_IMPL_RESTRICTION_H
