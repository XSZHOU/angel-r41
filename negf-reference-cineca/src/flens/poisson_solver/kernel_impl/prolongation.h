#ifndef POISSON_SOLVER_KERNEL_IMPL_PROLONGATION_H
#define POISSON_SOLVER_KERNEL_IMPL_PROLONGATION_H 1

namespace flens {

//-- dirichlet poisson 1d ------------------------------------------------------

template <typename VC, typename V>
    void
    dp1d_prolongation(const DenseVector<VC> &vc, DenseVector<V> &v);

//-- neumann poisson 1d ------------------------------------------------------

template <typename VC, typename V>
    void
    np1d_prolongation(const DenseVector<VC> &vc, DenseVector<V> &v);

//-- dirichlet poisson 2d ------------------------------------------------------

template <typename VC, typename V>
    void
    dp2d_prolongation(const GeMatrix<VC> &vc, GeMatrix<V> &v);

template <typename VC, typename V>
    void
    dp2d_prolongation_north(const GeMatrix<VC> &vc, GeMatrix<V> &v);

template <typename VC, typename V>
    void
    dp2d_prolongation_east(const GeMatrix<VC> &vc, GeMatrix<V> &v);

template <typename VC, typename V>
    void
    dp2d_prolongation_north_east(const GeMatrix<VC> &vc, GeMatrix<V> &v);

//-- neumann poisson 2d ------------------------------------------------------

template <typename VC, typename V>
    void
    np2d_prolongation(const GeMatrix<VC> &vc, GeMatrix<V> &v);


} // namespace flens

#include <poisson_solver/kernel_impl/prolongation.tcc>

#endif // POISSON_SOLVER_KERNEL_IMPL_PROLONGATION_H
