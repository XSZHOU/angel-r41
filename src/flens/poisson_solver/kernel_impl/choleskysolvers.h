#ifndef POISSON_SOLVER_KERNEL_IMPL_CHOLESKYSOLVERS_H
#define POISSON_SOLVER_KERNEL_IMPL_CHOLESKYSOLVERS_H 1

namespace flens {

//-- dirichlet poisson 1d ------------------------------------------------------

template <typename F, typename U>
    void
    dp1d_cholesky(int rh, const DenseVector<F> &f, DenseVector<U> &u);

} // namespace flens

#include <poisson_solver/kernel_impl/choleskysolvers.tcc>

#endif // POISSON_SOLVER_KERNEL_IMPL_CHOLESKYSOLVERS_H
