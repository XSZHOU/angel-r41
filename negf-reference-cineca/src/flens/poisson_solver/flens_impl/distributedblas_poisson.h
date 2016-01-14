#ifndef POISSON_SOLVER_FLENS_IMPL_DISTRIBUTEDBLAS_POISSON_H
#define POISSON_SOLVER_FLENS_IMPL_DISTRIBUTEDBLAS_POISSON_H 1

#include <poisson_solver/flens_impl/gauss_seidel.h>
#include <poisson_solver/flens_impl/distributedgridvector.h>
#include <poisson_solver/flens_impl/poissonmatrix.h>
#include <poisson_solver/flens_impl/prolongation.h>
#include <poisson_solver/flens_impl/restriction.h>

namespace flens {

//-- BLAS for DistributedGridVector2D ------------------------------------------

void
scal(double alpha, DistributedGridVector2D &x);

double
dot(const DistributedGridVector2D &x, const DistributedGridVector2D &y);

void
axpy(double alpha, const DistributedGridVector2D &x,
     DistributedGridVector2D &y);

void
copy(const DistributedGridVector2D &x, DistributedGridVector2D &y);

double
normInf(const DistributedGridVector2D &v);

double
normL2(const DistributedGridVector2D &v);

void
mv(double alpha, const DirichletPoisson2D &A, const DistributedGridVector2D &x,
   double beta, DistributedGridVector2D &y);

//-- Residual: r = f - A*u -----------------------------------------------------

void
residual(const DistributedGridVector2D &f,
         const DirichletPoisson2D &A, const DistributedGridVector2D &u,
         DistributedGridVector2D &r);

//-- Gauss-Seidel Red-Black: u = S(A,f)*u_1  -----------------------------------

void
mv(Transpose trans, double alpha,
   const GaussSeidelRedBlack<DirichletPoisson2D, DistributedGridVector2D> &GS,
   const DistributedGridVector2D &u_1, double beta, DistributedGridVector2D &u);

//-- Restriction: vc = R*v  ----------------------------------------------------
//-- dirichlet poisson 1d ------------------------------------------------------

void
mv(Transpose trans, double alpha, const Restriction &R,
   const DistributedGridVector2D &v, double beta, DistributedGridVector2D &vc);

//-- Prolongation: v += P*v_c  -------------------------------------------------
//-- dirichlet poisson 1d ------------------------------------------------------

void
mv(Transpose trans, double alpha, const Prolongation &P,
   const DistributedGridVector2D &vc, double beta, DistributedGridVector2D &v);

} // namespace flens

#endif // POISSON_SOLVER_FLENS_IMPL_DISTRIBUTEDBLAS_POISSON_H
