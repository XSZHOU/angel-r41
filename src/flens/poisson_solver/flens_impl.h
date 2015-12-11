#ifndef POISSON_SOLVER_FLENS_IMPL_H
#define POISSON_SOLVER_FLENS_IMPL_H 1

#include <poisson_solver/kernel_impl.h>

#include <poisson_solver/flens_impl/blas_poisson.h>
#include <poisson_solver/flens_impl/directsolvers.h>
#include <poisson_solver/flens_impl/gauss_seidel.h>
#include <poisson_solver/flens_impl/gridvector.h>
#include <poisson_solver/flens_impl/prolongation.h>
#include <poisson_solver/flens_impl/poissonmatrix.h>
#include <poisson_solver/flens_impl/restriction.h>

#ifdef USE_MPI
#   include <poisson_solver/flens_impl/distributedblas_poisson.h>
#   include <poisson_solver/flens_impl/distributedgridvector.h>
#   include <poisson_solver/flens_impl/mpicart.h>
#endif // USE_MPI

#endif // POISSON_SOLVER_FLENS_IMPL_H
