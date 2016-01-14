#ifndef POISSON_SOLVER_FLENS_IMPL_BLAS_POISSON_H
#define POISSON_SOLVER_FLENS_IMPL_BLAS_POISSON_H 1

#include <poisson_solver/flens_impl/gauss_seidel.h>
#include <poisson_solver/flens_impl/gridvector.h>
#include <poisson_solver/flens_impl/poissonmatrix.h>
#include <poisson_solver/flens_impl/prolongation.h>
#include <poisson_solver/flens_impl/restriction.h>

namespace flens {

//-- BLAS for GridVector1D -----------------------------------------------------

void
axpy(double alpha, const GridVector1D &x, GridVector1D &y);

void
copy(const GridVector1D &x, GridVector1D &y);

double
normInf(const GridVector1D &v);

double
normL2(const GridVector1D &v);

//-- BLAS for GridVector2D -----------------------------------------------------

void
copy(const GridVector2D &x, GridVector2D &y);

void
scal(double alpha, GridVector2D &x);

double
dot(const GridVector2D &x, const GridVector2D &y);

void
axpy(double alpha, const GridVector2D &x, GridVector2D &y);

double
normInf(const GridVector2D &v);

double
normL2(const GridVector2D &v);

void
mv(double alpha, const DirichletPoisson2D &A, const GridVector2D &x,
   double beta, GridVector2D &y);

//-- BLAS for StaggeredGridVector2D --------------------------------------------

template <bool DX, bool DY>
void
axpy(double alpha, const StaggeredGridVector2D<DX,DY> &x,
     StaggeredGridVector2D<DX,DY> &y);

template <bool DX, bool DY>
void
copy(const StaggeredGridVector2D<DX,DY> &x, StaggeredGridVector2D<DX,DY> &y);

template <bool DX, bool DY>
double
normInf(const StaggeredGridVector2D<DX,DY> &v);

template <bool DX, bool DY>
double
normL2(const StaggeredGridVector2D<DX,DY> &v);

//-- Residual: r = f - A*u -----------------------------------------------------

void
residual(const GridVector1D &f,
         const DirichletPoisson1D &A, const GridVector1D &u,
         GridVector1D &r);

void
residual(const GridVector1D &f,
         const NeumannPoisson1D &A, const GridVector1D &u,
         GridVector1D &r);

void
residual(const GridVector2D &f,
         const DirichletPoisson2D &A, const GridVector2D &u,
         GridVector2D &r);

template <bool DX, bool DY>
void
residual(const StaggeredGridVector2D<DX,DY> &f,
         const DirichletPoisson2D &A,
         const StaggeredGridVector2D<DX,DY> &u,
         StaggeredGridVector2D<DX,DY> &r);

template <bool DX, bool DY>
void
residual(const StaggeredGridVector2D<DX,DY> &f,
         const NeumannPoisson2D &A,
         const StaggeredGridVector2D<DX,DY> &u,
         StaggeredGridVector2D<DX,DY> &r);

//-- Gauss-Seidel: u = S(A,f)*u_1  ---------------------------------------------

void
mv(Transpose trans, double alpha,
   const GaussSeidel<DirichletPoisson1D, GridVector1D> &GS,
   const GridVector1D &u_1, double beta, GridVector1D &u);

void
mv(Transpose trans, double alpha,
   const GaussSeidel<NeumannPoisson1D, GridVector1D> &GS,
   const GridVector1D &u_1, double beta, GridVector1D &u);

//-- Gauss-Seidel Red-Black: u = S(A,f)*u_1  -----------------------------------

void
mv(Transpose trans, double alpha,
   const GaussSeidelRedBlack<DirichletPoisson2D, GridVector2D> &GS,
   const GridVector2D &u_1, double beta, GridVector2D &u);

template <bool DX, bool DY>
void
mv(Transpose trans, double alpha,
   const GaussSeidelRedBlack<DirichletPoisson2D,
                             StaggeredGridVector2D<DX,DY> > &GS,
   const StaggeredGridVector2D<DX,DY> &u_1, double beta,
   StaggeredGridVector2D<DX,DY> &u);

template <bool DX, bool DY>
void
mv(Transpose trans, double alpha,
   const GaussSeidelRedBlack<NeumannPoisson2D,
                             StaggeredGridVector2D<DX,DY> > &GS,
   const StaggeredGridVector2D<DX,DY> &u_1, double beta,
   StaggeredGridVector2D<DX,DY> &u);

//-- Restriction: vc = R*v  ----------------------------------------------------
//---- dirichlet poisson 1d ----------------------------------------------------
void
mv(Transpose trans, double alpha, const Restriction &R,
   const GridVector1D &v, double beta, GridVector1D &vc);

//---- dirichlet poisson 2d ----------------------------------------------------
void
mv(Transpose trans, double alpha, const Restriction &R,
   const GridVector2D &v, double beta, GridVector2D &vc);

//-- neumann poisson 1d --------------------------------------------------------
void
mv(Transpose trans, double alpha, const Restriction_NBC &R,
   const GridVector1D &v, double beta, GridVector1D &vc);

//-- neumann poisson 2d --------------------------------------------------------
template <bool DX, bool DY>
void
mv(Transpose trans, double alpha, const Restriction &R,
   const StaggeredGridVector2D<DX,DY> &v, double beta,
   StaggeredGridVector2D<DX,DY> &vc);

//-- Prolongation: v += P*v_c  -------------------------------------------------
//-- dirichlet poisson 1d ------------------------------------------------------
void
mv(Transpose trans, double alpha, const Prolongation &P,
   const GridVector1D &vc, double beta, GridVector1D &v);

void
mv(Transpose trans, double alpha, const Prolongation &P,
   const GridVector2D &vc, double beta, GridVector2D &v);

//-- neumann poisson 1d --------------------------------------------------------
void
mv(Transpose trans, double alpha, const Prolongation_NBC &P,
   const GridVector1D &vc, double beta, GridVector1D &v);

//-- neumann poisson 2d --------------------------------------------------------
template <bool DX, bool DY>
void
mv(Transpose trans, double alpha, const Prolongation &R,
   const StaggeredGridVector2D<DX,DY> &v, double beta,
   StaggeredGridVector2D<DX,DY> &vc);

} // namespace flens

#include <poisson_solver/flens_impl/blas_poisson.tcc>

#endif // POISSON_SOLVER_FLENS_IMPL_BLAS_POISSON_H
