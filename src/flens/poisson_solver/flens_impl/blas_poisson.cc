#include <poisson_solver/flens_impl/blas_poisson.h>
#include <poisson_solver/kernel_impl.h>

namespace flens {

//-- BLAS for GridVector1D -----------------------------------------------------

void
axpy(double alpha, const GridVector1D &x, GridVector1D &y)
{
    axpy(alpha, x.grid, y.grid);
}

void
copy(const GridVector1D &x, GridVector1D &y)
{
    y = x;
}

double
normInf(const GridVector1D &v)
{
    return dp1d_normInf(v.grid);
}

double
normL2(const GridVector1D &v)
{
    double h = 1./v.rh;
    return sqrt(h*dp1d_norm2sqr(v.grid));
}

//-- BLAS for GridVector2D -----------------------------------------------------

void
copy(const GridVector2D &x, GridVector2D &y)
{
    dp2d_copy(x.grid, y.grid);
    y.rh = x.rh;
}

void
scal(double alpha, GridVector2D &x)
{
    dp2d_scal(alpha, x.grid);
}

double
dot(const GridVector2D &x, const GridVector2D &y)
{
    return x.rh*x.rh*dp2d_dot(x.grid, y.grid);
}

void
axpy(double alpha, const GridVector2D &x, GridVector2D &y)
{
    axpy(alpha, x.grid, y.grid);
}

double
normInf(const GridVector2D &v)
{
    return dp2d_normInf(v.grid);
}

double
normL2(const GridVector2D &v)
{
    double h = 1./v.rh;
    return h*sqrt(dp2d_norm2sqr(v.grid));
}

void
mv(double alpha, const DirichletPoisson2D &A, const GridVector2D &x,
   double beta, GridVector2D &y)
{
    if ((y.grid.numRows()!=x.grid.numRows())
     || (y.grid.numCols()!=x.grid.numCols())) {
        assert(beta==0);
        y.grid.resize(x.grid.numRows(), x.grid.numCols(),
                      x.grid.firstRow(), x.grid.firstCol());
        y.rh = x.rh;
    }
    dp2d_mv(A.rh, alpha, x.grid, beta, y.grid);
}

//-- Residual: r = f - A*u -----------------------------------------------------
//-- dirichlet poisson 1d ------------------------------------------------------

void
residual(const GridVector1D &f,
         const DirichletPoisson1D &A, const GridVector1D &u,
         GridVector1D &r)
{
    dp1d_residual(A.rh, f.grid, u.grid, r.grid);
}

//-- dirichlet poisson 2d ------------------------------------------------------
void
residual(const GridVector2D &f,
         const DirichletPoisson2D &A, const GridVector2D &u,
         GridVector2D &r)
{
    if ((r.grid.numRows()!=u.grid.numRows())
     || (r.grid.numCols()!=u.grid.numCols())) {
        r.grid.resize(u.grid.numRows(), u.grid.numCols(),
                      u.grid.firstRow(), u.grid.firstCol());
        r.rh = u.rh;
    }
    dp2d_residual(A.rh, f.grid, u.grid, r.grid);
}

//-- neumann poisson 1d --------------------------------------------------------

void
residual(const GridVector1D &f,
         const NeumannPoisson1D &A, const GridVector1D &u,
         GridVector1D &r)
{
    dp1d_residual(A.rh, f.grid, u.grid, r.grid);
}

//-- Gauss-Seidel: u = S(A,f)*u_1 ----------------------------------------------

void
mv(Transpose trans, double alpha,
   const GaussSeidel<DirichletPoisson1D, GridVector1D> &GS,
   const GridVector1D &u_1, double beta, GridVector1D &u)
{
    assert(trans==NoTrans);
    assert(alpha==1.);
    assert(beta==0.);
    assert(ADDRESS(u)==ADDRESS(u_1));

    dp1d_gauss_seidel(GS.A.rh, GS.f.grid, u.grid);
}

void
mv(Transpose trans, double alpha,
   const GaussSeidel<NeumannPoisson1D, GridVector1D> &GS,
   const GridVector1D &u_1, double beta, GridVector1D &u)
{
    assert(trans==NoTrans);
    assert(alpha==1.);
    assert(beta==0.);
    assert(ADDRESS(u)==ADDRESS(u_1));

    dp1d_gauss_seidel(GS.A.rh, GS.f.grid, u.grid);
    np1d_normalize(u.grid);
}

//-- Gauss-Seidel Red-Black: u = S(A,f)*u_1 ------------------------------------

void
mv(Transpose trans, double alpha,
   const GaussSeidelRedBlack<DirichletPoisson2D, GridVector2D> &GS,
   const GridVector2D &u_1, double beta, GridVector2D &u)
{
    assert(trans==NoTrans);
    assert(alpha==1.);
    assert(beta==0.);
    assert(ADDRESS(u)==ADDRESS(u_1));

    dp2d_gauss_seidel_red(GS.A.rh, GS.f.grid, u.grid);
    dp2d_gauss_seidel_black(GS.A.rh, GS.f.grid, u.grid);
}

//-- Restriction: vc = R*v -----------------------------------------------------
//-- dirichlet poisson 1d ------------------------------------------------------
void
mv(Transpose trans, double alpha, const Restriction &R,
   const GridVector1D &v, double beta, GridVector1D &vc)
{
    assert(trans==NoTrans);
    assert(alpha==1.);
    assert(beta==0.);

    dp1d_restriction(v.grid, vc.grid);
}

//-- dirichlet poisson 2d ------------------------------------------------------
void
mv(Transpose trans, double alpha, const Restriction &R,
   const GridVector2D &v, double beta, GridVector2D &vc)
{
    assert(trans==NoTrans);
    assert(alpha==1.);
    assert(beta==0.);

    dp2d_restriction_hw(v.grid, vc.grid);
}

//-- neumann poisson 1d --------------------------------------------------------
void
mv(Transpose trans, double alpha, const Restriction_NBC &R,
   const GridVector1D &v, double beta, GridVector1D &vc)
{
    assert(trans==NoTrans);
    assert(alpha==1.);
    assert(beta==0.);

    np1d_restriction(v.grid, vc.grid);
}

//-- Prolongation: v += P*v_c --------------------------------------------------
//-- dirichlet poisson 1d ------------------------------------------------------
void
mv(Transpose trans, double alpha, const Prolongation &P,
   const GridVector1D &vc, double beta, GridVector1D &v)
{
    assert(trans==NoTrans);
    assert(alpha==1.);
    assert(beta==1.);

    dp1d_prolongation(vc.grid, v.grid);
}

void
mv(Transpose trans, double alpha, const Prolongation &P,
   const GridVector2D &vc, double beta, GridVector2D &v)
{
    assert(trans==NoTrans);
    assert(alpha==1.);
    assert(beta==1.);

    dp2d_prolongation(vc.grid, v.grid);
}

//-- neumann poisson 1d --------------------------------------------------------
void
mv(Transpose trans, double alpha, const Prolongation_NBC &P,
   const GridVector1D &vc, double beta, GridVector1D &v)
{
    assert(trans==NoTrans);
    assert(alpha==1.);
    assert(beta==1.);

    np1d_prolongation(vc.grid, v.grid);
}

} // namespace flens
