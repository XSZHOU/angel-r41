#include <poisson_solver/flens_impl/distributedblas_poisson.h>
#include <poisson_solver/flens_impl.h>

namespace flens {

//-- BLAS for DistributedGridVector2D ------------------------------------------

void
scal(double alpha, DistributedGridVector2D &x)
{
    dp2d_scal(alpha, x.grid);
}

double
dot(const DistributedGridVector2D &x, const DistributedGridVector2D &y)
{
    double r = dp2d_dot(x.localGrid(), y.localGrid());
    double R;
    y.mpiCart.comm.Allreduce(&r, &R, 1, MPI::DOUBLE, MPI::SUM);
    return x.rh*x.rh*R;
}

void
axpy(double alpha, const DistributedGridVector2D &x, DistributedGridVector2D &y)
{
    axpy(alpha, x.grid, y.grid);
}

void
copy(const DistributedGridVector2D &x, DistributedGridVector2D &y)
{
    y.grid.resize(x.grid.numRows(), x.grid.numCols(),
                  x.grid.firstRow(), x.grid.firstCol());
    y = x;
}

double
normInf(const DistributedGridVector2D &v)
{
    double r = dp2d_normInf(v.grid);

    double R;
    v.mpiCart.comm.Allreduce(&r, &R, 1, MPI::DOUBLE, MPI::MAX);
    return R;
}

double
normL2(const DistributedGridVector2D &v)
{
    double h = 1./v.rh;
    double r = dp2d_norm2sqr(v.localGrid());

    double R;
    v.mpiCart.comm.Allreduce(&r, &R, 1, MPI::DOUBLE, MPI::SUM);
    return h*sqrt(R);
}

void
mv(double alpha, const DirichletPoisson2D &A, const DistributedGridVector2D &x,
   double beta, DistributedGridVector2D &y)
{
    if ((y.grid.numRows()!=x.grid.numRows())
     || (y.grid.numCols()!=x.grid.numCols())) {
        assert(beta==0);
        y.grid.resize(x.grid.numRows(), x.grid.numCols(),
                      x.grid.firstRow(), x.grid.firstCol());
        y = x;
    }
    dp2d_mv(A.rh, alpha, x.grid, beta, y.grid);
    y.setGhostNodes();
}

//-- Residual: r = f - A*u -----------------------------------------------------
//---- dirichlet poisson 2d ----------------------------------------------------
void
residual(const DistributedGridVector2D &f,
         const DirichletPoisson2D &A, const DistributedGridVector2D &u,
         DistributedGridVector2D &r)
{
    if ((r.grid.numRows()!=f.grid.numRows())
     || (r.grid.numCols()!=f.grid.numCols())) {
         r.grid.resize(f.grid.numRows(), f.grid.numCols(),
                       f.grid.firstRow(), f.grid.firstCol());
         r = f;
    }

    DistributedGridVector2D::LocalGrid R = r.localGrid();
    dp2d_residual(A.rh, f.localGrid(), u.localGrid(), R);
    r.setGhostNodes();
}

//-- Gauss-Seidel Red-Black: u = S(A,f)*u_1 ------------------------------------

void
mv(Transpose trans, double alpha,
   const GaussSeidelRedBlack<DirichletPoisson2D, DistributedGridVector2D> &GS,
   const DistributedGridVector2D &u_1, double beta, DistributedGridVector2D &u)
{
    assert(trans==NoTrans);
    assert(alpha==1.);
    assert(beta==0.);
    assert(ADDRESS(u)==ADDRESS(u_1));

    DistributedGridVector2D::LocalGrid U = u.localGrid();
    dp2d_gauss_seidel_red(GS.A.rh, GS.f.localGrid(), U);
    u.setGhostNodes();
    dp2d_gauss_seidel_black(GS.A.rh, GS.f.localGrid(), U);
    u.setGhostNodes();
}

//-- Restriction: vc = R*v -----------------------------------------------------
//---- dirichlet poisson 2d ----------------------------------------------------
void
mv(Transpose trans, double alpha, const Restriction &R,
   const DistributedGridVector2D &v, double beta, DistributedGridVector2D &vc)
{
    assert(trans==NoTrans);
    assert(alpha==1.);
    assert(beta==0.);

    DistributedGridVector2D::LocalGrid VC = vc.localGrid();
    dp2d_restriction_hw(v.localGrid(), VC);
}

//-- Prolongation: v += P*v_c --------------------------------------------------
//-- dirichlet poisson 1d ------------------------------------------------------
void
mv(Transpose trans, double alpha, const Prolongation &P,
   const DistributedGridVector2D &vc, double beta, DistributedGridVector2D &v)
{
    DistributedGridVector2D::LocalGrid V = v.localGrid();
    dp2d_prolongation(vc.localGrid(), V);
    dp2d_prolongation_north(vc.localGrid(), V);
    dp2d_prolongation_east(vc.localGrid(), V);
    dp2d_prolongation_north_east(vc.localGrid(), V);
    v.setGhostNodes();
}

} // namespace flens
