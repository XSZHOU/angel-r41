namespace flens {

//-- BLAS for StaggeredGridVector2D --------------------------------------------

template <bool DX, bool DY>
void
axpy(double alpha, const StaggeredGridVector2D<DX,DY> &x,
     StaggeredGridVector2D<DX,DY> &y)
{
    axpy(alpha, x.grid, y.grid);
}

template <bool DX, bool DY>
void
copy(const StaggeredGridVector2D<DX,DY> &x, StaggeredGridVector2D<DX,DY> &y)
{
    y = x;
}

template <bool DX, bool DY>
double
normInf(const StaggeredGridVector2D<DX,DY> &v)
{
    return dp2d_normInf(v.grid);
}

template <bool DX, bool DY>
double
normL2(const StaggeredGridVector2D<DX,DY> &v)
{
    double h = 1./v.rh;
    return h*sqrt(dp2d_norm2sqr(v.grid));
}

//-- Residual: r = f - A*u -----------------------------------------------------

template <bool DX, bool DY>
void
residual(const StaggeredGridVector2D<DX,DY> &f,
         const DirichletPoisson2D &A,
         const StaggeredGridVector2D<DX,DY> &u,
         StaggeredGridVector2D<DX,DY> &r)
{
    dp2d_residual(A.rh, f.grid, u.grid, r.grid);
}

template <bool DX, bool DY>
void
residual(const StaggeredGridVector2D<DX,DY> &f,
         const NeumannPoisson2D &A,
         const StaggeredGridVector2D<DX,DY> &u,
         StaggeredGridVector2D<DX,DY> &r)
{
    dp2d_residual(A.rh, f.grid, u.grid, r.grid);
}

//-- Gauss-Seidel Red-Black: u = S(A,f)*u_1  -----------------------------------

template <bool DX, bool DY>
void
mv(Transpose trans, double alpha,
   const GaussSeidelRedBlack<DirichletPoisson2D,
                             StaggeredGridVector2D<DX,DY> > &GS,
   const StaggeredGridVector2D<DX,DY> &u_1, double beta,
   StaggeredGridVector2D<DX,DY> &u)
{
    assert(trans==NoTrans);
    assert(alpha==1.);
    assert(beta==0.);
    assert(ADDRESS(u)==ADDRESS(u_1));

    p2d_gauss_seidel_red(GS.A.rh, GS.f.grid, u.grid);
    p2d_gauss_seidel_black(GS.A.rh, GS.f.grid, u.grid);
}

template <bool DX, bool DY>
void
mv(Transpose trans, double alpha,
   const GaussSeidelRedBlack<NeumannPoisson2D,
                             StaggeredGridVector2D<DX,DY> > &GS,
   const StaggeredGridVector2D<DX,DY> &u_1, double beta,
   StaggeredGridVector2D<DX,DY> &u)
{
    assert(trans==NoTrans);
    assert(alpha==1.);
    assert(beta==0.);
    assert(ADDRESS(u)==ADDRESS(u_1));

    p2d_gauss_seidel_red(GS.A.rh, GS.f.grid, u.grid);
    p2d_gauss_seidel_black(GS.A.rh, GS.f.grid, u.grid);
    np2d_normalize(u.grid);
}

//-- neumann poisson 2d --------------------------------------------------------
template <bool DX, bool DY>
void
mv(Transpose trans, double alpha, const Restriction &R,
   const StaggeredGridVector2D<DX,DY> &v, double beta,
   StaggeredGridVector2D<DX,DY> &vc)
{
    np2d_restriction_hw(v.grid, vc.grid);
}

//-- neumann poisson 2d --------------------------------------------------------
template <bool DX, bool DY>
void
mv(Transpose trans, double alpha, const Prolongation &R,
   const StaggeredGridVector2D<DX,DY> &vc, double beta,
   StaggeredGridVector2D<DX,DY> &v)
{
    np2d_prolongation(vc.grid, v.grid);
}

} // namespace flens
