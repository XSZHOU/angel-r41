#include <flens/flens.h>
#include <poisson_solver/flens_impl.h>
#include <fftw3.h>

namespace flens {

//-- optimization --------------------------------------------------------------

template <typename VB, typename MA, typename VX, typename VR>
void
residual(const Vector<VB> &b, const Matrix<MA> &A, const Vector<VX> &x,
         Vector<VR> &r)
{
    typedef typename VR::ElementType T;

    copy(b.impl(), r.impl());
    mv(NoTrans, T(-1), A.impl(), x.impl(), T(1), r.impl());
}

// r = b -A*x
template <typename VB, typename MA, typename VX, typename VR>
void
copy(const VectorClosure<OpSub, VB, VectorClosure<OpMult, MA, VX> > &b_Ax,
     Vector<VR> &r)
{
    residual(b_Ax.left(), b_Ax.right().left(), b_Ax.right().right(), r.impl());
}

//== problem set ===============================================================

//-- Dirichlet Poisson 1D (low level) ------------------------------------------

// problem 1: -u''(x) = pi^2*sin(pi*x), u(0) = 0, u(1) = 0
// solution:   u(x) = sin(pi*x)
template <typename F, typename U, typename SOL>
void
problem1(int rh, DenseVector<F> &f, DenseVector<U> &u, DenseVector<SOL> &sol)
{
    int N = rh-1;
    double h = 1./rh;

    for (int i=1; i<=N; ++i) {
        double x = h*i;
        f(i) = sin(M_PI*x)*M_PI*M_PI;
        sol(i) = sin(M_PI*x);
    }
}

//-- Dirichlet Poisson 1D ------------------------------------------------------

// problem 1: -u''(x) = pi^2*sin(pi*x), u(0) = 0, u(1) = 0
// solution:   u(x) = sin(pi*x)
void
problem1(GridVector1D &f, GridVector1D &u, GridVector1D &sol)
{
    double h = 1./f.rh;
    GridVector1D::Grid &F = f.grid;
    GridVector1D::Grid &SOL = sol.grid;

    int i0 = F.firstIndex(),
        i1 = F.lastIndex();

    for (int i=i0; i<=i1; ++i) {
        double x = h*i;
        F(i) = sin(M_PI*x)*M_PI*M_PI;
        SOL(i) = sin(M_PI*x);
    }
}

// problem 2: -u''(x) = 0, u(0) = 1, u(1) = 1
// solution:   u(x) = 1
void
problem2(GridVector1D &f, GridVector1D &u, GridVector1D &sol)
{
    GridVector1D::Grid &F = f.grid;
    GridVector1D::Grid &SOL = sol.grid;
    GridVector1D::Grid &U = u.grid;

    int i0 = F.firstIndex(),
        i1 = F.lastIndex();

    F = 0;
    SOL = 1;
    U(i0) = U(i1) = 1;
}

// problem 3: -u''(x) = 0, u(0) = 0, u(1) = 1
// solution:   u(x) = x
void
problem3(GridVector1D &f, GridVector1D &u, GridVector1D &sol)
{
    double h = 1./f.rh;
    GridVector1D::Grid &F = f.grid;
    GridVector1D::Grid &SOL = sol.grid;
    GridVector1D::Grid &U = u.grid;

    int i0 = U.firstIndex(),
        i1 = U.lastIndex();

    U(i0) = U(i1) = 1;
    F = 0;
    U(i0) = 0;
    U(i1) = 1;
    for (int i=i0; i<=i1; ++i) {
        double x = h*i;
        SOL(i) = x;
    }
}

// problem 4: -u''(x) = 0, u(0) = 0, u(1) = 0
// solution:   u(x) = 0
void
problem4(GridVector1D &f, GridVector1D &u, GridVector1D &sol)
{
    GridVector1D::Grid &F = f.grid;
    GridVector1D::Grid &SOL = sol.grid;
    GridVector1D::Grid &U = u.grid;

    int i0 = U.firstIndex(),
        i1 = U.lastIndex();

    F = 0;
    U(i0) = 0;
    U(i1) = 0;
    SOL = 0;
}

//-- Dirichlet Poisson 2D ------------------------------------------------------

// problem 1: -u_xx -u_yy = 5*pi^2*sin(pi*x)*pi^2*sin(2*pi*y), BC: u = 0
// solution:   u(x) = sin(pi*x) * sin(2*pi*y)
void
problem1(GridVector2D &f, GridVector2D &u, GridVector2D &sol)
{
    double h = 1./f.rh;
    GridVector2D::Grid &F = f.grid;
    GridVector2D::Grid &SOL = sol.grid;
    for (int i=F.firstRow(); i<=F.lastRow(); ++i) {
        for (int j=F.firstCol(); j<=F.lastCol(); ++j) {
            double x = h*i;
            double y = h*j;
            F(i,j) = 5*M_PI*M_PI * sin(M_PI*x) * sin(2*M_PI*y);
            SOL(i,j) = sin(M_PI*x) * sin(2*M_PI*y);
        }
    }
}

#ifdef USE_MPI
void
problem1(DistributedGridVector2D &f, DistributedGridVector2D &u,
         DistributedGridVector2D &sol)
{
    double h = 1./f.rh;
    DistributedGridVector2D::LocalGrid F = f.localGrid();
    DistributedGridVector2D::LocalGrid SOL = sol.localGrid();

    for (int i=F.firstRow(); i<=F.lastRow(); ++i) {
        for (int j=F.firstCol(); j<=F.lastCol(); ++j) {
            double x = (i+f.i0)*h;
            double y = (j+f.j0)*h;
            F(i,j) = 5*M_PI*M_PI * sin(M_PI*x) * sin(2*M_PI*y);
            SOL(i,j) = sin(M_PI*x) * sin(2*M_PI*y);
        }
    }
}
#endif // USE_MPI

// problem 2: -u''(x) = 0, BC: u = 1
// solution:   u(x) = 1
void
problem2(GridVector2D &f, GridVector2D &u, GridVector2D &sol)
{
    int N = f.rh-1;

    f.grid = 0;
    u.grid(0,_) = 1;
    u.grid(N+1,_) = 1;
    u.grid(_,0) = 1;
    u.grid(_,N+1) = 1;
    sol.grid = 1;
}

void
problem2(StaggeredGridVector2D<false, false> &f,
         StaggeredGridVector2D<false, false> &u,
         StaggeredGridVector2D<false, false> &sol)
{
    int N = f.rh-1;

    f.grid = 0;
    u.grid(0,_) = 1;
    u.grid(N+1,_) = 1;
    u.grid(_,0) = 1;
    u.grid(_,N+1) = 1;
    sol.grid = 1;
}

// problem 2: -u''(x) = 2x + 2y - 2,  BC: u' = 0
// solution:   u(x) = x^2/2 - x^3/3 + y^2/2 - y^3/3 - 1/6
void
problem2(StaggeredGridVector2D<true, true> &f,
         StaggeredGridVector2D<true, true> &u,
         StaggeredGridVector2D<true, true> &sol)
{
    double h = 1./f.rh;
    double dx = 0.5;
    double dy = 0.5;

    StaggeredGridVector2D<true, true>::Grid &F = f.grid;
    StaggeredGridVector2D<true, true>::Grid &SOL = sol.grid;
    for (int i=F.firstRow(); i<=F.lastRow(); ++i) {
        for (int j=F.firstCol(); j<=F.lastCol(); ++j) {
            double x = (i+dx)*h;
            double y = (j+dy)*h;
            F(i,j) = 2*x + 2*y - 2;
            SOL(i,j) = x*x/2 - x*x*x/3 + y*y/2 - y*y*y/3 - 1./6;
        }
    }
}


#ifdef USE_MPI
void
problem2(DistributedGridVector2D &f, DistributedGridVector2D &u,
         DistributedGridVector2D &sol)
{
    f.grid = 0;
    sol.grid = 1;

    MpiCart mpiCart = f.mpiCart;
    int m = f.m,
        n = f.n;

    if (mpiCart.row==0) {
        u.grid(0,_) = 1;
    }
    if (mpiCart.row==mpiCart.numRows-1) {
        u.grid(m+1,_) = 1;
    }
    if (mpiCart.col==0) {
        u.grid(_,0) = 1;
    }
    if (mpiCart.col==mpiCart.numCols-1) {
        u.grid(_,n+1) = 1;
    }
}
#endif // USE_MPI

// problem 3: -u_xx -u_yy = 2*pi^2*sin(pi*x) * pi^2*sin(pi*y), BC: u = 0
// solution:   u(x) = sin(pi*x) * sin(pi*y)
void
problem3(GridVector2D &f, GridVector2D &u, GridVector2D &sol)
{
    double h = 1./f.rh;
    GridVector2D::Grid &F = f.grid;
    GridVector2D::Grid &SOL = sol.grid;
    for (int i=F.firstRow(); i<=F.lastRow(); ++i) {
        for (int j=F.firstCol(); j<=F.lastCol(); ++j) {
            double x = h*i;
            double y = h*j;
            F(i,j) = 2*M_PI*M_PI * sin(M_PI*x) * sin(M_PI*y);
            SOL(i,j) = sin(M_PI*x) * sin(M_PI*y);
        }
    }
}

#ifdef USE_MPI
void
problem3(DistributedGridVector2D &f, DistributedGridVector2D &u,
         DistributedGridVector2D &sol)
{
    double h = 1./f.rh;
    DistributedGridVector2D::LocalGrid F = f.localGrid();
    DistributedGridVector2D::LocalGrid SOL = sol.localGrid();

    for (int i=F.firstRow(); i<=F.lastRow(); ++i) {
        for (int j=F.firstCol(); j<=F.lastCol(); ++j) {
            double x = (i+f.i0)*h;
            double y = (j+f.j0)*h;
            F(i,j) = 2*M_PI*M_PI * sin(M_PI*x) * sin(M_PI*y);
            SOL(i,j) = sin(M_PI*x) * sin(M_PI*y);
        }
    }
}
#endif // USE_MPI

//== error statistic ===========================================================

void
errorStat(int it, const DirichletPoisson1D &A, const GridVector1D &f,
          const GridVector1D &u, const GridVector1D &solution)
{
    if (it>=0) {
        std::cout.width(3);
        std::cout << it << ") | ";
    }

    GridVector1D r(f.rh), error(f.rh);
    r = f - A*u;
    error = solution - u;

    double rNormInf = normInf(r);
    double rNormL2 = normL2(r);
    double errorNormInf = normInf(error);
    double errorNormL2 = normL2(error);

    std::cout.precision(12);
    std::cout.setf(std::ios::fixed);
    std::cout.width(20);
    std::cout << rNormL2 << " | ";

    std::cout.precision(12);
    std::cout.setf(std::ios::fixed);
    std::cout.width(20);
    std::cout << rNormInf << " | ";

    std::cout.precision(12);
    std::cout.setf(std::ios::fixed);
    std::cout.width(20);
    std::cout << errorNormL2 << " | ";

    std::cout.precision(12);
    std::cout.setf(std::ios::fixed);
    std::cout.width(20);
    std::cout << errorNormInf << " | " << std::endl;
}

void
errorStat(int it, const DirichletPoisson2D &A, const GridVector2D &f,
          const GridVector2D &u, const GridVector2D &solution)
{
    if (it>=0) {
        std::cout.width(3);
        std::cout << it << ") | ";
    }

    GridVector2D r(f.rh), error(f.rh);
    r = f - A*u;
    error = solution - u;

    double rNormInf = normInf(r);
    double rNormL2 = normL2(r);
    double errorNormInf = normInf(error);
    double errorNormL2 = normL2(error);

    std::cout.precision(12);
    std::cout.setf(std::ios::fixed);
    std::cout.width(20);
    std::cout << rNormL2 << " | ";

    std::cout.precision(12);
    std::cout.setf(std::ios::fixed);
    std::cout.width(20);
    std::cout << rNormInf << " | ";

    std::cout.precision(12);
    std::cout.setf(std::ios::fixed);
    std::cout.width(18);
    std::cout << errorNormL2 << " | ";

    std::cout.precision(12);
    std::cout.setf(std::ios::fixed);
    std::cout.width(18);
    std::cout << errorNormInf << " | " << std::endl;
}

template <typename MatType, bool DirectionX, bool DirectionY>
void
errorStat(int it, const MatType &A,
          const StaggeredGridVector2D<DirectionX, DirectionY> &f,
          const StaggeredGridVector2D<DirectionX, DirectionY> &u,
          const StaggeredGridVector2D<DirectionX, DirectionY> &solution)
{
    if (it>=0) {
        std::cout.width(3);
        std::cout << it << ") | ";
    }

    StaggeredGridVector2D<DirectionX, DirectionY> r(f.rh), error(f.rh);
    r = f - A*u;
    error = solution - u;

    double rNormInf = normInf(r);
    double rNormL2 = normL2(r);
    double errorNormInf = normInf(error);
    double errorNormL2 = normL2(error);

    std::cout.precision(12);
    std::cout.setf(std::ios::fixed);
    std::cout.width(20);
    std::cout << rNormL2 << " | ";

    std::cout.precision(12);
    std::cout.setf(std::ios::fixed);
    std::cout.width(20);
    std::cout << rNormInf << " | ";

    std::cout.precision(12);
    std::cout.setf(std::ios::fixed);
    std::cout.width(20);
    std::cout << errorNormL2 << " | ";

    std::cout.precision(12);
    std::cout.setf(std::ios::fixed);
    std::cout.width(20);
    std::cout << errorNormInf << " | " << std::endl;
}

#ifdef USE_MPI
void
errorStat(int it, const DirichletPoisson2D &A,
          const DistributedGridVector2D &f, const DistributedGridVector2D &u,
          const DistributedGridVector2D &solution)
{
    MpiCart mpiCart = f.mpiCart;
    DistributedGridVector2D r(mpiCart, f.rh), error(mpiCart, f.rh);
    r = f - A*u;
    error = solution - u;

    double rNormInf = normInf(r);
    double rNormL2 = normL2(r);
    double errorNormInf = normInf(error);
    double errorNormL2 = normL2(error);

    if ((mpiCart.row==0) && (mpiCart.col==0)) {
        if (it>=0) {
            std::cout.width(3);
            std::cout << it << ") | ";
        }

        std::cout.precision(12);
        std::cout.setf(std::ios::fixed);
        std::cout.width(20);
        std::cout << rNormL2 << " | ";

        std::cout.precision(12);
        std::cout.setf(std::ios::fixed);
        std::cout.width(20);
        std::cout << rNormInf << " | ";

        std::cout.precision(12);
        std::cout.setf(std::ios::fixed);
        std::cout.width(20);
        std::cout << errorNormL2 << " | ";

        std::cout.precision(12);
        std::cout.setf(std::ios::fixed);
        std::cout.width(20);
        std::cout << errorNormInf << " | " << std::endl;
    }
}
#endif // USE_MPI

} // namespace flens
