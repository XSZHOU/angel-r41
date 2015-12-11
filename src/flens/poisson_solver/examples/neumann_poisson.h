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

//-- Dirichlet Poisson 1D ------------------------------------------------------

// problem 1: -u''(x) = pi^2*cos(pi*x), u'(0) = 0, u'(1) = 0
// solution:   u(x) = cos(pi*x)
void
problem1(GridVector1D &f, GridVector1D &u, GridVector1D &sol)
{
    double h = 1./f.rh;
    GridVector1D::Grid &F = f.grid;
    GridVector1D::Grid &SOL = sol.grid;

    int i0 = F.firstIndex(),
        i1 = F.lastIndex();

    for (int i=i0; i<=i1; ++i) {
        double x = h*(i+0.5);
        F(i) = cos(M_PI*x)*M_PI*M_PI;
        SOL(i) = cos(M_PI*x);
    }
    np1d_normalize(SOL);
}

// problem 2: -u''(x) = pi^2*cos(pi*x), u'(0) = 0, u'(1) = 0
// solution:   u(x) = cos(pi*x)
void
problem2(GridVector1D &f, GridVector1D &u, GridVector1D &sol)
{
    double h = 1./f.rh;
    GridVector1D::Grid &F = f.grid;
    GridVector1D::Grid &U = u.grid;
    GridVector1D::Grid &SOL = sol.grid;

    int i0 = F.firstIndex(),
        i1 = F.lastIndex();

    for (int i=i0; i<=i1; ++i) {
        double x = h*(i+0.5);
        U(i) = x*x*x; //cos(500*M_PI*x)*500*500*M_PI*M_PI;
        SOL(i) = 0;
    }
}

//== error statistic ===========================================================

void
errorStat(int it, const NeumannPoisson1D &A, const GridVector1D &f,
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
    std::cout.width(17);
    std::cout << rNormL2 << " | ";

    std::cout.precision(12);
    std::cout.setf(std::ios::fixed);
    std::cout.width(17);
    std::cout << rNormInf << " | ";

    std::cout.precision(12);
    std::cout.setf(std::ios::fixed);
    std::cout.width(16);
    std::cout << errorNormL2 << " | ";

    std::cout.precision(12);
    std::cout.setf(std::ios::fixed);
    std::cout.width(16);
    std::cout << errorNormInf << " | " << std::endl;
}

} // namespace flens
