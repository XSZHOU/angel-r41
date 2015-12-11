#include <flens/flens.h>
#include <fftw3.h>

namespace flens {

//== problem set ===============================================================

int N;
GridVector2D u, f, sol, r, error;

void
setProblemSize(int n)
{
    N = n;
    u = f = sol = r = GridVector2D(N+1);
}

// problem 1: -u_xx -u_yy = 5*pi^2*sin(pi*x)*pi^2*sin(2*pi*y), BC: u = 0
// solution:   u(x) = sin(pi*x) * sin(2*pi*y)
void
problem1()
{
    double h = 1./(N+1);
    for (int i=0; i<=N+1; ++i) {
        for (int j=0; j<=N+1; ++j) {
            double x = h*i;
            double y = h*j;
            f.grid(i,j) = 5*M_PI*M_PI * sin(M_PI*x) * sin(2*M_PI*y);
            sol.grid(i,j) = sin(M_PI*x) * sin(2*M_PI*y);
        }
    }
}

// problem 2: -u''(x) = 0, BC: u = 1
// solution:   u(x) = 1
void
problem2()
{
    f.grid = 0;
    u.grid(0,_) = 1;
    u.grid(N+1,_) = 1;
    u.grid(_,0) = 1;
    u.grid(_,N+1) = 1;
    sol.grid = 1;
}

// problem 3: -u_xx -u_yy = 2*pi^2*sin(pi*x) * pi^2*sin(pi*y), BC: u = 0
// solution:   u(x) = sin(pi*x) * sin(pi*y)
void
problem3()
{
    double h = 1./(N+1);
    for (int i=0; i<=N+1; ++i) {
        for (int j=0; j<=N+1; ++j) {
            double x = h*i;
            double y = h*j;
            f.grid(i,j) = 2*M_PI*M_PI * sin(M_PI*x) * sin(M_PI*y);
            sol.grid(i,j) = sin(M_PI*x) * sin(M_PI*y);
        }
    }
}

//== error statistic ===========================================================

void
errorStat(int it=-1)
{
    if (it>=0) {
        std::cout.width(3);
        std::cout << it << ") | ";
    }

    dp2d_residual(N+1, f.grid, u.grid, r.grid);

    std::cout.precision(12);
    std::cout.setf(std::ios::fixed);
    std::cout.width(17);
    std::cout << normL2(r) << " | ";

    std::cout.precision(12);
    std::cout.setf(std::ios::fixed);
    std::cout.width(17);
    std::cout << normInf(r) << " | ";

    error = sol - u;

    std::cout.precision(12);
    std::cout.setf(std::ios::fixed);
    std::cout.width(16);
    std::cout << normL2(error) << " | ";

    std::cout.precision(12);
    std::cout.setf(std::ios::fixed);
    std::cout.width(16);
    std::cout << normInf(error) << " | " << std::endl;
}

//== output routines for gnuplot ===============================================

template <typename U>
void
write(int it, const char *file, const GeMatrix<U> &u)
{
    int N = u.length()-2;
    double h = 1./(N+1);

    std::ostringstream s;
    s << file << std::setw(3) << std::setfill('0') << it << ".dat";

    std::ofstream out(s.str().c_str());
    for (int i=0; i<=N+1; ++i) {
        for (int j=0; j<=N+1; ++j) {
            double x = i*h;
            double y = j*h;
            out << x << " " << y << " " << u(i,j) << std::endl;
        }
    }
}

} // namespace flens
