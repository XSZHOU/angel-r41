#include <poisson_solver/flens_impl.h>

#include <neumann_poisson.h>

using namespace flens;
using namespace std;

typedef NeumannPoisson1D   MatType;
typedef GridVector1D       VecType;

typedef GaussSeidel<MatType, VecType>       Smoother;
typedef StationaryIterativeSolver<Smoother> SVR;
typedef Restriction_NBC                     Res;
typedef Prolongation_NBC                    Pro;

typedef MultiGrid<MatType, VecType, Res, Pro, Smoother, SVR>  MG;

void
allocate(int rh, int p, NeumannPoisson1D *A,
         GridVector1D *f, GridVector1D *r, GridVector1D *u)
{
    for (int l=p; l>=0; --l, rh/=2) {
        A[l] = NeumannPoisson1D(rh);
        f[l] = r[l] = u[l] = GridVector1D(rh, -1);
    }
}

int
main()
{
    const int lMax = 3, lMin = 2;
    const int rh = (1<<lMax);
    const int p = lMax-lMin;
    MatType A[p+1];
    VecType f[p+1], r[p+1], u[p+1], solution(rh,-1);

    allocate(rh, p, A, f, r, u);
    //problem1(f[p], u[p], solution);

    for (int i=0; i<rh; ++i) {
        double x = i/double(rh);
        u[p].grid(i) = x;
    }
    cerr << "u[p] = " << u[p] << endl;

    Res R;
    Pro P;

    u[p-1] = R*u[p];
    cerr << "u[p-1] = " << u[p-1] << endl;

    u[p] = 0;
    u[p] += P*u[p-1];
    cerr << "u[p] = " << u[p] << endl;
}
