#include <poisson_solver/flens_impl.h>

#include <neumann_poisson.h>

using namespace flens;
using namespace std;

typedef NeumannPoisson1D   MatType;
typedef GridVector1D       VecType;

typedef GaussSeidel<MatType, VecType>       Smoother;
typedef StationaryIterativeSolver<Smoother> SVR;
typedef Restriction_NBC                     R;
typedef Prolongation_NBC                    P;

typedef MultiGrid<MatType, VecType, R, P, Smoother, SVR>  MG;

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
    const int lMax = 11, lMin = 2;
    const int rh = (1<<lMax);
    const int p = lMax-lMin;
    MatType A[p+1];
    VecType f[p+1], r[p+1], u[p+1], solution(rh,-1);

    allocate(rh, p, A, f, r, u);
    problem1(f[p], u[p], solution);

    SVR svr(A[0], f[0], u[0]);
    MG mg(A, f, r, u, svr);

    errorStat(0, A[p], f[p], u[p], solution);
    for (int it=1; it<=40; ++it) {
        mg.vCycle(p, 1, 1);
        errorStat(it, A[p], f[p], u[p], solution);
    }
}
