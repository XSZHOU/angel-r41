#include <poisson_solver/flens_impl.h>

#include <dirichlet_poisson.h>

using namespace flens;
using namespace std;

typedef DirichletPoisson2D MatType;
typedef GridVector2D       VecType;

typedef GaussSeidelRedBlack<MatType, VecType>  Smoother;
typedef FastPoissonSolver<MatType, VecType>    DS;

typedef MultiGrid<MatType, VecType, Restriction, Prolongation, Smoother, DS>  MG;

void
allocate(int rh, int p, DirichletPoisson2D *A,
         GridVector2D *f, GridVector2D *r, GridVector2D *u)
{
    for (int l=p; l>=0; --l, rh/=2) {
        A[l] = DirichletPoisson2D(rh);
        f[l] = GridVector2D(rh);
        r[l] = GridVector2D(rh);
        u[l] = GridVector2D(rh);
    }
}

int
main()
{
    int l = 10, lMin = 2;
    int rh = (1<<l);
    int p = l-lMin;
    MatType A[p+1];
    VecType f[p+1], r[p+1], u[p+1], solution(rh);

    allocate(rh, p, A, f, r, u);
    problem2(f[p], u[p], solution);

    DS ds(A[0], f[0], u[0]);
    MG mg(A, f, r, u, ds);

    errorStat(0, A[p], f[p], u[p], solution);
    for (int it=1; it<=20; ++it) {
        mg.vCycle(p, 1, 1);
        errorStat(it, A[p], f[p], u[p], solution);
    }
}
