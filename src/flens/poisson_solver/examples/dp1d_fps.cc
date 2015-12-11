#include <poisson_solver/flens_impl.h>

#include <dirichlet_poisson.h>

using namespace flens;
using namespace std;

typedef DirichletPoisson1D MatType;
typedef GridVector1D       VecType;

typedef GaussSeidel<MatType, VecType>       Smoother;
typedef FastPoissonSolver<MatType, VecType> DS;

typedef MultiGrid<MatType, VecType, Restriction, Prolongation, Smoother, DS>  MG;

int
main()
{
    int l = 10;
    int rh = (1<<l);

    MatType A(rh);
    VecType f(rh), r(rh), u(rh), solution(rh);

    problem3(f, u, solution);

    DS ds(A, f, u);

    errorStat(0, A, f, u, solution);
    ds.solve();
    errorStat(1, A, f, u, solution);
}
