#include <poisson_solver/flens_impl.h>
#include <dirichlet_poisson.h>
#include <timer.h>

using namespace flens;
using namespace std;

typedef DirichletPoisson2D       MatType;
typedef DistributedGridVector2D  VecType;

typedef GaussSeidelRedBlack<MatType, VecType>  Smoother;
typedef StationaryIterativeSolver<Smoother>    DS;

typedef MultiGrid<MatType, VecType, Restriction, Prolongation, Smoother, DS>  MG;

void
allocate(MpiCart mpiCart, int rh, int p, DirichletPoisson2D *A,
         DistributedGridVector2D *f, DistributedGridVector2D *r,
         DistributedGridVector2D *u)
{
    for (int l=p; l>=0; --l, rh/=2) {
        A[l] = DirichletPoisson2D(rh);
        f[l] = DistributedGridVector2D(mpiCart, rh);
        r[l] = DistributedGridVector2D(mpiCart, rh);
        u[l] = DistributedGridVector2D(mpiCart, rh);
    }
}

int
main(int argc, char **argv)
{
    MpiCart mpiCart(argc, argv);

    int l = 12, lMin = 8;
    if (argc==5) {
        l    = std::atoi(argv[3]);
        lMin = std::atoi(argv[4]);
    }

    int rh = (1<<l);
    int p = l-lMin;
    MatType A[p+1];
    VecType f[p+1], r[p+1], u[p+1], solution(mpiCart, rh);

    allocate(mpiCart, rh, p, A, f, r, u);
    problem2(f[p], u[p], solution);

    DS ds(A[0], f[0], u[0]);
    MG mg(A, f, r, u, ds);

    timer t;
    
    t.tic();
    errorStat(0, A[p], f[p], u[p], solution);
    for (int it=1; it<=100; ++it) {
        mg.vCycle(p, 1, 1);
        errorStat(it, A[p], f[p], u[p], solution);
    }
    if (mpiCart.comm.Get_rank()==0) {
        cerr << mpiCart.numRows << " "
             << mpiCart.numCols << " "
             << l << " " << lMin << " "
             << t.toc() << endl;
    }
}
