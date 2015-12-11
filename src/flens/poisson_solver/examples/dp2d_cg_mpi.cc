#include <poisson_solver/flens_impl.h>
#include <dirichlet_poisson.h>

using namespace flens;
using namespace std;

typedef DirichletPoisson2D       MatType;
typedef DistributedGridVector2D  VecType;

int
main(int argc, char **argv)
{
    MpiCart mpiCart(argc, argv);

    int l = 10;
    int rh = (1<<l);
    MatType A(rh);
    VecType f(mpiCart, rh), u(mpiCart, rh), solution(mpiCart, rh);

    problem2(f, u, solution);

    errorStat(0, A, f, u, solution);
    cout << "iterations = " << cg(A, u, f) << endl;
    errorStat(1, A, f, u, solution);
}
