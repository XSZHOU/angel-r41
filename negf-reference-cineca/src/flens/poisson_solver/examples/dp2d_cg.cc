#include <poisson_solver/flens_impl.h>

#include <dirichlet_poisson.h>

using namespace flens;
using namespace std;

typedef DirichletPoisson2D MatType;
typedef GridVector2D       VecType;

int
main()
{
    int l = 8;
    int rh = (1<<l);
    MatType A(rh);
    VecType f(rh), u(rh), solution(rh);

    problem2(f, u, solution);

    errorStat(0, A, f, u, solution);
    cout << "iterations = " << cg(A, u, f) << endl;
    errorStat(1, A, f, u, solution);
}
