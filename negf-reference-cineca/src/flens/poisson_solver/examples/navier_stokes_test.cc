#include <poisson_solver/flens_impl.h>
#include <poisson_solver/examples/navier_stokes.h>

#include <dirichlet_poisson.h>

using namespace flens;
using namespace std;

typedef StaggeredGridVector2D<false, true>  GridTypeU;
typedef StaggeredGridVector2D<true, false>  GridTypeV;
typedef StaggeredGridVector2D<true, true>   GridTypeP;

typedef GaussSeidelRedBlack<NeumannPoisson2D, GridTypeP>  Smoother;
typedef StationaryIterativeSolver<Smoother>               DS;


int
main()
{
    int l = 2;

    int rh = (1<<l);

    GridTypeU u(rh), fx(rh);
    GridTypeV v(rh), fy(rh);
    GridTypeP p(rh), rhs(rh);

    setBoundaryCondition(rh, u.grid, v.grid);

    cout << "u = " << u << endl;
    cout << "v = " << v << endl;

    double re = 1000;
    double pr = 7;
    double delta = 0.5;

    double dt = computeTimeStep(rh, u.grid, v.grid, re, pr, delta);

    cout << "dt = " << dt << endl;

    double gamma = 0.9;
    double gx = 0;
    double gy = -9.81;

    computeForce(rh, u.grid, v.grid, dt, gamma, re, gx, gy, fx.grid, fy.grid);

    cout << "fx = " << fx << endl;
    cout << "fy = " << fy << endl;

    computePoissonRhs(rh, fx.grid, fy.grid, dt, rhs.grid);

    cout << "rhs = " << rhs << endl;

    NeumannPoisson2D NP(rh);
    DS ds(NP, rhs, p);

    ds.solve();

    cout << "p = " << p << endl;

    computeVelocity(rh, fx.grid, fy.grid, p.grid, dt, u.grid, v.grid);

    cout << "u = " << u << endl;
    cout << "v = " << v << endl;
}
