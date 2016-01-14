#define NDEBUG

#include <poisson_solver/flens_impl.h>
#include <navier_stokes.h>

#include <dirichlet_poisson.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <iostream>

using namespace flens;
using namespace std;

typedef StaggeredGridVector2D<false, true>  GridTypeU;
typedef StaggeredGridVector2D<true, false>  GridTypeV;
typedef StaggeredGridVector2D<true, true>   GridTypeP;

typedef GaussSeidelRedBlack<NeumannPoisson2D, GridTypeP>  Smoother;
typedef StationaryIterativeSolver<Smoother>               DS;

typedef MultiGrid<NeumannPoisson2D, GridTypeP,
                  Restriction, Prolongation,
                  Smoother, DS>                           MG;


//-- allocator for multigrid ---------------------------------------------------

template <typename MatType, typename VecType>
void
allocate(int rh, int p, MatType *A, VecType *f, VecType *r, VecType *u)
{
    for (int l=p; l>=0; --l, rh/=2) {
        A[l] = MatType(rh);
        f[l] = VecType(rh);
        r[l] = VecType(rh);
        u[l] = VecType(rh);
    }
}

int
main()
{
    int l = 7, lMin = 2;
    int rh = (1<<l);
    int p_mg = l-lMin;

    //-- setup multigrid method ------------------------------------------------
    NeumannPoisson2D A_mg[p_mg+1];
    GridTypeP f_mg[p_mg+1], r_mg[p_mg+1], u_mg[p_mg+1];

    allocate(rh, p_mg, A_mg, f_mg, r_mg, u_mg);
    DS ds(A_mg[0], f_mg[0], u_mg[0]);
    MG mg(A_mg, f_mg, r_mg, u_mg, ds);
    //--------------------------------------------------------------------------

    GridTypeU u(rh), fx(rh);
    GridTypeV v(rh), fy(rh);
    GridTypeP &p = u_mg[p_mg], &rhs = f_mg[p_mg];

    ParticleField<double> particleField;

    double re = 1000;
    double pr = 7;
    double delta = 0.15;
    double gamma = 0.9;
    double gx = 0;
    double gy = -9.81;
    double dt = 0;
    double t=0, tSnapshot=0;
    double eps = 0.0000001;
    int    snapshot=0;


    while (t<40) {
        setBoundaryCondition(rh, u.grid, v.grid, 1);

        cout << "time = " << t << endl;
        if (t>=tSnapshot) {
            tSnapshot += 0.1;
            ++snapshot;
            writeVelocity("data/velocity", snapshot, rh, u.grid, v.grid);
            particleField.writeToFile("data/particle", snapshot);
        }

        dt = computeTimeStep(rh, u.grid, v.grid, re, pr, delta);
        computeForce(rh, u.grid, v.grid, dt, gamma, re, gx, gy, fx.grid, fy.grid);
        computePoissonRhs(rh, fx.grid, fy.grid, dt, rhs.grid);

        for (int it=1; it<=600; ++it) {
            mg.vCycle(p_mg, 1, 1);
            if (normL2(r_mg[p_mg])<eps) {
                cout << "it = " << it << endl;
                break;
            }
        }

        computeVelocity(rh, fx.grid, fy.grid, p.grid, dt, u.grid, v.grid);
        t += dt;

        particleField.update(rh, u.grid, v.grid, dt);
    }
}
