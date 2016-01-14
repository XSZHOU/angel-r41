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
typedef StaggeredGridVector2D<true, true>   GridTypeT;

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
    GridTypeT T(rh), T2(rh);

    ParticleField<double> particleField;

    double re = 1000;
    double pr = 7;
    double delta = 0.5;
    double beta = 2.1e-4;
    double gamma = 0.9;
    double gx = 0;
    double gy = -9.81;
    double dt = 0;
    double t=0, tSnapshot=0;
    double eps = 0.0000001;
    int    snapshot=0;

    u.grid = 0;
    v.grid = 0;
    T.grid = 20;
    while (t<=600) {
        setBoundaryCondition(rh, u.grid, v.grid, 0);
        setTemperatureBoundaryCondition(T.grid, 20);

        cout << "time = " << t;
        if (t>=tSnapshot) {
            cout << ", create snapshot " << snapshot;
            tSnapshot += 0.1;
            ++snapshot;
            writeVelocity("data/velocity", snapshot, rh, u.grid, v.grid);
            particleField.writeToFile("data/particle", snapshot);
            writeTemperature("data/temperature", snapshot, rh, T.grid);
        }
        cout << endl;

        dt = computeTimeStep(rh, u.grid, v.grid, re, pr, delta);
        if (dt>0.1) {
            dt =0.1;
        }

        computeForce(rh, u.grid, v.grid, dt, gamma, re, gx, gy, fx.grid, fy.grid);
        computeBouyantForce(rh, T.grid, beta, dt, gy, fy.grid);
        computePoissonRhs(rh, fx.grid, fy.grid, dt, rhs.grid);
        //cout << "fx =  " << fx << endl;
        //cout << "fy =  " << fy << endl;
        //cout << "rhs = " << rhs << endl;

        for (int it=1; it<=600; ++it) {
            mg.vCycle(p_mg, 1, 1);
            if (normL2(r_mg[p_mg])<eps) {
                cout << "it = " << it << endl;
                break;
            }
        }

        //cout << "p = " << p << endl;

        computeVelocity(rh, fx.grid, fy.grid, p.grid, dt, u.grid, v.grid);

        //cout << "u = " << u << endl;
        //cout << "v = " << v << endl;
        cout << "dt = " << dt << endl;

        t += dt;
        computeTemperature(rh, u.grid, v.grid, dt, re, pr, T.grid, T2.grid);
        T.grid = T2.grid;

        particleField.update(rh, u.grid, v.grid, dt);
    }
}
