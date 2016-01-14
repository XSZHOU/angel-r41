#include <poisson_solver/flens_impl.h>

#include <dirichlet_poisson.h>

using namespace flens;
using namespace std;

const int l = 10, lMin = 2;
const int rh = (1<<l);
const int p = l-lMin;

DenseVector<Array<double > >  f[p+1], r[p+1], u[p+1], solution(rh);

fftw_plan plan;

void
mg_vCycle(int l, int rh, int v1, int v2)
{
    if (l==0) {
        dp1d_fastpoissonsolver(rh, plan, f[0], u[0]);
    } else {
        for (int v=1; v<=v1; ++v) {
            dp1d_gauss_seidel(rh, f[l], u[l]);
        }
        dp1d_residual(rh, f[l], u[l], r[l]);
        dp1d_restriction(r[l], f[l-1]);

        u[l-1] = 0;
        mg_vCycle(l-1, rh/2, v1, v2);

        dp1d_prolongation(u[l-1], u[l]);
        for (int v=1; v<=v2; ++v) {
            dp1d_gauss_seidel(rh, f[l], u[l]);
        }
    }
}

void
init()
{
    int rh = (1<<l);
    for (int l=p; l>=0; --l, rh/=2) {
        f[l] = r[l] = u[l] = DenseVector<Array<double> >(_(0, rh));
        if (l==0) {
            plan = dp1d_fastpoissonsolver_init(rh, u[l], FFTW_MEASURE);
        }
    }
}

int
main()
{
    init();

    problem1(rh, f[p], u[p], solution);

    std::cerr << "f[0] = " << f[0] << std::endl;
    std::cerr << "u[0] = " << u[0] << std::endl;

    dp1d_residual(rh, f[p], u[p], r[p]);
    cerr << dp1d_normInf(r[p]) << endl;
    for (int it=1; it<=20; ++it) {
        mg_vCycle(p, rh, 1, 1);
        dp1d_residual(rh, f[p], u[p], r[p]);
        cerr << dp1d_normInf(r[p]) << endl;
    }
}
