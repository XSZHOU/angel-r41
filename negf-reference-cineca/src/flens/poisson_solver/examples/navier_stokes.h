#ifndef POISSON_FLENS_NAVIER_STOKES_H
#define POISSON_FLENS_NAVIER_STOKES_H 1

namespace flens {

template <typename U, typename V>
    void
    setBoundaryCondition(int rh, GeMatrix<U> &u, GeMatrix<V> &v, double u0=1);

//------------------------------------------------------------------------------

template <typename TMP>
    void
    setTemperatureBoundaryCondition(GeMatrix<TMP> &T, double t0=1);

//------------------------------------------------------------------------------

template <typename U, typename V>
    double
    computeTimeStep(int rh, const GeMatrix<U> &u, const GeMatrix<V> &v,
                    double re, double pr, double delta);

//------------------------------------------------------------------------------

template <typename U, typename V, typename FX, typename FY>
    void
    computeForce(int rh, const GeMatrix<U> &u, const GeMatrix<V> &v,
                 double dt, double gamma, double re, double gx, double gy,
                 GeMatrix<FX> &fx, GeMatrix<FY> &fy);

//------------------------------------------------------------------------------

template <typename TMP, typename FY>
    void
    computeBouyantForce(int rh, const GeMatrix<TMP> &T, double beta, double dt,
                        double gy, GeMatrix<FY> &fy);

//------------------------------------------------------------------------------

template <typename FX, typename FY, typename RHS>
    void
    computePoissonRhs(int rh, const GeMatrix<FX> &fx, const GeMatrix<FY> &fy,
                      double dt, GeMatrix<RHS> &rhs);

//------------------------------------------------------------------------------

template <typename FX, typename FY, typename P, typename U, typename V>
    void
    computeVelocity(int rh, const GeMatrix<FX> &fx, const GeMatrix<FY> &fy,
                    const GeMatrix<P> &p, double dt,
                    GeMatrix<U> &u, GeMatrix<V> &v);

//------------------------------------------------------------------------------

template <typename U, typename V, typename TMP>
    void
    computeTemperature(int rh, const GeMatrix<U> &u, const GeMatrix<V> &v,
                       double dt, double re, double pr, const GeMatrix<TMP> &T,
                       GeMatrix<TMP> &T2);

//------------------------------------------------------------------------------

template <typename U, typename V>
    void
    writeVelocity(const char *filename, int counter, int rh,
                  const GeMatrix<U> &u, const GeMatrix<V> &v);

//------------------------------------------------------------------------------

template <typename TMP>
    void
    writeTemperature(const char *filename, int counter, int rh,
                     const GeMatrix<TMP> &T);

//------------------------------------------------------------------------------

struct Particle {
    double x, y, u, v;
};

template <typename T>
class ParticleField
{
    public:
        ParticleField();

        ~ParticleField();

        template <typename U, typename V>
            void
            update(int rh, const GeMatrix<U> &u, const GeMatrix<V> &v, double dt);

        void
        writeToFile(const char *filename, int counter);

    private:
        int numParticles;
        Particle *particle;
};

} // namespace flens

#include <poisson_solver/examples/navier_stokes.tcc>

#endif // POISSON_FLENS_NAVIER_STOKES_H
