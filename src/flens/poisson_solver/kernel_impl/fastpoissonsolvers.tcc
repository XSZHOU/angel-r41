namespace flens {

//-- dirichlet poisson 1d ------------------------------------------------------

template <typename U>
fftw_plan
dp1d_fastpoissonsolver_init(int rh, DenseVector<U> &u, unsigned fftw_flags)
{
    int N = rh - 1;
    return fftw_plan_r2r_1d(N, &u(1), &u(1), FFTW_RODFT00, fftw_flags);
}

template <typename F, typename U>
void
dp1d_fastpoissonsolver(int rh, fftw_plan &plan,
                       const DenseVector<F> &f, DenseVector<U> &u)
{
    double hh = 1./(rh*rh);
    int N = rh - 1;

    // setup q from f and boundary nodes in u
    u(_(1,N)) = hh*f(_(1,N));
    u(1) += u(0);
    u(N) += u(N+1);

    //- step 1 -
    fftw_execute(plan);

    //- step 2 -
    for (int k=1; k<=N; ++k) {
        u(k) /= 2*(1-cos(k*M_PI/(N+1)));
    }

    //- step 3 -
    fftw_execute(plan);
    u(_(1,N)) /= 2*(N+1);
}

//-- dirichlet poisson 2d ------------------------------------------------------

template <typename U>
fftw_plan
dp2d_fastpoissonsolver_init(int rh, GeMatrix<U> &u, unsigned fftw_flags)
{
    int N = rh - 1;

    // - setup ffwt for 2D -
    int n[2], nembed[2];
    n[0] = N;
    n[1] = N;

    nembed[0] = N+2;
    nembed[1] = N+2;

    fftw_r2r_kind kind[2];
    kind[0] = FFTW_RODFT00;
    kind[1] = FFTW_RODFT00;

    fftw_plan plan = fftw_plan_many_r2r(2, n, 1,
                                        &u(1,1), nembed, 1, 0,
                                        &u(1,1), nembed, 1, 0,
                                        kind,
                                        fftw_flags);

    return plan;
}


template <typename F, typename U>
void
dp2d_fastpoissonsolver(int rh, fftw_plan &plan,
                       const GeMatrix<F> &f, GeMatrix<U> &u)
{
    double hh = 1./(rh*rh);
    int N = rh - 1;

    // setup q from f and boundary nodes in u
    u(_(1,N),_(1,N)) = hh*f(_(1,N),_(1,N));
    u(1,_(1,N)) += u(0,_(1,N));
    u(N,_(1,N)) += u(N+1,_(1,N));
    u(_(1,N),1) += u(_(1,N),0);
    u(_(1,N),N) += u(_(1,N),N+1);

    //- step 1 -
    fftw_execute(plan);

    //- step 2 -
    for (int k=1; k<=N; ++k) {
        for (int l=1; l<=N; ++l) {
            u(k,l) /= 2*(1-cos(k*M_PI/(N+1)))+2*(1-cos(l*M_PI/(N+1)));
        }
    }

    //- step 3 -
    fftw_execute(plan);
    u(_(1,N),_(1,N)) /= 4*(N+1)*(N+1);
}

} // namespace flens
