namespace flens {

template <typename F, typename U>
void
dp1d_cholesky(int rh, const DenseVector<F> &f, DenseVector<U> &u)
{
    int N = rh - 1;
    double hh = 1./(rh*rh);

    // forward substitution
    u(1) = hh*f(1)+u(0);
    for (int k=2; k<N; ++k) {
        u(k) = hh*f(k) + u(k-1)*(k-1)/k;
    }
    u(N) = hh*f(N)+u(N+1) + u(N-1)*(N-1)/N;

    // backward substitution
    u(N) = u(N)*N/(N+1);
    for (int k=N-1; k>=1; --k) {
        u(k) = (u(k)+u(k+1))*k/(k+1);
    }
}

} // namespace flens
