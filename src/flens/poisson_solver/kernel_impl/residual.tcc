namespace flens {

//-- dirichlet poisson 1d ------------------------------------------------------

template <typename U, typename F, typename R>
void
dp1d_residual(int rh, const DenseVector<F> &f, const DenseVector<U> &u,
              DenseVector<R> &r)
{
    int rhh = rh*rh;
    int i0 = u.firstIndex()+1,
        i1 = u.lastIndex()-1;

    for (int i=i0; i<=i1; ++i) {
        r(i) = f(i) + rhh*(u(i-1)-2*u(i)+u(i+1));
    }

}

//-- dirichlet poisson 2d ------------------------------------------------------

template <typename U, typename F, typename R>
void
dp2d_residual(int rh, const GeMatrix<F> &f, const GeMatrix<U> &u,
              GeMatrix<R> &r)
{
    int rhh = rh*rh;
    int i0 = u.firstRow()+1,
        j0 = u.firstCol()+1;
    int i1 = u.lastRow()-1,
        j1 = u.lastCol()-1;

    for (int i=i0; i<=i1; ++i) {
        for (int j=j0; j<=j1; ++j) {
            r(i,j) = f(i,j)-rhh*(4*u(i,j)-u(i-1,j)-u(i+1,j)-u(i,j-1)-u(i,j+1));
        }
    }
}

} // namespace flens
