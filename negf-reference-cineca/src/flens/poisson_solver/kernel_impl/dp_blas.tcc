namespace flens {

template <typename X>
void
dp2d_copy(const GeMatrix<X> &x,  GeMatrix<X> &y)
{
    int i0 = x.firstRow()+1;
    int j0 = x.firstCol()+1;

    int i1 = x.lastRow()-1;
    int j1 = x.lastCol()-1;

    for (int i=i0; i<=i1; ++i) {
        for (int j=j0; j<=j1; ++j) {
            y(i,j) = x(i,j);
        }
    }
}

template <typename X>
void
dp2d_scal(double alpha, GeMatrix<X> &x)
{
    int i0 = x.firstRow()+1;
    int j0 = x.firstCol()+1;

    int i1 = x.lastRow()-1;
    int j1 = x.lastCol()-1;

    for (int i=i0; i<=i1; ++i) {
        for (int j=j0; j<=j1; ++j) {
            x(i,j) *= alpha;
        }
    }
}

template <typename X, typename Y>
double
dp2d_dot(const GeMatrix<X> &x, const GeMatrix<Y> &y)
{
    int i0 = x.firstRow()+1;
    int j0 = x.firstCol()+1;

    int i1 = x.lastRow()-1;
    int j1 = x.lastCol()-1;

    double result = 0;
    for (int i=i0; i<=i1; ++i) {
        for (int j=j0; j<=j1; ++j) {
            result += x(i,j)*y(i,j);
        }
    }
    return result;
}

template <typename T, typename X, typename Y>
void
dp2d_mv(int rh, T alpha, const GeMatrix<X> &x, T beta, GeMatrix<Y> &y)
{
    int rhh = rh*rh;

    int i0 = x.firstRow()+1;
    int j0 = x.firstCol()+1;

    int i1 = x.lastRow()-1;
    int j1 = x.lastCol()-1;

    for (int i=i0; i<=i1; ++i) {
        for (int j=j0; j<=j1; ++j) {
            y(i,j) = alpha*rhh*(4*x(i,j)-x(i-1,j)-x(i+1,j)-x(i,j-1)-x(i,j+1)) + beta*y(i,j);
        }
    }
}

} // namespace flens
