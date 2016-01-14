#include <cmath>

namespace flens {

//-- dirichlet poisson 1d ------------------------------------------------------

template <typename U>
double
dp1d_normInf(const DenseVector<U> &u)
{
    int i0 = u.firstIndex(),
        i1 = u.lastIndex();
    int I = i0;

    for (int i=i0; i<=i1; ++i) {
        if (std::abs(u(i))> std::abs(u(I))) {
            I = i;
        }
    }
    return std::abs(u(I));
}

template <typename U>
double
dp1d_norm2sqr(const DenseVector<U> &u)
{
    int i0 = u.firstIndex()+1,
        i1 = u.lastIndex()-1;

    double result = 0;
    for (int i=i0; i<=i1; ++i) {
        result += u(i)*u(i);
    }
    return result;
}

//-- dirichlet poisson 2d ------------------------------------------------------

template <typename U>
double
dp2d_normInf(const GeMatrix<U> &u)
{
    int I = u.firstRow()+1,
        J = u.firstCol()+1;

    for (int i=u.firstRow()+1; i<=u.lastRow()-1; ++i) {
        for (int j=u.firstCol()+1; j<=u.lastCol()-1; ++j) {
            if (std::abs(u(i,j))> std::abs(u(I,J))) {
                I = i;
                J = j;
            }
        }
    }
    return std::abs(u(I,J));
}

template <typename U>
double
dp2d_norm2sqr(const GeMatrix<U> &u)
{
    double result = 0;
    for (int i=u.firstRow()+1; i<=u.lastRow()-1; ++i) {
        for (int j=u.firstCol()+1; j<=u.lastCol()-1; ++j) {
            result += u(i,j)*u(i,j);
        }
    }
    return result;
}

} // namespace flens
