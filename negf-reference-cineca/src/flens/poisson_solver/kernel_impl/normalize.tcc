#include <cmath>

namespace flens {

//-- neumann poisson 1d --------------------------------------------------------

template <typename U>
void
np1d_normalize(DenseVector<U> &u)
{
    int i0 = u.firstIndex()+1,
        i1 = u.lastIndex()-1;

    double sum = 0;
    for (int i=i0; i<=i1; ++i) {
        sum += u(i);
    }
    sum /= (i1-i0+1);
    for (int i=i0; i<=i1; ++i) {
        u(i) -= sum;
    }

    // set neumann boundary conditions
    u(i0-1) = u(i0);
    u(i1+1) = u(i1);
}

//-- neumann poisson 2d --------------------------------------------------------

template <typename U>
void
np2d_bc(GeMatrix<U> &u)
{
    int i0 = u.firstRow()+1;
    int j0 = u.firstCol()+1;

    int i1 = u.lastRow()-1;
    int j1 = u.lastCol()-1;

    // set neumann boundary conditions
    u(i0-1,_(j0,j1)) = u(i0,_(j0,j1));
    u(i1+1,_(j0,j1)) = u(i1,_(j0,j1));

    u(_(i0,i1),j0-1) = u(_(i0,i1),j0);
    u(_(i0,i1),j1+1) = u(_(i0,i1),j1);
}

template <typename U>
void
np2d_normalize(GeMatrix<U> &u)
{
    int i0 = u.firstRow()+1;
    int j0 = u.firstCol()+1;

    int i1 = u.lastRow()-1;
    int j1 = u.lastCol()-1;

    double sum = 0;
    for (int i=i0; i<=i1; ++i) {
        for (int j=j0; j<=j1; ++j) {
            sum += u(i,j);
        }
    }
    sum /= (i1-i0+1)*(j1-j0+1);
    for (int i=i0; i<=i1; ++i) {
        for (int j=j0; j<=j1; ++j) {
            u(i,j) -= sum;
        }
    }

    // set neumann boundary conditions
    np2d_bc(u);
}

} // namespace flens
