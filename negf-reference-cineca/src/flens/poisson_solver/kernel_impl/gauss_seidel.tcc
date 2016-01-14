namespace flens {

//-- dirichlet poisson 1d ------------------------------------------------------

template <typename F, typename U>
void
dp1d_gauss_seidel(int rh, const DenseVector<F> &f, DenseVector<U> &u)
{
    double hh = 1./(rh*rh);
    int i0 = u.firstIndex()+1,
        i1 = u.lastIndex()-1;

    for (int i=i0; i<=i1; ++i) {
        u(i) = 0.5*(u(i-1)+u(i+1)+hh*f(i));
    }
}

//-- dirichlet poisson 2d ------------------------------------------------------

template <typename F, typename U>
void
dp2d_gauss_seidel_red(int rh, const GeMatrix<F> &f, GeMatrix<U> &u)
{
    int rhh = rh*rh;
    double hh = 1./rhh;

    int M = u.lastRow()-1;
    int N = u.lastCol()-1;

    for (int i=1; i<=M; i+=2) {
        for (int j=1; j<=N; j+=2) {
            u(i,j) = 0.25*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)+hh*f(i,j));
        }
    }
    int i0 = (u.firstRow()==0) ? 2 : 0;
    int j0 = (u.firstCol()==0) ? 2 : 0;
    for (int i=i0; i<=M; i+=2) {
        for (int j=j0; j<=N; j+=2) {
            u(i,j) = 0.25*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)+hh*f(i,j));
        }
    }
}

template <typename F, typename U>
void
dp2d_gauss_seidel_black(int rh, const GeMatrix<F> &f, GeMatrix<U> &u)
{
    int rhh = rh*rh;
    double hh = 1./rhh;

    int M = u.lastRow()-1;
    int N = u.lastCol()-1;

    int j0 = (u.firstCol()==0) ? 2 : 0;
    for (int i=1; i<=M; i+=2) {
        for (int j=j0; j<=N; j+=2) {
            u(i,j) = 0.25*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)+hh*f(i,j));
        }
    }
    int i0 = (u.firstRow()==0) ? 2 : 0;
    for (int i=i0; i<=M; i+=2) {
        for (int j=1; j<=N; j+=2) {
            u(i,j) = 0.25*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)+hh*f(i,j));
        }
    }
}

//-- poisson 2d ----------------------------------------------------------------

template <typename T>
void
firstRedPoint(T i0, T j0, T &iR, T &jR)
{
    if ((i0+j0)%2==0) {
        iR = i0;
        jR = j0;
    } else {
        iR = i0+1;
        jR = j0;
    }
}

template <typename T>
void
secondRedPoint(T i0, T j0, T &iR, T &jR)
{
    if ((i0+j0)%2==0) {
        iR = i0+1;
        jR = j0+1;
    } else {
        iR = i0;
        jR = j0+1;
    }
}

template <typename T>
void
firstBlackPoint(T i0, T j0, T &iR, T &jR)
{
    if ((i0+j0)%2==0) {
        iR = i0+1;
        jR = j0;
    } else {
        iR = i0;
        jR = j0;
    }
}

template <typename T>
void
secondBlackPoint(T i0, T j0, T &iR, T &jR)
{
    if ((i0+j0)%2==0) {
        iR = i0;
        jR = j0+1;
    } else {
        iR = i0+1;
        jR = j0+1;
    }
}

template <typename F, typename U>
void
p2d_gauss_seidel_red(int rh, const GeMatrix<F> &f, GeMatrix<U> &u)
{
    int rhh = rh*rh;
    double hh = 1./rhh;

    int i0 = u.firstRow()+1;
    int j0 = u.firstCol()+1;

    int i1 = u.lastRow()-1;
    int j1 = u.lastCol()-1;

    int iR, jR;

    firstRedPoint(i0, j0, iR, jR);
    for (int i=iR; i<=i1; i+=2) {
        for (int j=jR; j<=j1; j+=2) {
            u(i,j) = 0.25*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)+hh*f(i,j));
        }
    }

    secondRedPoint(i0, j0, iR, jR);
    for (int i=iR; i<=i1; i+=2) {
        for (int j=jR; j<=j1; j+=2) {
            u(i,j) = 0.25*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)+hh*f(i,j));
        }
    }
}

template <typename F, typename U>
void
p2d_gauss_seidel_black(int rh, const GeMatrix<F> &f, GeMatrix<U> &u)
{
    int rhh = rh*rh;
    double hh = 1./rhh;

    int i0 = u.firstRow()+1;
    int j0 = u.firstCol()+1;

    int i1 = u.lastRow()-1;
    int j1 = u.lastCol()-1;

    int iB, jB;

    firstBlackPoint(i0, j0, iB, jB);
    for (int i=iB; i<=i1; i+=2) {
        for (int j=jB; j<=j1; j+=2) {
            u(i,j) = 0.25*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)+hh*f(i,j));
        }
    }

    secondBlackPoint(i0, j0, iB, jB);
    for (int i=iB; i<=i1; i+=2) {
        for (int j=jB; j<=j1; j+=2) {
            u(i,j) = 0.25*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)+hh*f(i,j));
        }
    }
}

} // namespace flens
