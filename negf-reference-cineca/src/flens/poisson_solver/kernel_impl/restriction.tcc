namespace flens {

//-- dirichlet poisson 1d ------------------------------------------------------

template <typename V, typename VC>
void
dp1d_restriction(const DenseVector<V> &v, DenseVector<VC> &vc)
{
    int i0 = vc.firstIndex()+1,
        i1 = vc.lastIndex()-1;

    for (int i=i0, I=2*i; i<=i1; ++i, I+=2) {
        vc(i) = 0.25*v(I-1) + 0.5*v(I) + 0.25*v(I+1);
    }
}

//-- neumann poisson 1d --------------------------------------------------------

template <typename V, typename VC>
void
np1d_restriction(const DenseVector<V> &v, DenseVector<VC> &vc)
{
    int i0 = vc.firstIndex()+1,
        i1 = vc.lastIndex()-1;

    for (int i=i0, I=2*i; i<=i1; ++i, I+=2) {
        vc(i) = 0.5*(v(I) + v(I+1));
    }
    vc(i0-1) = vc(i0);
    vc(i1+1) = vc(i1);
}

//-- dirichlet poisson 2d ------------------------------------------------------

template <typename V, typename VC>
void
dp2d_restriction_hw(const GeMatrix<V> &v, GeMatrix<VC> &vc)
{
    int i0 = vc.firstRow()+1,
        j0 = vc.firstCol()+1;
    int m = vc.lastRow()-1,
        n = vc.lastCol()-1;

    for (int i=i0, I=2*i0; i<=m; ++i, I+=2) {
        for (int j=j0, J=2*j0; j<=n; ++j, J+=2) {
            vc(i,j) = (                 v(I-1,J)
                     +   v(I  ,J-1) + 4*v(I  ,J) + v(I  ,J+1)
                                    +   v(I+1,J))/8;
        }
    }
}

//-- neumann poisson 2d --------------------------------------------------------

template <typename V, typename VC>
void
np2d_restriction_hw(const GeMatrix<V> &v, GeMatrix<VC> &vc)
{
    int i0 = vc.firstRow()+1;
    int i1 = vc.lastRow()-1;
    int j0 = vc.firstCol()+1;
    int j1 = vc.lastCol()-1;

    for (int i=i0, I=2*i0; i<=i1; ++i, I+=2) {
        for (int j=j0, J=2*j0; j<=j1; ++j, J+=2) {
            vc(i,j) = 0.25*(v(I,J) + v(I+1,J) + v(I,J+1) + v(I+1,J+1));
        }
    }
//    np2d_bc(vc);
}

} // namespace flens
