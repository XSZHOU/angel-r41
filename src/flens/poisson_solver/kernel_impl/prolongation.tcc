namespace flens {

//-- dirichlet poisson 1d ------------------------------------------------------

template <typename VC, typename V>
void
dp1d_prolongation(const DenseVector<VC> &vc, DenseVector<V> &v)
{
    int i0 = vc.firstIndex()+1,
        i1 = vc.lastIndex()-1;

    for (int i=i0, I=2*i; i<=i1; ++i, I+=2) {
        v(I-1) += 0.5*vc(i);
        v(I)   += vc(i);
        v(I+1) += 0.5*vc(i);
    }
}

//-- neumann poisson 1d ------------------------------------------------------

template <typename VC, typename V>
void
np1d_prolongation(const DenseVector<VC> &vc, DenseVector<V> &v)
{
    int i0 = vc.firstIndex()+1,
        i1 = vc.lastIndex();

    for (int i=i0, I=2*i; i<=i1; ++i, I+=2) {
        v(I-1)   += 0.75*vc(i-1) + 0.25*vc(i);
        v(I) += 0.25*vc(i-1) + 0.75*vc(i);
    }
    int I0 = v.firstIndex(),
        I1 = v.lastIndex();

    v(I0) = v(I0+1);
    v(I1) = v(I1-1);
}

//-- dirichlet poisson 2d ------------------------------------------------------

template <typename VC, typename V>
void
dp2d_prolongation(const GeMatrix<VC> &vc, GeMatrix<V> &v)
{
    int i0 = vc.firstRow()+1,
        j0 = vc.firstCol()+1;
    int i1 = vc.lastRow()-1,
        j1 = vc.lastCol()-1;

    for (int i=i0, I=2*i0; i<=i1; ++i, I+=2) {
        for (int j=j0, J=2*j0; j<=j1; ++j, J+=2) {
            v(I-1,J-1)+=.25*vc(i,j); v(I-1,J)+=.5*vc(i,j); v(I-1,J+1)+=.25*vc(i,j);
            v(I  ,J-1)+=.50*vc(i,j); v(I  ,J)+=   vc(i,j); v(I  ,J+1)+=.50*vc(i,j);
            v(I+1,J-1)+=.25*vc(i,j); v(I+1,J)+=.5*vc(i,j); v(I+1,J+1)+=.25*vc(i,j);
        }
    }
}

template <typename VC, typename V>
void
dp2d_prolongation_north(const GeMatrix<VC> &vc, GeMatrix<V> &v)
{
    int j0 = vc.firstCol()+1;
    int i1 = vc.lastRow(),
        j1 = vc.lastCol()-1;

    int i=i1, I=2*i1;
    for (int j=j0, J=2*j0; j<=j1; ++j, J+=2) {
        v(I-1,J-1)+=.25*vc(i,j); v(I-1,J)+=.5*vc(i,j); v(I-1,J+1)+=.25*vc(i,j);
        v(I  ,J-1)+=.50*vc(i,j); v(I  ,J)+=   vc(i,j); v(I  ,J+1)+=.50*vc(i,j);
    }
}

template <typename VC, typename V>
void
dp2d_prolongation_east(const GeMatrix<VC> &vc, GeMatrix<V> &v)
{
    int i0 = vc.firstRow()+1;
    int i1 = vc.lastRow()-1,
        j1 = vc.lastCol();

    int j=j1, J=2*j1;
    for (int i=i0, I=2*i0; i<=i1; ++i, I+=2) {
        v(I-1,J-1)+=.25*vc(i,j); v(I-1,J)+=.5*vc(i,j);
        v(I  ,J-1)+=.50*vc(i,j); v(I  ,J)+=   vc(i,j);
        v(I+1,J-1)+=.25*vc(i,j); v(I+1,J)+=.5*vc(i,j);
    }
}

template <typename VC, typename V>
void
dp2d_prolongation_north_east(const GeMatrix<VC> &vc, GeMatrix<V> &v)
{
    int i1 = vc.lastRow(),
        j1 = vc.lastCol();

    int i=i1, I=2*i1;
    int j=j1, J=2*j1;
    v(I-1,J-1)+=.25*vc(i,j); v(I-1,J)+=.5*vc(i,j);
    v(I  ,J-1)+=.50*vc(i,j); v(I  ,J)+=   vc(i,j);
}

//-- neumann poisson 2d --------------------------------------------------------

template <typename VC, typename V>
void
np2d_prolongation(const GeMatrix<VC> &vc, GeMatrix<V> &v)
{
    int i0 = vc.firstRow()+1;
    int j0 = vc.firstCol()+1;
    int i1 = vc.lastRow()-1;
    int j1 = vc.lastCol()-1;

    for (int i=i0, I=2*i0; i<=i1; ++i, I+=2) {
        for (int j=j0, J=2*j0; j<=j1; ++j, J+=2) {
            v(I,J)+=     0.75*(0.25*vc(i-1,  j) + 0.75*vc(  i,  j))
                       + 0.25*(0.25*vc(i-1,j-1) + 0.75*vc(i-1,j-1));

            v(I,J+1)+=   0.25*(0.25*vc(i-1,j+1) + 0.75*vc(i,j+1))
                       + 0.75*(0.25*vc(i-1,  j) + 0.75*vc(i,  j));

            v(I+1,J)+=   0.75*(0.75*vc(i,j)   + 0.25*vc(i+1,j))
                       + 0.25*(0.75*vc(i,j-1) + 0.25*vc(i+1,j-1));


            v(I+1,J+1)+= 0.75*(0.75*vc(i,  j) + 0.25*vc(i+1,  j))
                       + 0.25*(0.75*vc(i,j+1) + 0.25*vc(i+1,j+1));
        }
    }
    np2d_bc(v);
}

} // namespace flens
