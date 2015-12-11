#include <flens/flens.h>

using namespace flens;
using namespace std;

typedef GeMatrix<FullStorage<double, ColMajor> > GEMatrix;
typedef DenseVector<Array<double> >              DEVector;

int
main()
{
    int m = 5,
        n = 3;
    GEMatrix A(m, n), U(m,m), VT(n,n);
    DEVector s(min(m,n));

    for (int i=1; i<=m; ++i) {
        for (int j=1; j<=n; ++j) {
            A(i,j) = i + j;
        }
    }

    cout << "A = " << A << endl;

    // note that LAPACK will overwrite matrix A
    svd(A, s, U, VT);

    cout << "s = " << s << endl;
    cout << "U = " << U << endl;
    cout << "V^T = " << VT << endl;
    return 0;
}












