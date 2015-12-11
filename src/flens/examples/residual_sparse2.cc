#include <flens/flens.h>
#include <timer.h>
#include <fstream>

using namespace flens;
using namespace std;

typedef SparseGeMatrix<CRS<double> > Mat;
typedef DenseVector<Array<double> >  Vec;

int
main()
{
    ifstream in("sparse_pattern.dat");
    int n, k;
    in >> n >> k;

    Mat A(n, n, k);
    Vec x(n), b1(n), b2(n), r(n);

    {
        timer t;

        t.tic();
        for (int i=1; i<=k*n; ++i) {
            int row, col, incr;

            in >> row >> col >> incr;
            A(row, col) += incr;
        }
        A.finalize();
        cout << "fill: " << t.toc() << endl;

        for (int i=1; i<=n; ++i) {
            int x_, b1_, b2_;

            in >> x_ >> b1_ >> b2_;
            x(i) = x_;
            b1(i) = b1_;
            b2(i) = b2_;
        }
    }

    {
        timer t;

        t.tic();
        r = b1 - A*x;
        cout << "r = b - A*x: " << t.toc()
             << " (" << nrm2(r) << ")" << endl;
    }

    {
        timer t;

        t.tic();
        r = b2 - x*A;
        cout << "r = b - x*A: " << t.toc()
             << " (" << nrm2(r) << ")" << endl;
    }
}
