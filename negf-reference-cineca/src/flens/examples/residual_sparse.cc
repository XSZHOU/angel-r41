#include <flens/flens.h>
#include <timer.h>

using namespace flens;
using namespace std;

typedef SparseGeMatrix<CRS<double> > Mat;
typedef DenseVector<Array<double> >  Vec;

void
matrixFill(Mat &A)
{
    int n = A.numRows();
    for (int i=1; i<=n; ++i) {
        if (i>1) {
            A(i,i-1) += -1;
        }
        A(i,i) += 4;
        if (i<n) {
            A(i,i+1) += -1;
        }
    }
    A.finalize();
}

int
main(int argc, char **argv)
{
    int n = 10000000;
    if (argc>1) {
        n = atoi(argv[1]);
    }
    cout << "n = " << n << endl;

    Mat A(n, n, 3);
    Vec b(n), x(n), r(n);
    b = 2;
    b(1) = b(n) = 3;
    x = 0.99;

    {
        timer t;

        t.tic();
        matrixFill(A);
        cout << "fill: " << t.toc() << endl;
    }

    {
        timer t;

        t.tic();
        for (int k=1; k<=20; ++k) {
            r = b - A*x;
        }
        cout << "r = b - A*x: " << t.toc()
             << " (" << nrm2(r) << ")" << endl;
    }

    {
        timer t;

        t.tic();
        for (int k=1; k<=20; ++k) {
            r = A*x - b;
        }
        cout << "r = A*x - b: " << t.toc()
             << " (" << nrm2(r) << ")" << endl;
    }

    {
        timer t;

        t.tic();
        for (int k=1; k<=20; ++k) {
            r = -A*x;
        }
        cout << "r = -A*x: " << t.toc()
             << " (" << nrm2(r) << ")" << endl;
    }

}
