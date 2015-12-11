#include <flens/flens.h>
#include <fstream>
#include <sstream>

using namespace flens;
using namespace std;

typedef DenseVector<Array<double> > Vec;

int
main(int argc, char **argv)
{
    ofstream out("sparse_pattern.dat");
    int n = 10000000;
    int k = 3;
    if (argc>1) {
        n = atoi(argv[1]);
    }
    if (argc>2) {
        k = atoi(argv[2]);
    }
    out << n << endl;
    out << k << endl;

    Vec x(n);
    for (int i=1; i<=n; ++i) {
        int pos = rand() % n +1;
        int incr = rand() % 10 +1;
        x(pos) += incr;
    }

    // b1 =  A*x
    // b2 = x'*A
    Vec b1(n), b2(n);
    for (int i=1; i<=k*n; ++i) {
        int row = rand() % n +1;
        int col = rand() % n +1;
        int incr = rand() % 10 +1;
        out << row << " " << col << " " << incr << endl;
        b1(row) += incr*x(col);
        b2(col) += incr*x(row);
    }

    for (int i=1; i<=n; ++i) {
        out << x(i) << " " << b1(i) << " " << b2(i) << endl;
    }
}
