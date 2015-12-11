#include <flens/flens.h>
#include <timer.h>

using namespace flens;
using namespace std;

int
main(int argc, char **argv)
{
    {
        TinyVector<double,10> x, y, z;

        timer t;

        t.tic();
        for (int i=1; i<=8000000; ++i) {
            y = x;
            z = x + y;
        }
        cout << t.toc() << endl;
    }

    {
        DenseVector<Array<double> > x(10), y(10), z(10);

        timer t;

        t.tic();
        for (int i=1; i<=8000000; ++i) {
            y = x;
            z = x + y;
        }
        cout << t.toc() << endl;
    }

    {
        TinyVector<double,10> x, y;
        TinyMatrix<double,10,10> A;

        timer t;

        for (int i=0; i<10; ++i) {
            x(i) = i;
            for (int j=0; j<10; ++j) {
                A(i,j) = i+j;
            }
        }

        t.tic();
        for (int i=1; i<=8000000; ++i) {
            y = A*x;
            //y = transpose(A)*x;
        }
        cout << t.toc() << endl;

        cout << "y = " << y;
    }

    {
        DenseVector<Array<double> >               x(10), y(10);
        GeMatrix<FullStorage<double, RowMajor> >  A(10,10);

        timer t;

        for (int i=1; i<=10; ++i) {
            x(i) = i-1;
            for (int j=1; j<=10; ++j) {
                A(i,j) = i+j-2;
            }
        }

        t.tic();
        for (int i=1; i<=8000000; ++i) {
            y = A*x;
            //y = transpose(A)*x;
        }
        cout << t.toc() << endl;

        cout << "y = " << y;
    }

}
