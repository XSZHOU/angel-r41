#include <flens/flens.h>

using namespace flens;
using namespace std;

/*
 *  test copying between sparse and general matrices
 */

typedef GeMatrix<FullStorage<float, ColMajor> > DEMatrix;
typedef SparseGeMatrix<CRS<float> > SPMatrix;

int
main()
{
    /*
    TODO: implement the initWith method

    SPMatrix A(3,3);
    A(1,2) = 0.2;
    A(2,1) = 0.12;
    A(2,3) = 0.24;
    A(3,2) = 0.75;
    A.finalize();

    cout << "sparse matrix A: " << A;

    // convert sparse to dense matrix
    DEMatrix B;
    B = A;
    cout << "copied to dense matrix B: " << B;


    // first way to convert dense to sparse matrix
    // ATTENTION: there is no need to call C.finalize() !
    SPMatrix C;
    C = B;
    cout << "B to sparse matrix C: " << C;

    // second way to convert dense to sparse matrix:
    // one can specify an epsillon under that the values
    // are treated as zeros (to avoid numerical problems)
    // ATTENTION: there is no need to call C.finalize() !
    SPMatrix D;
    D.initWith(B, 1E-12);
    cout << "B to sparse matrix D with eps: " << C;

    */

    return 0;
}
