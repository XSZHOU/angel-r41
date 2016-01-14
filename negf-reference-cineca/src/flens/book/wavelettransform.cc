#include <flens/flens.h>
#include <wavelets/bspline.h>
#include <wavelets/dwt.h>
#include <wavelets/wavelets.h>

using namespace std;
using namespace flens;
using namespace wavelets;

typedef DenseVector<Array<double> > DDeVector;

int
main()
{
    DDeVector a = N(2),                    a_ = N_(2,2),
              b = scalCoeff2waveCoeff(a_), b_ = scalCoeff2waveCoeff(a);
              
    DDeVector signal(_(0,7));
    signal = _(0,7);
    DDeVector single(_(0,7));
    single(_(0,3)) = dwt<PeriodicPadding>(signal, a_);
    single(_(4,7)) = dwt<PeriodicPadding>(signal, b_);
    cout << idwt<PeriodicPadding>(single, a, b) << endl;
    return 0;
}

