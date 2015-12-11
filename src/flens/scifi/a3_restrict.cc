#include <cassert>
#include <cmath>
#include <iostream>

#include <flens/flens.h>
#include <scifi/restriction.h>

using namespace flens;
using namespace std;

//
// Assignment 3:  edit file restriction.h
//

int
main()
{
    Restriction R;
    DenseVector<Array<double> > x(_(0,8)),
                                y(_(0,4));
                                
    x(_(1,7)) = 1, 3, 3, 5, 5, 7, 7;
    y = R*x;
    
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;

}
