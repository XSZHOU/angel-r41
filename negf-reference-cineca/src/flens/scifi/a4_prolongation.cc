#include <cassert>
#include <cmath>
#include <iostream>

#include <flens/flens.h>
#include <scifi/prolongation.h>

using namespace flens;
using namespace std;

//
// Assignment 4:  edit file prolongation.h
//

int
main()
{
    Prolongation P;
    DenseVector<Array<double> > x(_(0,4)),
                                y(_(0,8));
                                
    x(_(1,3)) = 1, 2, 3;
    y += P*x;
    
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;

}
