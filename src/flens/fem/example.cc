#include <ctime>
#include <fem/fem.h>

using namespace flens;
using namespace std;


double
cputime()
{
    return double(clock())/CLOCKS_PER_SEC;
}

double
g(const FEM::Node &node)
{
    /*
    if (node(2)==0) {
        return node(1)*(1-node(1));
    }
    */
    if (node(1)==1) {
        return 0.25*node(2)*(3-node(2));
    }
    return 0;
}

double
f(const FEM::Node &node)
{
    double v = 0;
    if (node(1)<0.25) {
        v = node(1);
    } else {
        v = 0.5 - node(1);
    }
    return 20*v;
}

int
main()
{
	cerr << "initialize mesh...";
    FEM fem("nodes", "triangles", f, g);
    cerr << "DONE" << endl;

	cerr << "assembling...";
    double time = cputime();
    fem.assemble();
    time = cputime() - time;
    cerr << "DONE (" << time << "s)" << endl;

// solve benoetigt eigentlich den pardiso loeser
//    fem.solve();
//    fem.writeSolution("sol");    

    //cout << "A = " << fem.A << endl;
    //cout << "b = " << fem.b << endl;


    return 0;
}
