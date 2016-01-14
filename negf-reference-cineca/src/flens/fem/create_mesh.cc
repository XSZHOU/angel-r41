#include <fstream>
#include <iostream>

using namespace std;

int
_kellerNodeNumber(int i, int j, int m)
{
    return i+(m+1)*j+1;
}

/*
 * a, b : size of rectangular domain
 * m, n : number of slices
 */
void
writeKellerTriangulation(double a, double b, int m, int n,
                         const char *nodesFile, const char *trianglesFile)
{
    std::ofstream nodes(nodesFile), triangles(trianglesFile);
    double dx = a/m, dy = b/n;

    nodes << (m+1)*(n+1) << std::endl;
    for (int j=0; j<=n; ++j) {
        for (int i=0; i<=m; ++i) {
            nodes << dx*i << " " << dy*j << std::endl;
        }
    } 

    triangles << std::endl << m*n*2 << std::endl;
    for (int j=1; j<=n; ++j) {
        for (int i=1; i<=m; ++i) {
            triangles << _kellerNodeNumber(i-1,j-1,m) << " "
                      << _kellerNodeNumber(i,j-1,m) << " "
                      << _kellerNodeNumber(i-1,j,m) << std::endl;
            triangles << _kellerNodeNumber(i,j,m) << " "
                      << _kellerNodeNumber(i-1,j,m) << " "
                      << _kellerNodeNumber(i,j-1,m) << std::endl;
        }
    }
}

int
main()
{
    writeKellerTriangulation(1, 1, 100, 100, "nodes", "triangles");
    return 0;
}

