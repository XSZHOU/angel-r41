#include <cmath>
#include <fem/fem.h>
#include <set>
#include <sstream>

namespace flens {
    
//------------------------------------------------------------------------------

typedef std::pair<int,int> NodePair;

bool
operator<(const NodePair &p1, const NodePair &p2)
{
    if (p1.first<p2.first) {
        return true;
    }
    if ((p1.first==p2.first) || (p1.second<p2.second)) {
        return true;
    }
    return false;
}

NodePair
makeNodePair(int n1, int n2)
{
    return (n1<n2) ? NodePair(n1, n2) : NodePair(n2,n1);
}

//------------------------------------------------------------------------------

Mesh::Mesh(const char *nodesFile, const char *trianglesFile)
{
    std::ifstream inNodes(nodesFile), inTriangles(trianglesFile);
    int _numNodes, _numTriangles;

    inNodes >> _numNodes;
    _nodes.resize(_numNodes, 2);
    for (int i=1; i<=numNodes(); ++i) {
        inNodes >> _nodes(i,1) >> _nodes(i,2);
    }

    inTriangles >> _numTriangles;
    _triangles.resize(_numTriangles, 3);
    for (int i=1; i<=numTriangles(); ++i) {
        inTriangles >> _triangles(i,1) >> _triangles(i,2) >> _triangles(i,3);
    }
    
    _nu.resize(_numNodes);
    _eta.resize(_numNodes);
    
    // -------------------------------------------------------------------------
    typedef std::map<NodePair,int> EdgeCount;
    EdgeCount edgeCount;

    for (int i=1; i<=numTriangles(); ++i) {
        Triangle t = triangle(i);
        for (int i=1; i<=3; ++i) {
            NodePair e = makeNodePair(t(i), t((i<3) ? i+1 : 1));
            ++edgeCount[e];
        }
    }

    // create set of boundary nodes
    std::set<int> bn;
    typedef EdgeCount::const_iterator EdgeCountIt;
    for (EdgeCountIt it=edgeCount.begin(); it!=edgeCount.end(); ++it) {
        if (it->second==1) {
            bn.insert(it->first.first);
            bn.insert(it->first.second);
        }
    }
    
    // boundary nodes have indices nu(M+1),...,nu(N)
    // where M = number of inner nodes, N = number of nodes
    int counter = numNodes();
    for (std::set<int>::const_iterator it=bn.begin(); it!=bn.end(); ++it) {
        _nu(counter) = *it;
        _eta(*it) = counter--;
    }
    _numInteriorNodes = counter;

    // innern nodes have indices nu(1),...,nu(M)
    counter = 1;
    for (int i=1; i<=numNodes(); ++i) {
        if (_eta(i)==0) {
            _nu(counter) = i;
            _eta(i) = counter++;
        }
    }
}
    
int
Mesh::numTriangles() const
{
    return _triangles.numRows();
}

const Mesh::Triangle
Mesh::triangle(int index) const
{
    return _triangles(index,_);
}

int
Mesh::numNodes() const
{
    return _nodes.numRows();
}
 
Mesh::Node
Mesh::node(int index) const
{
    return _nodes(index,_);
}

int
Mesh::numInteriorNodes() const
{
    return _numInteriorNodes;
}

int
Mesh::nu(int index) const
{
    return _nu(index);
}

int
Mesh::eta(int index) const
{
    return _eta(index);
}

//------------------------------------------------------------------------------

double
det(const DenseVector<Array<double> > &a, const DenseVector<Array<double> > &b)
{
    return a(1)*b(2) - a(2)*b(1);
}

double
dot(const DenseVector<Array<double> > &a, const DenseVector<Array<double> > &b)
{
    return a(1)*b(1) + a(2)*b(2);
}

//------------------------------------------------------------------------------

extern "C" {

void
pardiso_(int *pt, int *maxfct, int *mnum, int *mtype, int *phase,
         int *n, double *a, int *ia, int *ja, int *perm, int *nrhs,
         int *iparm, int *msglvl, double *b, double *sol, int *error);
}

template <typename E1, typename E2>
void
pardiso(DenseVector<E1> &x,
        SparseSymmetricMatrix<CRS<double> > &A, DenseVector<E2> &b)
{
    int pt[2*64*10];
    memset(pt, 0, 2*sizeof(int)*64*10);

    int maxfct = 1;
    int mnum = 1;
    int mtype = 2;
    int phase = 13;
    int n = A.numRows();
    double *a = A.engine().values().data();
    int *ia = A.engine().rows().data();
    int *ja = A.engine().columns().data();
    int perm[n];
    int nrhs = 1;
    int iparm[64];
    iparm[0] = 0;
    iparm[2] = 1;
    int msglvl = 1;
    double *y = b.data();
    double *sol = x.data();
    int error;

    pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs,
             iparm, &msglvl, y, sol, &error);
}

//------------------------------------------------------------------------------

FEM::FEM(const char *nodesFile, const char *trianglesFile,
         const Func &_f, const Func &_g)
    : mesh(nodesFile, trianglesFile), f(_f), g(_g),
      u(mesh.numNodes()), b(mesh.numInteriorNodes()),
      A(mesh.numInteriorNodes(), mesh.numInteriorNodes()),
      Am(3,3), bm(3), S(3, DGeMatrix(3,3))
{
    for (int i=mesh.numInteriorNodes()+1; i<=mesh.numNodes(); ++i) {
        u(i) = g(mesh.node(mesh.nu(i)));
    }
    S[0] =  0.5, -0.5,    0,
           -0.5,  0.5,    0,
              0,    0,    0;

    S[1] =    1, -0.5, -0.5,
           -0.5,    0,  0.5,
           -0.5,  0.5,    0;

    S[2] =  0.5,    0, -0.5,
              0,    0,    0,
           -0.5,    0,  0.5;
}

void
FEM::assemble()
{
    int M = mesh.numInteriorNodes();
    for (int m=1; m<=mesh.numTriangles(); ++m) {
        Triangle t = mesh.triangle(m);
        
        compute_Am_bm(t);
        
        for (int r=1; r<=3; ++r) {
            int k=mesh.eta(t(r));
            if (k<=M) {
                b(k) += bm(r);
            }
            for (int s=r; s<=3; ++s) {
                int l=mesh.eta(t(s));

				if ((k>M) && (l>M)) {
					continue;
				}
                
                if ((k<=M) && (l<=M)) {
                    A(k,l) += Am(r,s);
                } else {
                    int k_ = std::min(k,l);
                    int l_ = std::max(k,l);
                    b(k_) -= Am(r,s)*u(l_);
                }
            }
        }
    }
    A.finalize();
}

void
FEM::compute_Am_bm(const Triangle &t)
{
    const Node &a1 = mesh.node(t(1)),
               &a2 = mesh.node(t(2)),
               &a3 = mesh.node(t(3));

/*
    DenseVector<Array<double> > b1 = a2-a1,
                                b2 = a3-a1;
*/
    DenseVector<Array<double> > b1(2), b2(2);
    b1 = a2(1) - a1(1), a2(2) - a1(2);
    b2 = a3(1) - a1(1), a3(2) - a1(2);

    double absDetB = std::abs(det(b1,b2));

    double gamma[3];
    gamma[0] =  dot(b2, b2)/absDetB,
    gamma[1] = -dot(b1, b2)/absDetB,
    gamma[2] =  dot(b1, b1)/absDetB;

    for (int i=1; i<=3; ++i) {
        for (int j=i; j<=3; ++j) {
            Am(i,j) = gamma[0]*S[0](i,j)
                     +gamma[1]*S[1](i,j)
                     +gamma[2]*S[2](i,j);
        }
    }
    
    bm(1) = f(a1)*absDetB/6;
    bm(2) = f(a2)*absDetB/6;
    bm(3) = f(a3)*absDetB/6;
}

void
FEM::solve()
{
	std::cerr << "install pardiso or write your own solver :-)" << std::endl;
	/*
    int M = mesh.numInteriorNodes();
    DenseVector<ArrayView<double> > u_ = u(_(1,M));
    pardiso(u_, A, b);
    */
}

void
FEM::writeSolution(const char *filename)
{
    std::ostringstream pointsFile, facesFile;
    pointsFile << filename << ".points";
    facesFile << filename << ".faces";
    
    std::ofstream pointsOut(pointsFile.str().c_str());
    std::ofstream facesOut(facesFile.str().c_str());
    
    pointsOut << mesh.numNodes() << std::endl;
    for (int i=1; i<=mesh.numNodes(); ++i) {
        Node node = mesh.node(i);
        
        pointsOut << node(1) << " " << node(2) << " ";
        pointsOut << u(mesh.eta(i)) << " 0" << std::endl;
    }

    facesOut << mesh.numTriangles() << std::endl;
    for (int i=1; i<=mesh.numTriangles(); ++i) {
        Triangle t = mesh.triangle(i);
        
        facesOut << "3 " << t(1) << " " << t(2) << " " << t(3) << " ";
        facesOut << "0 0 0 0" << std::endl;
    }

}


} // namespace flens
