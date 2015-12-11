#include <flens/flens.h>
#include <fstream>
#include <iostream>

namespace flens {


class Mesh
{
    public:
        typedef GeMatrix<FullStorage<double, RowMajor> >    NodeTable;
        typedef GeMatrix<FullStorage<int, RowMajor> >       TriangleTable;

        typedef const DenseVector<ConstArrayView<int> >     Triangle;
        typedef const DenseVector<ConstArrayView<double> >  Node;

        Mesh(const char *nodesFile, const char *trianglesFile);
    
        int numTriangles() const;
        Triangle triangle(int index) const;
        
        int numNodes() const;
        Node node(int index) const;
        
        int numInteriorNodes() const;
        int nu(int index) const;
        int eta(int index) const;
    
    private:
        NodeTable                 _nodes;
        TriangleTable             _triangles;
    
        int                       _numInteriorNodes;
        DenseVector<Array<int> >  _nu, _eta;
};

class FEM
{
    public:
        typedef Mesh::Triangle   Triangle;
        typedef Mesh::Node       Node;
        typedef double (*Func)(const Node &);
        
        FEM(const char *nodesFile, const char *trianglesFile,
            const Func &f, const Func &g);
        
        void
        assemble();

        void
        compute_Am_bm(const Triangle &t);

        void
        solve();

        void
        writeSolution(const char *filename);

        Mesh                                             mesh;
        Func                                             f, g;
        DenseVector<Array<double> >                      u, b;
        SparseSymmetricMatrix<CRS<double> >              A;

        typedef GeMatrix<FullStorage<double, RowMajor> > DGeMatrix;
        DGeMatrix                                        Am;
        DenseVector<Array<double> >                      bm;
        std::vector<DGeMatrix>                           S;

};

} // namespace flens

