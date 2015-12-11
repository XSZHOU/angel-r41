namespace flens {

template <typename M>
class GaussSeidel;

template <typename M>
struct TypeInfo<GaussSeidel<M> >
{
    typedef GaussSeidel<M> Impl;
    typedef double         ElementType;
};

template <typename MA>
class GaussSeidel
    : public SymmetricMatrix<GaussSeidel<MA> >
{
    public:
        typedef DenseVector<Array<double> > VectorType;
        
        GaussSeidel(const MA &_A, const VectorType &_b, double _omega)
            : A(_A), b(_b), omega(_omega)
        {
        }

        const MA         &A;
        const VectorType &b;
        double           omega;
}; 

template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
void
mv(ALPHA alpha, const GaussSeidel<MA> &A, const DenseVector<VX> &x,
   BETA beta, DenseVector<VY> &y)
{
    assert(alpha==ALPHA(1));
    assert(beta==BETA(0));
    assert(&x==&y);
    mv(A, y);
}

//----------------------------------------------------------------------------------------

template <typename VX>
void
mv(const GaussSeidel<Poisson1D> &GS, DenseVector<VX> &x)
{
    assert(x.length()==GS.b.length());

    //-- BEGIN
    // your code here
    
    /* NOTE:
        x.firstIndex()    -> first index of vector x
        x.lastIndex()     -> last index of vector x
        x.length()        -> vector length
        x(i)              -> i-th element of x
        
        GS.b              -> is the right-hand side vector of A*x = b
        GS.A.rh           -> reciprocal of h, i.e. 1/h = N+1
        GS.omega          -> relaxation parameter for gauss seidel

    */

    //-- END
}

} // namespace flens
