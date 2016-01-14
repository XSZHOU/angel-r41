namespace flens {

class Poisson1D;

template <>
struct TypeInfo<Poisson1D>
{
    typedef Poisson1D   Impl;
    typedef double      ElementType;
};

class Poisson1D
    : public SymmetricMatrix<Poisson1D>
{
    public:
        Poisson1D() {}
    
        Poisson1D(int _rh) : rh(_rh) {}
    
        int rh;
}; 

//----------------------------------------------------------------------------------------

template <typename ALPHA, typename VX, typename BETA, typename VY>
void
mv(ALPHA alpha, const Poisson1D &A, const DenseVector<VX> &x,
   BETA beta, DenseVector<VY> &y)
{
    assert(&x!=&y);
    
    if (y.length()!=x.length()) {
        y.resize(x.length(), x.firstIndex());
    }

    int N = A.rh - 1;
    double h = 1./(N+1);
    for (int i=1; i<=N; ++i) {
        y(i) = alpha*(-x(i-1) + 2*x(i) - x(i+1))/(h*h) + beta*y(i);
    }

    //-- BEGIN
   // your code here

    
    /* NOTE:
        x.firstIndex()    -> first index of vector x
        x.lastIndex()     -> last index of vector x
        x.length()        -> vector length
        x(i)              -> i-th element of x
        
        A.rh              -> reciprocal of h, i.e. 1/h = N+1
    */

    //-- END
}

} // namespace flens
