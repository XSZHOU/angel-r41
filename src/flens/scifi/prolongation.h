namespace flens {

class Prolongation;

template <>
struct TypeInfo<Prolongation>
{
    typedef Prolongation Impl;
    typedef double       ElementType;
};

class Prolongation
    : public GeneralMatrix<Prolongation>
{
};

template <typename ALPHA, typename VX, typename BETA, typename VY>
void
mv(Transpose trans, ALPHA alpha, const Prolongation &A, const DenseVector<VX> &x,
   BETA beta, DenseVector<VY> &y)
{
    assert(trans==NoTrans);
    assert(alpha==ALPHA(1));
    assert(beta==BETA(1));
    
    mv(A, x, 1, y);
}

//------------------------------------------------------------------------------

template <typename VX, typename BETA, typename VY>
void
mv(const Prolongation &R, const DenseVector<VX> &x, BETA beta, DenseVector<VY> &y)
{
    assert(beta==BETA(1));

    int n = x.lastIndex()-1;

    assert(y.firstIndex()<=1);
    assert(y.lastIndex()>=2*n+1);
    
    //-- BEGIN
    // your code here
    
    //-- END
}

} // namespace flens
