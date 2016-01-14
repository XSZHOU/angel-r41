namespace flens {

//-- GaussSeidel ---------------------------------------------------------------

template <typename M, typename V>
GaussSeidel<M,V>::GaussSeidel(const M &_A, const V &_f)
    : A(_A), f(_f)
{
}

//-- GaussSeidelRedBlack -------------------------------------------------------

template <typename M, typename V>
GaussSeidelRedBlack<M,V>::GaussSeidelRedBlack(const M &_A, const V &_f)
    : A(_A), f(_f)
{
}

} // namespace flens
