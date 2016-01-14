#include <poisson_solver/flens_impl/gridvector.h>

namespace flens {

//== GridVector1D ==============================================================

GridVector1D::GridVector1D()
    : rh(0)
{
}

GridVector1D::GridVector1D(int _rh, int firstIndex)
    : rh(_rh), grid(_(firstIndex, rh))
{
}

GridVector1D &
GridVector1D::operator=(double value)
{
    grid = value;
    return *this;
}

//------------------------------------------------------------------------------

std::ostream &
operator<<(std::ostream &out, const GridVector1D &v)
{
    out << v.grid;
    return out;
}

//== GridVector2D ==============================================================

GridVector2D::GridVector2D()
    : rh(0)
{
}

GridVector2D::GridVector2D(int _rh)
    : rh(_rh), grid(_(0,rh),_(0,rh))
{
}

GridVector2D &
GridVector2D::operator=(double value)
{
    grid = value;
    return *this;
}

//------------------------------------------------------------------------------

std::ostream &
operator<<(std::ostream &out, const GridVector2D &v)
{
    out << v.grid;
    return out;
}

} // namespace flens
