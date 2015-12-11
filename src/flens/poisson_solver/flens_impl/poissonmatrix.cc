#include <poisson_solver/flens_impl/poissonmatrix.h>

namespace flens {

//-- DirichletPoisson2D --------------------------------------------------------

DirichletPoisson1D::DirichletPoisson1D()
{
}

DirichletPoisson1D::DirichletPoisson1D(int _rh)
    : rh(_rh)
{
}

//-- NeumannPoisson1D --------------------------------------------------------

NeumannPoisson1D::NeumannPoisson1D()
{
}

NeumannPoisson1D::NeumannPoisson1D(int _rh)
    : rh(_rh)
{
}

//-- DirichletPoisson2D --------------------------------------------------------

DirichletPoisson2D::DirichletPoisson2D()
{
}

DirichletPoisson2D::DirichletPoisson2D(int _rh)
    : rh(_rh)
{
}

//-- NeumannPoisson2D ----------------------------------------------------------

NeumannPoisson2D::NeumannPoisson2D()
{
}

NeumannPoisson2D::NeumannPoisson2D(int _rh)
    : rh(_rh)
{
}

} // namespace flens
