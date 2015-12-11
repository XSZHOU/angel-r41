namespace flens {


//== DistributedGridVector2D ===================================================

template <typename RHS>
DistributedGridVector2D &
DistributedGridVector2D::operator=(const Vector<RHS> &rhs)
{
    copy(rhs.impl(), *this);
    return *this;
}

template <typename RHS>
DistributedGridVector2D &
DistributedGridVector2D::operator+=(const Vector<RHS> &rhs)
{
    axpy(1., rhs.impl(), *this);
    return *this;
}

template <typename RHS>
DistributedGridVector2D &
DistributedGridVector2D::operator-=(const Vector<RHS> &rhs)
{
    axpy(-1., rhs.impl(), *this);
    return *this;
}

} // namespace flens
