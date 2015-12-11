#include <cassert>

#include <transitions.h>

namespace flens {
	
void
mv(Transpose trans, double alpha, const Prolongation &P, 
   const DistributedVector &x, double beta, DistributedVector &y)
{
    assert(trans==NoTrans);
    assert(alpha==1.);
    assert(beta==1.);
    
    mv(P, x, y);
}

void
mv(const Prolongation &P, const DistributedVector &u_c, DistributedVector &u)
{
    const DistributedVector::Grid &U_C   = u_c.grid;
    DistributedVector::Grid       &U     = u.grid;

    int I = 2*u_c.firstRow;
    for (int i=u_c.firstRow; i<=u_c.m; ++i, I+=2) {
        int J = 2*u_c.firstCol;
        for (int j=u_c.firstCol; j<=u_c.n; ++j, J+=2) {
            U(I-1,J-1)+=.25*U_C(i,j);  U(I-1,J)+=.5*U_C(i,j); U(I-1,J+1)+=.25*U_C(i,j);
            U(I  ,J-1)+=.50*U_C(i,j);  U(I  ,J)+=   U_C(i,j); U(I  ,J+1)+=.50*U_C(i,j);
            U(I+1,J-1)+=.25*U_C(i,j);  U(I+1,J)+=.5*U_C(i,j); U(I+1,J+1)+=.25*U_C(i,j);
        }
    }
    
    int i = u_c.m+1; 
    I = 2*i;
    int J = 2*u_c.firstCol;
    for (int j=u_c.firstCol; j<=u_c.n; ++j, J+=2) {
        U(I-1,J-1)+=.25*U_C(i,j);  U(I-1,J)+=.5*U_C(i,j); U(I-1,J+1)+=.25*U_C(i,j);
        U(I  ,J-1)+=.50*U_C(i,j);  U(I  ,J)+=   U_C(i,j); U(I  ,J+1)+=.50*U_C(i,j);
    }


    int j = u_c.n+1; 
    J = 2*j;
    I = 2*u_c.firstRow;
    for (int i=u_c.firstRow; i<=u_c.m; ++i, I+=2) {
        U(I-1,J-1)+=.25*U_C(i,j);  U(I-1,J)+=.5*U_C(i,j);
        U(I  ,J-1)+=.50*U_C(i,j);  U(I  ,J)+=   U_C(i,j);
        U(I+1,J-1)+=.25*U_C(i,j);  U(I+1,J)+=.5*U_C(i,j);
    }

    i = u_c.m+1; 
    j = u_c.n+1; 
    J = 2*j;
    I = 2*i;

    U(I-1,J-1)+=.25*U_C(i,j);  U(I-1,J)+=.5*U_C(i,j);
    U(I  ,J-1)+=.50*U_C(i,j);  U(I  ,J)+=   U_C(i,j);
    
    u.updateBoundary();    
}

//------------------------------------------------------------------------------

void
mv(Transpose trans, double alpha, const Restriction &R, 
   const DistributedVector &x, double beta, DistributedVector &y)
{
    assert(trans==NoTrans);
    assert(alpha==1.);
    assert(beta==0.);
    
    mv(R, x, y);
}

void
mv(const Restriction &R, const DistributedVector &f, DistributedVector &f_c)
{
    const DistributedVector::Grid &F   = f.grid;
    DistributedVector::Grid       &F_C = f_c.grid;
    int I = 2*f_c.firstRow;
    for (int i=f_c.firstRow; i<=f_c.m; ++i, I+=2) {
        int J = 2*f_c.firstCol;
        for (int j=f_c.firstCol; j<=f_c.n; ++j, J+=2) {
            F_C(i,j) = (                 F(I-1,J)
                        + F(I  ,J-1) + 4*F(I  ,J) + F(I  ,J+1)
                                        +F(I+1,J))/8; 
        }
    }

}
	
} // namespace flens
