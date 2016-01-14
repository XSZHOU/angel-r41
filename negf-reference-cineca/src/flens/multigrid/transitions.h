#ifndef MULTIGRID_TRANSITIONS_H
#define MULTIGRID_TRANSITIONS_H 1

#include <flens/flens.h>
#include <distributedvector.h>

namespace flens {
	
//--- Prolongation -------------------------------------------------------------

DefaultTypeInfo(Prolongation, double);

struct Prolongation 
    : public GeneralMatrix<Prolongation>
{    
};

// common 'flens format': P*x
void
mv(Transpose trans, double alpha, const Prolongation &P, 
   const DistributedVector &u_c, double beta, DistributedVector &u);

// trans, alpha and beta are fixed for Prolongation.
void
mv(const Prolongation &P, const DistributedVector &u_c, DistributedVector &u);




//--- Restriction --------------------------------------------------------------

DefaultTypeInfo(Restriction, double);

struct Restriction 
    : public GeneralMatrix<Restriction>
{    
};

// common 'flens format': R*x
void
mv(Transpose trans, double alpha, const Restriction &R, 
   const DistributedVector &u_c, double beta, DistributedVector &);

// trans, alpha and beta are fixed for Restriction. 
void
mv(const Restriction &R, const DistributedVector &f, DistributedVector &f_c);

} // namespace flens

#endif // MULTIGRID_TRANSITIONS_H
