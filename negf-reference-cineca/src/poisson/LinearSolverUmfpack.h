#ifndef LINEARSOLVERUMFPACK_H_
#define LINEARSOLVERUMFPACK_H_

#include "all.h"

#include "CSRMatrix.h"

using namespace std;

namespace negf {

    /** Interface to the sequential direct linear solver UMFPACK - I'm a fan of this one. <BR>
     *  Note that UMFPACK uses CSC and not CSR matrix format! */
    class LinearSolverUmfpack {

    public:

        LinearSolverUmfpack(CSRMatrix<double>* matrix_, double * rhs_, double * solution_);//!< also sets up CSR to CSC reordering
        ~LinearSolverUmfpack(); //!< does nothing

        void solve();   //!< solve Ax=b

    protected:
        CSRMatrix<double> * matrix;   //!< pointer to the matrix A
        double            * rhs;      //!< pointer to the right-hand side b of Ax=b
        double            * solution; //!< pointer to the address where the solution x of Ax=b is stored

        int                 fidx;     //!< stores whether the initial pcol, irow arrays are 0-based (C indices) or 1-based (fortran indices)

        // irow, pcol, nonzeros for compressed sparse column (csc) format
        // re-ordering is done in constructor
        vector<int>    csc_irow;            //!< irow[ii] gives the row number of some entry (NOT column ii!)
        vector<int>    csc_pcol;            //!< pcol[ii] gives the position of the first nonzero entry of column ii in the irow and nonzeros arrays
        vector<double> csc_nonzeros;        //!< stores the matrix entries
    };

} // end namespace

#endif /*LINEARSOLVERUMFPACK_H_*/
