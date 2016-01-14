/*
Copyright (c) 2009 Sebastian Steiger, Integrated Systems Laboratory, ETH Zurich.
Comments, suggestions, criticism or bug reports are welcome: steiger@purdue.edu. 

This file is part of ANGEL, a simulator for LEDs based on the NEGF formalism.
The software is distributed under the Lesser GNU General Public License (LGPL).
ANGEL is free software: you can redistribute it and/or modify it under the terms 
of the Lesser GNU General Public License v3 or later. ANGEL is distributed
without any warranty; without even the implied warranty of merchantability or 
fitness for a particular purpose. See also <http://www.gnu.org/licenses/>.
*/
#ifndef LINEARSOLVERPARDISO_H_NEGF
#define LINEARSOLVERPARDISO_H_NEGF

#include "all.h"

#include "CSRMatrix.h"
#include "LinearSolver.h"

// external fortran function
extern "C" {
	int pardiso_(void*  pt[64],    /* pardiso internal memory pointer's */
				 int*  maxfct ,    /* max. number of factorized matrices having the SAME structure that should be kept in mem*/
				 int*  mnum,       /* number of actual matrix for the solution phase (there could be several) */
				 int*  mtype,      /* matrix type */
				 int*  phase,      /* control of the solver (fill in reduction, factorization, fwd/backwrd subst, term) */
				 int*  n,          /* size of A */
				 void* A,          /* the matrix ... in sparse format */
				 int*  prow,       /* sparse pointer to rows */
				 int*  icol,       /* colum indices of sparse entries in A */
				 int*  perm,       /* user fill in reducing permutation ... */
				 int*  nrhs,       /* number of right hand sides that need to be solved for */
				 int   iparam[64], /* pardiso control array */
				 int*  msglvl, 	   /* message level information */
				 void* b,          /* b(n,nrhs) rhs vector/matrix */
				 void* x,          /* x(n,nrhs) */
				 int*  error);     /* error */

	int pardisoinit_(
				 void* pt[64],
				 int*  mtype,
				 int   iparam[64]
	);
}

using namespace std;

namespace negf {

	/** Interface to the direct linear solver PARDISO, www.pardiso-project.org */
	class LinearSolverPardiso: public LinearSolver {

	public:
		
		// swig has a problem with default arguments, therefore we better have 2 constructors
		LinearSolverPardiso(CSRMatrix<double>* matrix_, double * rhs_, double * solution_);
		~LinearSolverPardiso();

		void solve();
	
	protected:
	
		void * pardiso_handle[64];	//!< different problem instances?
		int    pardiso_mtype;		//!< what kind of matrix
		int    pardiso_iparam[64];	//!< configuration options
		int *  pardiso_perm;		//!< permutation array
		int *  pardiso_prow;		//!< 1-based CSR array
		int *  pardiso_icol;		//!< 1-based CSR array
	};


} // end namespace

#endif /*LINEARSOLVERPARDISO_H_*/
