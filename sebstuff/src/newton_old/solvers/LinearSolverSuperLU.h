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
#ifndef LINEARSOLVERSUPERLU_H_
#define LINEARSOLVERSUPERLU_H_

#ifndef NOSUPERLU

#include "all.h"

#include "CSRMatrix.h"
#include "LinearSolver.h"
#include "SRC/pdsp_defs.h"

 
using namespace std;

namespace negf {

	/** Interface to the linear solver SuperLU - discontinued */
	class LinearSolverSuperLU: public LinearSolver {

	public:
		
		LinearSolverSuperLU(CSRMatrix<double>* matrix_, double * rhs_, double * solution_);
		~LinearSolverSuperLU();

		void solve();
	
	protected: 
	
		vector<double> B_nnz;		//!< will be RHS before solution, result after solution
		vector<int>    A_row_ptr_cache;
		
		NCformat* Astore;			//!< some storage format
		DNformat* Bstore;			//!< some storage format
		SuperMatrix A, AC, L, U, B;	//!< A - original matrix, AC - permuted, B - RHS,
		vector<int> perm_c;			//!< permutation vector (columns)
		vector<int> perm_r;			//!< permutation vector (rows)
	
		pdgstrf_options_t pdgstrf_options;
		//superlumt_options_t superlumt_options;
		int      nprocs;			//!< number of pthreads for factorization
		yes_no_t refact;			//!< refactorize?
		yes_no_t usepr;				
		int      panel_size;
		int      relax;
		double   partial_pivoting;
		double   drop_tol;   
		double*  work;
		int      lwork;
		int      info;	
		trans_t  trans;
		Gstat_t  Gstat;				//!< numerical statistics
	
		// irow, pcol, nonzeros for compressed sparse column (csc) format
		// re-ordering is done in constructor
		vector<int>    csc_irow;	//!< CSC format array
		vector<int>    csc_pcol;	//!< CSC format array
		vector<double> csc_nonzeros;//!< CSC format array
	};

} // end namespace

#endif // NOSUPERLU

#endif /*LINEARSOLVERSUPERLU_H_*/
