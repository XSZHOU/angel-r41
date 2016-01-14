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
#ifndef LINEARSOLVERUMFPACK_H_
#define LINEARSOLVERUMFPACK_H_

#include "all.h"

#include "CSRMatrix.h"
#include "LinearSolver.h"

using namespace std;

namespace negf {

	/** Interface to the sequential direct linear solver UMFPACK - I'm a fan of this one. <BR>
	 *  Note that UMFPACK uses CSC and not CSR matrix format! */
	class LinearSolverUmfpack: public LinearSolver {

	public:
		
		LinearSolverUmfpack(CSRMatrix<double>* matrix_, double * rhs_, double * solution_);//!< also sets up CSR to CSC reordering
		~LinearSolverUmfpack(); //!< does nothing

		void solve();	//!< solve Ax=b
	
	protected:
	
		// irow, pcol, nonzeros for compressed sparse column (csc) format
		// re-ordering is done in constructor
		vector<int>    csc_irow;			//!< irow[ii] gives the row number of some entry (NOT column ii!)
		vector<int>    csc_pcol;			//!< pcol[ii] gives the position of the first nonzero entry of column ii in the irow and nonzeros arrays
		vector<double> csc_nonzeros;		//!< stores the matrix entries
	};

} // end namespace

#endif /*LINEARSOLVERUMFPACK_H_*/
