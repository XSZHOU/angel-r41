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
#ifndef LINEARSOLVER_H_
#define LINEARSOLVER_H_

#include "all.h"

#include "CSRMatrix.h"

using namespace std;

namespace negf {

	/** Base class for various types of linear solvers.
	 *  - ILS     - parallel (OpenMP) iterative sparse solver, proprietary, very stable, does MPS reordering and incomplete LU factoring, very fast as soon as solution is close
	 *  - PARDISO - parallel (OpenMP) direct sparse solver, free, reliable, MPS reordering
	 *  - UMFPACK - direct sparse solver, open source, very reliable
	 *  - PETSC   - parallel (MPI) iterative sparse solver, open source, many options but I haven't found a stable set so far, restricted reordering options
	 *  - MUMPS   - parallel (MPI) direct sparse solver, open source w/ license
	 *  - ILUPACK - iterative sparse solver, did not work so far
	 *  - PASTIX  - direct or iterative sparse solver, did not really work so far
	 */
	class LinearSolver {

	public:
		
		LinearSolver(CSRMatrix<double>* matrix_, double * rhs_, double * solution_); //!< constructor common to all classes, just copies pointers
		
		virtual ~LinearSolver() {}

		virtual void solve() = 0;	  //!< to be implemented by derived classes
	
	protected:
				
		CSRMatrix<double> * matrix;	  //!< pointer to the matrix A
		double			  * rhs;      //!< pointer to the right-hand side b of Ax=b
		double 			  * solution; //!< pointer to the address where the solution x of Ax=b is stored
		
		int					nthreads; //!< number of threads used when doing OpenMP computations (else: 1)
		int					fidx;	  //!< stores whether the initial pcol, irow arrays are 0-based (C indices) or 1-based (fortran indices)
	};


} // end namespace

#endif /*LINEARSOLVER_H_*/

