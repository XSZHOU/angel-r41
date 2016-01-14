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
#include "LinearSolverUmfpack.h"

#include "Include/umfpack.h"

using namespace negf;
/*
umfpack_*_symbolic:
Pre-orders the columns of A to reduce fill-in. Returns an opaque Symbolic object as a
void * pointer. The object contains the symbolic analysis and is needed for the numeric
factorization. This routine requires only O(|A|) space, where |A| is the number of nonzero
entries in the matrix. It computes upper bounds on the nonzeros in L and U, the floatingpoint
operations required, and the memory usage of umfpack * numeric. The Symbolic
object is small; it contains just the column pre-ordering, the supernodal column elimination
tree, and information about each frontal matrix. It is no larger than about 13n integers if A
is n-by-n.

umfpack_*_numeric:
Numerically scales and then factorizes a sparse matrix into PAQ, PRAQ, or PR-1AQ into
the product LU, where P and Q are permutation matrices, R is a diagonal matrix of scale
factors, L is lower triangular with unit diagonal, and U is upper triangular. Requires the
symbolic ordering and analysis computed by umfpack_*_symbolic or umfpack_*_qsymbolic.
Returns an opaque Numeric object as a void * pointer. The object contains the numerical
factorization and is used by umfpack_*_solve. You can factorize a new matrix with a
different values (but identical pattern) as the matrix analyzed by umfpack * symbolic or
umfpack * qsymbolic by re-using the Symbolic object (this feature is available when using
UMFPACK in a C or Fortran program, but not in MATLAB). The matrix U will have zeros
on the diagonal if A is singular; this produces a warning, but the factorization is still valid.

umfpack_*_solve:
Solves a sparse linear system (Ax = b, ATx = b, or systems involving just L or U), using
the numeric factorization computed by umfpack_*_numeric. Iterative refinement with sparse
backward error [3] is used by default. The matrix A must be square. If it is singular,
then a divide-by-zero will occur, and your solution with contain IEEE Infs or Nans in the
appropriate places.

umfpack_*_free_symbolic:
Frees the Symbolic object created by umfpack * symbolic or umfpack * qsymbolic.

umfpack_*_free_numeric:
Frees the Numeric object created by umfpack * numeric.

Be careful not to free a Symbolic object with umfpack_*_free numeric. Nor should you attempt
to free a Numeric object with umfpack_*_free symbolic. Failure to free these objects will lead to
memory leaks.


di, dl: take doubles and return doubles
zi, zl: first arrays hold Re, second arrays hold Im
i --> int, l --> long

Alternatives:
-------------
defaults - ets default parameters
qsymbolic - user specifies own column preordering
wsolve - does not dynamically allocate memory

Control Parameters:
-------------------					default
Control[UMFPACK_PRL] 				1 			printing level
Control[UMFPACK_DENSE_ROW] 			0.2 		dense row parameter
Control[UMFPACK_DENSE_COL] 			0.2 		dense column parameter
Control[UMFPACK_PIVOT_TOLERANCE] 	0.1 		partial pivoting tolerance
Control[UMFPACK_BLOCK_SIZE] 		32 			BLAS block size
Control[UMFPACK_STRATEGY] 			0 (auto) 	select strategy
Control[UMFPACK_ALLOC_INIT] 		0.7 		initial memory allocation
Control[UMFPACK_IRSTEP] 			2 			max iter. refinement steps
Control[UMFPACK_2BY2_TOLERANCE] 	0.01 		defines 'large' entries
Control[UMFPACK_FIXQ] 				0 (auto) 	fix or modify Q
Control[UMFPACK_AMD_DENSE] 			10 			AMD dense row/column parameter
Control[UMFPACK_SYM_PIVOT_TOLERANCE] 0.001 		for diagonal entries
Control[UMFPACK_SCALE] 				1 (sum) 	row scaling (none, sum, or max)
Control[UMFPACK_FRONT_ALLOC_INIT] 	0.5			frontal matrix allocation ratio
Control[UMFPACK_DROPTOL] 			0 			drop tolerance
Control[UMFPACK_AGGRESSIVE] 		1 (yes) 	aggressive absorption in AMD and COLAMD
*/



LinearSolverUmfpack::LinearSolverUmfpack(CSRMatrix<double>* matrix_, double * rhs_, double * solution_):
	LinearSolver(matrix_, rhs_, solution_)
{STACK_TRACE(
	this->csc_irow.clear();
	this->csc_pcol.clear();
	this->csc_nonzeros.clear();
);}


LinearSolverUmfpack::~LinearSolverUmfpack()
{STACK_TRACE(
);}

void LinearSolverUmfpack::solve()
{//STACK_TRACE(
	double * null = (double *) NULL ;
	void * Symbolic, * Numeric ;
	
	int *    prow = matrix->get_prow();
	int *    icol = matrix->get_icol();
	double * nonzeros = matrix->get_nonzeros();
	uint      n = matrix->get_size();
	uint     nnz = matrix->get_num_nonzeros();
	double * b  = this->rhs;
	double * x  = this->solution;
	//int fidx = prow[0]; // 0 --> C ordering, 1 --> fortran ordering
	// this->fidx !!!
	
	
	// --------------------------------------------------------------------------------
	// re-sort csr to csc
	// new pcol, irow will always have C indices
	// --------------------------------------------------------------------------------
	// calculate elements in each colum
	vector<int> col_count;
	col_count.assign(n, 0);
	for(uint ii = 0; ii < nnz; ii++) {
		int tmp = icol[ii] - fidx;
		NEGF_ASSERT(tmp>=0 && tmp < (int)col_count.size(), "something went wrong.");
		col_count[tmp]++;
	}
	// init new space
	this->csc_irow.assign(nnz, 0);
	this->csc_pcol.assign(n + 1, 0);
	this->csc_nonzeros.assign(nnz, 0.0);
 
	// build pcol from col_count (sum up col counts, gives me number of elements in each col)
	csc_pcol[0] = /*fidx*/0;
	for(uint ii = 0; ii < n; ii++) {
		csc_pcol[ii + 1] = csc_pcol[ii] + col_count[ii] /*+ fidx*/;	
	}	
	
	vector<int> offsets(n);	
	// reorder arrays
	int csc_index = 0;		 
	for(uint rr = 0; rr < n; rr++) {
		for(int ii = prow[rr]; ii < prow[rr + 1]; ii++) {					
			// (rr,icol[ii]) is the element I have now
			// to save this in csc, I do the following:
			
			// 1. get starting index of column icol[ii] in csc	
			int tmp = icol[ii-fidx]-fidx;
			csc_index = csc_pcol[tmp]/*-fidx*/ + offsets[tmp];	
					
			// 2. store values at starting index + offset
			csc_nonzeros[csc_index] = nonzeros[ii-fidx];
			
			// 3. store current row in irow at index + offset
			csc_irow[csc_index] = rr /*+ fidx*/;
			
			// 4. increase offset
			offsets[tmp]++;		
		} 
	}
	
	/*
	// umfpack needs C indices
	int fidx = Ap[0];
	if (fidx==1) {
		for (uint ii=0; ii < n+1; ii++) {
			Ap[ii]--;
		}
		for (uint ii=0; ii < nz; ii++) {
			Ai[ii]--;
		}
	}*/
	
	int * Ap = &csc_pcol[0];
	int * Ai = &csc_irow[0];
	double * Ax = &csc_nonzeros[0];
	
	(void) umfpack_di_symbolic (n, n, Ap, Ai, Ax, &Symbolic, null, null) ;
	(void) umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, null, null) ;
	umfpack_di_free_symbolic (&Symbolic) ;
	(void) umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null) ;
	umfpack_di_free_numeric (&Numeric) ;
	
	/*
	// restore original indices
	if (fidx==1) {
		for (uint ii=0; ii < n+1; ii++) {
			Ap[ii]++;
		}
		for (uint ii=0; ii < nz; ii++) {
			Ai[ii]++;
		}
	}
	*/
//);}
	}
