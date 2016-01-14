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
#include "LinearSolverSuperLU.h"
#ifndef NOSUPERLU
using namespace negf;

LinearSolverSuperLU::LinearSolverSuperLU(CSRMatrix<double>* matrix_, double * rhs_, double * solution_):
	LinearSolver(matrix_, rhs_, solution_)
{STACK_TRACE(
	
	const uint 		size 	= matrix->get_size();
	const uint 		nnz  	= matrix->get_num_nonzeros();
	const int*    	icol    = matrix->get_icol();		// size n+1
	const int*    	prow    = matrix->get_prow();		// size nnz
	const double* 	nonzeros= matrix->get_nonzeros();	// size nnz
	const int       fidx    = prow[0]; // 0 --> C ordering, 1 --> fortran ordering
	
	this->perm_c.resize(size, 0);
	this->perm_r.resize(size, 0);
	
	this->nprocs           = 1;
	this->usepr            = NO;
	this->panel_size       = sp_ienv(1);
	this->relax            = sp_ienv(2);
	this->partial_pivoting = 1.0;
	this->drop_tol         = 0.0;   
	this->work             = 0;
	this->lwork            = 0;
	this->info             = 0;	
	this->trans            = NOTRANS;
	
	char* tmp = getenv("OMP_NUM_THREADS");
	if(tmp != 0) {
		nprocs = atoi(tmp);
		ostringstream sout;
		sout << "SuperLU: using " << nprocs << " processors for factorization";
		logmsg->emit(LOG_INFO_L2, sout.str().c_str()); 		
	}
	
	StatAlloc(size, nprocs, panel_size, relax, &Gstat);
	StatInit(size, nprocs, &Gstat);
	
	// --------------------------------------------------------------------------------
	// re-sort csr to csc
	// new pcol, irow will always have C indices
	// --------------------------------------------------------------------------------
	// calculate elements in each colum
	vector<int> col_count;
	col_count.assign(size, 0);
	for(uint ii = 0; ii < nnz; ii++) {
		int tmp = icol[ii] - fidx;
		NEGF_ASSERT(tmp>=0 && tmp < (int)col_count.size(), "something went wrong.");
		col_count[tmp]++;
	}
	// init new space
	this->csc_irow.assign(nnz, 0);
	this->csc_pcol.assign(size + 1, 0);
	this->csc_nonzeros.assign(nnz, 0.0);
 
	// build pcol from col_count (sum up col counts, gives me number of elements in each col)
	csc_pcol[0] = /*fidx*/0;
	for(uint ii = 0; ii < size; ii++) {
		csc_pcol[ii + 1] = csc_pcol[ii] + col_count[ii] /*+ fidx*/;	
	}	
	
	vector<int> offsets(size);	
	// reorder arrays
	int csc_index = 0;		 
	for(uint rr = 0; rr < size; rr++) {
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
	
	
	// -----------------------------------------------
	// setup matrices 
	// -----------------------------------------------
	
	// set up Astore, A and A_row_ptr_cache
	Astore  = new NCformat;
	Astore->nnz    = nnz;
	Astore->nzval  = &csc_nonzeros[0];
	Astore->colptr = &csc_pcol[0];
	Astore->rowind = &csc_irow[0];
	
	A.Stype = SLU_NC;
	A.Dtype = SLU_D;
	A.Mtype = SLU_GE;
	A.ncol  = size;
	A.nrow  = size;
	A.Store = Astore; 
	
	A_row_ptr_cache.assign(A.ncol + 1, 0);
	
	// set up B_nnz, Bstore and B
	B_nnz.assign(size, 0.0);
	
	Bstore = new DNformat;	
	Bstore->lda  = A.nrow;
	Bstore->nzval = &B_nnz[0];
	
	B.Stype = SLU_DN;
	B.Dtype = SLU_D;
	B.Mtype = SLU_GE;
	B.nrow  = A.nrow;
	B.ncol  = 1;
	B.Store = Bstore;
	
	// set up L, U and AC
	L.Store = U.Store = AC.Store = 0;
);}

LinearSolverSuperLU::~LinearSolverSuperLU()
{STACK_TRACE(
    if ( lwork >= 0 ) {
        Destroy_SuperNode_SCP(&L);
        Destroy_CompCol_NCP(&U);
    }    
    pdgstrf_finalize(&pdgstrf_options, &AC);
	
	StatFree(&Gstat);
	delete Bstore; Bstore = 0;
	delete Astore; Astore = 0;
);}



void LinearSolverSuperLU::solve()
{STACK_TRACE(	
	// ------------------------
	// factorization
	// ------------------------
	// calculate permutation
	int permc_spec = 1;
	get_perm_c(permc_spec, &A, &perm_c[0]);
	
	this->refact= NO;
	
	logmsg->emit(LOG_INFO/*_L2*/, "SuperLU: factorizing matrix");
	// initialize
	/*
	fact_t fact = DOFACT; // specifies whether or not the factored form of the matrix is supplied.
	trans_t trans = NOTRANS;  // specifies the form of the system of equations: 
	// NOTRANS: A * X = B;      TRANS: A**T * X = B;       CONJ: A**H * X = B
	pdgstrf_init(nprocs, fact, trans, refact, panel_size, relax,
			partial_pivoting, usepr, drop_tol, &perm_c[0], &perm_r[0],
	 		work, lwork, &A, &AC, &superlumt_options, &Gstat);	
	 */
	pdgstrf_init(nprocs, refact, panel_size, relax,
			partial_pivoting, usepr, drop_tol, &perm_c[0], &perm_r[0],
	 		work, lwork, &A, &AC, &pdgstrf_options, &Gstat);	
	 
	// factorize	 
	//pdgstrf(&superlumt_options, &AC, &perm_r[0], &L, &U, &Gstat, &info);
	pdgstrf(&pdgstrf_options, &AC, &perm_r[0], &L, &U, &Gstat, &info);
	NEGF_FASSERT(info == 0, "info of pdgstrf returned %d", info);
	logmsg->emit(LOG_INFO_L2, "SuperLU: successfully factorized matrix");
	
	// ------------------------
	// rest
	// ------------------------
	int size = matrix->get_size();
	
	//#pragma omp parallel for schedule(static, 5000)
	for(int ii = 0; ii < size; ii++) {
		B_nnz[ii] = this->rhs[ii];
	}
	
	logmsg->emit(LOG_INFO/*_L2*/, "SuperLU: rest");
	dgstrs(trans, &L, &U, &perm_r[0], &perm_c[0], &B, &Gstat, &info);
	NEGF_FASSERT(info == 0, "info of dgstrs returned %d", info);
	
	//#pragma omp parallel for schedule(static, 5000)
	for(int ii = 0; ii < size; ii++) {
		this->solution[ii] = B_nnz[ii];		// B_nnz was overwritten with result
	}		 	
);}

#endif // NOSUPERLU

