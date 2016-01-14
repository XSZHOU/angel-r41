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
#include "LinearSolverPardiso.h"
using namespace negf;

/** - set up permutation arrays
 *  - initialize PARDISO
 *  - set PARDISO options/configuration
 *  - set PARDISO matrix structure */
LinearSolverPardiso::LinearSolverPardiso(CSRMatrix<double>* matrix_, double * rhs_, double * solution_):
	LinearSolver(matrix_, rhs_, solution_)
{STACK_TRACE(
	if(matrix->symmetric()) {
		logmsg->emit(LOG_WARN, "It is a terribly BAD IDEA to use pardiso2 for symmetric matrices! it seems broken to me!");
	}
	
	// init permutations array to zero
	this->pardiso_perm = new int[this->matrix->get_size()];
	for(uint ii = 0; ii < matrix->get_size(); ii++) {
		this->pardiso_perm[ii] = 0;
	}
	// init handle and iparam
	for(int ii = 0; ii < 64; ii++) {
		this->pardiso_handle[ii] = 0;
		this->pardiso_iparam[ii] = 0;
	}
	this->pardiso_mtype = 11;	// real, unsymmetric
	
	pardisoinit_(this->pardiso_handle, &(this->pardiso_mtype), this->pardiso_iparam);
	
	this->pardiso_iparam[0]  = 1;  			// use custom iparams
	this->pardiso_iparam[1]  = /*2*/0;  	// fill-in reduction reordering. 2-->metis DOES NOT WORK with qw-3d_tensor example!
	this->pardiso_iparam[2]  = nthreads;  	// number of threads
	this->pardiso_iparam[3]  = 0;  			// don't do preconditioned cgs
	this->pardiso_iparam[4]  = 0;  			// no user permutation
	this->pardiso_iparam[5]  = 0;  			// write solution to X
	//this->pardiso_iparam[6]   			// (output) number of iterative refinement steps used
	this->pardiso_iparam[7]  = /*1*/100;  	// max numbers of iterative refinement steps
	// this->pardiso_iparam[8]				// not used, must be 0
	this->pardiso_iparam[9]  = 13;			// eps pivot (nonsymmetric->13, symmetric->8)
	this->pardiso_iparam[10] = 1/*0*/;			// use (non-)symmetric scaling vectors
	//this->pardiso_iparam[11] = 0;
	this->pardiso_iparam[12] = 1/*0*/;			// improved accuracy using (non-)symmetric matchings
	// this->pardiso_iparam[13]				// number of perturbed pivots
	// this->pardiso_iparam[14]				// 
	// this->pardiso_iparam[15]				// 
	// this->pardiso_iparam[16]				// 
	this->pardiso_iparam[17] = -1;			// number of nonzeros will only be determined when given as -1
	this->pardiso_iparam[18] = -1;			// mflops for LU factorization only determined when given as -1
	// this->pardiso_iparam[19]				// number of CG iterations
	this->pardiso_iparam[20] = 1;			// pivoting for symmetric indefinite matrices
	// this->pardiso_iparam[21]				// number of positive eigenvalus for symmetric indefinite systems
	// this->pardiso_iparam[22]				// number of negative eigenvalus for symmetric indefinite systems
	this->pardiso_iparam[23] = 1;			
	this->pardiso_iparam[24] = 1;
);}


LinearSolverPardiso::~LinearSolverPardiso()
{STACK_TRACE(
	int phase = -1;
		
	int maxfct    = 1;  // keep only single matrix --> only one numerical factorization in memory
	int mnum      = 1;	// only one matrix
	int n 		  = matrix->get_size();
	int nrhs 	  = 1;  // only 1 RHS to be solved
	int msg_level = 0;	// 0 - nothing, 1 - statistical info
	int error     = 0;	// just init with something
	
	double * dummy_soln = new double[n];
	double * dummy_rhs  = new double[n];
	
	pardiso_(this->pardiso_handle, &maxfct, &mnum, &(this->pardiso_mtype), &phase, &n,
			 (void*)matrix->get_nonzeros(), this->pardiso_prow, this->pardiso_icol, this->pardiso_perm, &nrhs,
			 this->pardiso_iparam, &msg_level, (void*)dummy_rhs, (void*)dummy_soln, &error);
		 
	delete [] dummy_soln;
	delete [] dummy_rhs;
	delete [] this->pardiso_perm;
	if (this->pardiso_prow != matrix->get_prow()) {
		delete [] pardiso_prow;
	}
	if (this->pardiso_icol != matrix->get_icol()) {
		delete [] pardiso_icol;
	}
);}

void LinearSolverPardiso::solve()
{STACK_TRACE(
	// pardiso_perm, pardiso_iparam, pardiso_handle --> see constructor
	
	int maxfct    = 1;  // keep only single matrix --> only one numerical factorization in memory
	int mnum      = 1;	// only one matrix
	int n 		  = matrix->get_size();
	int nnz 	  = matrix->get_num_nonzeros();
	int nrhs 	  = 1;  // only 1 RHS to be solved
	int msg_level = 0;	// 0 - nothing, 1 - statistical info
	int error     = 0;	// just init with something

	// set up pardiso_prow, pardiso_icol
	// check if matrix is built using fortran ordering
	if(matrix->get_prow()[0] == 0) { //C ordering
		logmsg->emit(LOG_INFO_L2, "Pardiso: mapping to fortran indices");
		// create new arrays containing the modified indices
		int * matrix_icol =  matrix->get_icol();
		this->pardiso_icol = new int[nnz];
		for(int ii = 0; ii < nnz; ii++) {
			this->pardiso_icol[ii] = matrix_icol[ii] + 1;
		}
		int * matrix_prow =  matrix->get_prow();
		this->pardiso_prow = new int[n+1];
		for(int ii = 0; ii < n + 1; ii++) {
			this->pardiso_prow[ii] = matrix_prow[ii] + 1;
		}
	} else {
		// we already have the right ordering, do not copy
		this->pardiso_prow = matrix->get_prow();
		this->pardiso_icol = matrix->get_icol();
	}
	
	bool continue_after_error = true;
	
	bool altogether = true;
	if (altogether) {			
		int phase = 13; // analysis + numerical factorization + solve + iterative refinement
		
		bool scale_lines_individually = true;
		if (scale_lines_individually)
		{
			// ------------------------------------------------------------------------
			// instead of Ax=b solve DAx=Db, D a diagonal matrix s.th. (DA)_ii=1.0
			// therefore D_ii = 1.0/A_ii, (Db)_i = b_i / A_ii
			// ------------------------------------------------------------------------
			double * vals = matrix->get_nonzeros();
			// 1. create a vector with scalings for each line
			double * scalings = new double[n];
			for (int ii=0; ii<n; ii++) {
				int jdiag = -1;	// will store position of diagonal element of row i in icol and vals vectors
				double maxelem = 0.0;
				for (int jj=pardiso_prow[ii]; jj < pardiso_prow[ii+1]; jj++) {
					/*
					if (pardiso_icol[jj-1]-1==ii) {	// note that pardiso_prow and pardiso_icol are normal C vectors and therefore have C indices
						jdiag = jj-1;
						break;
					}*/
					if (fabs(vals[jj-1]) > maxelem) {
						jdiag = jj-1;
						maxelem = fabs(vals[jj-1]);
					}
				}
				if (jdiag==-1) {
					/*
					logmsg->emit(LOG_INFO,"ii=%d: diagonal does not appear in sparsity pattern. setting scaling to 1.0",ii);
					scalings[ii] = 1.0;
					*/
					NEGF_EXCEPTION("Seems like an entire matrix row was zero.");
				} else {
					/*
					if (vals[jdiag]==0.0) {
						logmsg->emit(LOG_INFO,"ii=%d: zero diagonal element encountered. setting scaling to 1.0",ii);
						scalings[ii] = 1.0;
					} else {
						scalings[ii] = vals[jdiag];
						if (scalings[ii]>1e8) scalings[ii] = 1e8;
						if (scalings[ii]<-1e8) scalings[ii] = -1e8;
						logmsg->emit(LOG_INFO,"ii=%d: setting scaling to 1.0/%e", ii, scalings[ii]);
					}
					*/
					scalings[ii] = maxelem;
					if (scalings[ii]>1e10) scalings[ii] = 1e10;
					//logmsg->emit(LOG_INFO,"ii=%d: setting scaling to 1.0/%e", ii, scalings[ii]);
				}
			}
			// 2. scale jacobian and RHS
			for (int ii=0; ii<n; ii++) {
				rhs[ii] /= scalings[ii];
				for (int jj=pardiso_prow[ii]; jj < pardiso_prow[ii+1]; jj++) {
					vals[jj-1] /= scalings[ii];
				}
			}
			// save matrix to file
			//matrix->save_to_file("scaled_matrix.csr");
			// 3. solve
			pardiso_(this->pardiso_handle, &maxfct, &mnum, &(this->pardiso_mtype), &phase, &n,
						 (void*)vals, this->pardiso_prow, this->pardiso_icol, this->pardiso_perm, &nrhs,
						 this->pardiso_iparam, &msg_level, (void*)this->rhs, (void*)(this->solution), &error);
			// 4. scale back jacobian and RHS
			for (int ii=0; ii<n; ii++) {
				rhs[ii] *= scalings[ii];
				for (int jj=pardiso_prow[ii]; jj < pardiso_prow[ii+1]; jj++) {
					vals[jj-1] *= scalings[ii];
				}
			}
			
			delete [] scalings;
		} else 
		{
			// ------------------------------------------------------------------------
			// try scaling the whole RHS
			// ------------------------------------------------------------------------
			double max_rhs = 0.0;
			for (int ii=0; ii<n; ii++) {
				max_rhs = max(max_rhs, std::fabs(rhs[ii]));
			}
			NEGF_ASSERT(max_rhs > 0.0, "Entire RHS was 0.0!!!!");
			
			double scale = 1.0;
			double * new_rhs = new double[n];
			
			int iter = 0;
			do {
				scale = negf_math::pow(10.0, iter);
				logmsg->emit(LOG_INFO_L2,  "solving with RHS scaling=%e", scale);
				for (int ii = 0; ii < n; ii++) {
					new_rhs[ii] = rhs[ii] * scale;
				}
				pardiso_(this->pardiso_handle, &maxfct, &mnum, &(this->pardiso_mtype), &phase, &n,
						 (void*)matrix->get_nonzeros(), this->pardiso_prow, this->pardiso_icol, this->pardiso_perm, &nrhs,
						 this->pardiso_iparam, &msg_level, (void*)new_rhs, (void*)(this->solution), &error);
				if (error!=0) logmsg->emit(LOG_INFO_L1,"Pardiso reported error %d with RHS scaling %e.",error,scale);
				iter = iter + 2;
			} while (error!=0 && iter <= 20);
			if (error!=0) {
				iter = 1;
				do {
					scale = negf_math::pow(10.0, -iter);
					logmsg->emit(LOG_INFO_L2,  "solving with RHS scaling=%e", scale);
					for (int ii = 0; ii < n; ii++) {
						new_rhs[ii] = rhs[ii] * scale;
					}
					pardiso_(this->pardiso_handle, &maxfct, &mnum, &(this->pardiso_mtype), &phase, &n,
							 (void*)matrix->get_nonzeros(), this->pardiso_prow, this->pardiso_icol, this->pardiso_perm, &nrhs,
							 this->pardiso_iparam, &msg_level, (void*)new_rhs, (void*)(this->solution), &error);
					if (error!=0) logmsg->emit(LOG_INFO_L1,"Pardiso reported error %d with RHS scaling %e.",error,scale);
					iter = iter + 2;
				} while (error!=0 && iter <= 10);
				if (error!=0) {
					logmsg->emit(LOG_INFO,"saving matrix to pardiso_error.csr...");
					this->matrix->save_to_file("pardiso_error.csr");
					int old_error = error;
					cout << "solving dummy RHS for testing..." << endl;
					for (int ii = 0; ii < n; ii++) {
						new_rhs[ii] = 1.0;
					}
					pardiso_(this->pardiso_handle, &maxfct, &mnum, &(this->pardiso_mtype), &phase, &n,
							 (void*)matrix->get_nonzeros(), this->pardiso_prow, this->pardiso_icol, this->pardiso_perm, &nrhs,
							 this->pardiso_iparam, &msg_level, (void*)new_rhs, (void*)(this->solution), &error);
					logmsg->emit(LOG_INFO,"Dummy RHS [1,1,...] gave error code %d",error);
					error = old_error;
					if (continue_after_error) {
						cout << "solving w/ scaling 1 and continuing." << endl;
						scale = 1.0;
						for (int ii = 0; ii < n; ii++) {
							new_rhs[ii] = rhs[ii];
						}
						pardiso_(this->pardiso_handle, &maxfct, &mnum, &(this->pardiso_mtype), &phase, &n,
								 (void*)matrix->get_nonzeros(), this->pardiso_prow, this->pardiso_icol, this->pardiso_perm, &nrhs,
								 this->pardiso_iparam, &msg_level, (void*)new_rhs, (void*)(this->solution), &error);
					}
				}
			}
			logmsg->emit(LOG_INFO_L2,  "scaling back solution");
			for (int ii = 0; ii < n; ii++) {
				solution[ii] = solution[ii] / scale;
			}
			
			delete [] new_rhs;
		}
		std::ostringstream sout;
		sout << "pardiso statistics:\n"
		     << "error:                        " << error << "\n"
		     << "number iterative refin steps: " << this->pardiso_iparam[7] << "\n"
			 << "number pertubed pivots:       " << this->pardiso_iparam[13] << "\n"
			 << "peak memory:                  " << this->pardiso_iparam[14] << "\n"
			 << "permanent memory:             " << this->pardiso_iparam[15] << "\n"
			 << "fact memory:                  " << this->pardiso_iparam[16] << "\n"
			 << "nonzeros:                     " << this->pardiso_iparam[17] << "\n";
		logmsg->emit(LOG_INFO_L2, sout.str().c_str());
		if(error != 0) 
		{
			string err("pardiso error: ");
			switch(error) {
				case -1: err.append("input inconsistent");	break;
				case -2: err.append("not enough memory");	break;
				case -3: err.append("reordering problem");	break;
				case -4: err.append("zero pivot numerical fact. or iterative refinement problem");	break;
				case -5: err.append("unclassified (internal) error"); break;
				case -6: err.append("preordering failed");	break;
				case -7: err.append("diagonal matrix problem");	break;
				default: err.append("unknown error code");
			}
			if (!continue_after_error) NEGF_EXCEPTION(err.c_str());
		}
	} else {
		// --------------------------------------
		// LU factorization
		// --------------------------------------
		int phase = 12; // analysis + numerical factorization
		
		pardiso_(this->pardiso_handle, &maxfct, &mnum, &(this->pardiso_mtype), &phase, &n,
				 (void*)matrix->get_nonzeros(), this->pardiso_prow, this->pardiso_icol, this->pardiso_perm, &nrhs,
				 this->pardiso_iparam, &msg_level, (void*)(this->rhs), (void*)(this->solution), &error);	// solution --> meaningless
		std::ostringstream sout;
		sout << "pardiso statistics:\n"
		     << "error:                        " << error << "\n"
		     << "number iterative refin steps: " << this->pardiso_iparam[7] << "\n"
			 << "number pertubed pivots:       " << this->pardiso_iparam[13] << "\n"
			 << "peak memory:                  " << this->pardiso_iparam[14] << "\n"
			 << "permanent memory:             " << this->pardiso_iparam[15] << "\n"
			 << "fact memory:                  " << this->pardiso_iparam[16] << "\n"
			 << "nonzeros:                     " << this->pardiso_iparam[17] << "\n";
		logmsg->emit(LOG_INFO_L2, sout.str().c_str());
		if(error != 0) 
		{
			string err("pardiso says: ");
			switch(error) {
				case -1: err.append("input inconsistent");	break;
				case -2: err.append("not enough memory");	break;
				case -3: err.append("reordering problem");	break;
				case -4: err.append("zero pivot numerical fact. or iterative refinement problem");	break;
				case -5: err.append("unclassified (internal) error"); break;
				case -6: err.append("preordering failed");	break;
				case -7: err.append("diagonal matrix problem");	break;
				default: err.append("unknown error code");
			}
			NEGF_EXCEPTION(err.c_str());
		}
		
		// --------------------------------------
		// backward substitution
		// --------------------------------------
		phase = 33;
		pardiso_(this->pardiso_handle, &maxfct, &mnum, &(this->pardiso_mtype), &phase, &n,
				 (void*)matrix->get_nonzeros(), this->pardiso_prow, this->pardiso_icol, this->pardiso_perm, &nrhs,
				 this->pardiso_iparam, &msg_level, (void*)(this->rhs), (void*)(this->solution), &error);
		//cout << "iparam[6]=" << this->pardiso_iparam[6] << endl;
		if(error != 0) 
		{
			string err("pardiso says: ");
			switch(error) {
				case -1: err.append("input inconsistent");	break;
				case -2: err.append("not enough memory");	break;
				case -3: err.append("reordering problem");	break;
				case -4: err.append("zero pivot numerical fact. or iterative refinement problem");	break;
				case -5: err.append("unclassified (internal) error"); break;
				case -6: err.append("preordering failed");	break;
				case -7: err.append("diagonal matrix problem");	break;
				default: err.append("unknown error code");
			}
			NEGF_EXCEPTION(err.c_str());
		}
	}
	logmsg->emit(LOG_INFO_L2,"pardiso finished.");
);}
