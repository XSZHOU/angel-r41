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
#include "LinearSolverMumps.h"

#ifndef NOMUMPS
using namespace negf;

#define MUMPS_USE_COMM_WORLD -987654 // when id.comm_fortran is set to this value. MPI_COMM_WORLD is used as MPI communicator
#define ICNTL(I) icntl[(I)-1]        // macro s.th. indices match documentation
#define INFO(I) info[(I)-1]          // macro s.th. indices match documentation
 
MumpsSlave::MumpsSlave() throw (Exception *):
	root(0)
{STACK_TRACE(	
	logmsg->emit_all(LOG_INFO_L3,"p%d MumpsSlave()", mpi->get_rank());

	// clear
	this->n    = 0; 	
	this->nnz  = 0; 
	this->irn  = NULL;
	this->jcn  = NULL;
	this->a    = NULL;

	if (mpi->get_rank()!=root) {
		this->prepare();
	} else {
		// prepare() for master is called by LinearSolverMumps constructor
		// because n, nnz, irn, jcn, a, mumps_rhs memory locations need to be assigned first
	}
	logmsg->emit_all(LOG_INFO_L3,"p%d MumpsSlave() is done", mpi->get_rank());
);}


MumpsSlave::~MumpsSlave()
{STACK_TRACE(
	logmsg->emit_all(LOG_INFO,"\n\np%d ~MUMPSSLAVE\n\n", mpi->get_rank());
	
	// terminate MUMPS instance
	this->terminate();
);}


void MumpsSlave::terminate()
{STACK_TRACE(
	logmsg->emit(LOG_INFO_L3,"   termination");
	id.job          = -2;	// finalize
	dmumps_c(&id);	
	this->check_error(id.INFO(1));
);}


void MumpsSlave::prepare()
{STACK_TRACE(
	// -----------------------------------------------------------
	// Initialize MUMPS, even without knowledge of the matrix
	// -----------------------------------------------------------
	logmsg->emit(LOG_INFO_L3,"   initialization");
	id.par          = 1;	// host also takes part in factorization and solve phases
	id.sym          = 0;	// unsymmetric matrix
	id.comm_fortran = MUMPS_USE_COMM_WORLD;
	id.job          = -1;	// initialize
	dmumps_c(&id);

	// --------------------------
	// Define problem on host
	// --------------------------
	if (mpi->get_rank()==root) {
		id.n      = this->n;
		id.nz     = this->nnz;
		id.irn    = this->irn;
		id.jcn    = this->jcn;
		//id.a      = this->a;
		//id.rhs    = this->mumps_rhs;
	} else {
		id.n      = 0;
		id.nz     = 0;
		id.irn    = NULL;
		id.jcn    = NULL;
		//id.a      = NULL;
		//id.rhs    = NULL;
	}

	// --------------------------
	// communicate n, nnz
	// --------------------------
	//mpi->broadcast(  this->n, root);
	//mpi->broadcast(this->nnz, root);

	// -------------------------
	// MUMPS configuration
	// -------------------------
	logmsg->emit(LOG_INFO_L3,"   configuration");

	id.ICNTL(1)     =  6;	// output stream for error messages. negative or 0 --> suppressed
	id.ICNTL(2)     =  /*6*/-1;	// output stream for diagnostic printing, statistics, and warnings
	id.ICNTL(3)     =  /*6*/-1;	// output stream for global information on host
	id.ICNTL(4)     =  2;	// level of printing error. <=0: no msgs. 1: only errors. 2: errors, warnings and main statistics. 3, 4, 5
	id.ICNTL(5)     =  0;   // (default 0) 0 - matrix is given in assembled format. 1 - different
	id.ICNTL(6)     =  7;	// (default 7 = automatic choice) option for permuting and scaling the matrix. 0 - no column permutation. 1-6: different variants
	id.ICNTL(7)     =  7;	// (default 7 = automatic choice) determines pivot order. 0...6
	id.ICNTL(8)     = 77;	// (default 77) scaling strategy
	id.ICNTL(9)     =  1;   // (default 1) 1 = solve Ax=b. otherwise solve A^T x = b
	id.ICNTL(10)    =  0;   // (default 0) max. # iterative refinement steps (if NRHS=1)
	id.ICNTL(11)    =  0;   // (default 0) a positive value returns statistics on the linear system
	id.ICNTL(12)    =  0;   // (default 0) only relevant for symmetric matrices
	id.ICNTL(13)    =  0;   // (deflaut 0) usage of scalapack
	id.ICNTL(14)    = 20;   // (default 20) percentage increase in estimated working phase
	id.ICNTL(18)    =  0;   // (default 0) strategy for the distributed input matrix
	id.ICNTL(19)    =  0;   // (default 0) Schur complement matrix
	id.ICNTL(20)    =  0;   // (default 0) 0 - RHS is in dense forn in id.rhs. 1 - RHS is in sparse form in other fields
	id.ICNTL(21)    =  0;   // (default 0) 0 - solution is stored in RHS of host process. 1 - solution is kept distributed
	id.ICNTL(22)    =  0;   // (default 0) in-core / out-of-core
	id.ICNTL(23)    =  0;   // (default 0) maximum size of working memory in MB
	id.ICNTL(24)    =  0;   // (default 0) detection of null pivot rows
);}


void MumpsSlave::slave_solve() throw (Exception *)
{STACK_TRACE(
	LoggerLevel level = logmsg->get_level();
	//logmsg->set_level(LOG_INFO_L3);
	logmsg->emit_all(LOG_INFO_L3,"p%d MumpsSlave::slave_solve()", mpi->get_rank());
	//cout << "p" << mpi->get_rank() << " slave_solve()." << endl;
	//this->test_operation();

	if (mpi->get_rank()==root) {
		id.a      = this->a;
		id.rhs    = this->mumps_rhs;
	} else {
		id.a      = NULL;
		id.rhs    = NULL;
	}

	// solve!
	//id.job          = 6;	// performs analysis (job=1), factorization (job=2) and solution (job=3) in a row
	//dmumps_c(&id);
	logmsg->emit(LOG_INFO_L3,"   analysis");
	id.job          = 1;	// performs analysis
	dmumps_c(&id);
	this->check_error(id.INFO(1));
	logmsg->emit(LOG_INFO_L3,"   factorization");
	id.job          = 2;	// performs factorization
	dmumps_c(&id);
	this->check_error(id.INFO(1));
	logmsg->emit(LOG_INFO_L3,"   solution");
	id.job          = 3;	// performs solution
	dmumps_c(&id);
	this->check_error(id.INFO(1));
	
	logmsg->set_level(level);

	mpi->synchronize_processes();
	logmsg->emit_all(LOG_INFO_L3,"p%d slave_solve() is done.", mpi->get_rank());
);}


/** Here an endless loop is entered in which the object waits for MPI signals. <BR>
 *  Depending on the signal slaves are created, deleted or solve() operations are called. <BR>
 *  With signal -1 the constructor is exited. */
MumpsSlaveCollector::MumpsSlaveCollector() throw (Exception *):
	max_num_slaves(100)
{STACK_TRACE(
	// Initialize Mumps
	
	this->slaves.resize(max_num_slaves, NULL);
	
	int root = 0;
	if (mpi->get_rank()!=root) {
		// endless loop
		while(true) {
			// receive command from master thread
			int command = 0;
			mpi->broadcast(command, root);
			
			if (command == -1) {
				// terminate
				break;
			}
			
			int num = 0;
			mpi->broadcast(num, root);
			
			if (command == 1) {	// set up new MumpsSlave with number given by master process
				this->create_slave(num);
			}			
			if (command == 2) {	// delete MumpsSlave with number given by master process
				this->delete_slave(num);
			}			
			if (command == 3) {	// solve linear problem, number given by master process
				this->solve(num);
			}
		}
	}
);}


MumpsSlaveCollector::~MumpsSlaveCollector()
{STACK_TRACE(
	for (uint ii=0; ii<slaves.size(); ii++) {
		if (slaves[ii]!=NULL) {
			delete slaves[ii];
		}
	}
	
	// Finalize Mumps
);}

void MumpsSlaveCollector::create_slave(uint num)
{STACK_TRACE(
	NEGF_ASSERT(num < slaves.size() && slaves[num]==NULL, "could not create MumpsSlave.");
	slaves[num] = new MumpsSlave();
	constants::mumps_counter++;
);}

void MumpsSlaveCollector::delete_slave(uint num)
{STACK_TRACE(
	NEGF_ASSERT(num < slaves.size() && slaves[num]!=NULL, "could not delete MumpsSlave.");
	delete slaves[num]; 
	slaves[num] = NULL;
	// mumps_counter is NOT decreased
);}

MumpsSlave * MumpsSlaveCollector::get_slave(uint num)
{STACK_TRACE(
	NEGF_ASSERT(num < slaves.size() && slaves[num]!=NULL, "invalid MumpsSlave access.");
	return slaves[num];
);}

void MumpsSlaveCollector::solve(uint num)
{STACK_TRACE(
	NEGF_ASSERT(num < slaves.size() && slaves[num]!=NULL, "could not solve MumpsSlave.");
	slaves[num]->slave_solve();
);}


LinearSolverMumps::LinearSolverMumps(CSRMatrix<double>* matrix_, double * rhs_, double * solution_) throw (Exception *):
	MumpsSlave(),
	LinearSolver(matrix_, rhs_, solution_)
{STACK_TRACE(
	this->n   = matrix->get_size();
	this->nnz = matrix->get_num_nonzeros();
	logmsg->emit(LOG_INFO,"LinearSolverMumps() problem %d: matrix size %d, %d nonzeros", constants::mumps_counter, n, nnz);
	
	int * prow = matrix->get_prow();
	int * icol = matrix->get_icol();
	this->fidx = prow[0];

	// MUMPS needs indices in XY format, NOT CSR Format!
	this->irn = new int[nnz];
	int counter=0;
	for (int ii=0; ii<n; ii++) {
		for (int jj=prow[ii]-fidx; jj<prow[ii+1]-fidx; jj++) {
			NEGF_FASSERT(counter<nnz, "counter<nnz failed: counter=%d, nnz=%d, ii=%d", counter, nnz, ii);
			this->irn[counter] = ii+1;
			counter++;
		}
	}
	NEGF_ASSERT(counter==nnz, "counter==nnz failed.");
	if (fidx == 1) {	// matrix already has fortran-based indices
		this->jcn = icol;
	} else {
		logmsg->emit(LOG_INFO,"Creating 1-based irn, jcn");
		this->jcn = new int[nnz];
		for (int ii=0; ii<nnz; ii++) {
			jcn[ii] = icol[ii] + 1;
		}
	}
	this->a = matrix->get_nonzeros();
	
	// mumps problem number
	this->my_mumps_number = constants::mumps_counter;
	constants::mumps_counter++;
	
	// ---------------------------------------------------
	// send sign to slaves to set up new mumps problem
	// ---------------------------------------------------
	logmsg->emit(LOG_INFO_L3,"Getting slaves ready for next problem (%d)", my_mumps_number);
	int create = 1;	// will call a MumpsSlave constructor in which prepare() is called 
	mpi->broadcast(create, root);
	mpi->broadcast(this->my_mumps_number, root);

	// prepare() was called for slaves in MumpsSlave() constructor
	// prepare() initializes the MUMPS instance and assigns the matrix in the host process
	this->prepare();
);}


LinearSolverMumps::~LinearSolverMumps() throw (Exception *)
{STACK_TRACE(
	// send sign to slaves that the problem of this LinearSolver can be deleted
	int command = 2;	// will call a MumpsSlave destructor
	mpi->broadcast(command, root);
	mpi->broadcast(this->my_mumps_number, root);
);}


// called by host only!
void LinearSolverMumps::solve() throw (Exception *)
{STACK_TRACE(
	// copy RHS to solution vector - will be overwritten with solution
	for (int ii=0; ii<n; ii++) {
		this->solution[ii] = this->rhs[ii];
	}
	this->mumps_rhs = this->solution;
	
	logmsg->emit(LOG_INFO_L3,"LinearSolverMumps::solve()");
	// tell the slaves that a problem is about to be solved
	int command = 3;	// will call slave_solve() in MumpsSlave of a given number
	mpi->broadcast(command, root);
	mpi->broadcast(this->my_mumps_number, root);
	
	this->slave_solve();
);}


void MumpsSlave::check_error(int err)
{STACK_TRACE(
	switch (err) {
	case   0: return;
	case  -1: NEGF_EXCEPTION("An error occurred on a slave processor.");
	case  -2: NEGF_EXCEPTION("nnz is out of range.");
	case  -3: NEGF_EXCEPTION("invalid value for id.job");
	case  -4: NEGF_EXCEPTION("Error in user-provided permutation array.");
	case  -5: NEGF_EXCEPTION("Problem in REAL workspace allocation.");
	case  -6: NEGF_EXCEPTION("Matrix is singular in structure.");
	case  -7: NEGF_EXCEPTION("Problem in INT workspace allocation.");
	case  -8: NEGF_EXCEPTION("Main internal integer workarray IS is too small for factorization.");
	case  -9: NEGF_EXCEPTION("Main internal complex workassay S too small.");
	case -10: NEGF_EXCEPTION("Numerically singular matrix.");
	case -11: NEGF_EXCEPTION("Internal real/complex workarray S too small for solution.");
	case -12: NEGF_EXCEPTION("Internal real/complex workarray S too small for iterative refinement.");
	case -13: NEGF_EXCEPTION("Error in a fortran allocate statement.");
	case -14: NEGF_EXCEPTION("Internal integer workarray IS too small for solution.");
	case -15: NEGF_EXCEPTION("Integer workarray IS too small for iterative refinement and/or error analysis.");
	case -16: NEGF_EXCEPTION("N is out of range.");
	case -17: NEGF_EXCEPTION("Internal send buffer allocated dynamically is too small.");
	case -20: NEGF_EXCEPTION("Internal reception buffer allocated dynamically is too small.");
	case -21: NEGF_EXCEPTION("PAR value not allowed because one one processor available.");
	case -22: NEGF_EXCEPTION("A pointer array provided by the user is not associated, has insufficient size or should not be associated.");
	case -23: NEGF_EXCEPTION("MPI was not initialized.");
	case   1: NEGF_EXCEPTION("IRN of JCN index out of range.");
	case   2: NEGF_EXCEPTION("During error analysis the max-norm of the computed solution was zero.");
	case   4: NEGF_EXCEPTION("JCN was modified internally.");
	case   8: NEGF_EXCEPTION("Warning from iterative refiniement routine.");
	default: NEGF_EXCEPTION("Some other error occurred.");
	}
);}


void MumpsSlave::test_operation()
{STACK_TRACE(
	// copied from the MUMPS manual! 
	DMUMPS_STRUC_C id;
	int      n   = 2;			// matrix size
	int     nz   = 2;			// num. nonzeros
	int    irn[] = {1,2};		// row indices, 1-based
	int    jcn[] = {1,2};		// column indices, 1-based
	double   a[] = {1.0, 2.0};	// nonzeros of matrix
	double rhs[] = {1.0, 4.0};	// RHS
	
	// Initialize a MUMPS instance. Use MPI_COMM_WORLD. 
	id.job = -1; 
	id.par =  1; 
	id.sym =  0;
	id.comm_fortran = MUMPS_USE_COMM_WORLD;
	dmumps_c(&id);
	// Define the problem on the host 
	if (mpi->get_rank() == 0) {
		id.n   =   n; 
		id.nz  =  nz; 
		id.irn = irn; 
		id.jcn = jcn;
		id.a   =   a; 
		id.rhs = rhs;
	}
	id.ICNTL(1) = -1; // no outputs
	id.ICNTL(2) = -1; 
	id.ICNTL(3) = -1; 
	id.ICNTL(4) =  0;
	// Call the MUMPS package. 
	id.job = 6;
	dmumps_c(&id);
	// Terminate instance
	id.job = -2; 
	dmumps_c(&id); 
	if (mpi->get_rank() == 0) {	printf("Solution is : (%8.2f %8.2f)\n", rhs[0],rhs[1]);	}
);}

#endif // NOMUMPS
