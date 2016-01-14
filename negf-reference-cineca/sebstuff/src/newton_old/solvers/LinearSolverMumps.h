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
#ifndef LINEARSOLVERMUMPS_H_NEGF
#define LINEARSOLVERMUMPS_H_NEGF

#ifndef NOMUMPS

#include "all.h"

#include "CSRMatrix.h"
#include "LinearSolver.h"

extern "C" {
#include "dmumps_c.h" // double precision mumps
}

using namespace std;

namespace negf {

	/** A single MumpsSlave instance will be set up for all MPI processes started with mpiexec */
	class MumpsSlave {
	public:
		MumpsSlave() throw (Exception *);	//!< calls prepare() if it is not root
		virtual ~MumpsSlave();
		
		void slave_solve() throw (Exception *); //!< performs a MUMPS solution
	protected:

		void prepare();				//!< set up matrix arrays (for root), configuration parameters
		void terminate();			//!< tell MUMPS to terminate
		void check_error(int err);	//!< check MUMPS error code
		void test_operation();		//!< calculates an example from the MUMPS manual

		int  root;					//!< fixed to 0 (for MPI)
		
		DMUMPS_STRUC_C id;

		// things only used by master process (root, host)
		int        n; 		//!< matrix size
		int 	 nnz; 		//!< number of nonzeros
		int    * irn;		//!< 1-based prow
		int    * jcn;		//!< 1-based icol
		double *   a;		//!< nonzero entries
		double * mumps_rhs;	//!< null pointer for slaves, same address as solution pointer for host
	};
	
	/** a slave process constructs this class
	 * it then permanently waits for commands from the master thread by mpi->broadcast()
	 * -1 --> exit
	 *  1 --> create new slave with a certain number
	 *  2 --> delete a slave with a certain number
	 *  3 --> solve a problem in slave with a certain number
	 *  for choices 1,2,3 the number must be mpi->broadcasted as well.
	 * 
	 * The master process does NOT construct this class!
	 */
	class MumpsSlaveCollector {
	public:
		MumpsSlaveCollector() throw (Exception *);
		virtual ~MumpsSlaveCollector();			//!< deletes all remaining slaves
		
	private:
		void 			create_slave(uint num); //!< create a MumpsSlave object and insert it into the slaves vector at a specific position
		void 			delete_slave(uint num);	//!< delete a MumpsSlave with a specific number
		MumpsSlave * 	get_slave(uint num);	//!< get a slave with a specific number
		void 			solve(uint num);		//!< called by root only! broadcasts a signal s.th. all slaves perform slave_solve()
		
		uint 					mumps_counter;	//!< keeps track of the number of slaves (MUMPS problems) created
		const uint 				max_num_slaves; //!< maximum number of simultaneous problems that can be handled
		vector<MumpsSlave *> 	slaves;			//!< the slaves
	};
	
	
	/** Interface to the direct sparse solver MUMPS */
	class LinearSolverMumps: public MumpsSlave, public LinearSolver {

	public:
		
		LinearSolverMumps(CSRMatrix<double>* matrix_, double * rhs_, double * solution_) throw (Exception *);
		~LinearSolverMumps() throw (Exception *);

		void solve() throw (Exception *); //!< solve Ax=b
	
	protected:
		
		int my_mumps_number; //!< which problem is it - there might be several LinearSolverMumps instances present at the time!
	};

} // end namespace

#endif // NOMUMPS
#endif /*LINEARSOLVERMUMPS_H_NEGF*/
