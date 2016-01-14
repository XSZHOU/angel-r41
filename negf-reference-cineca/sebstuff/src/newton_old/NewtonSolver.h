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
#ifndef NEWTONSOLVER_H_NEGF
#define NEWTONSOLVER_H_NEGF

#include "all.h"

#include "Equation.h"
#include "CSRMatrix.h"
#include "LinearSolver.h"
#include "LinearSolverUmfpack.h"
#ifndef NEWTON_WITHOUT_OUTPUTDATA
	#include "OutputData.h"
#endif

using namespace std;

namespace negf { class NewtonSolverTest; class Selfconsistator; }

namespace negf {

/** This class implements the Newton iteration scheme to solve systems of nonlinear equations
 *
 *  Remember: we want to solve \f$ \vec{F}(\vec{x}) = 0 \f$.
 *  To do this we iterate \f$ \vec{x}^{(n+1)} = \vec{x}^{(n)} 
          - \left( \frac{\partial F_i}{\partial x_j} \right)^{-1} \vec{F}(\vec{x}^{(n)}) \f$.
 *  \f$ \vec{F}(\vec{x}^{(n)}) \f$ is called right-hand side and \f$ \left( \frac{\partial F_i}{\partial x_j} \right) \f$ Jacobian.
 *  \f$ - \left( \frac{\partial F_i}{\partial x_j} \right)^{-1} \vec{F}(\vec{x}^{(n)}) \f$ is the (Newton) update.
*/
class NewtonSolver
{
	friend class negf::NewtonSolverTest;
	friend class negf::Selfconsistator;
	friend class OuterLoop;
	// syntax might differ vor compiler versions earlier than g++-4.0.x

	public:

		NewtonSolver(vector<Equation *> & eqnhandler,
					 unsigned int max_iterations_ 	= 30,		//!< suggestion: 30
					 double relative_err_crit_		= 1e-10,	//!< suggestion: 1e-10
					 double absolute_err_crit_ 		= 1e-10,	//!< suggestion: 1e-10
					 bool bank_rose_ 				= false,	//!< suggestion: false
					 const string & name_			= "noname",	//!< suggestion: "yourname"
					 bool verbose_ 					= false,	//!< suggestion: true
					 bool debug_ 					= false		//!< suggestion: false
					 ) throw(Exception *);
		~NewtonSolver();	//!< de-allocates the Jacobian matrix, linear solver, update and RHS arrays

		// --------------------------------
		// essential functions for solving
		// --------------------------------
		bool solve() throw(Exception *);	//!< solve nonlinear Equation system be Newton iteration until convergence is achieved or a max. number of iterations is reached.
		void reset() throw(Exception *);	//!< set iteration counter, norms to zero

		// -----------------------------
		// setup functions
		// ----------------------------

		void set_relative_err_crit(double eps) throw(Exception *);	//!< set relative error criterion (w/ default convergence crit. this is a mixture rel/abs)
		void set_absolute_err_crit(double eps) throw(Exception *);  //!< set absolute error criterion (w/ default convergence crit. it doesn't make a difference)
		void set_verbose(bool yesorno) { this->verbose = yesorno; }	//!< set whether the NewtonSolver produces output
		void set_max_iterations(uint max_it) throw(Exception *);	//!< set the maximum number of iterations performed
		void set_min_iterations(uint min_it) throw(Exception *);	//!< set the minimum number of iterations performed, regardless of whether convergence is achieved
		void set_name(string name_) { NEGF_ASSERT(name_.length()>0, "invalid name"); this->name = name_; }	//!< set the name of the NewtonSolver
		void set_debug(bool flag = true) { debug = flag; }			//!< set a flag if debugging is done
		void set_bank_rose(bool flag = true) { bank_rose = flag; }	//!< set a flag whether the Bank-Rose update scheme is used (scales the update vector)
	
#ifndef NEWTON_WITHOUT_OUTPUTDATA	
		// for debugging
		void add_outputdata(OutputData * outputdata) throw(Exception *); 	//!< assign the OutputData object for writing files after each iteration
		void clear_outputdata() 			{ this->debug_data.clear(); }	//!< remove all OutputData objects
#endif
		
		// use these with caution!
		void set_timestamp(uint new_timestamp) { this->timestamp = new_timestamp; }	//!< set the NewtonSolver timestamp
		void set_last_converged_timestamp(uint tstamp_) { this->last_converged_timestamp = tstamp_; }	//!< set the timestamp where convergence was last reached

		// -----------------------------
		// access functions
		// -----------------------------

		uint 		get_max_iterstions() 	const { return this->max_iterations; }	//!< get the maximum number of iterations performed before aborting
		string 		get_name() 				const { return name; }					//!< get the name of the NewtonSolver
		bool 		get_debug() 			const { return debug; }					//!< check whether debugging is done
		bool 		get_bank_rose() 		const { return bank_rose; }				//!< check whether the Bank-Rose algorithm is used (scales update vectors)
		double		get_relative_err_crit() const { return relative_err_crit; }		//!< get the relative error criterion
		double		get_absolute_err_crit() const { return absolute_err_crit; }		//!< get the absolute error criterion
		uint 		get_timestamp() 		const { return this->timestamp;}		//!< get the current timestamp of the NewtonSolver
		
		const vector<Equation *> & get_handler() 	const { return this->handler; }	//!< get the vector storing all involved Equation objects
		
		double 		get_jac_assembly_time() const { return this->jac_assembly_time; }	//!< get the aggregate time spent in assembling the Jacobian
		double 		get_rhs_assembly_time() const { return this->rhs_assembly_time; }	//!< get the aggregate time spent in assembling the RHS
		double 		get_linear_solver_time() const { return this->linear_solver_time; }	//!< get the aggregate time spent in solving the linear system
		
		void 		display_dependency_info() const throw(Exception *);	//!< display some debug information about the mutual dependencies of the Equation objects
		
		
	protected:

		// ------------------------
		// methods
		// ------------------------
		//public:
		void step() throw(Exception *);								//!< perform a Newton step
		//protected:
		void update_dependencies() 				throw(Exception *);	//!< update the values of all needed Equation objects which are not Newton solution variables
		void calculate_jacobian() 				throw(Exception *);	//!< assemble the Jacobian
		void calculate_rhs() 					throw(Exception *);	//!< assemble the RHS
		void calculate_update() 				throw(Exception *);	//!< solve the linear system
		void update_newton_vars(uint tstamp) 	throw(Exception *);	//!< called when Bank-Rose is turned off.
		void bank_rose_update()								throw(Exception *);	//!< perform Bank-Rose update algorithm
		void update_newton_vars(double alpha, uint tstamp)  throw(Exception *);	//!< update using scaled update vector. used for Bank-Rose update

		// calculation of norms and residuals
		void calculate_rhs_norm()				throw(Exception *);	//!< calculate 2-norm of RHS (total and for each individual Equation)
		void calculate_update_norm()			throw(Exception *);	//!< calculate 2-norm of update (total and for each individual Equation)
		void calculate_new_variable_norm()		throw(Exception *);	//!< calculate 2-norm of Equation variables
		void calculate_residuum()				throw(Exception *);	//!< calculate residuum |Ax-b| to check correctness of linear solver solution

		bool converged() 				const throw(Exception *);	//!< test convergence
		bool global_convergence() 		const throw(Exception *);	//!< helper. must only be called from converged()!
		bool eqn_convergence() 			const throw(Exception *);	//!< helper. must only be called from converged()!
		bool mixed_convergence() 		const throw(Exception *);	//!< helper. must only be called from converged()!
		bool individual_convergence() 	const throw(Exception *);	//!< helper. must only be called from converged()!


		// -----------------------
		// class variables
		// -----------------------
		vector<Equation *> handler;			//!< list of all Equation objects (implicit and explicit)
		vector<Equation *> newton_list; 	//!< list of Equations which are Newton solution variables (typically all implicit equations)

		CSRMatrix<double> * jacobian;		//!< the Jacobian
		double * 			rhs;			//!< the right-hand side (Newton function)
		double *			update;			//!< the update

		LinearSolver * linear_solver;		//!< solves Ax=b

		double old_variable_norm;			//!< 2-norm of total variable vector
		double new_variable_norm;			//!< 2-norm of total variable vector
		double rhs_norm;					//!< 2-norm of total RHS vector
		double update_norm;					//!< 2-norm of total update vector
		vector<double> old_variable_norm_eqn; //!< 2-norm of the variables for each Equation
		vector<double> new_variable_norm_eqn; //!< 2-norm of the variables for each Equation
		vector<double> rhs_norm_eqn;		//!< 2-norm of the RHS, split into each Equation
		vector<double> update_norm_eqn;		//!< 2-norm of the update, split into each Equation

		double relative_err_crit;			//!< relative error criterion to check convergence (default: 1e-10)
		double absolute_err_crit;			//!< absolute error criterion to check convergence (default: 1e-10)

		uint timestamp;						//!< timestamp the current values correspond to
		uint last_converged_timestamp;		//!< timestamp when convergence was last reached
		uint nvars; 						//!< total number of variables in the Newton solution vector
		
		uint max_iterations;				//!< maximum number of iterations performed before aborting
		uint iteration;						//!< counter for the iteration number
		uint min_iterations;				//!< the least number of iterations before converged() can be true. default: 2 (set in constructor)
		
		string name;						//!< the name of the NewtonSolver instance
		bool verbose;						//!< controls amount of screen output
		double jac_assembly_time;			//!< time spent assembling the Jacobian
		double rhs_assembly_time;			//!< time spent assembling the RHS
		double linear_solver_time;			//!< time spent in the solution of the linear system
		
		bool debug;							//!< flag that determines whether debugging is done
#ifndef NEWTON_WITHOUT_OUTPUTDATA
		vector<OutputData *> debug_data;	//!< contains OutputData objects to write intermediate results to file (after each iteration)
#endif
		
		bool bank_rose;						//!< flag if the Bank-Rose iteration scheme is used or not
		double alpha;						//!< damping factor used by Bank-Rose
		
		bool share_nonzeros_equally;		//!< switch between line-based distribution and nonzeros-based distribution (OpenMP)

};


} // end of namespace negf

#endif /*NEWTONSOLVER_H_NEGF*/

