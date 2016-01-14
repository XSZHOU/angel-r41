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
#include "NewtonSolver.h"
#ifdef _OPENMP
	#include <omp.h>
#endif
using namespace negf;
typedef vector<Equation *>::const_iterator eqncit;


/** Constructor sets up the network between all the given Equation's and determines the Jacobian structure.
	It also creates the LinearSolver and allocates space for the RHS and update vectors. <BR>
	1. Let all Equation objects in the handler "handshake", i.e. determine (recursively) all the dependencies <BR>
	2. Determine which Equation objects are solved in the Newton an let the Equation's know that. <BR>
	3. Determine the nonzero entries in the Jacobian for each Equation from the sparsity structure.
       One of the goals is that the Jacobian assembly is distributed equally on all  OpenMP threads. When quantum regions are
       included in the simulation, the interpolation stuff has the effect that some lines in the Jacobian are almost full
       whereas others are practically diagonal. So we cannot just distribute the lines equally. <BR>
    4. Create the Jacobian matrix an set its sparsity structure by recursively looking at the Equation sparsities. <BR>
    5. Create RHS and update vectors, LinearSolver <BR>
    6. (experimental) do a reordering to improve matrix condition. */
NewtonSolver::NewtonSolver(vector<Equation *> & handler_,
					 uint max_iterations_,
					 double relative_err_crit_,
					 double absolute_err_crit_,
					 bool bank_rose_,
					 const string & name_,
					 bool verbose_,
					 bool debug_) throw(Exception *):
	handler(handler_),
	old_variable_norm(0),
	new_variable_norm(0),
	update_norm(0),
	relative_err_crit(relative_err_crit_),
	absolute_err_crit(absolute_err_crit_),
	timestamp(0),
	last_converged_timestamp(0),
	max_iterations(max_iterations_),
	iteration(0),
	min_iterations(constants::min_newton_iterations),
	name(name_),
	verbose(verbose_),
	jac_assembly_time(0.0),
	rhs_assembly_time(0.0),
	linear_solver_time(0.0),
	debug(debug_),
	bank_rose(bank_rose_),
	alpha(1.0),
	share_nonzeros_equally(true)
{STACK_TRACE(	
	NEGF_ASSERT(max_iterations<1000, "too many iterations");
	NEGF_ASSERT(relative_err_crit>1e-16 && absolute_err_crit > 0, "error criteria too small");
	NEGF_ASSERT(relative_err_crit<0.1, "relative error criterum too big");
	
/*  - handshake equations
	- set up list of all newton eqns
	- set offset for every newton equation
	- check that all dependencies are included in the handler list
	- set up newton_var_list of every equation
	- determine nvars, the total amount of variables in the Newton solver 
*/
	this->nvars = 0;
	ostringstream msg;
	msg << "Equations {";
	for (eqncit eq = handler.begin(); eq != handler.end(); eq++)
	{
		(*eq)->handshake();
		if ((*eq)->is_solved_in_newton()) {
			logmsg->emit(LOG_INFO_L1,"Added equation \"%s\" (%d vars, reference size %6.2e, offset %d) to newton list.",
							(*eq)->get_name().c_str(),(*eq)->get_num_variables(), (*eq)->get_reference_size(),nvars);
			newton_list.push_back(*eq);
			(*eq)->set_my_offset(nvars);
			//cout << "Equation " << (*eq)->get_name() << " (" << (*eq)->get_num_variables() << " variables) gets offset " << nvars << endl;
			NEGF_FASSERT((*eq)->get_num_variables() > 0, "Equation \"%s\" has no variables!", (*eq)->get_name().c_str());
			this->nvars += (*eq)->get_num_variables();
		} else {
			msg << "\"" << (*eq)->get_name() << "\" (" << (*eq)->get_num_variables() << " vars), ";
		}
	}
	msg << "} were not added.";
	logmsg->emit(LOG_INFO_L2, msg.str().c_str());
	this->old_variable_norm_eqn.resize(newton_list.size(), 0.0);
	this->new_variable_norm_eqn.resize(newton_list.size(), 0.0);
	this->rhs_norm_eqn.resize(newton_list.size(), 0.0);
	this->update_norm_eqn.resize(newton_list.size(), 0.0);
	if (verbose)
		logmsg->emit(LOG_INFO_L1,"Newton solver has %d variables.",nvars);
	
	// consistency checks
	logmsg->emit(LOG_INFO_L1, "Performing consistency checks...");
	for (uint ii=0; ii<handler.size(); ii++)
	{
		msg.str("");
		msg << "all_eqns list of an equation \"" << handler[ii]->get_name() << "\" was not yet set up.";
		NEGF_ASSERT( handler[ii]->is_all_eqns_ready(), msg.str().c_str() );
		vector<Equation *>	all_eqns = handler[ii]->get_all_eqns();
		for (uint jj = 0; jj < all_eqns.size(); jj++) {
			msg.str("");
			msg << "null pointer found in all_eqns of equation \"" << handler[ii]->get_name() << "\".";
			NEGF_ASSERT(all_eqns[jj]!=NULL, msg.str().c_str() );
			msg.str("");
			msg << "Equation \"" << all_eqns[jj]->get_name() 
				<< "\" was found as dependency of \"" << handler[ii]->get_name()
				<< "\" but is not in the equation list of the Newton solver.";
			NEGF_ASSERT( find(handler.begin(), handler.end(), all_eqns[jj]) != handler.end(), msg.str().c_str());
		}

		// if its an explicit equation, check if theres any newton equation depending on it
		if (!handler[ii]->is_solved_in_newton()) {
			bool useful = false;
			for (uint jj = 0; jj < handler.size(); jj++) {
				if (jj!=ii && handler[jj]->is_solved_in_newton()) {
					vector<Equation *>	jj_eqns = handler[jj]->get_all_eqns();
					if ( find(jj_eqns.begin(), jj_eqns.end(), handler[ii]) != jj_eqns.end() ) {
						useful = true;
						break;
					}
				}
			}
			msg.str("");
			msg << "   equation \"" << handler[ii]->get_name() << "\" is not needed in the Newton iteration.";
			if (!useful)
				logmsg->emit(LOG_INFO_L1, msg.str().c_str());
		}
	}
	
	/*  determine for each equation which dependencies are solved in the newton
	    note: set_newton_var_list needs to be called AFTER the offsets of 
	          all equations have been set, since the list is ordered w.r.t. the offset !*/
	logmsg->emit(LOG_INFO_L1, "Setting the list of newton vars each eqn depends on (also indirectly)...");
	for (eqncit eq = handler.begin(); eq != handler.end(); eq++) {
		(*eq)->set_newton_var_list();
	}
	
	/*  store the number of nonzeros per line in the equation object
        this can be used s.th. less memory will need to be allocated each time the jacobian is calculated
		this yields a big performance improvement */
	logmsg->emit(LOG_INFO_L1, "Determining for each eqn the number of nonzeros for each line...");
	for (eqncit eq = handler.begin(); eq != handler.end(); eq++) {
		(*eq)->set_nonzeros();
	}

	// allocate memory for rhs vector, update vector and Jacobian
	if (verbose) logmsg->emit(LOG_INFO_L1,"Allocating memory...");
	this->rhs    = new double[nvars];
	this->update = new double[nvars];
	for (uint idx = 0; idx < nvars; ++idx) {
		rhs[idx] = 0.0;
		update[idx] = 0.0;
	}
	this->jacobian = new CSRMatrix<double>(nvars, nonsymmetric_matrix);

	// set matrix structure of Jacobian
	for (eqncit eq = newton_list.begin(); eq != newton_list.end(); eq++)
	{
		if (verbose) logmsg->emit(LOG_INFO_L1,"   equation \"%s\"...",(*eq)->get_name().c_str());
		uint offset = (*eq)->get_offset();
		uint nonzeros = 0;
		uint idcs[constants::eqn_array_size];	// hard coded max. number of nonzeros per line, defined in Equation.cpp
		for (uint ii = 0; ii < (*eq)->get_num_variables(); ++ii)
		{
			(*eq)->get_sparsity(ii, nonzeros, idcs);
			NEGF_FASSERT(nonzeros > 0, "A line in the jacobian (var. %d of eqn %s) would be completely empty from its sparsity pattern.",
				ii, (*eq)->get_name().c_str());
			for (uint jj = 0; jj < nonzeros; ++jj) {
				NEGF_FASSERT(idcs[jj]<nvars, "an entry in the sparsity pattern of eqn %s, line %s was out of range (%d>%d)",
						(*eq)->get_name().c_str(), jj, idcs[jj],nvars);
				jacobian->announce(ii+offset,idcs[jj]);	// note that announce(ii,jj) expects 0-based indices even though prow,icol will be (by default) 1-based!
			}
		}
	}
	if (verbose) logmsg->emit(LOG_INFO_L2,"Setting matrix structure...");
	jacobian->set_structure(); // sets arrays prow, icol, nonzeros
	if (verbose) {
		logmsg->emit(LOG_INFO_L1,"Structure of Jacobian has been set and memory was allocated (%d nonzeros).",
			jacobian->get_num_nonzeros());
	}
	
	// set up linear solver
	// just needs the address of the matrix, rhs and update, not the values!
	if (verbose) logmsg->emit(LOG_INFO_L3, "creating UMFPACK linear solver");
	this->linear_solver = new LinearSolverUmfpack(jacobian, rhs, update);
);}


NewtonSolver::~NewtonSolver()
{STACK_TRACE(
	delete linear_solver;
	delete jacobian;
	delete [] update;
	delete [] rhs;
);}


void NewtonSolver::reset() throw(Exception *)
{STACK_TRACE(
	this->iteration = 0;
	this->new_variable_norm = 0.0; 
	this->old_variable_norm = 0.0; 
	this->update_norm = 0.0; 
	this->old_variable_norm_eqn.assign( old_variable_norm_eqn.size(), 0.0);
	this->new_variable_norm_eqn.assign( new_variable_norm_eqn.size(), 0.0);
	this->update_norm_eqn.assign(update_norm_eqn.size(), 0.0);
	this->rhs_norm_eqn.assign( rhs_norm_eqn.size(), 0.0);
	/*this->linear_solver_time = 0.0;*/
	// timestamp is NOT reset!!!
	// this->timestamp++;
);}


void NewtonSolver::set_max_iterations(uint max_it) throw(Exception *)
{STACK_TRACE(
	NEGF_ASSERT(max_it>1, "max_iterations must be >=2.");
	this->max_iterations = max_it;
);}

void NewtonSolver::set_min_iterations(uint min_it) throw(Exception *)
{STACK_TRACE(
	NEGF_ASSERT(min_it < this->max_iterations, "min_iterations must be < max_iterations");
	this->min_iterations = min_it;
);}

/** Call this method to solve the nonlinear Equation system! <BR>
    It does the following steps in order:
    1. Initially increment the timestamp by 1 and update all non-Newton variables to make sure they were computed at least once.
    2. Iterate the system using step() and perform convergence checks.
    @return true if converged, false if not converged. */
bool NewtonSolver::solve() throw(Exception *)
{STACK_TRACE(
	if (verbose) {
		logmsg->emit(LOG_INFO_L1,"|===================================================================================|");
		logmsg->emit(LOG_INFO_L1,"|                              starting Newton iteration                            |");
		logmsg->emit(LOG_INFO_L1,"|===================================================================================|");
	}
	
	NEGF_ASSERT(this->iteration < this->max_iterations, "Reset newton solver first.");
		
	// advance timestamp of Newton eqns and Newton solver by 1
	// like this one makes sure that the explicit eqns are freshly computed
	for (uint eq = 0; eq < this->newton_list.size(); eq++)
	{
		NEGF_ASSERT(newton_list[eq]!=NULL, "null pointer encountered.");
	//	NEGF_FASSERT(newton_list[eq]->get_timestamp()==this->timestamp, 
	//		"Timestamp of implicit eqn %s (%d) did not equal newton timestamp (%d)",
	//		newton_list[eq]->get_name().c_str(), newton_list[eq]->get_timestamp(), this->timestamp);
		newton_list[eq]->set_timestamp(this->timestamp + 1);
	}
	this->timestamp++;
		
	// -----------------------------------------------
	// initial norms and RHS
	// -----------------------------------------------
	// calculate RHS
	this->update_dependencies();
	this->calculate_rhs();
	
	// Calculate norm of variables and norm of RHS
	this->calculate_new_variable_norm();
	this->calculate_rhs_norm();

	// screen output
	if (verbose) {
		logmsg->emit(LOG_INFO_L1,"Initial situation (timestamp %d):",this->timestamp);
		logmsg->emit(LOG_INFO_L1,"            equation     eqn_norm     rhs_norm");
		logmsg->emit(LOG_INFO_L1,"-----------------------------------------------------");
		for (uint ee = 0; ee < newton_list.size(); ee++) 
		{
			logmsg->emit(LOG_INFO_L1,"%20s     %7.3e    %7.3e",
						newton_list[ee]->get_name().c_str(), new_variable_norm_eqn[ee], rhs_norm_eqn[ee]);
		}
		logmsg->emit(LOG_INFO_L1,"-----------------------------------------------------");	
		logmsg->emit(LOG_INFO_L1,"%20s     %7.3e    %7.3e", "total", this->new_variable_norm, this->rhs_norm);
		logmsg->emit(LOG_INFO_L2,"convergence criteria: absolute %6.2e OR relative %6.2e", this->absolute_err_crit, this->relative_err_crit);
	}
	
	
	// -----------------------------------------------
	// perform Newton iteration
	// -----------------------------------------------
	
	do {
		//if (iteration > 3) this->set_linear_solver_type(linear_solver_type::pardiso);
		this->step();
		this->iteration++;
	} while (!this->converged()  &&  this->iteration < this->max_iterations);
	
	if (this->converged()) {
		logmsg->emit(LOG_INFO,"Newton \"%s\" converged in %d iterations!",	name.c_str(),iteration);
		this->last_converged_timestamp = this->timestamp;
		return true;
	} else {
		logmsg->emit(LOG_INFO,"Newton \"%s\" did not converge in %d iterations!", name.c_str(), max_iterations);
		return false;
	}
);}
	
/** Perform a single Newton iteration step: \f$ \vec{x}^{(n+1)} = \vec{x}^{(n)} - \left( \frac{\partial F_i}{\partial x_j} \right)^{-1} \vec{F}(\vec{x}^{(n)}) \f$ */
void NewtonSolver::step() throw(Exception *)
{STACK_TRACE(
	if (verbose) {
		logmsg->emit(LOG_INFO_L1,"========================================================================================");
		logmsg->emit(LOG_INFO_L1,"Performing Newton step %d...",this->timestamp);
	}

	// DO IT!
	//update_dependencies(); // already performed in calculate_jacobian, calculate_rhs
	this->calculate_jacobian();
	//if (this->timestamp==7 && this->name=="device") {
	//	cout << "SAVING..." << endl;
	//	this->jacobian->save_to_file("superlu_test.mat");
	//}
	this->calculate_rhs();
	this->calculate_update();
	
	this->old_variable_norm     = this->new_variable_norm;
	this->old_variable_norm_eqn = this->new_variable_norm_eqn;
	logmsg->emit(LOG_INFO_L2,"   performing update...");
	// the first update after the voltage has been changed will usually involve a rise in the RHS
	// because initially only the vertices at the contact are wrong, which does not result in a big RHS
	if (verbose && bank_rose && this->iteration==0) 
		logmsg->emit(LOG_INFO, "First iteration. No Bank-Rose is performed.");
	if (bank_rose && this->iteration!=0) {
		bank_rose_update();	// also calculates new RHS, RHS norm and increases timestamp
	} else {
		this->alpha = 1.0;
		this->calculate_rhs_norm();
	
		// security check of timestamps
		for (uint ee = 0; ee < newton_list.size(); ee++) {
			vector<Equation *> alleqns = newton_list[ee]->get_all_eqns();
			for (uint ii = 0; ii < alleqns.size(); ii++) {
				NEGF_FASSERT(alleqns[ii]->get_timestamp() >= this->timestamp,
						"In equation \"%s\" the timestamp was not up to date (was %d, should have been at least %d).",
						alleqns[ii]->get_name().c_str(), alleqns[ii]->get_timestamp(), this->timestamp);
			}
		}
		
		// update stored variables
		this->update_newton_vars(this->timestamp + 1); // also increases timestamp within equations
		this->timestamp++;
		
		// calc RHS mainly for screen output - doesn't really affect simulation time because we'd have to do it next iteration anyway
		this->calculate_rhs();
		this->calculate_rhs_norm();
	}
	
	this->calculate_update_norm();
	this->calculate_new_variable_norm();

	// screen output
	if (verbose) {
		logmsg->emit(LOG_INFO_L1,"finished! New timestamp of equations and solver: %d.",this->timestamp);
		logmsg->emit(LOG_INFO_L1,"            equation     old_norm     new_norm     new_rhs      update_norm  rel_update");
		logmsg->emit(LOG_INFO_L1,"---------------------------------------------------------------------------------------");
		for (uint ee = 0; ee < newton_list.size(); ee++) 
		{
			Equation * eqn = newton_list[ee];
			double relative_update = update_norm_eqn[ee] / new_variable_norm_eqn[ee];
			logmsg->emit(LOG_INFO_L1,"%20s     %7.3e    %7.3e    %7.3e    %7.3e    %7.3e",
					eqn->get_name().c_str(), old_variable_norm_eqn[ee], new_variable_norm_eqn[ee],
					rhs_norm_eqn[ee],        update_norm_eqn[ee]/this->alpha,       relative_update/this->alpha);
		}
		logmsg->emit(LOG_INFO_L1,"---------------------------------------------------------------------------------------");	
		logmsg->emit(LOG_INFO_L1,"%20s     %7.3e    %7.3e    %7.3e    %7.3e",
					"total", this->old_variable_norm, this->new_variable_norm, this->rhs_norm, this->update_norm/this->alpha);
		logmsg->emit(LOG_INFO_L2,"convergence criteria: absolute %6.2e OR relative %6.2e", this->absolute_err_crit, this->relative_err_crit);
	}
	
#ifndef NEWTON_WITHOUT_OUTPUTDATA
	// write data to file for debug
	if (debug) 
	{
		this->update_dependencies();
		
		for (uint ii = 0; ii < this->debug_data.size(); ii++) {
			NEGF_ASSERT(debug_data[ii]!=NULL, "encountered null pointer.");
			
			// make sure all the debug data is up to date
			for (uint jj = 0; jj < debug_data[ii]->get_num_dat_equations(); jj++) 
			{
				if (!debug_data[ii]->get_dat_equation(jj)->is_solved_in_newton()
					&& debug_data[ii]->get_dat_equation(jj)->get_timestamp()>=this->timestamp) 
				{
					logmsg->emit(LOG_INFO_L2,  "Updating explicit eqn %s to timestamp %d for saving debug file.",
									debug_data[ii]->get_dat_equation(jj)->get_name().c_str(), this->timestamp);
					debug_data[ii]->get_dat_equation(jj)->compute_values(this->timestamp);
				} else {
					NEGF_ASSERT(debug_data[ii]->get_dat_equation(jj)->get_timestamp() >= this->timestamp, "weird timestamp.");
				}
			}
			
			// write to file
			string filename = debug_data[ii]->get_filename();
			ostringstream out_file;
			out_file << filename << ".step" << timestamp;
			debug_data[ii]->set_filename(out_file.str());
#ifndef NODFISE
			debug_data[ii]->write_dat();
#endif
			debug_data[ii]->set_filename(filename);
		}
	}
#endif
);}


void NewtonSolver::update_dependencies() throw(Exception *)
{STACK_TRACE(
	timeval tp1;
	gettimeofday(&tp1, NULL);
	
	ostringstream depmsg;
	depmsg << "   updating dependency values";
	#ifdef _OPENMP
		depmsg << " using " << omp_get_max_threads() << " threads";
	#endif
	depmsg << "...";
	if (verbose) {
		if (logmsg->get_level()==LOG_INFO_L2) 
			logmsg->emit_noendl(LOG_INFO_L2, depmsg.str().c_str());
		else 
			logmsg->emit(LOG_INFO_L2, depmsg.str().c_str());
	}
	for (uint eq = 0; eq < newton_list.size(); eq++)
	{	
		newton_list[eq]->check_dependency_timestamps(this->timestamp);
	}
	
	timeval tp2;
	gettimeofday(&tp2, NULL);
	if (verbose) logmsg->emit(LOG_INFO_L2,  "  finished (%.3g[s])",tp2.tv_sec - tp1.tv_sec + 1.0e-6 * (tp2.tv_usec - tp1.tv_usec));

);}


void NewtonSolver::calculate_jacobian() throw(Exception *)
{//STACK_TRACE(
 // OpenMP does not go together with our STACK_TRACE macro; exceptions must be handled within same thread
	
	timeval tp1;
	gettimeofday(&tp1, NULL);
	
	ostringstream jacmsg;
	jacmsg << "   calculating Jacobian";
	#ifdef _OPENMP
		jacmsg << " using " << omp_get_max_threads() << " threads";
	#endif
	jacmsg << "...";
	if (verbose) {
		if (logmsg->get_level()==LOG_INFO_L1) 
			logmsg->emit_noendl(LOG_INFO_L1, jacmsg.str().c_str());
		else 
			logmsg->emit(LOG_INFO_L1, jacmsg.str().c_str());
	}
	for (uint eq = 0; eq < newton_list.size(); eq++)
	{	
		logmsg->emit(LOG_INFO_L2,"      equation \"%s\"...",newton_list[eq]->get_name().c_str());
		uint offset = newton_list[eq]->get_offset();
		
		// make sure the dependencies are up to date
		// is now done in a separate method	
		logmsg->emit(LOG_INFO_L3,"      eqn->check_dependency_timestamps()...");
		newton_list[eq]->check_dependency_timestamps(this->timestamp);

		logmsg->emit(LOG_INFO_L3,"      get_newton_derivatives()...");
		#pragma omp parallel shared(eq)
		{
		try {
		
		uint     nonzeros = 0;
		//uint   * indices  = new uint  [constants::eqn_array_size];	// hard coded max. number of nonzeros per line, defined in Equation.cpp
		//double * coeffs   = new double[constants::eqn_array_size];
		uint   indices[constants::eqn_array_size];
		double coeffs[constants::eqn_array_size];
		
		if (this->share_nonzeros_equally) 	// share nonzeros equally
		{
			uint my_first_line = 0;
			uint my_last_line = newton_list[eq]->get_num_variables()-1;
			#ifdef _OPENMP
				// determine number of nonzeros that each CPU has to calculate when distributed equally
				uint total_eqn_nonzeros = 0;
				for (uint line=0; line < newton_list[eq]->get_num_variables(); line++) {
					total_eqn_nonzeros += newton_list[eq]->get_nonzeros(line);
				}
				uint nz_per_CPU = max(int(ceil(double(total_eqn_nonzeros) / double(omp_get_num_threads()))), 
										3*constants::newton_calc_jac_chunk);
				#pragma omp master
				logmsg->emit(LOG_INFO_L3,"   total_eqn_nonzeros=%d, nz_per_CPU=%d, total lines: 0...%d",
									total_eqn_nonzeros,nz_per_CPU, newton_list[eq]->get_num_variables()-1);
				
				// determine the lines the current CPU has to calculate
				my_first_line = 100000000;
				my_last_line = 0;
				uint nz = 0;
				for (uint line=0; line < newton_list[eq]->get_num_variables(); line++) {
					if (omp_get_thread_num() * nz_per_CPU <= nz && my_first_line==100000000) {
						my_first_line = line;
					}
					if (line==newton_list[eq]->get_num_variables()-1) {
						my_last_line = line;
						break;
					}
					if ((omp_get_thread_num()+1) * nz_per_CPU <= nz) {
						my_last_line = line - 1;
						break;
					}
					nz += newton_list[eq]->get_nonzeros(line);
				}
				if (verbose) {
					uint tot_nz = 0;
					for (uint line=my_first_line; line<=my_last_line; line++) {
						tot_nz += newton_list[eq]->get_nonzeros(line);
					}
					#pragma omp critical
					logmsg->emit(LOG_INFO_L3,"   thread %d calculates lines %d to %d with %d nonzeros", 
						omp_get_thread_num(), my_first_line, my_last_line, tot_nz);
				}
			#endif
			for (uint ii=my_first_line; ii<=my_last_line; ii++)
			{
				nonzeros = 0;
				
				newton_list[eq]->get_newton_derivatives(ii, nonzeros, indices, coeffs);
	
				for (uint idx = 0; idx < nonzeros; idx++) {
					jacobian->set(offset+ii, indices[idx], coeffs[idx]);	// note that set(ii,jj) expects 0-based indices even though prow,icol will be (by default) 1-based!
				}
			}
		} else 	// share Jacobian lines equally (roughly; will be slow when there are a few lines having much more nonzeros than the rest)
		{
			#pragma omp for //schedule(dynamic, constants::newton_calc_jac_chunk)
			for (int ii=0; ii < (int) newton_list[eq]->get_num_variables(); ii++)
			{
				nonzeros = 0;
				try {
				newton_list[eq]->get_newton_derivatives(ii, nonzeros, indices, coeffs);
				} catch (Exception * e) { 			
					e->append_info(__LINE__,__FILE__,__DATE__,__TIME__, __func__); 
					#ifdef _OPENMP 
						cout << e->get_reason() << "\n   equation:" <<this->get_name()  << "\nterminating from thread " << omp_get_thread_num() <<"\n";
						exit(1);
					#else
						NEGF_EXCEPTION(e->get_reason().c_str());
					#endif
				}
	
				for (uint idx = 0; idx < nonzeros; idx++) {
					jacobian->set(offset+ii, indices[idx], coeffs[idx]);
				}
			}
		}	// else
		
		} catch (Exception * e) { 
			e->append_info(__LINE__,__FILE__,__DATE__,__TIME__, __func__); 
			#ifdef _OPENMP 
				cout << e->get_reason() << "\n   equation:" <<this->get_name()  << "\nterminating from thread " << omp_get_thread_num() <<"\n";
				exit(1);
			#else
				NEGF_EXCEPTION(e->get_reason().c_str());
			#endif
		}
		} // pragma omp parallel
	}

	timeval tp2;
	gettimeofday(&tp2, NULL);
	this->jac_assembly_time += tp2.tv_sec - tp1.tv_sec + 1.0e-6 * (tp2.tv_usec - tp1.tv_usec);
	if (verbose) logmsg->emit(LOG_INFO_L1,  "  finished (%.3g[s])",tp2.tv_sec - tp1.tv_sec + 1.0e-6 * (tp2.tv_usec - tp1.tv_usec));

}
//);}


void NewtonSolver::calculate_rhs() throw(Exception *)
{//STACK_TRACE(
 // OpenMP does not go together with our STACK_TRACE macro; exceptions must be handled within same thread
	
//	try {
	timeval tp1;
	gettimeofday(&tp1, NULL);

	ostringstream rhsmsg;
	rhsmsg <<"   calculating RHS";
	#ifdef _OPENMP
		rhsmsg << " using " << omp_get_max_threads() << " threads";
	#endif
	rhsmsg << "...";
	if (verbose) {
		if (logmsg->get_level()==LOG_INFO_L1) 
			logmsg->emit_noendl(LOG_INFO_L1, rhsmsg.str().c_str());
		else 
			logmsg->emit(LOG_INFO_L1, rhsmsg.str().c_str());
	}
	for (eqncit eq = newton_list.begin(); eq != newton_list.end(); eq++)
	{
		logmsg->emit(LOG_INFO_L2,"      equation \"%s\"...",(*eq)->get_name().c_str());
		uint offset = (*eq)->get_offset();
		
		// make sure the dependencies are up to date and have the same timestamp (=this->timestamp)
		// is now done in a separate method	
		logmsg->emit(LOG_INFO_L3,"      eqn->check_dependency_timestamps()...");
		(*eq)->check_dependency_timestamps(this->timestamp);
		
		#pragma omp parallel for schedule(dynamic, constants::newton_calc_rhs_chunk)
		for (int ii=0; ii<(int) (*eq)->get_num_variables(); ii++)
		{
			try {
			rhs[offset+ii] = (*eq)->get_newton_function(ii);
			} catch (Exception * e) { 
				e->append_info(__LINE__,__FILE__,__DATE__,__TIME__, __func__); 
				#ifdef _OPENMP 
					cout << e->get_reason() << "\n   equation:" <<this->get_name()  << "\nterminating from thread " << omp_get_thread_num() <<"\n";
					exit(1);
				#else
					NEGF_EXCEPTION(e->get_reason().c_str());
				#endif
			}
		}	
	}

	timeval tp2;
	gettimeofday(&tp2, NULL);
	this->rhs_assembly_time += tp2.tv_sec - tp1.tv_sec + 1.0e-6 * (tp2.tv_usec - tp1.tv_usec);
	
	if (verbose) logmsg->emit(LOG_INFO_L1,  "  finished (%.3g[s])",tp2.tv_sec - tp1.tv_sec + 1.0e-6 * (tp2.tv_usec - tp1.tv_usec));
	
//	} catch (Exception * e) { e->append_info(__LINE__,__FILE__,__DATE__,__TIME__, __func__); throw e; }
}
//);}


void NewtonSolver::calculate_update() throw(Exception *)
{STACK_TRACE(
	// create linear solver for solving Ax=b
	// syntax: 1st argument - A, 2nd argument - b, 3rd argument - x
	// the arguments are passed as reference, so they only have to be specified once.
	ostringstream updmsg;
	updmsg <<"   calculating update";
	#ifdef _OPENMP
		updmsg << " using " << omp_get_max_threads() << " threads and UMFPACK...";
	#endif
	if (verbose) {
		if (logmsg->get_level()==LOG_INFO_L1) 
			logmsg->emit_noendl(LOG_INFO_L1, updmsg.str().c_str());
		else 
			logmsg->emit(LOG_INFO_L1, updmsg.str().c_str());
	}

	timeval tp1;
	gettimeofday(&tp1, NULL);

	// solve Jx=F for x. the result is stored in update.
	if (verbose) logmsg->emit(LOG_INFO_L2, "      solving...");
	this->linear_solver->solve();

	timeval tp2;
	gettimeofday(&tp2, NULL);
	this->linear_solver_time += tp2.tv_sec - tp1.tv_sec + 1.0e-6 * (tp2.tv_usec - tp1.tv_usec);
	if (verbose) logmsg->emit(LOG_INFO_L1,  "  finished (%.3g[s])",tp2.tv_sec - tp1.tv_sec + 1.0e-6 * (tp2.tv_usec - tp1.tv_usec));
		
	// debug
	this->calculate_residuum();
);}


void NewtonSolver::calculate_rhs_norm() throw(Exception *)
{STACK_TRACE(
	if (verbose) logmsg->emit(LOG_INFO_L2,"   calculating RHS norms...");
	this->rhs_norm     = 0.0;
	for (uint ee = 0; ee < newton_list.size(); ee++)
	{
		Equation * eqn = newton_list[ee];
		const uint offset = eqn->get_offset();
		this->rhs_norm_eqn[ee]    = 0.0; 
		for (uint ii = 0; ii < eqn->get_num_variables(); ii++) {
			if (fabs(rhs[offset+ii])>1e50) {
				NEGF_FEXCEPTION( "%s: a RHS value was detected to be greater than 1e50.",eqn->get_name().c_str());
			}
			this->rhs_norm_eqn[ee] += rhs[offset+ii]*rhs[offset+ii];
		}
		NEGF_FASSERT( rhs_norm_eqn[ee]>=0.0, "Wrong calculation of RHS norm (%e is not >=0).",rhs_norm_eqn[ee]);
		
		this->rhs_norm    += rhs_norm_eqn[ee];
		this->rhs_norm_eqn[ee] = sqrt(rhs_norm_eqn[ee]);
	}
	this->rhs_norm    = sqrt(rhs_norm);
);}


void NewtonSolver::calculate_update_norm() throw(Exception *)
{STACK_TRACE(
	if (verbose) logmsg->emit(LOG_INFO_L2,"   calculating update norms...");
	this->update_norm  = 0.0;
	for (uint ee = 0; ee < newton_list.size(); ee++)
	{
		Equation * eqn = newton_list[ee];
		const uint offset = eqn->get_offset();
		this->update_norm_eqn[ee] = 0.0;
		uint warnings = 0;
		for (uint ii = 0; ii < eqn->get_num_variables(); ii++) {
			if (fabs(update[offset+ii])>1e100) {
				warnings++;
				if (warnings < 100) {
					logmsg->emit(LOG_WARN,"warning: when calculating the update norm, update value %d of eqn %s was more than 1e100 (%e).",
						ii, eqn->get_name().c_str(), update[offset+ii]);
				}
				if (warnings == 100) {
					logmsg->emit(LOG_WARN,"further warnings about huge update values are suppressed.");
				}
			}
			this->update_norm_eqn[ee] += update[offset+ii]*update[offset+ii];
		}
		if (warnings > 100) {
			logmsg->emit(LOG_WARN,"%d warnings in total about huge update values.",warnings);
		}
		NEGF_ASSERT(update_norm_eqn[ee]>=0.0 , "Wrong calculation of norms.");
		
		this->update_norm += update_norm_eqn[ee];
		this->update_norm_eqn[ee] = sqrt(update_norm_eqn[ee]);
	}
	this->update_norm = sqrt(update_norm);
);}


void NewtonSolver::calculate_new_variable_norm() throw(Exception *)
{STACK_TRACE(
	if (verbose) logmsg->emit(LOG_INFO_L2,"   calculating variable norms...");
	this->new_variable_norm = 0.0;
	for (uint ee = 0; ee < newton_list.size(); ee++)
	{
		Equation * eqn = newton_list[ee];
		NEGF_ASSERT(eqn->get_offset()+eqn->get_num_variables() <= nvars, "index overflow.");
		this->new_variable_norm_eqn[ee] = 0.0;
		for (uint ii = 0; ii < eqn->get_num_variables(); ii++) {
			this->new_variable_norm_eqn[ee] += eqn->get_value(ii) * eqn->get_value(ii);
		}
		this->new_variable_norm += new_variable_norm_eqn[ee];
		this->new_variable_norm_eqn[ee] = sqrt(new_variable_norm_eqn[ee]);
	}
	this->new_variable_norm = sqrt(new_variable_norm);
);}


void NewtonSolver::calculate_residuum() throw(Exception *)
{STACK_TRACE(
	if (verbose) logmsg->emit(LOG_INFO_L2,"   calculating residuum |Ax-b|...");
	uint    n		 = this->jacobian->get_size();
	int* 	icol 	 = this->jacobian->get_icol();		// length NNZ+1
	int* 	prow 	 = this->jacobian->get_prow();		// length n+1
	double* nonzeros = this->jacobian->get_nonzeros();	// length NNZ+1
	double* Ax       = new double[n];
	int fortran_index = prow[0];
	for (uint ii=0; ii < n; ii++)
	{
		Ax[ii] = 0.0;
		for (int jj=prow[ii]-fortran_index; jj < prow[ii+1]-fortran_index; jj++) {
			Ax[ii] += nonzeros[jj] * this->update[icol[jj]-fortran_index];
		}
	}
	double Axmb = 0.0;
	cout.precision(15);
	bool matrix_was_written = false;
	uint errors = 0;
	for (uint ii = 0; ii < n; ii++)
	{
		//cout << "Ax[" << ii << "]=" << Ax[ii] << ", b[" << ii << "]=" << rhs[ii] << endl;
		double tmp = Ax[ii] - this->rhs[ii];
		Axmb += tmp*tmp;
		
		// error checking
		if (tmp/(rhs[ii] + 1.0) > 1.0)
		{
			errors++;
			if (errors<100) {
				logmsg->emit(LOG_WARN,"Solution of linear system is severely incorrect: Ax[%d]=%e, rhs[%d]=%e",
									ii,Ax[ii],ii,rhs[ii]);
			}
			if (errors==100) {
				logmsg->emit(LOG_WARN,"Further warnings will be suppressed.");
			}
			if (errors>100) {
				continue;
			}
			
			bool eqn_found = false;
			for (uint jj = 0; jj < newton_list.size(); jj++) {
				if (ii-newton_list[jj]->get_offset() < newton_list[jj]->get_num_variables()) {
					logmsg->emit(LOG_WARN,"    variable has local index %d in equation \"%s\".",
						ii - newton_list[jj]->get_offset(), newton_list[jj]->get_name().c_str());
					eqn_found = true;
					break;
				}
			}
			NEGF_ASSERT(eqn_found, "Variable not found in newton list?!");
			if (!matrix_was_written)
			{
				logmsg->emit(LOG_ERROR,"writing jacobian and RHS to linear_solver_failure.[mat,rhs]*");
				this->jacobian->save_to_file("linear_solver_failure.mat");
				ofstream fout("linear_solver_failure.rhs");
				if(fout) {
					fout.precision(12);	
					for(uint jj = 0; jj < n; jj++) {
						fout << rhs[jj] << "\n"; 	
					}
					fout << "\n";
					fout.close();
				} else {
					NEGF_EXCEPTION("file linear_solver_failure.rhs could not be opened for write.");	
				}
				matrix_was_written = true;
			}
			//NEGF_EXCEPTION("Aborting.");
		}
	}
	Axmb = negf_math::sqrt(Axmb);
	if (errors!=0) logmsg->emit(LOG_WARN, "%d variables were solved incorrectly in total.", errors);
	if (verbose) logmsg->emit(LOG_INFO_L2,  "    |Ax-b|=%e or %e per variable", Axmb, Axmb/jacobian->get_size());
	delete [] Ax;
);}



/** this method is called when the Bank-Rose scheme is OFF */
void NewtonSolver::update_newton_vars(uint tstamp) throw(Exception *)
{STACK_TRACE(
	this->update_newton_vars(1.0, tstamp);
);}


/** used for Bank-Rose update; multiplies update by alpha before updating */
void NewtonSolver::update_newton_vars(double alphafac, uint tstamp) throw(Exception *)
{STACK_TRACE(
	if (alphafac!=1.0) { 
		for (uint ii = 0; ii < nvars; ii++)
			update[ii] = update[ii] * alphafac;
	}
		
	for (eqncit eq = newton_list.begin(); eq != newton_list.end(); eq++)
	{
		NEGF_ASSERT((*eq)->get_offset()+(*eq)->get_num_variables() <= nvars, "index overflow.");
		(*eq)->newton_update(&update[(*eq)->get_offset()], tstamp);  // minus sign is included in this routine
		//if (this->get_name()=="KPlevels") {
		//	cout << (*eq)->get_name() << " var 0 value=" << (*eq)->get_value(0) << endl;
		//}
	}
);}


/** Bank-Rose Newton iteration scheme <BR>
 *  see e.g. the PhD thesis of A. Scholze for the essentials <BR>
 *  The algorithm simply makes sure that the RHS norm is smaller after the update than before.
 *  This is done by multiplying the update vector with a constant alpha, 0<alpha<=1. <BR>
 *  Assumption: The RHS, Jacobian and Update of the step have already been calculated.
 */
void NewtonSolver::bank_rose_update() throw(Exception *)
{STACK_TRACE(
	const double eps       = 1e-14;
	const double alpha_min = 1e-6;
	double factor = 10.0; // factor by which kappa is increased at every iteration, see also below
	NEGF_ASSERT(this->timestamp >= this->last_converged_timestamp, 
				"timestamp < last_converged_timestamp encountered in Newton solver.");

	// store the initial RHS norm and calculate an up-to-date RHS norm
	double store_rhs = this->rhs_norm;
	this->calculate_rhs_norm();
	NEGF_ASSERT(rhs_norm >= 0.0, "invalid RHS norm.");
	if (this->rhs_norm==0.0) {
		logmsg->emit(LOG_WARN,"WARNING: RHS norm was 0.0!! Either you EXACTLY hit the solution or something was wrong!");
		return;
	}
	
	// find a good first increase of kappa
	// look at update scheme of kappa below --> like this one always starts with alpha=1.0
	double kappa_start = 1.0 / rhs_norm * 0.01;
	double kappa = 0.0;
	
	// IMPORTANT!!! revert the norms
	this->rhs_norm = store_rhs;	
	
	while (true)
	{
		// set new alpha (factor with which update vector is scaled)
		this->alpha = 1.0 / (1.0 + kappa*rhs_norm);
				
		// perform scaled update
		update_newton_vars(alpha, this->timestamp+1); // does NOT change the Newton solver's timestamp
		
		// increase Newton solver's timestamp by 1 --> in calculate_rhs(), explicit equations will be freshly computed
		this->timestamp++;
		
		// calculate new norm of RHS
		double old_rhs_norm = this->rhs_norm;
		this->calculate_rhs();
		this->calculate_rhs_norm();
		
		if (verbose) logmsg->emit(LOG_INFO_L1,"  Bank-Rose: trying alpha=%6.2g: 1-Fm+1/Fm=%7.3g (Fm+1=%.3e, kappa=%6.2e)",
										this->alpha, 1.0-rhs_norm/old_rhs_norm, rhs_norm, kappa);		
		
		// determine kappa for next iteration
		if (alpha < 0.1) factor = 1e4; // faster decrease
		kappa = max(kappa,kappa_start)*factor;
		
		// determine alpha in next iteration (if it was not converged)
		double alpha_next = 1.0 / (1.0 + kappa*old_rhs_norm);
		
		// check if update has improved solution
		if (   alpha*eps    <= 1.0 - rhs_norm/old_rhs_norm 
			|| alpha_next  <  alpha_min						
			/*|| this->timestamp - this->last_converged_timestamp <  2*/) {		
			logmsg->emit(LOG_INFO_L2, "Bank-Rose finished w/ alpha=%e and new RHS norm %e",this->alpha, this->rhs_norm);
			break;
		} else 	// update did not improve solution
		{					
			// decrease Newton solver timestamp
			this->timestamp--;
			
			// restore implicit equations (value and timetamp)
			update_newton_vars(-1.0/alpha, this->timestamp);
						
			// move back timestamp of explicit equations
			for (eqncit eq = handler.begin(); eq != handler.end(); eq++) {
				if ((*eq)->get_timestamp()==this->timestamp+1) {
					NEGF_ASSERT(!(*eq)->is_solved_in_newton(), "something went wrong.");
					logmsg->emit(LOG_INFO_L3,"      reverting %s to timestamp %d",(*eq)->get_name().c_str(),this->timestamp);
					(*eq)->set_timestamp(this->timestamp);
				}
			}
			
			// restore RHS norm
			this->rhs_norm = old_rhs_norm;
		}
	}
);}


void NewtonSolver::set_relative_err_crit(double eps) throw(Exception *)
{STACK_TRACE(
	NEGF_ASSERT(eps > 1e-15 || eps < 1, "You have chosen a bad relative error criterion.");
	this->relative_err_crit = eps;
);}


void NewtonSolver::set_absolute_err_crit(double eps) throw(Exception *)
{STACK_TRACE(
	NEGF_ASSERT(eps > 1e-15, "You have chosen a bad absolute error criterion.");
	this->absolute_err_crit = eps;
);}


/** Check convergence of the Newton iteration. Several criteria are implemented, see the descriptions. */
bool NewtonSolver::converged() const throw(Exception *)
{STACK_TRACE(
	// note that the initial RHS computation already increases the timestamp by 1
	if (timestamp - this->last_converged_timestamp < this->min_iterations + 1) {
		if (verbose)
			logmsg->emit(LOG_INFO_L1,
				"No convergence because less than %d iterations were performed.", this->min_iterations);
		return false;
	}
	
	// return this->global_convergence();  // uses Newton vector as a whole, very simple, not physical (units!)
	// return this->eqn_convergence();	   // uses norms of the variable vectors of the separate equations
	 return this->mixed_convergence();  // DESSIS convergence criterion, similar to eqn_convergence
	//return this->individual_convergence(); // uses individual variable values (no norms)
);}


/** Check convergence by testing relative and absolute error of norm of total update vector (includes all newton equations) */
bool NewtonSolver::global_convergence() const throw(Exception *)
{STACK_TRACE(
	if (   update_norm / new_variable_norm > this->relative_err_crit	&& update_norm > absolute_err_crit ) {
		logmsg->emit(LOG_INFO_L1,"No global convergence: update_norm=%7.3e, new_norm=%7.3e",
					update_norm, new_variable_norm);
		return false;
	} else {
		return true;
	}
);}


/** Check convergence by testing relative and absolute error of norm of update vector for each Newton Equation */
bool NewtonSolver::eqn_convergence() const throw(Exception *)
{STACK_TRACE(
	for (uint ee = 0; ee < newton_list.size(); ee++) 
	{
		if ( new_variable_norm_eqn[ee]==0.0 &&  update_norm_eqn[ee]!=0.0)
			return false;
		if (bank_rose) {
			if (   update_norm_eqn[ee]/alpha / new_variable_norm_eqn[ee] > this->relative_err_crit
				&& rhs_norm_eqn[ee] 			 						 > this->absolute_err_crit ) {
				if (verbose) {
					logmsg->emit(LOG_INFO_L1,"Neither absolute nor relative convergence for eqn \"%s\".",
								newton_list[ee]->get_name().c_str());
					logmsg->emit(LOG_INFO_L2,"Equation update_norm=%7.3e, new_norm=%7.3e, alpha=%7.3e",
								update_norm_eqn[ee], new_variable_norm_eqn[ee], alpha);
				}
				return false;
			}
		} else {
			if (   update_norm_eqn[ee] / new_variable_norm_eqn[ee] > this->relative_err_crit
				&& rhs_norm_eqn[ee] 							   > this->absolute_err_crit ){
				if (verbose) {
					logmsg->emit(LOG_INFO_L1,"Neither absolute nor relative convergence for eqn \"%s\".",
								newton_list[ee]->get_name().c_str());
					logmsg->emit(LOG_INFO_L2,"Equation update_norm=%7.3e, new_norm=%7.3e",
								update_norm_eqn[ee], new_variable_norm_eqn[ee]);
				}
				return false;
			}
		}
	}
	return true;
);}


/** Check convergence by a mixed relative/absolute error of norm of update vector for each Newton Equation. <BR>
 *  For the criterion see Steiger et al., J. Comp. Electron. 7, 509 (2008). It involves a reference size as a measure
 *  of the variable size for each Equation. Apart from the reference size, only the relative error criterion is used. */
bool NewtonSolver::mixed_convergence() const throw(Exception *)
{STACK_TRACE(
	for (uint ee = 0; ee < newton_list.size(); ee++) 
	{
		double ref = std::sqrt(newton_list[ee]->get_num_variables()) * newton_list[ee]->get_reference_size();
		if ( new_variable_norm_eqn[ee]==0.0 &&  update_norm_eqn[ee]!=0.0)
			return false;
		if (bank_rose) {
			if (   update_norm_eqn[ee]/alpha / (new_variable_norm_eqn[ee]+ref) > this->relative_err_crit ) {
				if (verbose) {
					logmsg->emit(LOG_INFO_L1,"SDEVICE convergence criterion for eqn \"%s\" failed.",
								newton_list[ee]->get_name().c_str());
					logmsg->emit(LOG_INFO_L1,"(|update|=%7.2e/alpha=%.2g)/(|new|=%7.2e +ref=%7.2e)=%7.3e > %.1e",
								newton_list[ee]->get_name().c_str(),
								update_norm_eqn[ee], alpha, new_variable_norm_eqn[ee], ref, 
								update_norm_eqn[ee]/alpha / (new_variable_norm_eqn[ee]+ref),
								this->relative_err_crit);
								
					
					const uint offset = newton_list[ee]->get_offset();
					double max_update = 0.0;
					int max_update_var = -1;
					for (uint ii = 0; ii < newton_list[ee]->get_num_variables(); ii++) {
						if (negf_math::abs(max_update) < negf_math::abs(update[offset+ii])) {
							max_update = update[offset+ii];
							max_update_var = ii;
						}
					}
					logmsg->emit(LOG_INFO_L1,"Maximum update was %e for variable %d",max_update, max_update_var);
				}
				return false;
			}
		} else {
			if ( update_norm_eqn[ee] / (new_variable_norm_eqn[ee]+ref) > this->relative_err_crit ){
				if (verbose) {
					logmsg->emit(LOG_INFO_L1,"SDEVICE convergence criterion for eqn \"%s\" failed.",
								newton_list[ee]->get_name().c_str());
					logmsg->emit(LOG_INFO_L1,"|update|=%7.2e/(|new|=%7.2e +ref=%7.2e)=%7.3e > %.1e",
								newton_list[ee]->get_name().c_str(),
								update_norm_eqn[ee], new_variable_norm_eqn[ee], ref, 
								update_norm_eqn[ee] / (new_variable_norm_eqn[ee]+ref),
								this->relative_err_crit);
				}
				return false;
			}
		}
	}
	return true;
);}


/** Check convergence by testing a mixed relative/absolute error of each and every individual variable. <BR>
 *  It involves a reference size as a measure  of the variable size for each Equation.
 *  Apart from the reference size, only the relative error criterion is used. */
bool NewtonSolver::individual_convergence() const throw(Exception *)
{STACK_TRACE(
	
	for (uint ee = 0; ee < newton_list.size(); ee++) 
	{
		Equation * eqn = newton_list[ee];
		double ref = eqn->get_reference_size();
		uint offset = eqn->get_offset();
		if ( new_variable_norm_eqn[ee]==0.0 &&  update_norm_eqn[ee]!=0.0)
			return false;
			
		for (uint ii = 0; ii < eqn->get_num_variables(); ii++) 
		{
			double new_variable_normm = fabs(eqn->get_value(ii));   // cheaper than sqrt(x*x)
			double update_normm = fabs(update[offset+ii]); // cheaper than sqrt(x*x)
			//double sqrt_N = std::sqrt(eqn->get_num_variables());
			if (bank_rose) {
				if (   update_normm/alpha / (new_variable_normm+ref) > this->relative_err_crit ) {
					if (verbose) {
						logmsg->emit(LOG_INFO_L1,"    convergence failed in eqn \"%s\" var %d (%6.2e>%6.2e, update=%6.2e, new_val=%6.2e, ref=%6.2e).", 
									eqn->get_name().c_str(), ii, update_normm/alpha / (new_variable_normm+ref), this->relative_err_crit,
									update_normm, new_variable_normm, ref);
					}
					return false;
				}
			} else {
				if ( update_normm / (new_variable_normm+ref) > this->relative_err_crit ){
					if (verbose) {
						logmsg->emit(LOG_INFO_L1,"    convergence failed in eqn \"%s\" var %d (%6.2e>%6.2e, update=%6.2e, new_val=%6.2e, ref=%6.2e).", 
									eqn->get_name().c_str(), ii, update_normm/alpha / (new_variable_normm+ref), this->relative_err_crit,
									update_normm, new_variable_normm, ref);
					}
					return false;
				}
			}
		}
	}
	return true;
);}


#ifndef NEWTON_WITHOUT_OUTPUTDATA	
/** The equations to be written into file after each iteration step are only .dat-equations */
void NewtonSolver::add_outputdata(OutputData * outputdata) throw(Exception *)
{STACK_TRACE(
	NEGF_ASSERT(outputdata!=NULL, "Tried to assign null pointer.");
	NEGF_ASSERT(find(debug_data.begin(), debug_data.end(), outputdata)
				==debug_data.end(), "ouput data object was already assigned.");
	for (uint ii = 0; ii < outputdata->get_num_dat_equations(); ii++) {
		const Equation * eqn = outputdata->get_dat_equation(ii);
		if (find(this->handler.begin(), this->handler.end(), eqn) ==this->handler.end()) 
		{
			if (eqn->is_solved_in_newton()) {
				NEGF_FEXCEPTION("Equation \"%s\" to be written as output was not found in the Newton solver.",
										eqn->get_name().c_str());
			} else {
			}
		}
	}
					
	debug_data.push_back(outputdata);
);}
#endif


/** Have some output about the dependencies */
void NewtonSolver::display_dependency_info() const throw(Exception *)
{STACK_TRACE(
	//ostringstream msg;
	for (uint ii = 0; ii < this->handler.size(); ii++)
	{
		NEGF_ASSERT(handler[ii]!=NULL, "null pointer encountered.");
		
		logmsg->emit(LOG_INFO,"%s has the following direct dependencies:", handler[ii]->get_name().c_str());
		const vector<Equation *> deps = handler[ii]->get_dependencies();
		for (uint jj = 0; jj < deps.size(); jj++) {
			NEGF_ASSERT(deps[jj]!=NULL, "null pointer encountered.");
			logmsg->emit(LOG_INFO,"    %s", deps[jj]->get_name().c_str());
		}
		
		logmsg->emit(LOG_INFO,"%s depends on the following Newton variables:", handler[ii]->get_name().c_str());
		NEGF_ASSERT(handler[ii]->is_newton_var_list_ready(), "set up newton var list first.");
		const vector<Equation *> newtons = handler[ii]->get_newton_var_list();
		for (uint jj = 0; jj < newtons.size(); jj++) {
			NEGF_ASSERT(newtons[jj]!=NULL, "null pointer encountered.");
			logmsg->emit(LOG_INFO,"    %s", newtons[jj]->get_name().c_str());
		}
		
		//msg.clear(); msg << "   ";
		logmsg->emit(LOG_INFO,"%s depends on all these equations:", handler[ii]->get_name().c_str());
		NEGF_ASSERT(handler[ii]->is_all_eqns_ready(), "set up all_eqns var list first.");
		const vector<Equation *> alleqs = handler[ii]->get_all_eqns();
		for (uint jj = 0; jj < alleqs.size(); jj++) {
			NEGF_ASSERT(alleqs[jj]!=NULL, "null pointer encountered.");
			//logmsg->emit(LOG_INFO,"    %s", alleqs[jj]->get_name().c_str());
			//msg << alleqs[jj]->get_name() << "     ";
			cout << alleqs[jj]->get_name() << "     ";
		}
		cout << endl;
		//logmsg->emit(LOG_INFO,"%s", msg.str().c_str());
	}
);}


