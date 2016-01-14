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
#include "Equation.h"
#ifdef _OPENMP
	#include <omp.h>
#endif
using namespace negf;


/** Common constructor for all derived classes */
Equation::Equation():
	offset_ready(false),
	all_eqns_ready(false),
	newton_var_list_ready(false),
	reference_size(1.0),
	max_update(1e100),
	timestamp(0),
	min_dep_timestamp(0),
	name("unknown"),
	max_nonzeros_per_line(constants::eqn_array_size),
	nonzeros_set(false)
{STACK_TRACE(    
);}


double Equation::get_value(uint idx) const throw(Exception *)
{STACK_TRACE(
	NEGF_ASSERT(idx < number_of_variables, "Invalid variable index.");
	NEGF_FASSERT(idx < current_variable_values.size(),
				"Variable %d of equation \"%s\" has not yet been computed.",idx,this->get_name().c_str());
	return current_variable_values[idx]; 
);}



void Equation::set_values(const vector<double> & vec, uint tstamp)
{STACK_TRACE(
	this->current_variable_values = vec;
	this->timestamp = tstamp;
);}


/** Compute values (only possible for explicit equations). The real work will be done by compute_value(uint line).
 *  @param tstamp the desired timestamp. Only if the Equation's timestamp is smaller will the values be computed.
 */
void Equation::compute_values(uint tstamp) throw(Exception *)
{// OpenMP does not go together with our STACK_TRACE macro
	try {
	
	if (this->get_timestamp()>=tstamp 
	   && this->current_variable_values.size() == this->get_num_variables()) { // nothing to do
		return;
	}
	logmsg->emit(LOG_INFO_L3,"      computing equation \"%s\", timestamp %d",
						this->get_name().c_str(), tstamp);

	if (this->current_variable_values.size() != this->get_num_variables()){
		this->current_variable_values.clear();
		this->current_variable_values.resize(this->get_num_variables(), 0.0);
	}
		
	// make sure time stamps of dependencies are correct; update if necessary
	this->check_dependency_timestamps(tstamp);
	
	} catch (Exception * e) { e->append_info(__LINE__,__FILE__,__DATE__,__TIME__, __func__); throw(e); }

	// now we're safe. note that compute_value(ii) is a const function and therefore easily parallelizable
	#pragma omp parallel for schedule(dynamic, constants::eqn_comp_vals_chunk)
	for (int ii = 0; ii < (int) this->get_num_variables(); ii++) {
		try {
		this->current_variable_values[ii] = this->compute_value(ii);
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
		
	this->set_timestamp(tstamp);
}


void Equation::newton_update(const double * update, uint tstamp)
{STACK_TRACE(
	NEGF_ASSERT(this->current_variable_values.size() == this->get_num_variables(), "Initialize variable vector first.");
	uint max_updates = 0;
	for (uint ii = 0; ii < this->get_num_variables(); ii++) {
		double upd = negf_math::sign(update[ii]) * min(negf_math::abs(update[ii]), this->max_update);
		if (upd!=update[ii]) {
			max_updates++;
		}
		this->current_variable_values[ii] -= upd;  // NOTE THE MINUS SIGN!
	}
	if (max_updates!=0) {
		logmsg->emit(LOG_WARN,"Warning: in eqn \"%s\" the update of %d variables had to be limited to +-%g.", 
				this->get_name().c_str(), max_updates, this->max_update);
	}
	timestamp = tstamp;
);}


void Equation::handshake() throw(Exception *)
{STACK_TRACE(
	all_eqns.clear();
	all_eqns.push_back(this);
	this->handshake(all_eqns);
	all_eqns_ready = true;
);}


void Equation::handshake(vector<Equation *> & the_list)
{STACK_TRACE(
	for (uint ii = 0; ii < this->dependencies.size(); ++ii)
	{
		Equation * eqn = dependencies[ii];
		
		NEGF_FASSERT(eqn != NULL, "Dependency %d of equation %s was not allocated so far.",
					ii, this->get_name().c_str());
		if (find(the_list.begin(), the_list.end(), eqn) == the_list.end())
		{
			the_list.push_back(eqn);
			eqn->handshake(the_list);
		}
	}
	
	/*
	for (eqcit dep_eqn = dependencies.begin(); dep_eqn != dependencies.end(); ++dep_eqn)
	{
		NEGF_FASSERT(*dep_eqn != 0, "A dependency was not allocated so far.");
		if (find(the_list.begin(), the_list.end(), *dep_eqn) == the_list.end())
		{
			the_list.push_back(*dep_eqn);
			(*dep_eqn)->handshake(the_list);
		}
	}
	*/
);}


uint Equation::get_offset() const throw(Exception *)
{STACK_TRACE(
	if (!this->is_solved_in_newton())
		NEGF_EXCEPTION("Cannot get global offset of equation because it is not solved in Newton.");
	if (!this->is_offset_ready())
		NEGF_EXCEPTION("Offset has not been set so far.");
	return offset;
);}


void Equation::set_newton_var_list() throw(Exception *)
{STACK_TRACE(
	if (!this->is_all_eqns_ready())
		NEGF_EXCEPTION("Please handshake first.");

	newton_var_list.clear();

	for (eqit eqn = all_eqns.begin(); eqn != all_eqns.end(); ++eqn)
	{
		if ((*eqn)->is_solved_in_newton()) {
			newton_var_list.push_back(*eqn);
		}
	}

	// sort the list ascending w.r.t. the equation's offset
	for (uint ii = 0; ii < newton_var_list.size(); ii++)
		if ( !(newton_var_list[ii]->is_offset_ready()) )
			NEGF_EXCEPTION("Please set the offset of all Newton equations first.");
	Equation * tmp_eqn;
	for (uint ii = 0; ii < newton_var_list.size(); ii++)
	{
		uint place = ii;
		while (place!=0)
		{
			if (newton_var_list[place]->get_offset() < newton_var_list[place-1]->get_offset())
			{ // swap places
				tmp_eqn = newton_var_list[place];
				newton_var_list[place] = newton_var_list[place-1];
				newton_var_list[place-1] = tmp_eqn;
				place--;
			}
			else
			{
				if  (newton_var_list[place]->get_offset() == newton_var_list[place-1]->get_offset()) {
					logmsg->emit( LOG_INFO,"equation %X (%s), place=%d, newton_var_list[place]=%X (%s), newton_var_list[place-1]=%X (%s)",
						this, this->get_name().c_str(), place, newton_var_list[place], 
						newton_var_list[place]->get_name().c_str(), newton_var_list[place-1], newton_var_list[place-1]->get_name().c_str() );
					NEGF_FEXCEPTION("Equations with the same offset (%d) were encountered!",
						newton_var_list[place]->get_offset());
				}
				else {
					break; // we're done
				}
			}
		}
	}
	this->newton_var_list_ready = true;
);}


const vector<Equation *> & Equation::get_newton_var_list() const throw(Exception *)
{STACK_TRACE(
	if (!this->is_newton_var_list_ready()) NEGF_EXCEPTION("newton_var_list was not set up yet!");
	return this->newton_var_list;
);}

const vector<Equation *> & Equation::get_all_eqns() const throw(Exception *)
{STACK_TRACE(
	if (!this->is_all_eqns_ready()) NEGF_EXCEPTION("all_eqns was not yet set up!");
	return this->all_eqns;
);}


void Equation::set_reference_size(double size) throw(Exception *)
{STACK_TRACE(
	NEGF_ASSERT(size>0.0, "Reference size must be strictly positive.");
	this->reference_size = size;
);}

void Equation::set_maximum_update(double max_update_) throw(Exception *)
{STACK_TRACE(
	NEGF_ASSERT(max_update_>0.0, "Maximum update must be strictly positive.");
	this->max_update = max_update_;
);}


/** Get the sparsity pattern of an entire line in the Jacobian. Methodology is similar to get_derivative.
 * @param line which line to calculate (which variable of the equation, not the global line)
 * @param nonzeros (result) stores the length of icols (number of nonzero entries)
 * @param icols (rsult) list of long integers containing the nonzero column indices of that line
 */
void Equation::get_sparsity(const uint & line, uint & nonzeros, uint icols[]) const throw(Exception *)
{STACK_TRACE(
	if (!this->is_solved_in_newton())
		NEGF_EXCEPTION("get_sparsity shall only be called for Newton variables.");
	if (!this->is_newton_var_list_ready())
		NEGF_EXCEPTION("Please set up the list of newton variables associated with this equation first.");

	nonzeros = 0;

	uint nonzeros_part;
	uint icol_part[constants::eqn_array_size];

	// assume that newton_var_list is ordered w.r.t. the offset of the equations
	for  (eqcit eqn = this->newton_var_list.begin(); eqn != this->newton_var_list.end(); ++eqn)
	{
		NEGF_FASSERT((*eqn)!=NULL, "null pointer encountered in newton_var_list of eqn %s!",this->get_name().c_str());
		nonzeros_part = 0;
		this->get_sparsity(line, *eqn, nonzeros_part, icol_part); // returns sorted icol_part
		uint other_eqn_offset = (*eqn)->get_offset();
		for (uint ii = 0; ii < nonzeros_part; ii++)
		{
			// skip checking duplicate entries of icol_part and icols because it shouldnt be possible
			NEGF_ASSERT(nonzeros < constants::eqn_array_size, "Overflow. Please make constants::eqn_array_size bigger.");
			icols[nonzeros] = icol_part[ii] + other_eqn_offset;
			nonzeros++;
		}
	}
	NEGF_FASSERT(nonzeros!=0, "Jacobian line %d of eqn %s is a priori entirely zero, which generates a noninvertible Jacobian.",
			line, this->get_name().c_str());
);}


/** Get the sparsity pattern of a line in the Jacobian w.r.t a certain newton variable
 *
 * This method is usually called from get_sparsity.
 * Note that the offset (global position in the newton vector) is NOT added to the indices.
 * @param line which line to calculate (which variable of the equation, not the global line)
 * @param nn which Newton variable to look at
 * @param icol_part a list of integers containing the nonzero column indices of that line w.r.t. Equation nn
 *        note that the range of these indices is within [0..nn->num_vars-1], i.e. local and NOT global
 */
void Equation::get_sparsity(const uint & line, const Equation * nn, uint & nonzeros, uint icol_part[]) const
{STACK_TRACE(
	uint nonzeros1;
	uint temp1 [constants::eqn_array_size];
	uint nonzeros2;
	uint temp2 [constants::eqn_array_size]; // nn->get_num_variables()

	nonzeros = 0;

	// -----------------------------------------------------
	// dependence of own variables on the newton variable nn
	// -----------------------------------------------------
	if (!this->is_solved_in_newton() && this==nn)
		NEGF_EXCEPTION("Inconsistency if this eqn is solved in the Newton or not.");
	if (this->is_solved_in_newton() && this!=nn)
		{ } // do nothing, only look at the dependencies
	if (!this->is_solved_in_newton() && this!=nn)
		{ }
		// in this case, we have an explicit equation which is not part of the newton and 
		// hence has no entry there --> do nothing, only look at dependencies
	if (this->is_solved_in_newton() && this==nn) {
		this->get_dependence_sparsity(line, nn, nonzeros, icol_part);
	}


	// ---------------------------------------------------------------------------
	// dependence of this equation's direct dependencies on the newton variable nn
	// ---------------------------------------------------------------------------
	for (eqcit dep = dependencies.begin(); dep != dependencies.end(); ++dep)
	{
		if (*dep == nn)
		{
			nonzeros2 = 0;
			this->get_dependence_sparsity(line, *dep, nonzeros2, temp2);
			this->add_idcs(nonzeros, icol_part, nonzeros2, temp2); // add temp2 to icol_part without duplicates
		}

		// do the recursion if:
		//  1. the dependency is NOT a Newton variable (i.e., there exists an explicit form!) and 
		//  2. the newton_var_list of the dependency contains nn (otherwise we'd have to work for nothing)
		if ( find(this->newton_var_list.begin(), this->newton_var_list.end(), *dep) == this->newton_var_list.end() )
		{
			if (!(*dep)->is_newton_var_list_ready())
				NEGF_FEXCEPTION("A non-newton equation (%s) was found where newton_var_list was not yet set up.",
						(*dep)->get_name().c_str());
			const vector<Equation *> & dep_newton_var_list = (*dep)->get_newton_var_list();
			if ( find(dep_newton_var_list.begin(), dep_newton_var_list.end(), nn) != dep_newton_var_list.end() )
			{
				nonzeros1 = 0;
				this->get_dependence_sparsity(line, *dep, nonzeros1, temp1);
				for (uint ii = 0; ii < nonzeros1; ii++)	// for all variables ii of the dependency that this line depends on...
				{
					nonzeros2 = 0;
					(*dep)->get_sparsity(temp1[ii], nn, nonzeros2, temp2);	// get the dependence of ii on the newton equation
					this->add_idcs(nonzeros, icol_part, nonzeros2, temp2); // add temp2 to icol_part without duplicates
				}
			}
		}
	}

	// sort!
	uint tmp_idx;
	for (uint ii=0; ii < nonzeros; ii++)
	{
		uint place = ii;
		while (place!=0)
		{
			if (icol_part[place] < icol_part[place-1]) // swap places
			{
				tmp_idx = icol_part[place];
				icol_part[place] = icol_part[place-1];
				icol_part[place-1] = tmp_idx;
				place--;
			}
			else
			{
				if  (icol_part[place] == icol_part[place-1]) {
					NEGF_FEXCEPTION("Index %d was encountered twice in icols! Equation:%s, Dependency:%s",
						icol_part[place],this->get_name().c_str(),nn->get_name().c_str());
				}
				else {
					break; // we're done
				}
			}
		}
	}
);}


/** Adds a list of indices to an existing one. The first elements (idcs1[0], idcs2[0]) contain the number of list items.
 * @param indices1 index array TO which indices2 is added
 * @param indices2 indices to add
 */
void Equation::add_idcs(uint & nonzeros1, uint idcs1[], uint & nonzeros2, uint idcs2[]) const
{STACK_TRACE(
	bool flag;
	for (uint ii = 0; ii < nonzeros2; ++ii)
	{
		// look if that index is already existing in array 1...
		flag = false;
		for (uint jj = 0; jj < nonzeros1; ++jj)
		{
			if (idcs1[jj]==idcs2[ii])
			{
				flag = true;
				break;
			}
		}
		// if not, add it at the end
		if (!flag)
		{
			NEGF_ASSERT(nonzeros1 < constants::eqn_array_size, "Overflow. Please make constants::eqn_array_size bigger.");
			idcs1[nonzeros1] = idcs2[ii];
			nonzeros1++;
		}
	}
);}


/** Adds a list of indices/coefficients to an existing one. The first elements (indices1[0], indices2[0]) contain the number of list items.
 * @param nonzeros1 (result) number of used entries in indices1, coeffs1
 * @param indices1 (result) index array TO which indices2 is added
 * @param coeffs1 (result) coefficient array TO which coeffs2 is added
 * @param nonzeros2 number of used entries in indices2, coeffs2
 * @param indices2 indices to add
 * @param coeffs2 coefficients to add
 */
void Equation::add_coeff(uint & nonzeros1, uint idcs1[], double cffs1[], uint & nonzeros2, uint idcs2[], double cffs2[]) const
{STACK_TRACE(
	bool flag;
	for (uint ii = 0; ii < nonzeros2; ++ii)
	{
		// look if that index is already existing in array 1...
		flag = false;
		for (uint jj = 0; jj < nonzeros1; ++jj)
		{
			if (idcs1[jj]==idcs2[ii])
			{
				flag = true;
				cffs1[jj] += cffs2[ii];	
				// not putting a break here allows multiple coefficients with the same index to be added
			}
		}
		if (!flag)
		{
			NEGF_ASSERT(nonzeros1 < constants::eqn_array_size, "Overflow. Please make constants::eqn_array_size bigger.");
			idcs1[nonzeros1] = idcs2[ii];
			cffs1[nonzeros1] = cffs2[ii];
			nonzeros1++;
		}
	}
);}


/** Compute a line of the Jacobian matrix \f$ \partial F^i / \partial x_j \f$ belonging to this Equation. <BR>
 *  The Jacobian is stored in a sparse format (indices, entries). <BR>
 *  The Equation MUST be a Newton solution variable, otherwise an error is thrown.
 *
 *  The method first makes a list of all solution variables that the equation depends on
 *  and calculates the ***total*** derivative w.r.t. these variables.
 *  If a dependency of the equation is not a Newton solution variable, the chain rule is applied.
 * 
 *  @param line     which variable (0...number_of_variables-1, i.e. without offset)
 *  @param nonzeros (result) stores the number of nonzeros in the Jacobian at that line.
 *  @param indices  (result) stores the column indices of the nonzero entries
 *  @param coeffs   (result) stores the values of the nonzero entries */
void Equation::get_newton_derivatives(uint line, uint & nonzeros, uint indices [], double coeffs []) const throw(Exception *)
{STACK_TRACE(
	NEGF_FASSERT(line<this->get_num_variables(), "invalid line (%d >= get_num_vars(%d)).",line,this->get_num_variables());
	if (!this->is_solved_in_newton()) {
		NEGF_FEXCEPTION("Attempting to calculate a line in the Jacobian corresponding to a quantity not solved for (\"%s\").",this->get_name().c_str());
	}
	if (!this->is_newton_var_list_ready()) {
		NEGF_FEXCEPTION("Please set up the newton_var_list of this equation (\"%s\")  first.",this->get_name().c_str());
	}

	// make sure time stamps are correct; update if necessary; afterwards, timetamp is not used anymore
	// OBSOLETE, is now done separately (better because this can become parallelized)
	// this->check_dependency_timetamps(tstamp);

	nonzeros = 0;

	// set up temporary index and coefficient arrays
	// used array size will be less or equal than the total number of nonzeros on this line of the jacobian
	// this gives huge performance improvement!
	uint arraysize = (this->nonzeros_was_set()) ? this->get_nonzeros(line) : constants::eqn_array_size;
	arraysize = arraysize + 1;	// must have minimum size 1!!
	uint   localnonzeros;
	uint   localindices[arraysize];
	double coeff[arraysize];
	/*
	uint *   localindices = new uint  [constants::eqn_array_size];
	double * coeff 		  = new double[constants::eqn_array_size];
	*/
	
	// note: it is assumed that the equations in newton_var_list are ascendingly ordered w.r.t. their offset
	for  (eqcit eqn = this->newton_var_list.begin(); eqn != this->newton_var_list.end(); ++eqn)
	{
		// compute derivatives w.r.t. newton variable *var; resulting index array has local indices
		localnonzeros = 0;
		if (constants::eqn_debug) {
			for (uint ii = 0; ii < min((uint)1000, arraysize); ii++) {
				localindices[ii] = 8888888;
				coeff[ii] = -888888.888;
			}
		}
		this->newton_derivative(line, *eqn, localnonzeros, localindices, coeff); // returns sorted localindices
		if (constants::eqn_debug) {
			for (uint ii = 0; ii < localnonzeros; ii++) {
				NEGF_FASSERT(localindices[ii]!=8888888, 
					"eqn %s: localindices array element %d was not assigned (der w.r.t. %s, localnonzeros=%d).",
					this->get_name().c_str(), ii, (*eqn)->get_name().c_str(), localnonzeros);
				NEGF_FASSERT(coeff[ii] != -888888.888, 
					"eqn %s: coeff array element %d was not assigned (der w.r.t. %s, localnonzeros=%d).",
					this->get_name().c_str(), ii, (*eqn)->get_name().c_str(), localnonzeros);
			}
		}
		uint other_eqn_offset = (*eqn)->get_offset();

		// add the stuff for this newton variable to the result arrays
		// we skip the duplicate checking because uniqueness should be guaranteed by offset
		for (uint ii = 0; ii < localnonzeros; ii++)
		{
			NEGF_ASSERT(nonzeros < arraysize, "Overflow of indices/coeffs array.");
			indices[nonzeros] = localindices[ii] + other_eqn_offset;
			
			// error checking
			if (std::isnan(coeff[ii])) {
				logmsg->emit( LOG_ERROR,"Not-A-Number Jacobian entry detected in \"%s\", derivative w.r.t. \"%s\"",
							this->get_name().c_str(), (*eqn)->get_name().c_str() );
				NEGF_EXCEPTION("Aborting.");
			}
			if (std::isinf(coeff[ii])) {
				logmsg->emit( LOG_ERROR,"Infinite Jacobian entry detected in \"%s\", derivative w.r.t. \"%s\"",
							this->get_name().c_str(), (*eqn)->get_name().c_str() );
				NEGF_EXCEPTION("Aborting.");
			}
			
			coeffs[nonzeros] = coeff[ii];
			nonzeros++;
		}
	}
	
	/*
	delete [] localindices;
	delete [] coeff;
	*/
	
	// check that the line is not entirely zero
	for (uint ii = 0; ii < nonzeros; ii++) {
		if (coeffs[ii]!=0.0) {
			return;
		}
	}
	NEGF_FEXCEPTION("The Jacobian line %d (%d nonzeros) of equation %s was found to be entirely zero.",
		line, nonzeros, this->get_name().c_str() );
);}


/** Calculate the derivative of an equation w.r.t a given Newton solution variable
 * 
 * This method is usually called by get_newton_derivatives.
 * If a dependency of the equation is not a Newton solution variable, the chain rule is applied.
 * If a dependency of the equation is a Newton solution variable but not the wanted, return zero
 * If a dependency of the equation is the wanted Newton solution variable, return the partial derivatives
 * @param line which line to calculate (which variable of the equation, not the global line)
 * @param nvar the Newton solution Equation
 * @param localnonzeros (result) stores the number of nonzero entries
 * @param localindices  (result) stores the indices with nonzero entries
 * @param localcoeffs   (result) stores the nonzero entries
 */
void Equation::newton_derivative(const uint & line, const Equation * nvar, uint & localnonzeros, 
								 uint localindices [], double localcoeffs []) const
{STACK_TRACE(
	NEGF_ASSERT(this->is_solved_in_newton(), "Method must only be called if eqn is included in Newton solver.");
	NEGF_ASSERT(this->is_newton_var_list_ready(), "Please prepare the list of newton variables first.");
	
	localnonzeros = 0;
	
	uint   nonzeros1;
	// used array size will be less or equal than the total number of nonzeros on this line of the jacobian
	// this gives huge performance improvement!
	uint arraysize = (this->nonzeros_was_set()) ? this->get_nonzeros(line) : constants::eqn_array_size;
	arraysize = arraysize + 1;	// must have minimum size 1!!
	uint   tmpidcs1[arraysize];
	double tmpcoeffs1[arraysize];

	// -------------------------------------------------------
	// dependence of own variables on the newton variable nvar
	// -------------------------------------------------------
	if (this==nvar) {
		this->get_newton_direct_derivatives(line, nvar, localnonzeros, localindices, localcoeffs);
		
		// error checking
		if (constants::eqn_debug) {
		for (uint ii = 0; ii < localnonzeros; ii++)
		{
			if (std::isnan(localcoeffs[ii])) {
				logmsg->emit( LOG_ERROR,"Not-A-Number Jacobian entry detected in \"%s\", direct derivative w.r.t. itself",
							this->get_name().c_str() );
				NEGF_EXCEPTION("Aborting.");
			}
			if (std::isinf(localcoeffs[ii])) {
				logmsg->emit( LOG_ERROR,"Infinite Jacobian entry detected in \"%s\", direct derivative w.r.t. itself",
							this->get_name().c_str() );
				NEGF_EXCEPTION("Aborting.");
			}
		}
		}
	} else {
		// do nothing, derivative w.r.t. nvar only through dependencies
	} // do nothing, derivative w.r.t. nvar only through dependencies

	// -----------------------------------------------------------------------------
	// dependence of this equation's direct dependencies on the newton variable nvar
	// -----------------------------------------------------------------------------
	for (eqcit dep = dependencies.begin(); dep != dependencies.end(); ++dep)
	{
		if (*dep==nvar)
		{
			nonzeros1 = 0;
			this->get_newton_direct_derivatives(line, *dep, nonzeros1, tmpidcs1, tmpcoeffs1);
			
			// error checking
			if (constants::eqn_debug) {
			for (uint ii = 0; ii < nonzeros1; ii++)
			{
				if (std::isnan(tmpcoeffs1[ii])) {
					logmsg->emit( LOG_ERROR,"Not-A-Number Jacobian entry detected in \"%s\", direct derivative w.r.t. \"%s\"",
								this->get_name().c_str(), (*dep)->get_name().c_str() );
					NEGF_EXCEPTION("Aborting.");
				}
				if (std::isinf(tmpcoeffs1[ii])) {
					logmsg->emit( LOG_ERROR,"Infinite Jacobian entry detected in \"%s\", direct derivative w.r.t. \"%s\"",
								this->get_name().c_str(), (*dep)->get_name().c_str() );
					NEGF_EXCEPTION("Aborting.");
				}
			}
			}
			
			add_coeff(localnonzeros, localindices, localcoeffs, nonzeros1, tmpidcs1, tmpcoeffs1);
		}

		// do the recursion (chain rule!) if:
		//  1. the dependency is NOT a Newton variable (i.e., there exists an explicit form!) and 
		//  2. the newton_var_list of the dependency contains nvar (otherwise we'd have to work for nothing)
		if ( find(this->newton_var_list.begin(), this->newton_var_list.end(), *dep) == this->newton_var_list.end() )
		{
			NEGF_FASSERT((*dep)->is_newton_var_list_ready(),
						"In equation %s newton_var_list was not yet set up.", (*dep)->get_name().c_str());
			const vector<Equation * > & dep_newton_var_list = (*dep)->get_newton_var_list();
			if ( find(dep_newton_var_list.begin(), dep_newton_var_list.end(), nvar) != dep_newton_var_list.end())
			{
				NEGF_FASSERT( !(*dep)->is_solved_in_newton(),
						"Dependency %s of eqn %s was marked as solved in newton even though it was not in the newton var list.",
						(*dep)->get_name().c_str(), this->get_name().c_str() );
	
				nonzeros1 = 0;
				this->get_newton_direct_derivatives(line, *dep, nonzeros1, tmpidcs1, tmpcoeffs1);
				
				// error checking
				if (constants::eqn_debug) {
				for (uint ii = 0; ii < nonzeros1; ii++)
				{
					if (std::isnan(tmpcoeffs1[ii])) {
						logmsg->emit( LOG_ERROR,"Not-A-Number Jacobian entry detected in \"%s\", direct derivative w.r.t. \"%s\"",
									this->get_name().c_str(), (*dep)->get_name().c_str() );
						NEGF_EXCEPTION("Aborting.");
					}
					if (std::isinf(tmpcoeffs1[ii])) {
						logmsg->emit( LOG_ERROR,"Infinite Jacobian entry detected in \"%s\", direct derivative w.r.t. \"%s\"",
									this->get_name().c_str(), (*dep)->get_name().c_str() );
						NEGF_EXCEPTION("Aborting.");
					}
				}
				}
				
				uint   nonzeros2 = 0;
				uint   tmpidcs2[(*dep)->get_max_nonzeros()+1];
				double tmpcoeffs2[(*dep)->get_max_nonzeros()+1];
			
				for (uint jj = 0; jj < nonzeros1; jj++)
				// for all vars of the dependency where the derivative of the equation is not zero...
				{
					 // compute da(b,c)/db; can only be done for explicit equations!!!
					nonzeros2 = 0;
					(*dep)->all_direct_derivatives(tmpidcs1[jj], nvar, nonzeros2, tmpidcs2, tmpcoeffs2);
					NEGF_FASSERT(nonzeros2 < (*dep)->get_max_nonzeros()+1, 
									"%s: overflow of arrays (nonzeros2<%s->get_max_nonzeros()+1)..",
									this->get_name().c_str(), (*dep)->get_name().c_str());
	
					// no error checking! this is performed in ExplicitEquation
									
					// apply the chain rule
					for (uint kk = 0; kk < nonzeros2; kk++)
						tmpcoeffs2[kk] = tmpcoeffs1[jj] * tmpcoeffs2[kk];

					// now tmpidcs2, tmpcoeffs2 store the stuff from var. jj of equation dep.
					// add this to the result arrays
					add_coeff(localnonzeros, localindices, localcoeffs, nonzeros2, tmpidcs2, tmpcoeffs2);
				}
			}
		}
	}

	// sort!
	uint   tmp_idx;
	double tmp_coeff;
	for (uint ii=0; ii < localnonzeros; ii++)
	{
		uint place = ii;
		while (place!=0)
		{
			if (localindices[place] < localindices[place-1]) // swap places
			{
				tmp_idx   = localindices[place];
				tmp_coeff = localcoeffs[place];
				localindices[place]   = localindices[place-1];
				localcoeffs[place]    = localcoeffs[place-1];
				localindices[place-1] = tmp_idx;
				localcoeffs[place-1]  = tmp_coeff;
				place--;
			}
			else {
				NEGF_ASSERT(localindices[place] != localindices[place-1],
							"The same index was encountered twice in localindices!");
				break; // we're done
			}
		}
	}
);}


/*  Make sure time stamps are correct; update to the given timestamp if necessary
 *  The Equation knows that all dependencies have at least the timestamp stored in "min_dep_timestamp"
 *  The method is not const because min_dep_timestamp is also updated.
 *  EDIT: min_dep_timestamp is currently out of use because "negative" timestamp updates cause trouble!
 *  @param tstamp the dependencies will be updated s.th. the have a timestamp >= tstamp.  */
void Equation::check_dependency_timestamps(uint tstamp) throw(Exception *)
{STACK_TRACE(
	if (tstamp > /*this->min_dep_timestamp*/ 0) {
		//logmsg->emit(LOG_INFO_L3,  
		//		"      eqn \"%s\" requested timestamp update of dependencies to %d (so far: %d)",
		//		this->get_name().c_str(), tstamp, this->min_dep_timestamp);
		for (uint nu = 0; nu < this->dependencies.size(); nu++) 
		{
			Equation * depcy = dependencies[nu];
			NEGF_ASSERT(depcy!=NULL, "a dependency was not allocated.");
			if (depcy->get_timestamp() < tstamp)
			{
				if (depcy->is_solved_in_newton()==false) {
					logmsg->emit(LOG_INFO_L3,  "      eqn \"%s\" updates dependency \"%s\" to timestamp %d (so far: %d)",
						this->get_name().c_str(), depcy->get_name().c_str(), tstamp, depcy->get_timestamp());
					depcy->compute_values(tstamp);
				} else {
					logmsg->emit(LOG_INFO,"Timestamp error when calling eqn \"%s\" from equation",
												depcy->get_name().c_str(), this->get_name().c_str() );
					logmsg->emit(LOG_INFO,"Requested time %d, other eqn has time %d.",
												tstamp, depcy->get_timestamp());
					NEGF_EXCEPTION("Trying to get value but wrong timestamp and qty is solved in newton.");
				}
			}
		}
		this->min_dep_timestamp = tstamp;
	}
);}


/** Call this right at the end of the preparation. <BR>
 *  The purpose is to have knowledge about the array sizes needed in <BR>
 *       Equation::get_newton_derivatives,                           <BR>
 *       Equation::newton_derivative and                             <BR>
 *       ExplicitEquation::all_direct_derivatives.                   <BR>
 *
 *  There the arrays are used in
 *  - (1)  this->newton_derivative(line, newton_var, localnonzeros, localindices, coeff)
 *  - (2)  this->get_newton_direct_derivatives(line, nvar, localnonzeros, localindices, localcoeffs)
 *  - (3)  this->get_newton_direct_derivatives(line, *dep, nonzeros1, tmpidcs1, tmpcoeffs1)
 *  - (4)  this->direct_derivatives(line, newton_var, nonzeros2, indices2, coeffs2)
 *  - (5)  this->direct_derivatives(line, *dep, nonzeros2, indices2, coeffs2)
 *  - (6)  (*dep)->all_direct_derivatives(tmpidcs1[jj], nvar, nonzeros2, tmpidcs2, tmpcoeffs2);
 *  For (6) the maximum of the other equation's nonzeros is taken. */
void Equation::set_nonzeros()
{//STACK_TRACE(
	logmsg->emit(LOG_INFO_L2, "Storing nonzeros per line for equation \"%s\"...",this->get_name().c_str());
	this->nonzeros_per_line.clear();
	this->nonzeros_per_line.resize(this->get_num_variables(), constants::eqn_array_size);
	
	#pragma omp parallel 
	{
	uint indices[constants::eqn_array_size]; // private to each thread
	#pragma omp for
	for (int line = 0; line < (int) this->get_num_variables(); line++)
	{
		try {
		uint nz_temp = 0;
		if (this->is_solved_in_newton())
		{
			this->get_sparsity(line, nz_temp, indices);	// will throw exception when called too early
			this->nonzeros_per_line[line] = nz_temp;
		} else {
			NEGF_ASSERT(this->is_newton_var_list_ready(), "set up newton_var_list first.");
			uint tmp = 0;
			for (uint ii = 0; ii < this->newton_var_list.size(); ii++) {
				this->all_direct_derivatives(line, newton_var_list[ii], nz_temp, indices, NULL);
				tmp += nz_temp+1;
			}
			for (uint ii = 0; ii < this->dependencies.size(); ii++) {
				if (!this->dependencies[ii]->is_solved_in_newton()) {
					this->all_direct_derivatives(line, dependencies[ii], nz_temp, indices, NULL);	// should really be direct_derivatives
					tmp = max(tmp+1, nz_temp);
				}
			}
			this->nonzeros_per_line[line] = tmp;
		}
		
		} catch(Exception * e) { 
			e->append_info(__LINE__,__FILE__,__DATE__,__TIME__, __func__); 		
			#ifdef _OPENMP 
				cout << e->get_reason() << "\n   equation:" <<this->get_name()  << "\nterminating from thread " << omp_get_thread_num() <<"\n";
				exit(1);
			#else
				NEGF_EXCEPTION(e->get_reason().c_str());
			#endif
		}
	} // pragma for
	} // pragma parallel
	
	
	// determine maximum and minimum number of nonzeros a line has
	this->max_nonzeros_per_line = 0;
	uint min_nonzeros = 100000000;
	for (uint ii=0; ii < this->get_num_variables(); ii++)
	{
		this->max_nonzeros_per_line = max(this->max_nonzeros_per_line, this->nonzeros_per_line[ii]);
		min_nonzeros = min(min_nonzeros, this->nonzeros_per_line[ii]);
	}
	
	logmsg->emit(LOG_INFO_L2, "Equation \"%s\": Nonzeros per line range between %d and %d.", 
		this->get_name().c_str(), min_nonzeros, this->max_nonzeros_per_line);
	this->nonzeros_set = true;
//);}
}

uint Equation::get_nonzeros(uint line) const
{STACK_TRACE(
	NEGF_ASSERT(this->nonzeros_was_set(), "nonzeros array has not yet been set up.");
	return this->nonzeros_per_line[line];
);}

uint Equation::get_max_nonzeros() const
{STACK_TRACE(
	if (this->nonzeros_was_set())
		return this->max_nonzeros_per_line;
	else
		return constants::eqn_array_size;
);}

