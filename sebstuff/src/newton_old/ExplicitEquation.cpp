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
#include "ExplicitEquation.h"
using namespace negf;


/** This function is used only when for some reason the ExplicitEquation is included in the Newton iteration,
 *  even though its variables could be substituted. <BR>
 *  The Newton dir. der. are derived from direct_derivative! <BR>
 *  \f$ a = a(b,c)  \Rightarrow  F(a,b,c) = a - a(b,c) = 0 \f$  */
void ExplicitEquation::get_newton_direct_derivatives(uint line, const Equation * eqn, 
							uint & nonzeros, uint indices[], double coeff[]) const
{STACK_TRACE(
	if (eqn==this) {
		nonzeros   = 1; 
		indices[0] = line;    
		coeff[0]   = 1.0;
		return;
	}

	if ( find(this->dependencies.begin(), this->dependencies.end(), eqn) == this->dependencies.end() )
		NEGF_EXCEPTION("Error, attempted to get direct derivative even though equation doesnt directly depend on it.");

	this->direct_derivatives(line, eqn, nonzeros, indices, coeff);
	for (uint ii = 0; ii < nonzeros; ii++)
		coeff[ii] = coeff[ii] * -1;

	return;
);}


/** This function collects the TOTAL derivative w.r.t. a newton variable from all the dependencies. <BR>
 *  When called with coeffs=NULL, the function only returns the sparse indices */
void ExplicitEquation::all_direct_derivatives(const uint & line, const Equation * newton_var, uint & nonzeros, 
												uint indices[], double coeffs[]) const
{STACK_TRACE(
	NEGF_ASSERT(newton_var!=NULL, "null pointer encountered.");
	nonzeros = 0;
	
	uint   nonzeros2;
	// used array size will be less or equal than the total number of nonzeros on this line of the jacobian
	// this gives huge performance improvement!
	uint arraysize = (this->nonzeros_was_set()) ? this->get_nonzeros(line) : constants::eqn_array_size;
	arraysize = arraysize + 1;	// arrays must have minimum size 1!
	uint   indices2[arraysize];
	double coeffs2[arraysize];		
	
	for (vector<Equation *>::const_iterator dep = dependencies.begin(); dep != dependencies.end(); ++dep)
	{
		NEGF_ASSERT((*dep)!=NULL, "null pointer encountered.");
		if (*dep==newton_var)  // no recursion
		{
			nonzeros2 = 0;
			if (coeffs!=NULL) {
				this->direct_derivatives(line, newton_var, nonzeros2, indices2, coeffs2);
			} else { 
				this->direct_derivatives(line, newton_var, nonzeros2, indices2, NULL);
			}
			NEGF_FASSERT(nonzeros2 < arraysize, "%s Overflow of arrays(#1, nonzeros2=%d < arraysize=%d), newton_var=%s",
					this->get_name().c_str(), nonzeros2, arraysize, newton_var->get_name().c_str());
			
			// error checking
			if (constants::eqn_debug && coeffs!=NULL) {
				for (uint ii = 0; ii < nonzeros2; ii++)
				{
					if (std::isnan(coeffs2[ii])) {
						logmsg->emit( LOG_ERROR,"Not-A-Number Jacobian entry detected derivative of \"%s\" w.r.t. \"%s\"",
									this->get_name().c_str(), (*dep)->get_name().c_str() );
						NEGF_EXCEPTION("Aborting.");
					}
					if (std::isinf(coeffs2[ii])) {
						logmsg->emit( LOG_ERROR,"Infinite Jacobian entry detected in derivative of \"%s\" w.r.t. \"%s\"",
									this->get_name().c_str(), (*dep)->get_name().c_str() );
						NEGF_EXCEPTION("Aborting.");
					}
				}
			}
			
			if (coeffs!=NULL) {
				this->add_coeff(nonzeros, indices, coeffs, nonzeros2, indices2, coeffs2);
			} else {
				this->add_idcs(nonzeros, indices, nonzeros2, indices2);
			}		
		}
		
		// do a recursion if 
		// 1. the dependency is not solved in the newton (i.e. its an explicit equation)
		// 2. the newton_var appears in the dependency's list of newton variables
		//if ( find(this->newton_var_list.begin(), this->newton_var_list.end(), *dep) == this->newton_var_list.end() )
		if (!(*dep)->is_solved_in_newton())	// faster
		{
			NEGF_ASSERT( (*dep)->is_newton_var_list_ready(),
						"An equation was encountered where newton_var_list was not yet set up.");
			const vector<Equation * > & dep_newton_var_list = (*dep)->get_newton_var_list();
			if ( find(dep_newton_var_list.begin(), dep_newton_var_list.end(), newton_var) != dep_newton_var_list.end())
			{
				NEGF_ASSERT( !(*dep)->is_solved_in_newton(),
						"A dependency was marked as solved by newton even though it was not in the newton var list.");

				//logmsg->emit(LOG_INFO,"%s: depcy %s depends on %s.",this->get_name().c_str(), (*dep)->get_name().c_str(), newton_var->get_name().c_str());

				nonzeros2 = 0;
				if (coeffs!=NULL) {
					this->direct_derivatives(line, *dep, nonzeros2, indices2, coeffs2);
				} else { 
					this->direct_derivatives(line, *dep, nonzeros2, indices2, NULL);
				}
				NEGF_FASSERT(nonzeros2 < arraysize, "%s: Overflow of arrays (#2, nonzeros2=%d < arraysize=%d). dep=%s",
					this->get_name().c_str(), nonzeros2, arraysize, (*dep)->get_name().c_str());
				
				// error checking
				if (constants::eqn_debug && coeffs!=NULL) {
					for (uint ii = 0; ii < nonzeros2; ii++)
					{
						if (std::isnan(coeffs2[ii])) {
							logmsg->emit( LOG_ERROR,"Not-A-Number Jacobian entry detected derivative of \"%s\" w.r.t. \"%s\"",
										this->get_name().c_str(), (*dep)->get_name().c_str() );
							NEGF_EXCEPTION("Aborting.");
						}
						if (std::isinf(coeffs2[ii])) {
							logmsg->emit( LOG_ERROR,"Infinite Jacobian entry detected in derivative of \"%s\" w.r.t. \"%s\"",
										this->get_name().c_str(), (*dep)->get_name().c_str() );
							NEGF_EXCEPTION("Aborting.");
						}
					}
				}
				
				uint   nonzeros3 = 0;
				uint   indices3[(*dep)->get_max_nonzeros()+1];
				double coeffs3[(*dep)->get_max_nonzeros()+1];

				for (uint ii = 0; ii < nonzeros2; ii++)
				{
					NEGF_ASSERT(ii < arraysize, "index overflow!");
										
					nonzeros3 = 0;
				
					// get the full derivative of variable indices2[ii] of the dependency w.r.t. newton_var
					if (coeffs!=NULL) {
						(*dep)->all_direct_derivatives(indices2[ii], newton_var, nonzeros3, indices3, coeffs3);
						NEGF_FASSERT(nonzeros3 < (*dep)->get_max_nonzeros()+1, 
							"(*dep)=%s: Overflow of arrays (nonzeros3=%d < (*dep)->get_max_nonzeros()+1=%d), newton_var=%s",
							(*dep)->get_name().c_str(), nonzeros3, (*dep)->get_max_nonzeros()+1, newton_var->get_name().c_str());

						// chain rule
						for (uint jj = 0; jj < nonzeros3; jj++)
							coeffs3[jj] = coeffs2[ii] * coeffs3[jj];
							
						this->add_coeff(nonzeros, indices, coeffs, nonzeros3, indices3, coeffs3);
					} else {
						(*dep)->all_direct_derivatives(indices2[ii], newton_var, nonzeros3, indices3, NULL);
						this->add_idcs(nonzeros, indices, nonzeros3, indices3);
					}
				}
			}
		}
	}
	
	//if (this->get_name()=="equasifermi_well0" && newton_var->get_name()=="potential_3d")
	//	cout << "equasifermi_well0::all_direct_derivatives("<<line<<","<<newton_var->get_name()<<"): nonzeros=" << nonzeros << endl;
		
	// sorting does not have to be done because these entries are added to the higher-level list from
	// Equation::newton_derivative
	
	//delete [] indices2;
	//delete [] coeffs2;
);}


/** Get the sparsity pattern w.r.t a direct dependence. <BR>
 *  The user must implement direct_derivatives in a way that it does not compute coeffs[] when a NULL pointer is handed over!
 *  @param line the line under investogation
 *  @param eqn the Equation with respect to which the sparsity is found. Must be in the dependency list!
 *  @param nonzeros (result) number of nonzeros on that line w.r.t. that dependency
 *  @param icols    (result) stores the indices */
void ExplicitEquation::get_dependence_sparsity(const uint & line, const Equation * eqn, 
										uint & nonzeros, uint icols[]) const
{STACK_TRACE(
	if (eqn==this) {
		nonzeros = 1;
		icols[0] = line;
		return;
	}
	this->direct_derivatives(line, eqn, nonzeros, icols, NULL);
);}


