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
#ifndef IMPLICITEQN_H_NEGF
#define IMPLICITEQN_H_NEGF

#include "all.h"

#include "Equation.h"

using namespace std;


namespace negf {
     

/** Base class for all equations given in the form \f$F(a, b, c, ...)=0\f$ where a is the Equation's variable
 *
 *  This is an ImplicitEquation for the a_i, and an explicit form might not exist.
 *  Equations of this type MUST be included in the Newton solver.
 *
 *  Equations using this base class have to implement the following things:
 *  - get_newton_function:           the actual function F(a,b,c)
 *  - get_newton_direct_derivatives: dF/da_j, dF_db_j, df/dc_j etc.	ATTENTION: ALWAYS MAKE DISTINCION COEFFS==NULL!
 * 
 *  Also the following class variables need to be initialized:
 *  - number_of_variables
 *  - its_type
 *  - dependencies
 *  - current_variable_values (if desired)
 *  - timestamp (if desired)
 *  The name is usually given by a separate routine set_name().
 */
class ImplicitEquation : public Equation

{
	public:
		ImplicitEquation() { solved_in_newton = true; }
		virtual ~ImplicitEquation() {}

		//! will be implemented by derived classes
		virtual double	get_newton_function(uint line) const = 0;

	protected:

		//! will be implemented by derived classes
		virtual void	get_newton_direct_derivatives(uint line, const Equation * eqn, uint & nonzeros, 
													  uint indices[], double coeff[]) const = 0;

		// methods which are not possible to do with an implicit equation
		double			compute_value(uint line) const
			{ NEGF_EXCEPTION("compute_value cannot be done for an implicit equation!"); return 0; }
		void 			all_direct_derivatives(const uint & line, const Equation * newton_var, uint & nonzeros, 
												uint indices[], double coeff[]) const 
			{ NEGF_EXCEPTION("direct_derivatives cannot be computed for an implicit equation!"); }
			
		// implemented from equation class
		void	get_dependence_sparsity (const uint & line, const Equation * eqn, 
												 uint & nonzeros, uint icol_parts[]) const;

};

} // end of namespace

#endif /*IMPLICITEQN_H_NEGF*/
