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
#ifndef EXPLICITEQN_H_NEGF
#define EXPLICITEQN_H_NEGF

#include "all.h"

#include "Equation.h"    

using namespace std;

namespace negf {


/** Base class for all equations given in the form \f$a_i=f(b, c, ...)\f$ where \f$a_i\f$ is the equation's variable i
 *
 *  Note also that no other own variables \f$a_j\f$ are necessary for \f$a_i\f$.
 *  Equations of this type need not be included in the Newton solver, but the user has the possibility to do so.
 *
 *  Equations using this base class have to implement the following things:
 *  - compute_value:           the actual function a(b,c)
 *  - direct_derivatives:      df/db_j, df/dc_j etc.		ATTENTION: ALWAYS MAKE DISTINCION COEFFS==NULL!
 * 
 *  Also the following class variables need to be initialized:
 *  - number_of_variables
 *  - its_type
 *  - dependencies
 *  - current_variable_values (if desired)
 *  - timestamp (if desired)
 *
 *  The name is usually given by a separate routine set_name() since every instance of the class might be named differently.
*/
class ExplicitEquation : public Equation
{
	friend class ExplicitEquationTester;
	
	public:
		ExplicitEquation() { solved_in_newton=false; }
		virtual ~ExplicitEquation() {}

		//! User has the chance to decide if this equation is solved in the Newton iteration
		void	set_solved_in_newton(bool solved_or_not) { solved_in_newton = solved_or_not; }

		//! Finding the Newton function from th explicit relationship is trivial.
		double	get_newton_function(uint line) const { return get_value(line) - compute_value(line); }

	protected:

		// will be implemented by derived classes
		virtual double	compute_value      (uint line) const = 0;
		virtual void	direct_derivatives (uint line, const Equation * newton_var, uint & nonzeros, 
											uint indices[], double coeff[]) const = 0;

		// for this kind of equation, the newton dir. der. are derived from direct_derivative!
		// a = a(b,c)  -->  F(a,b,c) = a - a(b,c) = 0
		void	get_newton_direct_derivatives  (uint line, const Equation * eqn, uint & nonzeros, 
												uint indices[], double coeff[]) const;

		// implemented from Equation class
		void 	get_dependence_sparsity(const uint & line, const Equation * eqn, 
											 uint & nonzeros, uint icols[]) const;
		void	all_direct_derivatives (const uint & line, const Equation * newton_var, uint & nonzeros, 
										uint indices[], double coeff[]) const;

};

} // end of namespace

#endif /*EXPLICITEQN_H_NEGF*/
