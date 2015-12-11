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
#include "ImplicitEquation.h"
using namespace negf;


/** Get the sparsity pattern w.r.t a direct dependence. <BR>
 *  The equationn must be in the dependency list. <BR>
 *  The user must implement direct_derivatives such that it does not compute coeffs[] when a NULL pointer is handed over! */
void ImplicitEquation::get_dependence_sparsity (const uint & line, const Equation * eqn, 
							   			 	    uint & nonzeros, uint icols[]) const
{STACK_TRACE(
	this->get_newton_direct_derivatives(line, eqn, nonzeros, icols, NULL);
);}
