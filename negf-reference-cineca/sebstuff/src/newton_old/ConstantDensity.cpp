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
#include "ConstantDensity.h"
using namespace negf;


/** @param type_ quantities::electron_density or quantities::hole_density
 *  @param grid_ the grid on which the density i s defined
 *  @param density_ the (constant) value. */
ConstantDensity::ConstantDensity(quantities::PhysicalQuantity type_, const Geometry * grid_, double density_):
	grid(grid_),
	density(density_)
{STACK_TRACE(
	NEGF_ASSERT(type_==quantities::electron_density || type_==quantities::hole_density, "Type must be a density");
	NEGF_ASSERT(grid_!=0 && density_ >= 0, "grid was null pointer or density <0.");
	this->its_type = type_;
	this->dependencies.clear();			// no dependencies
	this->number_of_variables = grid_->get_num_vertices();
	this->current_variable_values.clear();
	this->current_variable_values.reserve(grid_->get_num_vertices());
	for (uint ii = 0; ii < grid_->get_num_vertices(); ii++)
		this->current_variable_values.push_back(density_);
	this->timestamp = 99999;			// makes sure update is never needed
);}


void ConstantDensity::direct_derivatives(uint line, const Equation * eqn, 
							uint & nonzeros, uint icols [], double coeff []) const
{STACK_TRACE(
	NEGF_EXCEPTION("This should not be called.");
);}
