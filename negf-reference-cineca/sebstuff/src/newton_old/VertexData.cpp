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
#include "VertexData.h"
using namespace negf;


VertexData::VertexData(const Geometry * const grid_, quantities::PhysicalQuantity type_, 
					vector<double> & values):
	grid(grid_)
{STACK_TRACE(
	this->its_type = type_;
	NEGF_ASSERT(grid_!=0, "geometry was null pointer.");
	NEGF_ASSERT(values.size() == grid_->get_num_vertices(), 
					"size of the vector of values is not the number of vertices.");
	dependencies.clear();				// no dependencies
	this->number_of_variables = grid_->get_num_vertices();
	this->timestamp = 99999;			// makes sure update is never needed

	current_variable_values = values;
);}

VertexData::VertexData(const Geometry * const grid_, quantities::PhysicalQuantity type_):
	grid(grid_)
{STACK_TRACE(
	this->its_type = type_;
	NEGF_ASSERT(grid_!=0, "geometry was null pointer.");
	dependencies.clear();				// no dependencies
	this->number_of_variables = grid->get_num_vertices();
	this->timestamp = 99999;			// makes sure update is never needed

	current_variable_values.clear();
	current_variable_values.resize(this->get_num_variables(), 0.0);
);}



void VertexData::direct_derivatives(uint line, const Equation * eqn, uint & nonzeros, 
								 uint indices [], double coeff []) const
{STACK_TRACE(
	NEGF_ASSERT(line < grid->get_num_vertices(), "Invalid line.");
	if (eqn!=this)
		NEGF_EXCEPTION("Epsilon is not dependent on this equation.");
	nonzeros   = 1;
	indices[0] = line;
	if (coeff!=NULL)
		coeff[0]   = 1;
);}

