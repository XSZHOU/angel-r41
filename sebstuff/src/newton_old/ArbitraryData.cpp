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
#include "ArbitraryData.h"
using namespace negf;


/** @param type_ the physical type of the data
 *  @param values the values. */
ArbitraryData::ArbitraryData(quantities::PhysicalQuantity type_, const vector<double> & values)
{STACK_TRACE(
	this->its_type = type_;
	dependencies.clear();				// no dependencies
	this->number_of_variables = values.size();
	this->timestamp = 99999;			// makes sure update is never needed

	current_variable_values = values;
);}

