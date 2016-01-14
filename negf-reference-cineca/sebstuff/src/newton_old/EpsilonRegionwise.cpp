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
#include "EpsilonRegionwise.h"
using namespace negf;

EpsilonRegionwise::EpsilonRegionwise(const Geometry * grid_,
									 const MaterialDatabase * db_):
	grid(grid_),
	db(db_)
{STACK_TRACE(
	NEGF_ASSERT(grid_!=0 && db_!=0, "tries to assign null pointer.");
	
	this->its_type = quantities::dielectric_coeff;
	this->dependencies.clear();			// no dependencies
	this->number_of_variables = grid->get_num_elements();
	this->current_variable_values.clear();
	this->current_variable_values.resize(grid->get_num_elements(), 0.0);
	
	// -----------------------------
	// compute values on elements
	// -----------------------------
	//string fieldname = "electrostatic_permittivity";
	string fieldname = "static_dielectric_constant";
	for (uint ii = 0; ii < grid->get_num_elements(); ii++)
	{
		Region * region = grid->get_element(ii)->get_region();
		double value = constants::convert_from_SI(units::dielectric,
						region->get_material()->get(fieldname.c_str())
						* constants::SIeps0 );

		current_variable_values[ii] = value;
		
		map<Region *,double>::iterator iter = region_values.find(region);
		if (iter == region_values.end() ) {
			region_values[region] = value;
			logmsg->emit(LOG_INFO_L2,  "Region \"%10s\" has eps/eps0=%7.3g[C/Vm].",
							region->get_name().c_str(), 
							value / constants::convert_from_SI(units::dielectric,1.0)  / constants::SIeps0 );
		}
	}
	
	this->timestamp = 99999;			// makes sure update is never needed
);}


double EpsilonRegionwise::get_epsilon(Region * region) const
{STACK_TRACE(
	map<Region *,double>::const_iterator iter = region_values.find(region);
	if (iter != region_values.end() ) {
		return iter->second;
	} else {
		NEGF_FEXCEPTION("Region %s not found.", region->get_name().c_str());
		return -1.0;
	}
);}
