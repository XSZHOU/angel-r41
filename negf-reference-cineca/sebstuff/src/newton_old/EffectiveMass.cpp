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
#include "EffectiveMass.h"
using namespace negf;

EffectiveMass::EffectiveMass(const Geometry * grid_, const BoxMethod * boxmethod_,
				 const MaterialDatabase * db_, 
				 quantities::PhysicalQuantity electrons_or_holes_, string verts_or_elems_):
	grid(grid_),
	boxmethod(boxmethod_),
	db(db_),
	electrons_or_holes(electrons_or_holes_),
	verts_or_elems(verts_or_elems_)
{STACK_TRACE(
	NEGF_ASSERT(grid!=NULL && boxmethod!=NULL && db!=NULL, "encountered null pointer.");
	NEGF_ASSERT(electrons_or_holes==quantities::electron_density 
				|| electrons_or_holes==quantities::hole_density,
				"must specify electrons or holes, not something else.");
	NEGF_ASSERT(verts_or_elems=="vertex" || verts_or_elems=="element",
				"must specify vertex or element, not something else.");
	
	this->its_type = quantities::mass;
		
	// ----------------------------------------------
	// read in values defined on elements from file
	// ----------------------------------------------
	this->elem_values.assign(grid->get_num_elements(), 0.0);
	for (uint ii = 0; ii < grid->get_num_elements(); ii++)
	{
		double effmass;
		if (electrons_or_holes==quantities::electron_density) {
			effmass = grid->get_element(ii)->get_region()->get_material()->get("electron_effective_mass");
		} else if (electrons_or_holes==quantities::hole_density) {
			effmass = grid->get_element(ii)->get_region()->get_material()->get("hole_effective_mass");
		} else {
			NEGF_EXCEPTION("electrons or holes...");
		}
		elem_values[ii] = constants::convert_from_SI(units::mass, constants::SIm0 * effmass);
	}
	
	// -----------------------------------------------
	// set current_variable_values
	// -----------------------------------------------
	if (verts_or_elems=="element")
	{ 
		this->number_of_variables = this->grid->get_num_elements();
		this->current_variable_values.resize(this->number_of_variables, 0.0);
		for (uint ii = 0; ii < grid->get_num_elements(); ii++)
		{
			this->current_variable_values[ii] = elem_values[ii];
		}
	} else if (verts_or_elems=="vertex") 
	{
		this->number_of_variables = this->grid->get_num_vertices();
		this->current_variable_values.resize(this->number_of_variables, 0.0);
		
		const double * const * measure  = boxmethod->get_measure();
		for (uint ii = 0; ii < grid->get_num_vertices(); ii++)
		{
			const vector<Element *> & elems_near_vert = grid->get_elems_near(grid->get_vertex(ii));
			double total_voronoi = 0.0;
			for (uint jj = 0; jj < elems_near_vert.size(); jj++) 
			{
				uint elem_idx = elems_near_vert[jj]->get_index_global();
				uint local_vert_idx = elems_near_vert[jj]->get_local_index(grid->get_vertex(ii));
				
				this->current_variable_values[ii] += measure[elem_idx][local_vert_idx] * elem_values[elem_idx];
				
				total_voronoi += measure[elem_idx][local_vert_idx];
			}
			NEGF_ASSERT(total_voronoi > 0.0, "something went wrong.");
			this->current_variable_values[ii] = current_variable_values[ii] / total_voronoi;
		}
	} else {
		NEGF_EXCEPTION("location must be element or vertex.");
	}
	
	this->timestamp = 999999;	// makes sure eqn is never computed
);}


double EffectiveMass::get_element_value(uint elem_idx) const
{STACK_TRACE(
	NEGF_ASSERT(elem_values.size()>elem_idx, "invalid element index.");
	return elem_values[elem_idx];
);}

