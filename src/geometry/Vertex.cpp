/*
Copyright (c) 2009 Sebastian Steiger, Ratko Veprek, Integrated Systems Laboratory, ETH Zurich.
Comments, suggestions, criticism or bug reports are welcome: steiger@purdue.edu. 

This file is part of ANGEL, a simulator for LEDs based on the NEGF formalism.
The software is distributed under the Lesser GNU General Public License (LGPL).
ANGEL is free software: you can redistribute it and/or modify it under the terms 
of the Lesser GNU General Public License v3 or later. ANGEL is distributed
without any warranty; without even the implied warranty of merchantability or 
fitness for a particular purpose. See also <http://www.gnu.org/licenses/>.
*/
#include "Vertex.h"

using namespace negf;


/** constructor */ 
Vertex::Vertex(uint index_global_, double x_, double y_, double z_):
	dimension(3),
	index_internal(0),
	index_external(0),
	index_global(index_global_),
	contact(0),
	contact_vertex_index(0),
	heterointerface(false),
	at_quantized_region(false),
	minority_drain(-1)
{
	coord.push_back(x_);
	coord.push_back(y_);
	coord.push_back(z_);
}


Vertex::Vertex(uint index_global_, double x_, double y_):
	dimension(2),
	index_internal(0),
	index_global(index_global_),
	contact(0),
	contact_vertex_index(0),
	heterointerface(false),
	at_quantized_region(false),
	minority_drain(-1)
{
	coord.push_back(x_);
	coord.push_back(y_);
}


Vertex::Vertex(uint index_global_, double x_):
	dimension(1),
	index_internal(0),
	index_global(index_global_),
	contact(0),
	contact_vertex_index(0),
	heterointerface(false),
	at_quantized_region(false),
	minority_drain(-1)
{
	coord.push_back(x_);
}


void Vertex::set_contact_properties(Contact * contact_, uint idx)
{
	this->contact = contact_;
	this->contact_vertex_index = idx;	// index of the vertex within the contact list of vertices
}

/** NEVER change a coordinate during the simulation!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
void Vertex::set_coordinate(int xyz, double value)
{
	NEGF_ASSERT(xyz >=0 && xyz < this->dimension, "Invalid dimensionality.");
	this->coord[xyz] = value;
}

double Vertex::get_coordinate(unsigned short int ii) const
{STACK_TRACE(
	NEGF_ASSERT(ii < this->dimension, "Wrong dimensionality.");
	return coord[ii];
);}

double Vertex::get_distance_to(Vertex * other_vertex) const
{STACK_TRACE(
	NEGF_ASSERT(this->get_dimension()==other_vertex->get_dimension(), "incompatible dimensions.");
	double length2 = 0.0;
	for (uint ii = 0; ii < this->get_dimension(); ii++) {
		double tmp = this->get_coordinate(ii) - other_vertex->get_coordinate(ii);
		length2 += tmp * tmp;
	}
	return std::sqrt(length2);
);}


