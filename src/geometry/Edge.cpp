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
#include "Edge.h"
using namespace negf;

Edge::Edge(uint index_global_, Vertex * vertex1_, Vertex * vertex2_):
	vertex1(vertex1_),
	vertex2(vertex2_),
	index_global(index_global_),
	index_external(0),
	voronoi_face(0.0)
{STACK_TRACE(
	if (vertex1==0 || vertex2==0)
		NEGF_EXCEPTION("Tried to construct an edge with null vertex pointers!");
	if (vertex1->get_dimension() != vertex2->get_dimension())
		NEGF_EXCEPTION("Edge contains vertices with different dimensionality.");
	switch (vertex1->get_dimension())
	{
	case 1:
		length = fabs(vertex1->get_coordinate(0) - vertex2->get_coordinate(0));
		break;
	case 2:
		length= sqrt( (vertex1->get_coordinate(0) - vertex2->get_coordinate(0))
					* (vertex1->get_coordinate(0) - vertex2->get_coordinate(0))
					+ (vertex1->get_coordinate(1) - vertex2->get_coordinate(1))
					* (vertex1->get_coordinate(1) - vertex2->get_coordinate(1)) );
		break;
	case 3:
		length= sqrt( (vertex1->get_coordinate(0) - vertex2->get_coordinate(0))
					* (vertex1->get_coordinate(0) - vertex2->get_coordinate(0))
					+ (vertex1->get_coordinate(1) - vertex2->get_coordinate(1))
					* (vertex1->get_coordinate(1) - vertex2->get_coordinate(1))
					+ (vertex1->get_coordinate(2) - vertex2->get_coordinate(2))
					* (vertex1->get_coordinate(2) - vertex2->get_coordinate(2)) );
		break;
	}
);}


uint Edge::get_lower_vertex_index() const
{STACK_TRACE(
	NEGF_ASSERT(this->vertex1 != 0, "lower vertex is NULL!");
	return this->vertex1->get_index_global();
);}


uint Edge::get_upper_vertex_index() const
{STACK_TRACE(
	NEGF_ASSERT(this->vertex2 != 0, "upper vertex is NULL!");
	return this->vertex2->get_index_global();
);}


Vertex * Edge::get_lower_vertex() const
{STACK_TRACE(
	NEGF_ASSERT(this->vertex1 != 0, "lower vertex is NULL!");
	return this->vertex1;
);}


Vertex * Edge::get_upper_vertex() const
{STACK_TRACE(
	NEGF_ASSERT(this->vertex2 != 0, "upper vertex is NULL!");
	return this->vertex2;
);}


Vertex * Edge::get_other_vertex(Vertex * v) const
{STACK_TRACE(
	if (vertex1==v) {
		NEGF_ASSERT(this->vertex2 != 0, "upper vertex is NULL!");
		return this->vertex2;
	} else {
		NEGF_ASSERT(this->vertex2 == v, "vertex not found in vertex list of edge!");
		return this->vertex1;
	}
);}


bool Edge::verify() const
{STACK_TRACE(
	// check if all vertices have been set up
	if (vertex1 == NULL) {
		logmsg->emit(LOG_ERROR,  "lower vertex of edge %d is null", this->get_index_global());			
		return false;	
	}	
	if (vertex2 == NULL) {
		logmsg->emit(LOG_ERROR,  "upper vertex of edge %d is null", this->get_index_global());			
		return false;	
	}
	return true;
);}
