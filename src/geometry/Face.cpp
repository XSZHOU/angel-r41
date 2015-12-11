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
#include "Face.h"
using namespace negf;

Face::Face(uint index_global_, usint nvertex_, usint nedge_):
	nvertex(nvertex_),
	nedge(nedge_),
	index_global(index_global_),
	index_external(0),
	area(0.0)
{STACK_TRACE(
	// determine face_type
	this->type = face_type::FT_Other;
	switch (nvertex)
	{
	case 0:
		NEGF_EXCEPTION("You attempted to create a face with no vertices!");
	case 1:
		NEGF_FASSERT(nedge==0, "You attempted to create a face with one vertex and %d edges.",nedge);
		type = face_type::FT_Point;
		break;
	case 2:
		NEGF_FASSERT(nedge==1, "You attempted to create a face with 2 vertices and %d edges.",nedge);
		type = face_type::FT_Interval;
		break;
	case 3:
		NEGF_FASSERT(nedge==3, "You attempted to create a face with 3 vertices and %d edges.",nedge);
		type = face_type::FT_Triangle;
		break;
	case 4:
		NEGF_FASSERT(nedge==4, "You attempted to create a face with 4 vertices and %d edges.",nedge);
		type = face_type::FT_Rectangle;
		break;
	default:
		NEGF_FEXCEPTION("Currently only faces with 1, 2, 3 or 4 vertices are supported: you have nvertex=%d",nvertex);
		break;
	}
	
	vertices.clear();
	edges.clear();
	
	ready = false;
);}


void Face::add_vertex(Vertex * vert) 
{STACK_TRACE(
	NEGF_ASSERT(vertices.size() < this->nvertex, "vertices.size() < nvertex failed");
	NEGF_ASSERT(vert!=0, "null pointer encountered.");
	this->vertices.push_back(vert);
);}


void Face::add_edge(Edge * the_edge) 
{STACK_TRACE(
	NEGF_ASSERT(edges.size() < this->nedge, "edges.size() < nedge failed");
	NEGF_ASSERT(the_edge!=0, "null pointer encountered.");
	this->edges.push_back(the_edge);
);}


void Face::get_vertex_indices(long index_array[]) const
{STACK_TRACE(
	NEGF_ASSERT(vertices.size()==this->nvertex, "not all vertices were assigned yet.");
	for (usint lid = 0; lid < this->nvertex; lid++) {
		NEGF_ASSERT(vertices[lid] != NULL, "vertices[lid] = NULL!");
		index_array[lid] = vertices[lid]->get_index_global();
	}
);}


void Face::get_edge_indices(long index_array[]) const
{STACK_TRACE(
	NEGF_ASSERT(edges.size()==this->nedge, "not all edges were assigned yet.");
	for (usint lid = 0; lid < this->nedge; lid++) {
		NEGF_ASSERT(vertices[lid] != NULL, "vertices[lid] = NULL!");
		index_array[lid] = edges[lid]->get_index_global();
	}
);}


Vertex* Face::get_vertex(usint lid) const 
{STACK_TRACE(
	NEGF_ASSERT(vertices.size()==this->nvertex, "not all vertices were assigned yet.");
	NEGF_ASSERT(lid < this->nvertex, "lid < nvertex failed");
	return this->vertices[lid];	
);}


Edge* Face::get_edge(usint lid) const 
{STACK_TRACE(
	NEGF_ASSERT(edges.size()==this->nedge, "not all edges were assigned yet.");
	NEGF_ASSERT(lid < this->nedge, "lid < nedge failed");
	return this->edges[lid];
);}


void Face::prepare()
{STACK_TRACE(
	// -----------------------------
	// compute the face area
	// -----------------------------
	double a[3], b[3];	// connection vectors v0->v1, v0->v2
	double tmp1, tmp2, tmp3, rhomboeder; // helpers
	switch(this->type)
	{
	case face_type::FT_Point:
		this->area = 1.0;
		break;
	case face_type::FT_Interval:
		this->area = this->get_vertex(0)->get_distance_to(this->get_vertex(1));
		break;
	case face_type::FT_Triangle:
	case face_type::FT_Rectangle:
		NEGF_ASSERT(this->get_vertex(0)->get_dimension()==3, "how can the vertex not be 3D when the face is 2D?");

		a[0] = this->get_vertex(1)->get_coordinate(0) - this->get_vertex(0)->get_coordinate(0);
		a[1] = this->get_vertex(1)->get_coordinate(1) - this->get_vertex(0)->get_coordinate(1);
		a[2] = this->get_vertex(1)->get_coordinate(2) - this->get_vertex(0)->get_coordinate(2);
		b[0] = this->get_vertex(2)->get_coordinate(0) - this->get_vertex(0)->get_coordinate(0);
		b[1] = this->get_vertex(2)->get_coordinate(1) - this->get_vertex(0)->get_coordinate(1);
		b[2] = this->get_vertex(2)->get_coordinate(2) - this->get_vertex(0)->get_coordinate(2);
		
		tmp1 = a[1]*b[2]-a[2]*b[1];
		tmp2 = a[2]*b[0]-a[0]*b[2];
		tmp3 = a[0]*b[1]-a[1]*b[0];
		rhomboeder = std::sqrt(tmp1*tmp1+tmp2*tmp2+tmp3*tmp3);
		
		if (type==face_type::FT_Triangle) {
			this->area = 0.5 * rhomboeder;
		} else { // FT_Rectangle
			this->area = rhomboeder;
		}
		break;
	case face_type::FT_Other:
	default:
		NEGF_EXCEPTION("Don't know how to compute the face are for this type.");
	}
	
	this->ready = true;
);}


bool Face::verify() const
{STACK_TRACE(
	// ------------------------------------------------
	// check if all vertices and edges have been set up
	// ------------------------------------------------
	NEGF_ASSERT(vertices.size()==this->nvertex, "not all vertices were assigned yet.");
	for(usint ii = 0; ii < nvertex; ii++) 
	{
		if(vertices[ii] == NULL) {
			logmsg->emit(LOG_ERROR,  "vertex %d in element is null", ii);			
			return false;	
		}	
		for(unsigned int jj = ii + 1; jj < nvertex; jj++) {
			if(vertices[ii] == vertices[jj]) {
				logmsg->emit(LOG_ERROR,  "vertices %d and %d in element are equal", ii, jj);
				return false;	
			}
		}	
	}
	for(usint ii = 0; ii < nedge; ii++) 
	{
		if(edges[ii] == NULL) {
			logmsg->emit(LOG_ERROR,  "edge %d in element is null", ii);	
			return false;	
		}	
		for(unsigned int jj = ii + 1; jj < nvertex; jj++) {
			if(edges[ii] == edges[jj]) {
				logmsg->emit(LOG_ERROR,  "edges %d and %d in element are equal", ii, jj);
				return false;	
			}
		}	
	}

	// ----------------------------------------------------------------------------
	// check if the vertices of each edge are found in the vertex list of this face
	// ----------------------------------------------------------------------------
	NEGF_ASSERT(edges.size()==this->nedge, "not all edges were assigned yet.");
	bool ok1, ok2;
	for (usint ii = 0; ii < nedge; ii++)
	{
		ok1 = false;
		ok2 = false;
		Vertex * vert1 = edges[ii]->get_lower_vertex();
		Vertex * vert2 = edges[ii]->get_upper_vertex();
		for (usint jj = 0; jj < nvertex; jj++)
		{
			if (vert1 == vertices[jj]) ok1 = true;
			if (vert2 == vertices[jj]) ok2 = true;
			if (ok1 && ok2) break;
		}
		if (!ok1) {
			logmsg->emit(LOG_ERROR,  "lower vertex of edge %d (global index %d) in face %d was not found in vertex list.", 
					ii, edges[ii]->get_index_global(), index_global);
			return false;
		}
		if (!ok2) {
			logmsg->emit(LOG_ERROR,  "upper vertex of edge %d (global index %d) in face %d was not found in vertex list.", 
					ii, edges[ii]->get_index_global(), index_global);
			return false;
		}
	}

	return true;
);}

double Face::get_area() const
{STACK_TRACE(
	if (this->ready) {
		return this->area;
	} else {
		NEGF_EXCEPTION("Face area was not yet prepared.")
		return -1.0;
	}
);}

