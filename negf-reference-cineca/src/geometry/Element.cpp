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
#include "Element.h"
using namespace negf;

Element::Element(uint index_global_, element_type::ElementType type_) 
{STACK_TRACE(
	this->type = type_;
	switch (type)
	{
	case element_type::interval:
		num_verts = 2; num_edges = 1; num_faces = 0;
		break;
	case element_type::triangle:
		num_verts = 3; num_edges = 3; num_faces = 0;
		break;
	case element_type::rectangle:
		num_verts = 4; num_edges = 4; num_faces = 0;
		break;
	case element_type::tetrahedron:
		num_verts = 4; num_edges = 6; num_faces = 4;
		break;
	case element_type::pyramid:
		num_verts = 5; num_edges = 8; num_faces = 5;
		break;
	case element_type::prism:
		num_verts = 6; num_edges = 9; num_faces = 5;
		break;
	case element_type::cuboid:
		num_verts = 8; num_edges = 12; num_faces = 6; // is this right???
		break;
	case element_type::tetrabrick:
		num_verts = 7; num_edges = 12; num_faces = 7;
		break;
	case element_type::point:
	default:
		NEGF_EXCEPTION("Element type not supported.");
	}
	// vertices, edges and faces arrays will be left empty for now
	this->vertices.clear();
	this->edges.clear();
	this->faces.clear();
	this->region         = 0;
	this->element_volume = -1.0;
	this->diameter		 = -1.0;
	this->inv_jacobian   = NULL;
	this->index_global   = index_global_;
	this->index_external = 0;
	this->ready		     = false;
);}


Element::~Element() 
{STACK_TRACE(
	this->region = 0;
	if (this->inv_jacobian!=NULL)
		delete [] inv_jacobian;
);}


Vertex* Element::get_vertex(uint lid) const 
{STACK_TRACE(
	NEGF_ASSERT(lid < this->vertices.size(), "lid < vertices.size()");
	NEGF_ASSERT(vertices[lid] != 0, "vertices[lid] contains empty pointer!");
	return this->vertices[lid];	
);}


Edge* Element::get_edge(uint lid) const 
{STACK_TRACE(
	NEGF_ASSERT(lid < edges.size(), "lid < edges.size()");
	NEGF_ASSERT(edges[lid] != 0, "edges[lid] contains empty pointer!");
	return this->edges[lid];
);}


Face* Element::get_face(uint lid) const 
{STACK_TRACE(
	NEGF_ASSERT(lid < faces.size(), "lid < faces.size()");
	NEGF_ASSERT(faces.size()>lid, "faces was not yet set up");
	NEGF_ASSERT(faces[lid] != 0, "faces[lid] contains empty pointer!");
	return this->faces[lid];
);}


void Element::add_vertex(Vertex* vert) 
{STACK_TRACE(
	NEGF_ASSERT(vert!=0, "tried to assign null pointer");
	NEGF_ASSERT(vertices.size() < num_verts, "element already has enough vertices");
	this->vertices.push_back(vert);
);}


void Element::add_edge(Edge* new_edge) 
{STACK_TRACE(
	NEGF_ASSERT(new_edge!=0, "tried to assign null pointer");
	NEGF_ASSERT(edges.size() < num_edges, "element already has enough edges");
	this->edges.push_back(new_edge);
);}


void Element::add_face(Face* new_face) 
{STACK_TRACE(
	NEGF_ASSERT(new_face!=0, "tried to assign null pointer");
	NEGF_ASSERT(faces.size() < num_faces, "element already has enough faces");
	this->faces.push_back(new_face);
);}


void Element::set_region(Region* region_) 
{STACK_TRACE(
	this->region = region_;
);}


void Element::get_vertex_indices(uint index_array[]) const
{STACK_TRACE(
	for(uint ii = 0; ii < this->get_num_vertices(); ii++) {
		index_array[ii] = vertices[ii]->get_index_global();
	}
);}


void Element::get_edge_indices(uint index_array[]) const
{STACK_TRACE(
	for(uint ii = 0; ii < this->get_num_edges(); ii++) {
		index_array[ii] = edges[ii]->get_index_global();
	}
);}


void Element::get_face_indices(uint index_array[]) const
{STACK_TRACE(
	for(uint ii = 0; ii < this->get_num_faces(); ii++) {
		index_array[ii] = faces[ii]->get_index_global();
	}
);}


Region* Element::get_region() const
{STACK_TRACE(
	NEGF_ASSERT(this->region != 0, "Region was not set yet!");
	return region;
);}


double Element::get_volume() const
{STACK_TRACE(
	NEGF_ASSERT( this->ready && element_volume != -1.0, "Element has not been finalized properly!");
	return element_volume;
);}


double Element::get_diameter() const
{STACK_TRACE(
	NEGF_ASSERT( this->ready && diameter != -1.0, "Element has not been finalized properly!");
	return diameter;
);}


const double * Element::get_invjacobian() const
{STACK_TRACE(
	NEGF_ASSERT( this->ready, "Element has not been finalized properly!");
	return inv_jacobian;
);}


/** Check whether element is properly initialized:  <BR>
  * Check whether all vertices are set, no vertex two times set,
  * if element is assigned to a region, element volume was calculated and if Jacobian determinant is set
  */
bool Element::verify() const
{STACK_TRACE(
	// -------------------------------------------------------
	// check if all vertices, edges and faces have been set up
	// -------------------------------------------------------
	for(uint ii = 0; ii < this->num_verts; ii++) 
	{
		if(vertices[ii] == 0) {
			logmsg->emit(LOG_ERROR,  "vertex %d in element is null", ii);
			return false;	
		}	
		for(uint jj = ii+1; jj < this->num_verts; jj++) {
			if(vertices[ii] == vertices[jj]) {
				logmsg->emit(LOG_ERROR,  "vertices %d and %d in element are equal", ii, jj);
				return false;	
			}
		}	
	}
	for(uint ii = 0; ii < this->num_edges; ii++) 
	{
		if(edges[ii] == 0) {
			logmsg->emit(LOG_ERROR,  "edge %d in element is null", ii);
			return false;	
		}	
		for(uint jj = ii+1; jj < this->num_edges; jj++) {
			if(edges[ii] == edges[jj]) {
				logmsg->emit(LOG_ERROR,  "edges %d and %d in element are equal", ii, jj);
				return false;	
			}
		}	
	}
	for(uint ii = 0; ii < this->num_faces; ii++) 
	{
		if(faces[ii] == 0) {
			logmsg->emit(LOG_ERROR,  "face %d in element is null", ii);
			return false;	
		}	
		for(uint jj = ii+1; jj < this->num_faces; jj++) {
			if(faces[ii] == faces[jj]) {
				logmsg->emit(LOG_ERROR,  "faces %d and %d in element are equal", ii, jj);
				return false;	
			}
		}	
	}

	// ---------------------------------------------------------
	// check if element type, region and volume have been set up
	// ---------------------------------------------------------
	if(type == element_type::noshape){
		logmsg->emit(LOG_ERROR, "no type set to element");
		return false;	
	}
	if(region == 0) {
		logmsg->emit(LOG_ERROR, "no region set to element");
		return false;	
	}
	if(element_volume <= 0.0) {
		logmsg->emit(LOG_ERROR,  "volume of element is %5.5g and therefore incorrect", element_volume);
		return false;
	}

	// -------------------------------------------------------------------------------
	// check if the vertices of each edge are found in the vertex list of this element
	// -------------------------------------------------------------------------------
	bool ok1, ok2;
	for (uint ii = 0; ii < this->num_edges; ii++)
	{
		ok1 = false;
		ok2 = false;
		Vertex * vert1 = edges[ii]->get_lower_vertex();
		Vertex * vert2 = edges[ii]->get_upper_vertex();
		for (uint jj = 0; jj < this->num_verts; jj++)
		{
			if (vert1 == vertices[jj]) ok1 = true;
			if (vert2 == vertices[jj]) ok2 = true;
			if (ok1 && ok2) break;
		}
		if (!ok1) {
			logmsg->emit(LOG_ERROR,  "lower vertex of edge %d (global index %d) in element %d was not found in vertex list."
					,ii, edges[ii]->get_index_global(), index_global);
			return false;
		}
		if (!ok2) {
			logmsg->emit(LOG_ERROR,  "upper vertex of edge %d (global index %d) in element %d was not found in vertex list."
					, ii, edges[ii]->get_index_global(), index_global);
			return false;
		}
	}

	// --------------------------------------------------------------------------
	// check if the edges of each face are found in the edge list of this element
	// --------------------------------------------------------------------------
	bool ok;
	Edge * edg;
	for (uint ii = 0; ii < this->num_faces; ii++)
	{
		for (uint jj = 0; jj < (faces[ii])->get_num_edges(); jj++)
		{
			ok = false;
			edg = (faces[ii])->get_edge(jj);
			for (usint kk = 0; kk < this->get_num_edges(); kk++)
			{ 
				if (edg == edges[kk]) ok = true; 
			}
			if (!ok) {
				logmsg->emit(LOG_ERROR,  "edge %d of face %d was not found in the element's edge list.", jj, ii);
				return false;
			}
		}
	}
	
	// ------------------------------------------------
	// check for correct diameter
	// ------------------------------------------------
	if (this->get_type()==element_type::tetrahedron || this->get_type()==element_type::triangle 
		||this->get_type()==element_type::interval)
	{
		double max_edgelength = 0.0;
		for (uint ii=0; ii < this->get_num_edges(); ii++) {
			double edgelength = this->get_edge(ii)->get_length();
			if (edgelength > max_edgelength)
					max_edgelength = edgelength;
			}
			NEGF_ASSERT(max_edgelength > 0.0 && fabs(max_edgelength - this->get_diameter()) < 1e-14,
					"Diameter inconsistency.");
	}
	
	return true;
);}


void Element::prepare() 
{STACK_TRACE(
	NEGF_ASSERT(ready==false, "element has already beed prepared.");
	NEGF_ASSERT(num_verts==vertices.size() && num_edges==edges.size() && num_faces==faces.size(),
				"number of vx/ed/fac not consistent.");
	// check if vertices are all set
	for (uint ii = 0; ii < this->num_verts; ii++) {
		NEGF_ASSERT(vertices[ii]!=NULL, "Vertices are not completely set yet.");
	}
	this->compute_volume();
	this->compute_invjacobian();
	this->compute_diameter();
	this->ready = true;
);}


void Element::compute_volume() 
{STACK_TRACE(
	
	// check if vertices have all the same dimensionality
	NEGF_ASSERT(vertices.size()==this->num_verts, "wrong vertices array.");
	NEGF_ASSERT(vertices.size()>0, "no vertices!?!");
	NEGF_ASSERT(vertices[0]!=NULL, "vertex 0 was the null pointer.");
	usint dim = vertices[0]->get_dimension();
	for (uint ii = 1; ii < this->num_verts; ii++)
	{
		NEGF_ASSERT(vertices[ii]!=NULL, "a vertex is missing in an element.");
		if (vertices[ii]->get_dimension() != dim)
			NEGF_EXCEPTION("Not all vertices have the same dimension!");
	}

	// ---------------------------------------------------
	// calculate element volume
	// ---------------------------------------------------
	Edge * edge1 = 0;
	Edge * edge2 = 0;
	double vector1[2];
	double vector2[2];
	vector<Edge *> prism_heights;
	vector<Edge *> prism_other_edges;
	Vertex * vertex = 0;
	double prism_groundarea = 0.0;
	switch (dim)
	{
	case 1:
		this->element_volume = abs(vertices[0]->get_coordinate(0) - vertices[1]->get_coordinate(0));
		break;
	case 2:
		switch (this->type)
		{
		case element_type::triangle:
			NEGF_ASSERT(this->get_num_vertices()==3 && this->num_verts==3, "Is this really a triangle?");
			vector1[0] = vertices[1]->get_coordinate(0) - vertices[0]->get_coordinate(0);
			vector1[1] = vertices[1]->get_coordinate(1) - vertices[0]->get_coordinate(1);
			vector2[0] = vertices[2]->get_coordinate(0) - vertices[0]->get_coordinate(0);
			vector2[1] = vertices[2]->get_coordinate(1) - vertices[0]->get_coordinate(1);
			this->element_volume = 0.5 * fabs( vector1[0] * vector2[1] - vector1[1]*vector2[0] );
			break;
		case element_type::rectangle:
			// determine two edges with a common vertex
			// the rectangle are is then simply length1*length2
			NEGF_ASSERT(this->get_num_edges()==4, "Is this really a rectangle?");
			edge1 = edges[0];
			edge2 = 0;
			for (uint ii = 1; ii < 4; ii++)
			{
				if (   edges[ii]->get_lower_vertex()==edge1->get_lower_vertex()
					|| edges[ii]->get_upper_vertex()==edge1->get_lower_vertex()
					|| edges[ii]->get_lower_vertex()==edge1->get_upper_vertex()
					|| edges[ii]->get_upper_vertex()==edge1->get_upper_vertex() ) {
					edge2 = edges[ii];
					continue;
				}
			}
			NEGF_ASSERT(edge2!=0, "Could not find another edge with a common vertex.");
			this->element_volume = edge1->get_length() * edge2->get_length();
			break;
		default:
			NEGF_EXCEPTION("Can't compute element volume for this type.");
			break;
		}
		break;
	case 3:
		uint num_edges_connected_to_vertex0;
		switch (this->type)
		{
		case element_type::tetrahedron:
			NEGF_ASSERT(this->get_num_vertices()==4, "Is this really a tetrahedron?");
			// create jacobi matrix (jacobi is: Pi - P0, i = 1..3)
			double  m[9];
			for(uint ii = 0; ii < 3; ii++) {
				for(uint jj = 0; jj < 3; jj++)
					m[jj * 3 + ii] =  vertices[ii + 1]->get_coordinate(jj) - vertices[0]->get_coordinate(jj);
			}
			// calculate determinant
			double jacobi_det;
			jacobi_det =  m[0] * m[4] * m[8] 
						+ m[1] * m[5] * m[6]
						+ m[2] * m[3] * m[7]
						- m[6] * m[4] * m[2]
						- m[7] * m[5] * m[0]
						- m[8] * m[3] * m[1];
					
			// calculate element volume
			this->element_volume = 1.0 / 6.0 * (jacobi_det > 0.0 ? jacobi_det : -jacobi_det);
			break;
		case element_type::cuboid:
			NEGF_ASSERT(this->get_num_vertices()==8, "Is this really a cuboid?");
			// note: we assume that the edges are perpendicular to eachother!!! (we don't check this)
			num_edges_connected_to_vertex0 = 0;
			this->element_volume = 1.0;
			for (uint ii = 0; ii < this->get_num_edges(); ii++)
			{
				if (   this->get_edge(ii)->get_lower_vertex() == this->get_vertex(0)
					|| this->get_edge(ii)->get_upper_vertex() == this->get_vertex(0) )
				{
					this->element_volume = element_volume * this->get_edge(ii)->get_length();
					num_edges_connected_to_vertex0++;
				}
			}
			if (num_edges_connected_to_vertex0 != 3)
			{
				cout << num_edges_connected_to_vertex0 << endl;
				NEGF_EXCEPTION("Vertex 0 seemed to have not exactly 3 connected edges.");
				ready = false;
				break;
			}
			break;
		case element_type::prism:   // ground area: triangle
			NEGF_ASSERT(this->get_num_vertices()==6, "Is this really a prism?");
			
			// find the three edges which have the same length and point in the same direction
			prism_heights.clear();
			for (uint ii = 0; ii < this->get_num_edges(); ii++)
			{
				//cout <<"edge "<<ii<<endl;
				Edge * edge = this->get_edge(ii);
				double vec1[3];
				vec1[0] = edge->get_upper_vertex()->get_coordinate(0) - edge->get_lower_vertex()->get_coordinate(0);
				vec1[1] = edge->get_upper_vertex()->get_coordinate(1) - edge->get_lower_vertex()->get_coordinate(1);
				vec1[2] = edge->get_upper_vertex()->get_coordinate(2) - edge->get_lower_vertex()->get_coordinate(2);
				
				// find edges parallel to edge
				uint num_parallel_edges = 0;
				for (uint jj = 0; jj < this->get_num_edges(); jj++)
				{
					if (jj != ii)
					{
						// scalar product
						double vec2[3];
						vec2[0] = get_edge(jj)->get_upper_vertex()->get_coordinate(0) - get_edge(jj)->get_lower_vertex()->get_coordinate(0);
						vec2[1] = get_edge(jj)->get_upper_vertex()->get_coordinate(1) - get_edge(jj)->get_lower_vertex()->get_coordinate(1);
						vec2[2] = get_edge(jj)->get_upper_vertex()->get_coordinate(2) - get_edge(jj)->get_lower_vertex()->get_coordinate(2);
						double prod = negf_math::vector_scalar_product_3d(vec1, vec2);
						//cout << "   prod:" << prod <<", edge_length*edge_length:"<<edge->get_length()*edge->get_length()<<endl;
						if ( (fabs(fabs(prod) - edge->get_length()*edge->get_length())) / (edge->get_length()*edge->get_length()) < 1e-10 ) {
							num_parallel_edges++;
						}
					}
				}
				//cout << "num_parallel_edges=" << num_parallel_edges << endl;
				if (num_parallel_edges==2) {
					prism_heights.push_back(edge);
				}
			}
			NEGF_FASSERT(prism_heights.size()==3, "%d instead of 3 \"vertical\" edges were found.", prism_heights.size());
			
			// from the edge prism_heights[0] take lower vertex an find the other two edge which are connected to it
			// from those edges find the ground area and multiply it by prism_heights[0]->get_length()
			vertex = prism_heights[0]->get_lower_vertex();
			prism_other_edges.clear();
			for (uint ii = 0; ii < this->get_num_edges(); ii++)
			{
				if (this->get_edge(ii)!=prism_heights[0] && 
					(this->get_edge(ii)->get_lower_vertex()==vertex || this->get_edge(ii)->get_upper_vertex()==vertex)) {
					prism_other_edges.push_back(this->get_edge(ii));
				}
			}
			NEGF_ASSERT(prism_other_edges.size()==2, "not exactly 2 ground area edges were found.");
			double vec1[3];
			vec1[0] = prism_other_edges[0]->get_upper_vertex()->get_coordinate(0) - prism_other_edges[0]->get_lower_vertex()->get_coordinate(0);
			vec1[1] = prism_other_edges[0]->get_upper_vertex()->get_coordinate(1) - prism_other_edges[0]->get_lower_vertex()->get_coordinate(1);
			vec1[2] = prism_other_edges[0]->get_upper_vertex()->get_coordinate(2) - prism_other_edges[0]->get_lower_vertex()->get_coordinate(2);
			double vec2[3];
			vec2[0] = prism_other_edges[1]->get_upper_vertex()->get_coordinate(0) - prism_other_edges[1]->get_lower_vertex()->get_coordinate(0);
			vec2[1] = prism_other_edges[1]->get_upper_vertex()->get_coordinate(1) - prism_other_edges[1]->get_lower_vertex()->get_coordinate(1);
			vec2[2] = prism_other_edges[1]->get_upper_vertex()->get_coordinate(2) - prism_other_edges[1]->get_lower_vertex()->get_coordinate(2);
			double vec3[3];
			negf_math::vector_outer_product(vec1,vec2,vec3);
			
			prism_groundarea = negf_math::vector_norm_3d(vec3);
			this->element_volume = prism_groundarea * prism_heights[0]->get_length();
			break;
		default:
			NEGF_FEXCEPTION(
				"volume computation for element %d (dfise number %d, type %d) is not yet implemented.",
				this->get_index_global(), this->get_index_external(), this->get_type());
			break;
		}
		break;
	default:
		NEGF_EXCEPTION("Strange dimensionality.");
	}
	return;
);}

/** compute inverse of Jacobian of affine transformation reference element-->real element
 *  in the case of triangle/tetrahedron, the reference element is a 2D/3D simplex
 *  in the case of rectangle/cuboid, the reference element is a unit square/cube
 *  for other elements, it is not implemented. */
void Element::compute_invjacobian()
{STACK_TRACE(

	// check if vertices have all the same dimensionality
	usint dim = vertices[0]->get_dimension();
	for (uint ii = 1; ii < this->num_verts; ii++) {
		if (vertices[ii]->get_dimension() != dim)
			NEGF_EXCEPTION("Not all vertices have the same dimension!");
	}

	vector<Vertex *> adj_vertices;
	double m3[9]; double P03[3]; double Pi3[3]; double div;     // used in 3D case
	double m2[4]; double P02[2]; double Pi2[2];	double bc_m_ad; // used in 2D case
	double m1;													// used in 1D case
	double prod;	// helpers used in 2D/3D
	Vertex * swap;	
	double PxP0[3][3]; // helpers used in 3D
	double cross_prod[3];
	double rotP1P0[2]; // helpers used in 2D
	double P3P0[2];
	uint count, sideedge, num_parallel_edges; // helpers for prism
	switch (dim)
	{
	case 3:
		inv_jacobian = new double[9];
		for(uint ii = 0; ii < 3; ii++) {
			P03[ii] = this->get_vertex(0)->get_coordinate(ii);
		}
		switch(this->get_type())
		{
		case element_type::tetrahedron:			// (jacobi is: Pi - P0, i = 1..3)
			for(uint jj = 0; jj < 3; jj++) {
				Pi3[0] = this->get_vertex(jj+1)->get_coordinate(0);
				Pi3[1] = this->get_vertex(jj+1)->get_coordinate(1);
				Pi3[2] = this->get_vertex(jj+1)->get_coordinate(2);
				for(uint ii = 0; ii < 3; ii++) {
					m3[ii * 3 + jj] =  Pi3[ii] - P03[ii];
				}
			}
			break;
		case element_type::cuboid:					// (jacobi is: Pi - P0, i = 1..3)
			// find the three vertices P1, P2, P3 adjacent to P0
			for (uint ii=0; ii < this->get_num_edges(); ii++) {
				if (this->get_edge(ii)->get_lower_vertex()==this->get_vertex(0))
					adj_vertices.push_back(this->get_edge(ii)->get_upper_vertex());
				if (this->get_edge(ii)->get_upper_vertex()==this->get_vertex(0))
					adj_vertices.push_back(this->get_edge(ii)->get_lower_vertex());
			}
			NEGF_ASSERT( adj_vertices.size()==3, "Wrong number of vertices connected to vertex 0.");
			
			// check if adj_vertices[2]-P0 is on the same side as (adj_vertices[0]-P0 x adj_vertices[1]-P0)
			// if not, switch adj_vertices[0] and adj_vertices[1]
			//double PxP0[3][3];
			for (uint ii = 0; ii < 3; ii++)
				for (uint jj = 0; jj < 3; jj++)
					PxP0[ii][jj] = adj_vertices[ii]->get_coordinate(jj) - P03[jj];
			//double cross_prod[3];
			negf_math::vector_outer_product(PxP0[0],PxP0[1],cross_prod);
			prod = negf_math::vector_scalar_product_3d(cross_prod, PxP0[2]);
			if (prod < 0.0) {
				swap = adj_vertices[1];
				adj_vertices[1] = adj_vertices[0];
				adj_vertices[0] = swap;
			}
			
			// write Pi-P0 into m3
			for(uint jj = 0; jj < 3; jj++) {
				Pi3[0] = adj_vertices[jj]->get_coordinate(0);
				Pi3[1] = adj_vertices[jj]->get_coordinate(1);
				Pi3[2] = adj_vertices[jj]->get_coordinate(2);
				for(uint ii = 0; ii < 3; ii++) {
					m3[ii * 3 + jj] =  Pi3[ii] - P03[ii];
				}
			}
			break;
		case element_type::prism:	// (jacobi is: Pi - P0, i = 1..3); P3 is at the "top" of the prism (P0 at the bottom)
		
			// find the three edges which are connected to vertex 0
			count = 0;
			for (uint ii = 0; ii < this->get_num_edges(); ii++)
			{
				if (this->get_edge(ii)->get_lower_vertex()==this->get_vertex(0)) {
					PxP0[count][0] = this->get_edge(ii)->get_upper_vertex()->get_coordinate(0) - P03[0];
					PxP0[count][1] = this->get_edge(ii)->get_upper_vertex()->get_coordinate(1) - P03[1];
					PxP0[count][2] = this->get_edge(ii)->get_upper_vertex()->get_coordinate(2) - P03[2];
					count++;
				}
				if (this->get_edge(ii)->get_upper_vertex()==this->get_vertex(0)) {
					PxP0[count][0] = P03[0] - this->get_edge(ii)->get_lower_vertex()->get_coordinate(0);
					PxP0[count][1] = P03[1] - this->get_edge(ii)->get_lower_vertex()->get_coordinate(1);
					PxP0[count][2] = P03[2] - this->get_edge(ii)->get_lower_vertex()->get_coordinate(2);
					count++;					
				}
			}
			NEGF_ASSERT(count==3, "Not exactly 3 edges found.");
			
			// determine which edge goes to the top
			sideedge = 888888;
			for (uint ii = 0; ii < 3; ii++)
			{
				//cout << "edge "<<ii<<"..."<<endl;
				double PxP0_length = negf_math::vector_norm_3d(PxP0[ii]);
				
				// find edges parallel to edge
				num_parallel_edges = 0;
				for (uint jj = 0; jj < this->get_num_edges(); jj++)
				{
					// scalar product
					double vec2[3];
					vec2[0] = get_edge(jj)->get_upper_vertex()->get_coordinate(0) - get_edge(jj)->get_lower_vertex()->get_coordinate(0);
					vec2[1] = get_edge(jj)->get_upper_vertex()->get_coordinate(1) - get_edge(jj)->get_lower_vertex()->get_coordinate(1);
					vec2[2] = get_edge(jj)->get_upper_vertex()->get_coordinate(2) - get_edge(jj)->get_lower_vertex()->get_coordinate(2);
					double prod = negf_math::vector_scalar_product_3d(PxP0[ii], vec2);
					//cout << "   prod: " << prod << ", PxP0_length:" << PxP0_length << endl;
					if ( (fabs(fabs(prod) - PxP0_length*PxP0_length)) / (PxP0_length*PxP0_length) < 1e-8 ) {
						num_parallel_edges++;
					}
				}
				//cout << "num_parallel_edges: " << num_parallel_edges << endl;
				if (num_parallel_edges==3) {
					NEGF_ASSERT(sideedge==888888, "prism side edge was already found!");
					sideedge = ii;
				}
			}
			NEGF_ASSERT(sideedge!=888888, "prism side edge was not found.");
			
			// write Pi-P0 into m3
			count = 0;
			for (uint ii = 0; ii < 3; ii++)
			{	
				if (ii==sideedge) {
					for(uint jj = 0; jj < 3; jj++) {
						m3[jj * 3 + 2]     = PxP0[ii][jj];
					}
				} else {
					for(uint jj = 0; jj < 3; jj++) {
						m3[jj * 3 + count] = PxP0[ii][jj];
					}
					count++;
				}
			}
			NEGF_ASSERT(count==2, "something went wrong.");
			
			break;
		default:
			NEGF_EXCEPTION("Jacobian computation is not implemented for this element type.");
			break;
		}
		// inverse of matrix w/ entries { m[0..2], m[3..5], m[6..9] }
		div = - m3[2]*m3[4]*m3[6] + m3[1]*m3[5]*m3[6] + m3[2]*m3[3]*m3[7] 
			  - m3[0]*m3[5]*m3[7] - m3[1]*m3[3]*m3[8] + m3[0]*m3[4]*m3[8];
		inv_jacobian[0] = (- m3[5] * m3[7] + m3[4] * m3[8]) / div;
		inv_jacobian[1] = (  m3[2] * m3[7] - m3[1] * m3[8]) / div;
		inv_jacobian[2] = (- m3[2] * m3[4] + m3[1] * m3[5]) / div;
		inv_jacobian[3] = (  m3[5] * m3[6] - m3[3] * m3[8]) / div;
		inv_jacobian[4] = (- m3[2] * m3[6] + m3[0] * m3[8]) / div;
		inv_jacobian[5] = (  m3[2] * m3[3] - m3[0] * m3[5]) / div;
		inv_jacobian[6] = (- m3[4] * m3[6] + m3[3] * m3[7]) / div;
		inv_jacobian[7] = (  m3[1] * m3[6] - m3[0] * m3[7]) / div;
		inv_jacobian[8] = (- m3[1] * m3[3] + m3[0] * m3[4]) / div;
		break;
	case 2:
		inv_jacobian = new double[4];
		P02[0] = this->get_vertex(0)->get_coordinate(0);
		P02[1] = this->get_vertex(0)->get_coordinate(1);
		switch(this->get_type())
		{
		case element_type::triangle:
			for(uint jj = 0; jj < 2; jj++) {
				Pi2[0] = this->get_vertex(jj+1)->get_coordinate(0);
				Pi2[1] = this->get_vertex(jj+1)->get_coordinate(1);
				for(uint ii = 0; ii < 2; ii++) {
						m2[ii * 2 + jj] =  Pi2[ii] - P02[ii];
				}
			}
			break;
		case element_type::rectangle:
			// find the two vertices adjacent to P0
			for (uint ii=0; ii < this->get_num_edges(); ii++) {
				if (this->get_edge(ii)->get_lower_vertex()==this->get_vertex(0))
					adj_vertices.push_back(this->get_edge(ii)->get_upper_vertex());
				if (this->get_edge(ii)->get_upper_vertex()==this->get_vertex(0))
					adj_vertices.push_back(this->get_edge(ii)->get_lower_vertex());
			}
			NEGF_ASSERT( adj_vertices.size()==2, "Wrong number of vertices connected to vertex 0.");

			// check if adj_vertices[1] is on the left (counter-clockwise) side of adj_vertices[0]-P0
			// if not, switch adj_vertices[1] and adj_vertices[0]
			// if P1P0 = (x,y), then rotP1P0=(y,-x)
			rotP1P0[0] = adj_vertices[0]->get_coordinate(1) - P02[1];
			rotP1P0[1] = P02[0] - adj_vertices[0]->get_coordinate(0);
			P3P0[0] = adj_vertices[1]->get_coordinate(0) - P02[0];
			P3P0[1] = adj_vertices[1]->get_coordinate(1) - P02[1];
			prod = negf_math::vector_scalar_product_2d(P3P0, rotP1P0);
			if (prod < 0.0) {
				swap = adj_vertices[1];
				adj_vertices[1] = adj_vertices[0];
				adj_vertices[0] = swap;
			}
			
			for(uint jj = 0; jj < 2; jj++) {
				Pi2[0] = adj_vertices[jj]->get_coordinate(0);
				Pi2[1] = adj_vertices[jj]->get_coordinate(1);
				for(uint ii = 0; ii < 2; ii++) {
						m2[ii * 2 + jj] =  Pi2[ii] - P02[ii];
				}
			}
			break;
		default:
			NEGF_EXCEPTION("Jacobian computation is not implemented for this element type.");
			break;
		}
		bc_m_ad = m2[1]*m2[2] - m2[0]*m2[3];
		this->inv_jacobian[0] = - m2[3] / bc_m_ad;
		this->inv_jacobian[1] =   m2[1] / bc_m_ad;
		this->inv_jacobian[2] =   m2[2] / bc_m_ad;
		this->inv_jacobian[3] = - m2[0] / bc_m_ad;
		break;
	case 1:
		inv_jacobian = new double[1];
		switch(this->get_type())
		{
		case element_type::interval:
			NEGF_ASSERT(this->get_num_vertices()==2 && this->get_vertex(0)->get_dimension()==1,
							"something is fishy.");
			m1 = this->get_vertex(1)->get_coordinate(0) - this->get_vertex(0)->get_coordinate(0);
			NEGF_ASSERT(fabs(m1)>constants::min_vector_norm, "edge is too small.");
			this->inv_jacobian[0] = 1.0 / m1;
			break;
		default:
			NEGF_EXCEPTION("Jacobian computation is not implemented for this element type.");
			break;
		}
		break;
	default:
		NEGF_EXCEPTION("Dimensionality does not make sense for Jacobian computation.");
		break;
	}
);}

/** The diameter is defined as the biggest distance between two vertices */
void Element::compute_diameter()
{STACK_TRACE(
	double result = 0.0;
	uint dim = get_vertex(0)->get_dimension();
	for (uint ii = 0; ii < this->get_num_vertices(); ii++) {
		double vi[dim];
		for (uint idx = 0; idx < dim; idx++)
			vi[idx] = get_vertex(ii)->get_coordinate(idx);
			
		for (uint jj = ii+1; jj < this->get_num_vertices(); jj++) 
		{
			NEGF_ASSERT( get_vertex(jj)->get_dimension()==dim, "dimension mismatch.");
			double vj[dim];
			for (uint idx = 0; idx < dim; idx++)
				vj[idx] = get_vertex(jj)->get_coordinate(idx);
				
			double dist = 0.0;
			for (uint idx = 0; idx < dim; idx++)
				dist += (vi[idx]-vj[idx]) * (vi[idx]-vj[idx]);
			dist = sqrt(dist);
			
			if (dist > result)
				result = dist;
		}
	}
	this->diameter = result;
);}


/** Check if a vertex is inside this element.
	To do so, the (previously computed) Jacobian is used for the
	transformation to the reference element (simplex/unitcube/etc.)
	and the coordinates of the vertex are checked if they lie
	within the reference element.
*/
bool Element::contains(Vertex * vertex) const
{STACK_TRACE(
	NEGF_ASSERT(this->ready, "prepare() element first.");
	NEGF_ASSERT(this->inv_jacobian != NULL, "Compute element jacobian first!");
	const uint dim = this->get_vertex(0)->get_dimension();
	NEGF_ASSERT(vertex!=NULL && dim==vertex->get_dimension(), "dimension mismatch between element and vertex.");
	
	if (   this->get_type()==element_type::triangle || this->get_type()==element_type::tetrahedron  
		|| this->get_type()==element_type::interval)
	{
		// only entries up to dim are used, but init of variable length arrays is non-standard and supported in gcc, not in CC
		double tmp_coords[3];
		for(uint ii = 0; ii < dim; ii++) {
			tmp_coords[ii] = vertex->get_coordinate(ii) - this->get_vertex(0)->get_coordinate(ii);
			if (fabs(tmp_coords[ii]) > this->get_diameter())	// speeds up everything....
				return false;
		}
		for (uint ii = dim; ii < 3; ii++)
			tmp_coords[ii] = 0.0;			// not used
			
		//const uint nverts = this->get_num_vertices();
		double contrib[3+1]; 				//contrib[nverts]?
		for (uint ii=0; ii < 3+1; ii++) 	// dim+1 / nverts+1?
			contrib[ii] = 0.0;
		
		for(uint ii = 0; ii < dim; ii++) {
			for(uint jj = 0; jj < dim; jj++) {
				contrib[ii + 1] += inv_jacobian[ii * dim + jj] * tmp_coords[jj];
			}
			if(contrib[ii + 1] < (0.0 - constants::elem_contain_tol) || contrib[ii + 1] > (1.0 + constants::elem_contain_tol)) {
				return false;
			}
		}
		
		contrib[0] = 1.0;
		for(uint ii = 0; ii < dim; ii++) 	// nverts?
			contrib[0] -= contrib[ii+1];
			
		return (contrib[0] >= (0.0 - constants::elem_contain_tol) && contrib[0] <= (1.0 + constants::elem_contain_tol));
	}
	else if ( this->get_type()==element_type::rectangle || this->get_type()==element_type::cuboid )
	{
		// we assume that inv_jacobian is the inverse of the matrix (P1-P0,P3-P0) where 
		// P0 is vertex 0 of the element and P1,P3 are connected to P0 by edges, 
		// the numbering being counter-clockwise
		
		// compute P-P0
		// only entries up to dim are used, but init of variable length arrays is non-standard and supported in gcc, not in CC
		double tmp_coords[3];
		for(uint ii = 0; ii < dim; ii++) {
			tmp_coords[ii] = vertex->get_coordinate(ii) - this->get_vertex(0)->get_coordinate(ii);
			if (fabs(tmp_coords[ii]) > this->get_diameter())	// speeds up everything....
				return false;
		}
		for (uint ii = dim; ii < 3; ii++)
			tmp_coords[ii] = 0.0;			// not used

		// compute the relative coordinates: P = P0 + x*(P1-P0) + y*(P3-P0) (+ z*(P4-P0))
		double contrib[3+1];
		for (uint ii=0; ii < 3+1; ii++)
			contrib[ii] = 0.0;
		for(uint ii = 0; ii < dim; ii++) {
			for(uint jj = 0; jj < dim; jj++) {
				contrib[ii + 1] += inv_jacobian[ii * dim + jj] * tmp_coords[jj];
			}
			if(contrib[ii + 1] < (0.0 - constants::elem_contain_tol) || contrib[ii + 1] > (1.0 + constants::elem_contain_tol)) {
				return false;
			}
		}
		
		return true;
	}
	else {
		NEGF_EXCEPTION("Element type not implemented.");
		return false;
	}
);}


/** Given some vertex, find the coefficients of the linear combination of shape functions 
 *  of the element corners. The code is very duplicate to the routine Element::contains(...) */
void Element::get_linear_combination_coefficients(Vertex * vertex, vector<double> & contribs) const
{STACK_TRACE(

	NEGF_ASSERT(this->ready, "prepare() element first.");
	NEGF_ASSERT(vertex!=NULL, "faulty input vertex");
	const uint dim    = this->get_vertex(0)->get_dimension();
	NEGF_ASSERT(vertex!=NULL && dim==vertex->get_dimension(), "dimension mismatch between element and vertex.");
	
	if (   this->get_type()==element_type::triangle || this->get_type()==element_type::tetrahedron  
		|| this->get_type()==element_type::interval)
	{
		// only entries up to dim are used, but initialization of variable length arrays is non-standard and supported in gcc, not in CC
		double tmp_coords[3];
		for(uint ii = 0; ii < dim; ii++)
			tmp_coords[ii] = vertex->get_coordinate(ii) - this->get_vertex(0)->get_coordinate(ii);
		for (uint ii = dim; ii < 3; ii++)
			tmp_coords[ii] = 0.0;			// not used
	
		contribs.clear();
		contribs.resize(dim+1,0.0);	
		for(uint ii = 0; ii < dim; ii++) {
			for(uint jj = 0; jj < dim; jj++)
				contribs[ii + 1] += inv_jacobian[ii * dim + jj] * tmp_coords[jj];
		}
		contribs[0] = 1.0;
		for(uint ii = 0; ii < dim; ii++)
			contribs[0] -= contribs[ii+1];
		
		return;
	}
	else if ( this->get_type()==element_type::rectangle || this->get_type()==element_type::cuboid )
	{
		// we assume that inv_jacobian is the inverse of the matrix (P1-P0,P3-P0) where 
		// P0 is vertex 0 of the element and P1,P3 are connected to P0 by edges, 
		// the numbering being counter-clockwise
		
		// only entries up to dim are used, but initialization of variable length arrays is non-standard and supported in gcc, not in CC
		double tmp_coords[3];
		for(uint ii = 0; ii < dim; ii++)
			tmp_coords[ii] = vertex->get_coordinate(ii) - this->get_vertex(0)->get_coordinate(ii);
		for (uint ii = dim; ii < 3; ii++)
			tmp_coords[ii] = 0.0;			// not used
	
		// compute the relative coordinates: P = P0 + x*(P1-P0) + y*(P3-P0) (+ z*(P4-P0))
		double rel_coords[3];
		for (uint ii=0; ii < 3; ii++)
			rel_coords[ii] = 0.0;
		for(uint ii = 0; ii < dim; ii++) {
			for(uint jj = 0; jj < dim; jj++) {
				rel_coords[ii] += inv_jacobian[ii * dim + jj] * tmp_coords[jj];
			}
		}
		for (uint ii = dim; ii < 3; ii++)
			rel_coords[ii] = 0.0;			// not used
		
		contribs.clear();
		switch(dim) {
		case 2: // rectangle
			contribs.resize(4, 0.0);
			contribs[0] = (1.0 - rel_coords[0]) * (1.0 - rel_coords[1]);    // P0 has shapefunction (1-x)(1-y)
			contribs[1] = rel_coords[0] 		* (1.0 - rel_coords[1]);   	// P1 has shapefunction x(1-y)
			contribs[2] = rel_coords[0] 		* rel_coords[1];   			// P2 has shapefunction xy
			contribs[3] = (1.0 - rel_coords[0]) * rel_coords[1];   			// P3 has shapefunction (1-x)y
			break;
		case 3:
			contribs.resize(8, 0.0);
			contribs[0] = (1.0 - rel_coords[0]) * (1.0 - rel_coords[1]) * (1.0 - rel_coords[2]);    // P0 has shapefunction (1-x)(1-y)(1-z)
			contribs[1] = rel_coords[0] 		* (1.0 - rel_coords[1]) * (1.0 - rel_coords[2]);    // P1 has shapefunction x(1-y)(1-z)
			contribs[2] = rel_coords[0] 		* rel_coords[1] 		* (1.0 - rel_coords[2]);    // P2 has shapefunction xy(1-z)
			contribs[3] = (1.0 - rel_coords[0]) * rel_coords[1] 		* (1.0 - rel_coords[2]);    // P3 has shapefunction (1-x)y(1-z)
			contribs[4] = (1.0 - rel_coords[0]) * (1.0 - rel_coords[1]) * rel_coords[2];   			// P4 has shapefunction (1-x)(1-y)z
			contribs[5] = rel_coords[0] 		* (1.0 - rel_coords[1]) * rel_coords[2];  			// P5 has shapefunction x(1-y)z
			contribs[6] = rel_coords[0] 		* rel_coords[1] 		* rel_coords[2];    		// P6 has shapefunction xyz
			contribs[7] = (1.0 - rel_coords[0]) * rel_coords[1] 		* rel_coords[2];    		// P7 has shapefunction (1-x)yz
			break;
		default:
			NEGF_EXCEPTION("Wrong dimension.");
			break;
		}
		double check = 0.0;
		for (uint ii = 0; ii < contribs.size(); ii++)
			check += contribs[ii];
		NEGF_ASSERT(fabs(check - 1.0) < 1e-12, "something went wrong.");
		
		return;
	}
	else {
		NEGF_EXCEPTION("Element type not implemented.");
		return;
	}
	
);}


/** Get the index of a vertex within the element's list of vertices
 *  e.g. for a tetrahedron this will be 0, 1, 2 or 3.
 *  If the vertex cannot be found, an exception is thrown. */
uint Element::get_local_index(const Vertex * vertex) const
{STACK_TRACE(
	for (uint ii = 0; ii < this->get_num_vertices(); ii++) 
	{
		if (this->get_vertex(ii)==vertex)
			return ii;
	}
	NEGF_FEXCEPTION("Vertex %d was not found in element %d's list of vertices.",
				vertex->get_index_global(), this->get_index_global());
);}


/** Get the index of an edge within the element's list of edges
 *  e.g. for a triangle this will be 0, 1, 2 or 3.
 *  If the edge cannot be found, an exception is thrown. */
uint Element::get_local_index(const Edge * edge) const
{STACK_TRACE(
	for (uint ii = 0; ii < this->get_num_edges(); ii++) 
	{
		if (this->get_edge(ii)==edge)
			return ii;
	}
	for (uint ii=0; ii<this->get_num_vertices(); ii++) {
		logmsg->emit(LOG_INFO,"    vertex %d", this->get_vertex(ii)->get_index_global());
	}
	NEGF_FEXCEPTION("Edge %d (vertices %d, %d) was not found in the list of edges of element %d.",
				edge->get_index_global(), 
				edge->get_lower_vertex()->get_index_global(),
				edge->get_upper_vertex()->get_index_global(),				
				this->get_index_global());
);}

