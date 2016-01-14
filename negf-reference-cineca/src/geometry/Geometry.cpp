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
#include "Geometry.h"
using namespace negf;


/** Initialize geometry object with an expected number of vertex, edge etc. entities */
Geometry::Geometry(uint nvertices_, uint nedges_, uint nfaces_, uint nelements_) throw (Exception *):
	nvertices(nvertices_),
	nedges(nedges_),
	nfaces(nfaces_),
	nelems(nelements_),
	nregions(0),  // will only be determined during reading of grid file -->add_contact(.)
	ncontacts(0), // will only be determined during reading of grid file -->add_contact(.)
	max_connectivity(0),
	num_dfise_elems(0),
	num_dfise_regions(0),
	num_internal_vertices(0)
{STACK_TRACE(
	// ---------------------------------------------
	// reserve memory for vertex, edge etc. arrays
	// ---------------------------------------------
	vertices.reserve(nvertices);
	edges.   reserve(nedges);
	faces.   reserve(nfaces);
	elements.reserve(nelems);
	NEGF_ASSERT(vertices.capacity() >= nvertices, "not enough memory for vertex vector");
	NEGF_ASSERT(edges.capacity() >= nedges, "not enough memory for edge vector");
	NEGF_ASSERT(faces.capacity() >= nfaces, "not enough memory for face vector");
	NEGF_ASSERT(elements.capacity() >= nelems, "not enough memory for element vector");

	// -------------------------------------
	// reserve memory for X_near_Y arrays
	// -------------------------------------
	edges_near_vertex		.reserve(nvertices);	//needs even more space
	elems_near_face			.reserve(nfaces);
	elems_near_edge			.reserve(nedges);
	elems_near_vertex		.reserve(nvertices);
	faces_near_vertex		.reserve(nvertices);
	regions_near_vertex		.reserve(nvertices);
	NEGF_ASSERT(edges_near_vertex.capacity() 	   >= nvertices, "not enough memory for edges_near_vertex");
	NEGF_ASSERT(elems_near_face.capacity() 		   >= nfaces, 	 "not enough memory for elems_near_face");
	NEGF_ASSERT(elems_near_edge.capacity() 		   >= nedges,	 "not enough memory for elems_near_edge");
	NEGF_ASSERT(elems_near_vertex.capacity()	   >= nvertices, "not enough memory for elems_near_vertex");
	NEGF_ASSERT(faces_near_vertex.capacity() 	   >= nvertices, "not enough memory for faces_near_vertex");
	NEGF_ASSERT(regions_near_vertex.capacity() 	   >= nvertices, "not enough memory for regions_near_vertex");

	vector<Edge *> empty_edge_vector;
	for (uint ii = 0; ii < nvertices; ii++) {
		edges_near_vertex.push_back(empty_edge_vector);
		edges_near_vertex[ii].reserve(15);	// is this enough?
		NEGF_ASSERT((edges_near_vertex[ii]).capacity() >= 15, "not enough memory for edges_near_vertex (2)");
	}
	vector<Element *> empty_element_vector;
	for (uint ii = 0; ii < nfaces; ii++) {
		elems_near_face.push_back(empty_element_vector);
		elems_near_face[ii].reserve(5);	// is this enough?
		NEGF_ASSERT((elems_near_face[ii]).capacity() >= 5, "not enough memory for elems_near_face (2)");
	}
	for (uint ii = 0; ii < nedges; ii++) {
		elems_near_edge.push_back(empty_element_vector);
		elems_near_edge[ii].reserve(5);	// is this enough?
		NEGF_ASSERT((elems_near_edge[ii]).capacity() >= 5, "not enough memory for elems_near_edge (2)");
	}
	for (uint ii = 0; ii < nvertices; ii++) {
		elems_near_vertex.push_back(empty_element_vector);
		elems_near_vertex[ii].reserve(10);	// is this enough?
		NEGF_ASSERT((elems_near_vertex[ii]).capacity() >=10, "not enough memory for elems_near_vertex (2)");
	}
	vector<Face *> empty_face_vector;
	for (uint ii = 0; ii < nvertices; ii++) {
		faces_near_vertex.push_back(empty_face_vector);
		faces_near_vertex[ii].reserve(10);	// is this enough?
		NEGF_ASSERT((faces_near_vertex[ii]).capacity() >= 10, "not enough memory for faces_near_vertex (2)");
	}
	vector<Region *> empty_region_vector;
	for (uint ii = 0; ii < nvertices; ii++) {
		regions_near_vertex.push_back(empty_region_vector);
	}
	
	this->prepared = false;
);}

Geometry::~Geometry() 
{STACK_TRACE(	
	for(uint ii = 0; ii < elements.size(); ii++) {
		delete elements[ii]; elements[ii] = 0;
	}
	for(uint ii = 0; ii < faces.size(); ii++) {
		delete faces[ii]; faces[ii] = 0;	
	}	
	for(uint ii = 0; ii < edges.size(); ii++) {
		delete edges[ii]; edges[ii] = 0;	
	}
	for(uint ii = 0; ii < vertices.size(); ii++) {
		delete vertices[ii]; vertices[ii] = 0;
	}
	for(uint ii = 0; ii < regions.size(); ii++) {
		delete regions[ii];	regions[ii] = 0;	
	}
	for(uint ii = 0; ii < contacts.size(); ii++) {
		delete contacts[ii]; contacts[ii] = 0;	
	}
);}


/** Prepare geometry for calculation. Sets up the following quantities:
 * 1. the array marking if a vertex is "internal" or not
 * 2. the element volumes
 * 3. the maximum number of edges connected to a vertex;
 * 4. for each vertex a list of all edges connected to it
 * 5. for each edge a list of all elements adjacent to it
 * 6. for each contact determine the adjacent region
 */
void Geometry::prepare() throw (Exception *)
{STACK_TRACE(

	// check arrays
	NEGF_ASSERT(vertices.size()==nvertices, "vertices.size()!=nvertices");
	NEGF_ASSERT(edges.size()==nedges, 		"edges.size()!=nedges");
	NEGF_ASSERT(faces.size()==nfaces,		"faces.size()!=nfaces");
	NEGF_ASSERT(elements.size()==nelems, 	"elements.size()!=nelems");
	NEGF_ASSERT(regions.size()==nregions, 	"regions.size()!=nregions");
	NEGF_ASSERT(contacts.size()==ncontacts, "contacts.size()!=ncontacts");
	
	// set the dimension
	NEGF_ASSERT(vertices.size()!=0, "You attempted to finalize the geometry object with zero vertices!");
	this->dimension = vertices[0]->get_dimension();

	// set up internal indices
	int current_internal_index = 0;
	this->internal_to_global_vertex_indices.clear();
	for(uint ii = 0; ii < vertices.size(); ii++) {
		if(vertices[ii]->get_index_internal() != -1) {
			vertices[ii]->set_index_internal(current_internal_index++);
			this->internal_to_global_vertex_indices.push_back(/*ii*/vertices[ii]->get_index_global());
			NEGF_ASSERT(int(internal_to_global_vertex_indices.size()-1)==vertices[ii]->get_index_internal(), "something's wrong.");
		}	
	}
	this->num_internal_vertices = this->internal_to_global_vertex_indices.size();
	logmsg->emit(LOG_INFO,"There are %d internal vertices.", num_internal_vertices);
	
	// no prepare() and is_ready() exists for vertices and edges since this was all done in the constructor
	
	// compute the face areas
	logmsg->emit(LOG_INFO_L2, "   computing face areas");
	for(uint ii = 0; ii < this->nfaces; ii++) {
		NEGF_ASSERT(this->faces[ii] != NULL, "faces[ii] was NULL!");
		if (!faces[ii]->is_ready()) {
			faces[ii]->prepare();
		}
	}
	
	// compute the element volumes, jacobians and diameters
	logmsg->emit(LOG_INFO_L2, "   computing element volumes and inverse jacobians");
	for(int ii = 0; ii < (int) this->nelems; ii++) {
		NEGF_ASSERT(this->elements[ii] != NULL, "elements[ii] was NULL!");
		if (!this->elements[ii]->is_ready()) {
			this->elements[ii]->prepare();
		}
	}
	
	// compute edges_near_vertex
	logmsg->emit(LOG_INFO_L2, "   preparing array edges_near_vertex");
	this->edges_near_vertex.clear();
	this->edges_near_vertex.resize(nvertices);	// fills it with empty vectors
	uint vertex_idx = 0;
	for (uint ii = 0; ii < this->nedges; ii++)
	{
		NEGF_ASSERT(edges[ii]!=NULL, "encountered NULL edge.");
		vertex_idx = edges[ii]->get_lower_vertex()->get_index_global();
		NEGF_FASSERT(vertex_idx < nvertices, "vertex_idx=%d, nvertices=%d", vertex_idx, nvertices);
		(this->edges_near_vertex[vertex_idx]).push_back(edges[ii]);
		vertex_idx = edges[ii]->get_upper_vertex()->get_index_global();
		NEGF_FASSERT(vertex_idx < nvertices, "vertex_idx=%d, nvertices=%d", vertex_idx, nvertices);
		(this->edges_near_vertex[vertex_idx]).push_back(edges[ii]);
	}
	
	// compute max_connectivity
	this->max_connectivity = 0;
	int max_conn_vert = -1;
	for (uint ii = 0; ii < this->nvertices; ii++)
	{
		if (this->get_num_edges_near(get_vertex(ii)) > max_connectivity) {
			this->max_connectivity = get_num_edges_near(get_vertex(ii));
			max_conn_vert = ii;
		}
	}
	logmsg->emit(LOG_INFO_L2,  "   max_connectivity was determined to be %d at vertex %d",max_connectivity, max_conn_vert);
	
	// ------------------------------------------
	// compute X_near_Y
	// ------------------------------------------
	
	// compute array elements_near_edge
	logmsg->emit(LOG_INFO_L2,  "   preparing array elements_near_edge");
	this->elems_near_edge.clear();
	this->elems_near_edge.resize(nedges);	// fills it with empty vectors
	for (uint ii = 0; ii < this->nelems; ii++)
	{
		for (uint jj = 0; jj < elements[ii]->get_num_edges(); jj++)
		{
			int edge_idx = elements[ii]->get_edge(jj)->get_index_global();
			(this->elems_near_edge[edge_idx]).push_back(elements[ii]);
		}
	}
	for (uint ii = 0; ii < this->elems_near_edge.size(); ii++) {
		if (elems_near_edge[ii].size()==0) logmsg->emit(LOG_WARN,"warning: edge %d has 0 surrounding elements.",ii);
	}
	
	// compute array elems_near_face
	logmsg->emit(LOG_INFO_L2,  "   preparing array elements_near_face");
	this->elems_near_face.clear();
	this->elems_near_face.resize(nfaces);	// fills it with empty vectors
	for (int ii = 0; ii < (int) this->nfaces; ii++)
	{
		(this->elems_near_face[ii]).clear();
		uint vidx = this->get_face(ii)->get_vertex(0)->get_index_global();
		for (uint jj = 0; jj < this->get_num_edges_near(get_vertex(vidx)); jj++)
		{
			uint edge_idx = ((this->edges_near_vertex[vidx])[jj])->get_index_global();
			
			// add all the elements around this edge to the list (w/o duplicates)
			// ... if the face is part of the element (does not have to be the case)
			for (uint kk = 0; kk < this->get_num_elems_near(get_edge(edge_idx)); kk++)
			{
				Element * elem = this->elems_near_edge[edge_idx][kk];
			
				// check if that element has ii as a face - if not, do not add to list
				bool face_is_part_of_element = false;
				for (uint ll = 0; ll < elem->get_num_faces(); ll++) {
					if (elem->get_face(ll)==this->get_face(ii)) {
						face_is_part_of_element = true;
						break;
					}
				}
				if (!face_is_part_of_element) {
					continue;
				}
				
				// do not add if already in list
				bool found = false;
				for (uint ll = 0; ll < this->get_num_elems_near(get_face(ii)); ll++) {
					if (elem==this->elems_near_face[ii][ll]) {
						found = true;
						break;
					}
				}
				if (!found) {
					(this->elems_near_face[ii]).push_back(elem);
				}
			}
		}
	}	
	
	// compute array elements_near_vertex
	// even though its a lot of for-loops, they're small
	logmsg->emit(LOG_INFO_L2,  "   preparing array elements_near_vertex");
	this->elems_near_vertex.clear();
	this->elems_near_vertex.resize(nvertices);	// fills it with empty vectors
	for (int ii = 0; ii < (int) this->nvertices; ii++)
	{
		(this->elems_near_vertex[ii]).clear();
		for (uint jj = 0; jj < this->get_num_edges_near(get_vertex(ii)); jj++)
		{
			NEGF_FASSERT(this->edges_near_vertex[ii].size()>jj, "edges_near_vertex[%d].size()=%d > jj=%d failed.", ii, edges_near_vertex[ii].size(), jj);
			uint edge_idx = ((this->edges_near_vertex[ii])[jj])->get_index_global();
			// add all the elements around this edge to the list (w/o duplicates)
			for (uint kk = 0; kk < this->get_num_elems_near(get_edge(edge_idx)); kk++)
			{
				bool found = false;
				for (uint ll = 0; ll < this->get_num_elems_near(get_vertex(ii)); ll++) {
					NEGF_ASSERT(this->elems_near_edge[edge_idx].size()>kk, "elems_near_edge[edge_idx].size()>kk failed.");
					NEGF_ASSERT(this->elems_near_vertex[ii].size()>ll, "elems_near_vertex[ii].size()>ll failed.");
					if ((this->elems_near_edge[edge_idx])[kk]==(this->elems_near_vertex[ii])[ll]) {
						found = true;
						break;
					}
				}
				if (!found) {
					(this->elems_near_vertex[ii]).push_back((elems_near_edge[edge_idx])[kk]);
				}
			}
		}
	}
	
	// compute array faces_near_vertex (relies on elems_near_vertex array)
	logmsg->emit(LOG_INFO_L2,  "   preparing array faces_near_vertex");
	this->faces_near_vertex.clear();
	this->faces_near_vertex.resize(nvertices);	// fills it with empty vectors
	for (int ii = 0; ii < (int) this->nvertices; ii++)
	{
		for (uint jj = 0; jj < this->get_num_elems_near(get_vertex(ii)); jj++)
		{
			Element * elem = elems_near_vertex[ii][jj];
			
			// go through the faces of this element
			// add to array w/o duplicates if vertex is part of face 
			for (uint kk = 0; kk < elem->get_num_faces(); kk++)
			{
				Face * face  = elem->get_face(kk);
				bool found_vertex = false;
				for (uint ll = 0; ll < face->get_num_vertices(); ll++)
				{
					if (face->get_vertex(ll) == this->vertices[ii]) {
						found_vertex = true;
						break;
					}
				}
				if (!found_vertex)	// in that case vertex is not part of the face
					continue;
				
				bool found = false;
				for (uint ll = 0; ll < this->get_num_faces_near(get_vertex(ii)); ll++) {
					if (face == faces_near_vertex[ii][ll]) {
						found = true;
						break;
					}
				}
				if (!found) {
					(this->faces_near_vertex[ii]).push_back(face);
				}
			}
		}
	}
	
	// compute array regions_near_vertex
	logmsg->emit(LOG_INFO_L2,  "   preparing array regions_near_vertex");
	this->regions_near_vertex.clear();
	this->regions_near_vertex.resize(nvertices);	// fills it with empty vectors
	for (uint ii = 0; ii < this->nvertices; ii++)
	{
		const vector<Element *> & elems_near_vert = elems_near_vertex[ii];
		for (uint jj = 0; jj < this->get_num_elems_near(get_vertex(ii)); jj++)
		{
			Region * region = elems_near_vert[jj]->get_region();
			bool found = false;
			for (uint kk = 0; kk < this->get_num_regions_near(get_vertex(ii)); kk++) {
				if (regions_near_vertex[ii][kk]==region) {
					found = true;
					break;
				}
			}
			if (!found) {
				(this->regions_near_vertex[ii]).push_back(region);
			}
		}
	}
	
	// --------------------------------------------------------------------------------------------
	// NEGF: add all vertices lying entirely within a region "contact_*" to the nearest contact
	//       like this the Poisson equation will effectively only be solved on the other vertices
	// --------------------------------------------------------------------------------------------
	logmsg->emit(LOG_INFO_L2,  "   adding vertices in \"contact*\"-regions to contacts (NEGF special)");
	for (uint ii = 0; ii < this->nvertices; ii++)
	{
		Vertex * v = this->get_vertex(ii);
		const vector<Region *> & regs_near_vertex = this->get_regions_near(v);
		if (regs_near_vertex.size()==1 && regs_near_vertex[0]->get_name().substr(0,8)=="contact_") 
		{
			double mindist = 1e10;
			Contact * nearest_contact = 0;
			for (uint jj = 0; jj < this->get_num_contacts(); jj++) {
				double dist = this->get_distance(v, this->get_contact(jj));
				if (dist < mindist) {
					mindist = dist;
					nearest_contact = this->get_contact(jj);
				}
			}
			NEGF_ASSERT(nearest_contact!=0, "nearest contact could not be found.");
			nearest_contact->add_vertex(v); // also assigns the contact to the vertex
			NEGF_ASSERT(v->is_at_contact() && v->get_contact()==nearest_contact, "contact was not automatically asigned to vertex.");
			logmsg->emit(LOG_INFO,"vertex %d in region %s was now assigned to contact %s.",
					ii, regs_near_vertex[0]->get_name().c_str(), nearest_contact->get_name().c_str());
		}
	}
	
	// -------------------------------------------------------------------------------------
	// find the contact's adjacent region
	// strategy: for each contact vertex we assign the region of its adjacent elements,
	//           provided the element has at least 2 contact vertices
	// -------------------------------------------------------------------------------------
	logmsg->emit(LOG_INFO_L2,  "   preparing adjacent_material of contacts");
	for (uint ii = 0; ii < this->ncontacts; ii++)
	{
		NEGF_ASSERT( contacts[ii]->get_num_contact_vertices()>0, 
				"at least one contact vertex must exist for each contact.");
		const vector<Vertex *> & contact_vertices = contacts[ii]->get_contact_vertices();
		bool inconsistent_regions = false;
		for (uint jj=0; jj<contacts[ii]->get_num_contact_vertices(); jj++)
		{
			const Vertex * vert = contacts[ii]->get_contact_vertex(jj);
			const vector<Element*> & elems_near_vert = this->get_elems_near(vert);
			NEGF_ASSERT(elems_near_vert.size()>0, "A contact vertex had no adjacent element.");
			if (this->get_dimension()==1) {
				// NO ASSERT IN NEGF!
				//NEGF_ASSERT(contacts[ii]->get_num_contact_vertices()==1 && elems_near_vert.size()==1,
				//			"Unexpected 1D structure.");
				for (uint kk=0; kk < elems_near_vert.size(); kk++) {
					if (   (elems_near_vert[kk]->get_region()->get_name().length()>=7)
						&& (elems_near_vert[kk]->get_region()->get_name().substr(0,7)=="region_") ) {
						continue;
					} else {
						contacts[ii]->set_adjacent_region(elems_near_vert[kk]->get_region());
					}
				}
			} else {
				// rule out "boundary" vertices which are adjacent to several regions 
				const vector<Region *> & regions_near_vert = this->get_regions_near(vert);
				if (regions_near_vert.size()>1)
					continue;
				
				for (uint kk = 0; kk < elems_near_vert.size(); kk++) {
					// it is possible that the element has no edge in common with the contact 
					// (if the vertex is at the contact boundary). Therefore we must check if
					// there is a common edge, or equivalently if there is a second element vertex
					// in the list of contact vertices
					bool second_vertex = false;
					for (uint ll=0; ll < elems_near_vert[kk]->get_num_vertices(); ll++) {
						if ((elems_near_vert[kk])->get_vertex(ll) != vert) {
							if ( find(contact_vertices.begin(), contact_vertices.end(),elems_near_vert[kk]->get_vertex(ll))
								 !=contact_vertices.end() ) {
								second_vertex = true;
								break;
							}
						}
					}
					if (!second_vertex)
						continue;
					
					// now we assign the element's region to the contact and check for consistency
					if (contacts[ii]->get_adjacent_region()==NULL) { // default value in Contact constructor
						contacts[ii]->set_adjacent_region(elems_near_vert[kk]->get_region());
					} else {
						//NEGF_ASSERT(contacts[ii]->get_adjacent_region()==elems_near_vert[kk]->get_region(),
						//		"Inconsistency: elements adjacent to the same contact have different regions.");
						if (contacts[ii]->get_adjacent_region()!=elems_near_vert[kk]->get_region()) {
							inconsistent_regions = true; 
						}
					}
				}
			}
		}
		NEGF_FASSERT(contacts[ii]->get_adjacent_region()!=NULL,
				 "Adjacent region of contact %s could not be determined. Are there no internal vertices?",
				 contacts[ii]->get_name().c_str());
		if (inconsistent_regions) {
			logmsg->emit(LOG_WARN,  "WARNING: Ambiguity in adjacent region to contact %s. Taking region %s.",
							contacts[ii]->get_name().c_str(), contacts[ii]->get_adjacent_region()->get_name().c_str());
		}
	}
	
	// -------------------------------------------------------------------------------------------------
	// calculate for each contact its area (3D), length (2D); in 1D, area=1 and nothing has to be done
	// -------------------------------------------------------------------------------------------------
	logmsg->emit(LOG_INFO_L2, "   preparing area of contacts");
	for (uint ii = 0; ii < this->ncontacts; ii++)
	{		
		NEGF_ASSERT( contacts[ii]->get_num_contact_vertices()>0, 
				"at least one contact vertex must exist for each contact.");
		const vector<Vertex *> & contact_vertices = contacts[ii]->get_contact_vertices();
		
		if (this->get_dimension()==1)		// can't make a switch because of const vector<...> defs
		{
			contacts[ii]->set_area(1.0);
		} else if (this->get_dimension()==2)
		{
			// assemble list of all edges that constitute the contact
			vector<Edge *> contact_edges;
			for (uint jj = 0; jj < contact_vertices.size(); jj++)
			{
				const vector<Edge *> & edges_near_vert = this->get_edges_near(contact_vertices[jj]);
				for (uint kk = 0; kk < edges_near_vert.size(); kk++) 
				{
					// check if edge is within contact by looking at the second vertex of the edge
					Edge * edge = edges_near_vert[kk];
					Vertex * other_vertex;
					if (edge->get_lower_vertex()==contact_vertices[jj]) {
						other_vertex = edge->get_upper_vertex();
					} else if (edge->get_upper_vertex()==contact_vertices[jj]) {
						other_vertex = edge->get_lower_vertex();
					} else {
					NEGF_EXCEPTION("get_edges_near(vertx) was faulty.");
					}
					if (other_vertex->is_at_contact()) {
						NEGF_ASSERT(other_vertex->get_contact()==contacts[ii],
										"An edge seemed to be connected to two different contacts.");
					} else {
						continue;
					}
					
					// if edge is not already in the list, add it
					if (find(contact_edges.begin(), contact_edges.end(), edge) == contact_edges.end()) {
					 	contact_edges.push_back(edge);
					}
				}
			}
			
			// the length of the contact is simply the sum of all edge lengths
			double sum = 0.0;
			for (uint jj = 0; jj < contact_edges.size(); jj++)
				sum += contact_edges[jj]->get_length();
			
			contacts[ii]->set_area(sum);
		} else if (this->get_dimension()==3)
		{
			// assemble list of all faces that constitute the contact
			vector<Face *> contact_faces;
			for (uint jj = 0; jj < contact_vertices.size(); jj++)
			{
				const vector<Element *> & elems_near_vert = this->get_elems_near(contact_vertices[jj]);
				for (uint kk = 0; kk < elems_near_vert.size(); kk++) 
				{
					NEGF_ASSERT(elems_near_vert[kk]->get_num_faces()>0, "no faces but should be 3D element.");
					for (uint ll = 0; ll < elems_near_vert[kk]->get_num_faces(); ll++)
					{
						Face * face = elems_near_vert[kk]->get_face(ll);
						
						// check if all face vertices are within contact
						bool all_vertices_in_contact = true;
						for (uint mm = 0; mm < face->get_num_vertices(); mm++)
						{
							if (face->get_vertex(mm)->is_at_contact()) {
								NEGF_ASSERT(face->get_vertex(mm)->get_contact()==contacts[ii],
									"a face seemed to be connected to more than 1 contact.");
							} else {
								all_vertices_in_contact = false;
								break;
							}
						}
						
						// if not, skip face; if yes, add face to list if it is not already there
						if (!all_vertices_in_contact) {
							continue;
						} else {
							if (find(contact_faces.begin(), contact_faces.end(), face)==contact_faces.end()) {
								contact_faces.push_back(face);
							}
						}
					}
				}
			}
			
			// the area of the contact is simply the sum of all face areas
			double sum = 0.0;
			for (uint jj = 0; jj < contact_faces.size(); jj++)
				sum += contact_faces[jj]->get_area();
			
			contacts[ii]->set_area(sum);
		} else {
			NEGF_EXCEPTION("Strange dimensionality.");
		}
		logmsg->emit(LOG_INFO_L2,  "      contact %s: area is %e",
							contacts[ii]->get_name().c_str(), contacts[ii]->get_area());
	}

	this->prepared = true;
	logmsg->emit(LOG_INFO_L1,  "preparation of geometry object done.");
);}


uint Geometry::get_global_vertex_index(uint internal_index) const
{STACK_TRACE(
	NEGF_FASSERT(this->internal_to_global_vertex_indices.size() > internal_index, 
			"invalid internal index %d.",internal_index);
	return internal_to_global_vertex_indices[internal_index];
);}


uint Geometry::get_internal_vertex_index(uint global_index) const
{STACK_TRACE(
	return this->get_vertex(global_index)->get_index_internal();
);}


/** verify parsed geometry
 * 
 *  1. check if 'announced' number of vertices/edges/faces/elements/regions corresponds to 
 *     the created size of v/e/f/e/r.
 *  2. check if there are entries in the vertex/edge/face/element pointer lists pointing to NULL
 *  3. check if there are no vertices with same index
 *  4. check if boundary conditions are set
 *  5. check if the global vertex index corresponds to the position it was added
 *  6. check if no element has two times the same vertex
 *  7. perform elem->verify() and check if there are no elements with equal vertices/edges/faces
 *  8. perform face->verify() and check if there are no faces with equal vertices/edges
 *  9. perform edge->verify() and check if there are no edges with equal vertices
 * 10. check the arrays edges_near_vertex, elements_near_edge
 */
bool Geometry::verify() const throw (Exception *)
{STACK_TRACE(
	const bool check_equal_elements = false;
	const bool check_equal_faces    = false;
	const bool check_equal_edges    = false;
	const bool check_equal_vertices = false;
	const bool check_equal_internal_indices = false;

	logmsg->emit(LOG_INFO_L2, "********* checking geometry... ***********");

	logmsg->emit(LOG_INFO_L2, "   checking vector sizes");	
	// -----------
	// check sizes
	// -----------
	if(vertices.size() != nvertices) {
		logmsg->emit(LOG_WARN,  "inconsistent geometry. %d vertices have been announced, %d were defined", nvertices, vertices.size());			
		return false;
	}
	if(edges.size() != nedges) {
		logmsg->emit(LOG_WARN,  "inconsistent geometry. %d edges have been announced, %d were defined", nedges, edges.size());			
		return false;
	}
	if(faces.size() != nfaces) {
		logmsg->emit(LOG_WARN,  "inconsistent geometry. %d faces have been announced, %d were defined", nfaces, faces.size());			
		return false;
	}
	if(elements.size() != nelems) {
		logmsg->emit(LOG_WARN,  "inconsistent geometry. %d elements have been announced, %d were defined", nelems, elements.size());
		return false;
	}
	if(regions.size() != nregions) {
		logmsg->emit(LOG_WARN,  "inconsistent geometry. %d regions have been announced, %d were defined", nregions, regions.size());
		return false;
	}
	if(contacts.size() != ncontacts) {
		logmsg->emit(LOG_WARN,  "inconsistent geometry. %d contacts have been announced, %d were defined", ncontacts, contacts.size());
		return false;
	}
	// --------------------------------
	// check if there are null pointers
	// --------------------------------
	logmsg->emit(LOG_INFO_L2, "   checking for null pointers");
	for (uint ii = 0; ii < contacts.size(); ii++) {
		if (contacts[ii]==NULL) {
			logmsg->emit(LOG_ERROR,  "found NULL pointer in contact list.");
			return false;
		}
	}
	for (uint ii = 0; ii < regions.size(); ii++) {
		if (regions[ii]==NULL) {
			logmsg->emit(LOG_ERROR,  "found NULL pointer in region list.");
			return false;
		}
	}
	for (uint ii = 0; ii < elements.size(); ii++) {
		if (elements[ii]==NULL) {
			logmsg->emit(LOG_ERROR,  "found NULL pointer in element list.");
			return false;
		}
	}
	for (uint ii = 0; ii < faces.size(); ii++) {
		if (faces[ii]==NULL) {
			logmsg->emit(LOG_ERROR,  "found NULL pointer in face list.");
			return false;
		}
	}
	for (uint ii = 0; ii < edges.size(); ii++) {
		if (edges[ii]==NULL) {
			logmsg->emit(LOG_ERROR,  "found NULL pointer in edge list.");
			return false;
		}
	}
	for (uint ii = 0; ii < vertices.size(); ii++) {
		if (vertices[ii]==NULL) {
			logmsg->emit(LOG_ERROR,  "found NULL pointer in vertex list.");
			return false;
		}
	}

	// -----------------------------------------------------------------------------------
	// calculate number of boundary vertices and check whether all nodes have proper index
	// -----------------------------------------------------------------------------------
	logmsg->emit(LOG_INFO_L2, "   checking boundary vertices, internal indices");
	int current_index;  
	int bvert 		 = 0;
	uint global_index = 0;
	for(uint ii = 0; ii < vertices.size(); ii++) {
		// check global index
		if(global_index != vertices[ii]->get_index_global()) {
			logmsg->emit(LOG_ERROR,  "inconsistent geometry - vertices were not added sequentially. vertex %d was added as %d", 
								vertices[ii]->get_index_global(), global_index);
			return false;
		}
		// check internal index
		current_index = vertices[ii]->get_index_internal();
		if(current_index == -1) {
			bvert++;
		} else if(check_equal_internal_indices) {
			// check that no other node has my index
			for(uint jj = ii + 1; jj < vertices.size(); jj++) {
				if(current_index == vertices[jj]->get_index_internal()) {
					logmsg->emit(LOG_WARN,  "inconsistent geometry. vertex: %d has same internal index (%d) like vertex %d", 
							vertices[ii]->get_index_global(), vertices[ii]->get_index_internal(), vertices[jj]->get_index_global());
					return false;
				}
			}	
		}	
		global_index++;
	}
	logmsg->emit(LOG_INFO_L2,  "   there are %d vertices at boundary and %d vertices inside", bvert, nvertices - bvert);

	// -----------------------------------------------------
	// check regions
	// -----------------------------------------------------
	vector<bool> elem_found; elem_found.resize(this->get_num_elements(), false);
	for (uint ii = 0; ii < this->get_num_regions(); ii++) {
		const vector<Element *> elems = this->get_region(ii)->get_elements();
		for (uint jj=0; jj < elems.size(); jj++) {
			NEGF_ASSERT(!elem_found[elems[jj]->get_index_global()], "an element seemed to be contained in more than one region.");
			elem_found[elems[jj]->get_index_global()] = true;
		}
	}
	for (uint ii = 0; ii < elem_found.size(); ii++) {
		NEGF_ASSERT(elem_found[ii], "an element seemed to be not contained in any region.");
	}
	
	// -----------------------------------------------------
	// check elements
	// -----------------------------------------------------
	int next = 0; int count = 0;
	if (check_equal_elements)
		logmsg->init_progress_bar(LOG_INFO_L2,"   checking elements", nelems);
	uint idx       = 0;
	uint idx_other = 0;
	double min_volume = 1e10;
	int min_vol_element = -1;
	double total_volume = 0.0;
	for(uint ee = 0; ee < elements.size(); ee++) {
		Element * eit = elements[ee];
		// check if an element has two times the same vertex
		if(!eit->verify()) {
			logmsg->emit(LOG_ERROR,  "element %d claims to be invalid", idx);
			return false;
		}
		if (eit->get_volume() < min_volume)
		{
			min_volume = eit->get_volume();
			min_vol_element = eit->get_index_global();
		}
		total_volume += eit->get_volume();

		if (check_equal_elements)
		{
			// check if there is another element with the same vertices
			const uint eit_nvert = eit->get_num_vertices();
//			uint eit_vertex_indices[eit_nvert];		// stores the vertex indices of element eit
			uint eit_vertex_indices[30];
			eit->get_vertex_indices(eit_vertex_indices);
			sort(eit_vertex_indices, eit_vertex_indices + eit_nvert);
//			uint ecit_vertex_indices[eit_nvert];		// stores the vertex indices of element ecit
			uint ecit_vertex_indices[30];			// stores the vertex indices of element ecit
			idx_other = idx + 1;
			for(uint ff = ee + 1; ff < elements.size(); ff++) 
			{
				Element * ecit = elements[ff];
				if (eit->get_num_vertices() != ecit->get_num_vertices()) continue;
				ecit->get_vertex_indices(ecit_vertex_indices);
				sort(ecit_vertex_indices, ecit_vertex_indices + eit_nvert);
				bool ok = false;
				for (uint ii = 0; ii < eit_nvert; ii++)
				{
					if (ecit_vertex_indices[ii] != eit_vertex_indices[ii]) {
						ok = true;
						break;
					}
				}
				if (!ok) {
					logmsg->emit(LOG_ERROR,  "element %d has the same vertices as element %d", idx, idx_other);
					return false;
				}
	/*
				if((*eit)->get_vertices() == ecit->get_vertices()) {
					logmsg->emit(LOG_ERROR,  "element %d has the same vertices as element %d", idx, idx_other);
					return false;
				}
				if((*eit)->get_edges() == ecit->get_edges()) {
					logmsg->emit(LOG_ERROR,  "element %d has the same edges as element %d", idx, idx_other);
					return false;
				}
				if((*eit)->get_faces() == ecit->get_faces()) {
					logmsg->emit(LOG_ERROR,  "element %d has the same faces as element %d", idx, idx_other);
					return false;
				}
	*/
				idx_other++;
			}
			idx++;
			if(next == count++)  next = logmsg->set_progress_bar(count, nelems);
		}
		
	}
	if (check_equal_elements)
		logmsg->end_progress_bar();
	logmsg->emit(LOG_INFO_L2,"   element %d is smallest (volume %e)",min_vol_element,min_volume);
	logmsg->emit(LOG_INFO,"Total volume calculated from elements: %12.10e",total_volume);
	
	// ----------------
	// check faces
	// ----------------
	next = 0; count = 0;
	if (nfaces!=0)
	{
		if (check_equal_faces)
			logmsg->init_progress_bar(LOG_INFO_L2,"   checking faces", nfaces);
		idx = 0;
		idx_other = 0;
		for (uint ff = 0; ff < faces.size(); ff++) 
		{
			Face * fit = faces[ff];
			if(!fit->verify()) {
				logmsg->emit(LOG_ERROR,  "face %d claims to be invalid", idx);
				return false;
			}
			idx_other = idx + 1;
	
			if (check_equal_faces)
			{
				// check if there is another face with the same vertices
				int  fit_nvert = fit->get_num_vertices(); 
				// long fit_vertex_indices[fit_nvert];			// stores the vertex indices of face fit
				long fit_vertex_indices[20];
				fit->get_vertex_indices(fit_vertex_indices);
				sort(fit_vertex_indices, fit_vertex_indices + fit_nvert);
				// long fcit_vertex_indices[fit_nvert];		// stores the vertex indices of face fcit
				long fcit_vertex_indices[20];
				for(uint gg = ff + 1; gg < faces.size(); gg++) 
				{
					Face * fcit = faces[gg];
					if (fit->get_num_vertices() != fcit->get_num_vertices()) continue;
					fcit->get_vertex_indices(fcit_vertex_indices);
					sort(fcit_vertex_indices, fcit_vertex_indices + fit_nvert);
					bool ok = false;
					for (int ii = 0; ii < fit_nvert; ii++)
					{
						if (fcit_vertex_indices[ii] != fit_vertex_indices[ii]) {
							ok = true;
							break;
						}
					}
					if (!ok) {
						logmsg->emit(LOG_ERROR,  "face %d has the same vertices as face %d", idx, idx_other);
						return false;
					}
		/*
					if((*fit)->get_vertices() == fcit->get_vertices()) {
						logmsg->emit(LOG_ERROR,  "face %d has the same vertices as face %d", idx, idx_other);
						return false;
					}
					if((*fit)->get_edges() == fcit->get_edges()) {
						logmsg->emit(LOG_ERROR,  "face %d has the same edges as face %d", idx, idx_other);
						return false;
					}
		*/
					idx_other++;
				}
				idx++;
				if(next == count++)  next = logmsg->set_progress_bar(count, nfaces);
			}
		}
		if (check_equal_faces)
			logmsg->end_progress_bar();
	}

	// ----------------
	// check edges
	// ----------------
	next = 0; count = 0;
	if (check_equal_edges)
		logmsg->init_progress_bar(LOG_INFO_L2,"   checking edges", nedges);
	idx = 0;
	idx_other = 0;
	double min_length = 1e10;
	int min_edge = -1;
	for (uint ee = 0; ee != edges.size(); ee++) 
	{
		Edge * edit = edges[ee];
		if(!edit->verify()) {
			logmsg->emit(LOG_ERROR,  "edge %d claims to be invalid", idx);
			return false;
		}
		idx_other = idx + 1;
		if (edit->get_length() < min_length)
		{
			min_length = edit->get_length();
			min_edge = edit->get_index_global();
		}

		if (check_equal_edges)
		{
			// check if there is another edge with the same vertices
			for(uint ff = ee + 1; ff < edges.size(); ff++) 
			{
				Edge * edcit = edges[ff];
				if(   (   edit->get_lower_vertex() == edcit->get_lower_vertex()
					&& edit->get_upper_vertex() == edcit->get_upper_vertex() )
					||(   edit->get_lower_vertex() == edcit->get_upper_vertex()
					&& edit->get_upper_vertex() == edcit->get_lower_vertex() ) ) {
					logmsg->emit(LOG_ERROR,  "edge %d has the same vertices as edge %d", idx, idx_other);
					return false;
				}
				idx_other++;
			}
			idx++;
			if(next == count++)  next = logmsg->set_progress_bar(count, nedges);
		}
	}
	if (check_equal_edges)
		logmsg->end_progress_bar();
	logmsg->emit(LOG_INFO_L2,"   edge %d is shortest (length %e)",min_edge,min_length);

	// ---------------------------------------------------
	// check if there are (nearly) identical vertices
	// ---------------------------------------------------
	if (check_equal_vertices)
	{
		next = 0; count = 0;
		logmsg->init_progress_bar(LOG_INFO_L2, "   checking for near-equal vertices", nvertices);
		double epsilon = 1e-8;
		unsigned int dim = vertices[0]->get_dimension();
		for (idx = 0; idx != this->get_num_vertices(); idx++) {
			for (idx_other = idx + 1; idx_other != this->get_num_vertices(); idx_other++) {
				switch (dim)
				{
				case 1:
					if (fabs(vertices[idx]->get_coordinate(0) - vertices[idx_other]->get_coordinate(0)) < epsilon) 
					{
						logmsg->emit(LOG_ERROR,  "the distance between vertices %d and %d is smaller than %d", idx, idx_other, epsilon);
						return false;
					}
					break;
				case 2: // check if |x1-x2|+|y1-y2| < 2*eps
					if (fabs(vertices[idx]->get_coordinate(0) - vertices[idx_other]->get_coordinate(0)) < epsilon)
					{
						if (fabs(vertices[idx]->get_coordinate(1) - vertices[idx_other]->get_coordinate(1)) < epsilon)
						{
							logmsg->emit(LOG_ERROR,  "the distance between vertices %d and %d is smaller than %d", idx, idx_other, epsilon);
							return false;
						}
					}
					break;
				case 3: // check if |x1-x2|+|y1-y2|+|z1-z2| < 3*eps
					if (fabs(vertices[idx]->get_coordinate(0) - vertices[idx_other]->get_coordinate(0)) < epsilon)
					{
						if (fabs(vertices[idx]->get_coordinate(1) - vertices[idx_other]->get_coordinate(1)) < epsilon)
						{
							if (fabs(vertices[idx]->get_coordinate(2) - vertices[idx_other]->get_coordinate(2)) < epsilon)
							{
								logmsg->emit(LOG_ERROR,  "the distance between vertices %d and %d is smaller than %d", idx, idx_other, epsilon);
								return false;
							}
						}
					}
					break;
				}
			}
			if(next == count++)  next = logmsg->set_progress_bar(count, nvertices);
		}
		logmsg->end_progress_bar();
	}

	// -------------------------------------------------
	// check edges_near_vertex
	// -------------------------------------------------
	next = 0; count = 0;
	logmsg->init_progress_bar(LOG_INFO_L2,"   checking edges_near_vertex...",nvertices);
	for (uint ii = 0; ii < nvertices; ii++)
	{
		if (this->get_num_edges_near(get_vertex(ii)) <= 1 && this->get_dimension()>1 )
			logmsg->emit(LOG_ERROR,  "vertex %d seems to be connected to <=1 edges", ii);
		vector<Edge *> edge_array = get_edges_near(get_vertex(ii));
		//if (edge_array.size() != this->num_edges_near_vertex[ii])
		//	logmsg->emit(LOG_ERROR,  "stored number of connected edges to a vertex (%d) does not coincide with size of corresponding vector (%d).",
		//					this->num_edges_near_vertex[ii], edge_array.size());
		for (uint jj = 0; jj < edge_array.size(); jj++)
		{
			if (edge_array[jj]==NULL)
				logmsg->emit(LOG_ERROR,  "empty edge pointer found.");
			if (   !((edge_array[jj])->get_lower_vertex()==get_vertex(ii))
				&& !((edge_array[jj])->get_upper_vertex()==get_vertex(ii)) )
				logmsg->emit(LOG_ERROR,  "a connected edge does not contain the vertex.");
		}
		if(next == count++)  next = logmsg->set_progress_bar(count, nvertices);
	}
	logmsg->end_progress_bar();

	// -------------------------------------------------
	// check elems_near_edge
	// -------------------------------------------------
	next = 0; count = 0;
	logmsg->init_progress_bar(LOG_INFO_L2,"   checking elems_near_edge...",nedges);
	for (uint ii = 0; ii < nedges; ii++)
	{
		if (get_num_elems_near(get_edge(ii)) <= 0)
			logmsg->emit(LOG_ERROR,  "edge %d does not seem to be connected to any elements", ii);
		vector<Element *> elem_array = get_elems_near(get_edge(ii));
		//if (elem_array.size() != get_num_elems_near(get_edge(ii))) {
		//	logmsg->emit(LOG_ERROR,  "stored number of adj. elems near edge (%d) does not coincide with size of corresponding vector (%d).",
		//		get_num_elems_near(get_edge(ii)), elem_array.size());
		//}
		for (uint jj = 0; jj < get_num_elems_near(get_edge(ii)); jj++)
		{
			if (elem_array[jj]==NULL)
				logmsg->emit(LOG_ERROR,  "empty element pointer found.");
			bool ok = false;
			for (int kk = 0; kk < elem_array[jj]->get_num_edges(); kk++)
			{
				if (elem_array[jj]->get_edge(kk) == get_edge(ii)) ok = true;
			}
			if (!ok)
				logmsg->emit(LOG_ERROR,  "an adjacent element does not contain the edge.");
		}
		if(next == count++)  next = logmsg->set_progress_bar(count, nedges);
	}
	logmsg->end_progress_bar();
	
	// -------------------------------------------------
	// check adjacent_material of contacts
	// -------------------------------------------------
	logmsg->emit(LOG_INFO_L2,"   checking adjacent_region of contacts...");
	for (uint ii=0; ii < ncontacts; ii++)
	{
		NEGF_ASSERT(contacts[ii]->get_adjacent_region()!=NULL,
				"A contact's adjacent region was not found.");
	}


	logmsg->emit(LOG_INFO_L1,  "geometry looks good!");

	return true;
);}


/** add vertex to geometry 
 * the passed vertex will be deleteted when the geometry object is deleted
 * vertices have to be added in the sequence of their global indices
 * @param vert pointer to vertex object
 * @return index in vector vertices (== global index)
 */
uint Geometry::add_vertex(Vertex* vert) 
{
	vertices.push_back(vert);
	return vertices.size();
}


/** add edge to geometry
 * @param edge pointer to edge object
 * @return index in vector vertices (== global index)
 */
uint Geometry::add_edge(Edge* edge) 
{
	this->edges.push_back(edge);
	return edges.size();
}


/** add face to geometry */
uint Geometry::add_face(Face* face) 
{
	this->faces.push_back(face);
	return faces.size();
}


/** add element to geometry
 * the passed element will be deleted when the geometry object is deleted
 * @param elem pointer to element object
 * @return index in vector elements (== global element index)
 */
uint Geometry::add_element(Element* elem) 
{
	elements.push_back(elem);
	return elements.size();
}


/** add region to geometry (includes incrementing nregions)
 * the passed region will be deleted when the geometry object is deleted
 * @param reg pointer to region object
 * @return region index
 */
uint Geometry::add_region(Region* reg) 
{STACK_TRACE(
	reg->set_index(nregions);
	regions.push_back(reg);
	nregions++;
	return regions.size();
);}


/** add contact to geometry (includes incrementing ncontacts)
 * the passed contact will be deleted when the geometry object is deleted
 * @param cont pointer to contact object
 * @return contact index
 */
uint Geometry::add_contact(Contact* cont) 
{STACK_TRACE(
	cont->set_index(ncontacts);
	contacts.push_back(cont);
	ncontacts++;
	return contacts.size();
);}


/** return pointer to vertex with global idx vert_idx*/
Vertex* Geometry::get_vertex (uint vert_idx) const throw (Exception *)
{STACK_TRACE(
	NEGF_FASSERT(vert_idx < vertices.size(), "vert_idx(%d)<vertices.size()(%d), nvertices=%d",vert_idx,vertices.size(),nvertices);
	NEGF_FASSERT(vert_idx < nvertices, "vert_idx(%d)<nvertices(%d), vertices.size()=%d",vert_idx,nvertices,vertices.size());
	NEGF_ASSERT(vertices[vert_idx] != NULL, "vertices[idx] pointed to NULL!");
	return this->vertices[vert_idx];	
);}

/** return pointer to edge with global idx edge_idx*/
Edge* Geometry::get_edge (uint edge_idx) const throw (Exception *)
{STACK_TRACE(
	NEGF_FASSERT(edge_idx < nedges, "edge_idx(%d)<nedges(%d)",edge_idx,nedges);
	NEGF_ASSERT(edges[edge_idx] != NULL, "edges[idx] pointed to NULL!");
	return this->edges[edge_idx];	
);}

/** return pointer to face */
Face* Geometry::get_face(uint face_idx) const throw (Exception *)
{STACK_TRACE(
	NEGF_FASSERT(face_idx < this->faces.size(), "face_idx(%d) < this->faces.size()(%d)",face_idx,faces.size());
	NEGF_ASSERT(faces[face_idx] != NULL, "faces[idx] pointed to NULL!");
	return this->faces[face_idx];
);}


/** return pointer to element */
Element* Geometry::get_element(uint elem_idx) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(elem_idx >= 0 && elem_idx < nelems, "elem_idx >= 0 && elem_idx < nelems");
	NEGF_ASSERT(elements[elem_idx] != NULL, "elements[idx] pointed to NULL!");
	return this->elements[elem_idx];	
);}


/** return pointer to contact from index */
Contact* Geometry::get_contact(uint contact_idx) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(contacts.size()==ncontacts, "ncontacts does not coincide with the contact list.");
	NEGF_ASSERT(contact_idx >= 0 && contact_idx < ncontacts, "contact_idx >= 0 && contact_idx < ncontacts");
	NEGF_ASSERT(contacts[contact_idx] != NULL, "contacts[idx] pointed to NULL!");
	return this->contacts[contact_idx];	
);}


/** return pointer to contact from name */
Contact* Geometry::get_contact(string name) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(contacts.size()==ncontacts, "ncontacts does not coincide with the contact list.");

	int contact_idx = -1;
	for (uint ii = 0; ii <= this->get_num_contacts(); ii++)
	{
		NEGF_ASSERT(this->get_contact(ii) != NULL, "encountered null pointer.");
		if (this->get_contact(ii)->get_name()==name)
		{
			NEGF_ASSERT(contact_idx==-1, "Ambiguity in the contact names.");
			contact_idx = ii;
		}
	}
	NEGF_ASSERT(contact_idx=!-1, "contact not found.");
	NEGF_ASSERT(contacts[contact_idx] != NULL, "contacts[idx] pointed to NULL!");
	return this->contacts[contact_idx];	
);}


/** return pointer to region from index*/
Region* Geometry::get_region(uint region_idx) const throw (Exception *)
{STACK_TRACE(
	NEGF_FASSERT(region_idx < nregions, "region_idx=%d was >=  nregions=%d",region_idx,nregions);
	NEGF_ASSERT(regions[region_idx] != NULL, "regions[idx] pointed to NULL!");
	return this->regions[region_idx];	
);}


/** return pointer to region from name*/
Region* Geometry::get_region(string name) const throw (Exception *)
{STACK_TRACE(
	int region_idx = -1;
	for (uint ii = 0; ii <= this->get_num_regions(); ii++)
	{
		NEGF_ASSERT(this->get_region(ii) != NULL, "encountered null pointer.");
		if (this->get_region(ii)->get_name()==name)
		{
			NEGF_ASSERT(region_idx==-1, "Ambiguity in the region names.");
			region_idx = ii;
		}
	}
	NEGF_ASSERT(region_idx=!-1, "region not found.");
	NEGF_ASSERT(regions[region_idx] != NULL, "regions[idx] pointed to NULL!");
	return this->regions[region_idx];	
);}


Edge * Geometry::find_edge(Vertex * v1, Vertex * v2) const throw (Exception *)
{STACK_TRACE(
	if (v1==NULL || v2==NULL)
		NEGF_EXCEPTION("find_edge was called with NULL pointers.");
	if (v1==v2)
		NEGF_EXCEPTION("find_edge was called with identical vertices.");
	for (uint ii = 0; ii < nedges; ii++)
	{
		if (this->get_edge(ii)->get_lower_vertex()==v1)
			if (this->get_edge(ii)->get_upper_vertex()==v2)
				return (this->get_edge(ii));
		if (this->get_edge(ii)->get_upper_vertex()==v1)
			if (this->get_edge(ii)->get_lower_vertex()==v2)
				return (this->get_edge(ii));
	}
	NEGF_EXCEPTION("find_edge was not able to find an edge containing both vertices.");
	return NULL;
);}


double Geometry::get_area_from_edges(Edge * edge1, Edge * edge2) const throw (Exception *)
{STACK_TRACE(
	if (this->get_dimension() != 3)
		NEGF_EXCEPTION("This method only works for 3D.");

	if (edge1 == edge2)
		NEGF_EXCEPTION("get_area_from_edges was called with identical edges.");
	if (edge1==NULL || edge2==NULL)
		NEGF_EXCEPTION("get_area_from_edges was called with NULL pointers.");
	
	// ---------------------------------------------------------------------
	// find a face containing both edges by looking at its adjacent elements
	// ---------------------------------------------------------------------
	vector<Element *> edge1_elems = get_elems_near(edge1);
	vector<Element *> edge2_elems = get_elems_near(edge2);

	if (edge1_elems.size()==0 || edge2_elems.size()==0)
		NEGF_EXCEPTION("one of the edges does not seem to be part of any element!");
	
	// next follows a combination over 4 for-loops
	// say we have 4 adj. elems each, each elem contains 4 faces --> 4^4 = 256 loops
	// if we loop over all faces, we would have like 100'000 loops!
	Face * face = NULL;
	for (uint ii = 0; ii < edge1_elems.size(); ii++)
	{
		for (uint jj = 0; jj < edge1_elems[ii]->get_num_faces(); jj++)
		{
			for (uint kk = 0; kk < edge2_elems.size(); kk++)
			{
				for (uint ll = 0; ll < edge2_elems[kk]->get_num_faces(); ll++)
				{
					if (edge1_elems[ii]->get_face(jj) == edge2_elems[kk]->get_face(ll))
					{
						face = edge1_elems[ii]->get_face(jj);
						break;
					}
				}
				if (face) break;
			}
			if (face) break;
		}
		if (face) break;
	}
	
	if (!face)
		NEGF_EXCEPTION("No face was found containing both edges!");
	
		
	// --------------------------------------------------------------
	// depending on the number of vertices of that face, get the area
	// --------------------------------------------------------------
	// compute the cross-product of the edge vectors
	double edge1_vec[3];
	for (int ii = 0; ii < 3; ii++)
		edge1_vec[ii] = edge1->get_upper_vertex()->get_coordinate(ii) - edge1->get_lower_vertex()->get_coordinate(ii);
	double edge2_vec[3];
	for (int ii = 0; ii < 3; ii++)
		edge2_vec[ii] = edge2->get_upper_vertex()->get_coordinate(ii) - edge2->get_lower_vertex()->get_coordinate(ii);
	double cross_product[3];
	cross_product[0] = edge1_vec[1] * edge2_vec[2] - edge1_vec[2] * edge2_vec[1];
	cross_product[1] = edge1_vec[2] * edge2_vec[0] - edge1_vec[0] * edge2_vec[2];
	cross_product[2] = edge1_vec[0] * edge2_vec[1] - edge1_vec[1] * edge2_vec[0];
	// compute the area
	double area = sqrt(cross_product[0]*cross_product[0] + cross_product[1]*cross_product[1] + cross_product[2]*cross_product[2]);
	switch (face->get_num_vertices())
	{
	case 3:
		area = area / 2.;
		break;
	case 4:
		break;
	default:
		NEGF_EXCEPTION("Faces with other than 3 or 4 vertices not implemented.");
		break;
	}
	
	return area;
);}


double Geometry::get_distance(const Vertex * const v1, const Vertex * const v2) const throw (Exception *)
{STACK_TRACE(
	double a, b, c;
	switch (this->get_dimension())
	{
	case 1:
		return (abs (v1->get_coordinate(0) - v2->get_coordinate(0)));
		break;
	case 2:
		a = v1->get_coordinate(0) - v2->get_coordinate(0);
		b = v1->get_coordinate(1) - v2->get_coordinate(1);
		return sqrt(a*a + b*b);
		break;
	case 3:
		a = v1->get_coordinate(0) - v2->get_coordinate(0);
		b = v1->get_coordinate(1) - v2->get_coordinate(1);
		c = v1->get_coordinate(2) - v2->get_coordinate(2);
		return sqrt(a*a + b*b + c*c);
		break;
	default:
		NEGF_EXCEPTION("Strange dimensionality.");
	}
);}


/** Get the distance of a vertex to an edge 
 *  We first find lambda from P=v0+lambda*(v1-v0) where 
 *  P=nearest point to the vertex on the edge axis,
 *  v0=lower vertex, v1=upper vertex
 *  If lambda<0 or lambda>1, the nearest point is either v0 or v1
 */
double Geometry::get_distance(const Vertex * const vertex, const Edge * const edge) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(vertex!=0 && edge!=0 && vertex->get_dimension()==edge->get_lower_vertex()->get_dimension(),
				"faulty input.");
	double vector1[3], vector2[3], vector3[3], vector4[3];
	double lambda;
	switch (vertex->get_dimension())
	{
	case 3:
		// vector1 --> P1-P0, vector2 --> P-P0
		vector1[0] = edge->get_upper_vertex()->get_coordinate(0) - edge->get_lower_vertex()->get_coordinate(0);
		vector1[1] = edge->get_upper_vertex()->get_coordinate(1) - edge->get_lower_vertex()->get_coordinate(1);
		vector1[2] = edge->get_upper_vertex()->get_coordinate(2) - edge->get_lower_vertex()->get_coordinate(2);
		vector2[0] = vertex->get_coordinate(0) - edge->get_lower_vertex()->get_coordinate(0);
		vector2[1] = vertex->get_coordinate(1) - edge->get_lower_vertex()->get_coordinate(1);
		vector2[2] = vertex->get_coordinate(2) - edge->get_lower_vertex()->get_coordinate(2);
		NEGF_ASSERT(negf_math::vector_norm_3d(vector1) > constants::min_vector_norm, "edge is too small.");
		
		// P = P0 + lambda*(P1-P0) + nu*(vector perpendicular to P1-P0)
		lambda = negf_math::vector_scalar_product_3d(vector1,vector2) / negf_math::vector_scalar_product_3d(vector1,vector1);
		if (lambda <= 0.0) { // lower vertex is nearest
			//cout << "    get_distance_edge: lower vertex (" << edge->get_lower_vertex()->get_coordinate(0) 
			//	<< ", " << edge->get_lower_vertex()->get_coordinate(1) << ", " << edge->get_lower_vertex()->get_coordinate(2) 
			//	<< ") is nearest, returning " << negf_math::vector_norm_3d(vector2) << endl;
			return negf_math::vector_norm_3d(vector2);
		} else if (lambda >= 1.0) { // upper vertex is nearest
			vector3[0] = vertex->get_coordinate(0) - edge->get_upper_vertex()->get_coordinate(0);
			vector3[1] = vertex->get_coordinate(1) - edge->get_upper_vertex()->get_coordinate(1);
			vector3[2] = vertex->get_coordinate(2) - edge->get_upper_vertex()->get_coordinate(2);
			//cout << "    get_distance_edge: upper vertex (" << edge->get_upper_vertex()->get_coordinate(0) 
			//	<< ", " << edge->get_upper_vertex()->get_coordinate(1) << ", " << edge->get_upper_vertex()->get_coordinate(2) 
			//	<< ") is nearest, returning " << negf_math::vector_norm_3d(vector3) << endl;
			return negf_math::vector_norm_3d(vector3);
		} else {    // 0<lambda<1
			// P = P0 + lambda*(P1-P0) + vector4; vector4 is perpendicular to P1-P0
			vector4[0] = vector2[0] - lambda * vector1[0];
			vector4[1] = vector2[1] - lambda * vector1[1];
			vector4[2] = vector2[2] - lambda * vector1[2];
			//cout << "    get_distance_edge: projection is in between (lambda=" << lambda << "), returning " << negf_math::vector_norm_3d(vector4) << endl;
			return negf_math::vector_norm_3d(vector4);
		}
		break;
	case 2:
		vector1[0] = edge->get_upper_vertex()->get_coordinate(0) - edge->get_lower_vertex()->get_coordinate(0);
		vector1[1] = edge->get_upper_vertex()->get_coordinate(1) - edge->get_lower_vertex()->get_coordinate(1);
		vector2[0] = vertex->get_coordinate(0) - edge->get_lower_vertex()->get_coordinate(0);
		vector2[1] = vertex->get_coordinate(1) - edge->get_lower_vertex()->get_coordinate(1);
		NEGF_ASSERT(negf_math::vector_norm_2d(vector1) > constants::min_vector_norm, "edge is too small.");
		lambda = negf_math::vector_scalar_product_2d(vector1,vector2) / negf_math::vector_scalar_product_2d(vector1,vector1);
		if (lambda <= 0.0) { // lower vertex is nearest
			return negf_math::vector_norm_2d(vector2);
		} else if (lambda >= 1.0) { // upper vertex is nearest
			vector3[0] = vertex->get_coordinate(0) - edge->get_upper_vertex()->get_coordinate(0);
			vector3[1] = vertex->get_coordinate(1) - edge->get_upper_vertex()->get_coordinate(1);
			return negf_math::vector_norm_2d(vector3);
		} else {    // 0<lambda<1
			vector4[0] = vector2[0] - lambda * vector1[0];
			vector4[1] = vector2[1] - lambda * vector1[1];
			return negf_math::vector_norm_2d(vector4);
		}
		break;
	case 1:
		vector1[0] = edge->get_upper_vertex()->get_coordinate(0) - edge->get_lower_vertex()->get_coordinate(0);
		vector2[0] = vertex->get_coordinate(0)                   - edge->get_lower_vertex()->get_coordinate(0);
		NEGF_ASSERT(fabs(vector1[0]) > constants::min_vector_norm, "edge is too small.");
		lambda = vector1[0]*vector2[0] / vector1[0]*vector1[0];
		if (lambda <= 0.0) { // lower vertex is nearest
			return fabs(vector2[0]);
		} else if (lambda >= 1.0) { // upper vertex is nearest
			return fabs( vertex->get_coordinate(0) - edge->get_upper_vertex()->get_coordinate(0));
		} else {    // 0<lambda<1
			return fabs(vector2[0] - lambda * vector1[0]);
		}
		break;
	default:
		NEGF_EXCEPTION("Wrong dimensionality.");
		break;
	}
	NEGF_EXCEPTION("One never gets here.");
	return 0.0;
);}


/** Get distance of vertex to face
 *  The vertex can be written as v0+a*(v1-v0)+b*(v2-v0)+c*n
 *  where vi are face vertices, n is the face normal and a,b,c are unknown.
 *  we solve a linear system of equations to get a,b,c
 * depending on a<0, b>1 etc. we choose whats the distance.
 */
double Geometry::get_distance(const Vertex * const vertex, const Face * const face) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(vertex!=0 && face!=0, "faulty input.");
	NEGF_ASSERT(face->get_num_vertices()==3 || face->get_num_vertices()==4, "only triangles or rectangles are implemented.");
	NEGF_ASSERT(vertex->get_dimension()==face->get_vertex(0)->get_dimension(), "dimensionality mismatch.");
	NEGF_ASSERT(vertex->get_dimension()==3, "Faces exist only in 3D, no?");
	/*cout << "get_distance of face to vertex (";
	for (uint ii=0; ii<vertex->get_dimension(); ii++) cout << vertex->get_coordinate(ii) << ",";
	cout << ")" << endl << "   face vertices: " << endl;
	for (uint ii=0; ii<face->get_num_vertices(); ii++) {
		cout << "      " << ii << "=(";
		for (uint jj=0; jj<face->get_vertex(ii)->get_dimension(); jj++) cout << face->get_vertex(ii)->get_coordinate(jj) << ",";
		cout << ")" << endl;
	}*/
	
	double v0[3];
	v0[0] = face->get_vertex(0)->get_coordinate(0);
	v0[1] = face->get_vertex(0)->get_coordinate(1);
	v0[2] = face->get_vertex(0)->get_coordinate(2);
	double vector1[3], vector2[3], normal[3], vector3[3];
	vector1[0] = face->get_vertex(1)->get_coordinate(0) - v0[0];
	vector1[1] = face->get_vertex(1)->get_coordinate(1) - v0[1];
	vector1[2] = face->get_vertex(1)->get_coordinate(2) - v0[2];
	vector2[0] = face->get_vertex(2)->get_coordinate(0) - v0[0];
	vector2[1] = face->get_vertex(2)->get_coordinate(1) - v0[1];
	vector2[2] = face->get_vertex(2)->get_coordinate(2) - v0[2];
	
	negf_math::vector_outer_product(vector1,vector2,normal);
	bool other = false;
	if (face->get_num_vertices()==4 && negf_math::vector_norm_3d(normal)<=constants::min_vector_norm) {
		other = true;
		vector2[0] = face->get_vertex(3)->get_coordinate(0) - v0[0];
		vector2[1] = face->get_vertex(3)->get_coordinate(1) - v0[1];
		vector2[2] = face->get_vertex(3)->get_coordinate(2) - v0[2];
		negf_math::vector_outer_product(vector1,vector2,normal);
	}
	NEGF_ASSERT(negf_math::vector_norm_3d(normal)>constants::min_vector_norm, "too many edges of a face are almost parallel.");
	negf_math::vector_normalize_3d(normal);
	
	vector3[0] = vertex->get_coordinate(0) - v0[0];
	vector3[1] = vertex->get_coordinate(1) - v0[1];
	vector3[2] = vertex->get_coordinate(2) - v0[2];
	
	double a,b,c;
	negf_math::solve_linear_equation(vector1[0], vector2[0], normal[0], vector3[0],
					 				 vector1[1], vector2[1], normal[1], vector3[1],
					 				 vector1[2], vector2[2], normal[2], vector3[2],
					 				 a, b, c);
	//cout << "   v=v0+a*(v1-v0)+b*" << (other ? "(v3-v0)" : "(v2-v0)") << "+c*normal where a=" << a << ", b=" << b << ", c=" << c << endl;
	
	// v0 + a(v1-v0)+b(v2-v0) is now the projection of the vertex onto the face plane.
	// if the projection is outside the face, determine the minimum of its distance to the face edges
	bool projection_inside = false;
	switch(face->get_num_vertices())
	{
	case 3:
		if (a>=0.0 && b>=0.0 && a+b<=1.0)
			projection_inside = true;
		break;
	case 4:
		if (a>=0.0 && b>=0.0 && a<=1.0 && b<=1.0)
			projection_inside = true;
		break;
	default:
		NEGF_EXCEPTION("Implemented only for triangular and rectangular faces.");
		break;
	}
	double dist_projection = 0.0;
	if (!projection_inside) {
		Vertex projection(0, v0[0] + a*vector1[0]+b*vector2[0],
							 v0[1] + a*vector1[1]+b*vector2[1],
							 v0[2] + a*vector1[2]+b*vector2[2]);
		dist_projection = this->get_distance(&projection,face->get_edge(0));
		for (uint ii=1; ii<face->get_num_edges(); ii++) {
			dist_projection = min( dist_projection,this->get_distance(&projection,face->get_edge(ii)) );
		}
		//for (uint ii=0; ii<face->get_num_edges(); ii++) {
		//	cout << "   projection is outside. distance of it to edge " << ii << " is " << this->get_distance(&projection,face->get_edge(ii)) << endl;
		//}
	}
	
	// the result is simply pythagoras
	//cout << "   returning sqrt(" << dist_projection << "^2 + " << c << "^2)=" << sqrt(dist_projection*dist_projection + c*c) << endl; 
	return sqrt(dist_projection*dist_projection + c*c);
);}


/** Get the distance of a vertex to a contact
 *  the routine computes the distance to all faces (edges) the contact contains */
double Geometry::get_distance(const Vertex * const vertex, const Contact * const contact) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(vertex->get_dimension()==this->get_dimension(), "Wrong vertex dimension.");
	NEGF_ASSERT(find(this->contacts.begin(), this->contacts.end(), contact)!=this->contacts.end(),
				"Contact must be contained in geometry object.");
	vector<Vertex *> contact_verts = contact->get_contact_vertices();
	
	// assemble list of elements touching contact
	vector<Element *> adj_elems;
	for (uint ii = 0; ii < contact_verts.size(); ii++)
	{
		vector<Element *> vertex_elems = this->get_elems_near(contact_verts[ii]);
		for (uint jj = 0; jj < vertex_elems.size(); jj++) {
			if (find(adj_elems.begin(), adj_elems.end(), vertex_elems[jj])==adj_elems.end())
				adj_elems.push_back(vertex_elems[jj]);
		}
	}
	
	double distance;
	vector<Face *> contact_faces;
	vector<Edge *> contact_edges;
	switch(this->get_dimension())
	{
	case 3:
		// assemble list of contact faces
		for (uint ii=0; ii<adj_elems.size(); ii++)
		{
			for (uint jj=0; jj<adj_elems[ii]->get_num_faces(); jj++)
			{
				Face * face = adj_elems[ii]->get_face(jj);
				
				// determine for each face how many of its vertices are in the contact
				uint num_face_contactverts = 0;
				for (uint kk=0; kk<face->get_num_vertices(); kk++)
				{
					if (find(contact_verts.begin(), contact_verts.end(), face->get_vertex(kk)) != contact_verts.end())
						num_face_contactverts++;
				}
				// if there are >=3 verts in the contact, the face is a contact face
				if (num_face_contactverts >=3) {
					NEGF_ASSERT(find(contact_faces.begin(),contact_faces.end(),face)==contact_faces.end(),
						"a contact face belonged to more than one element (how is this possible?).");
					contact_faces.push_back(face);
				}
			}
		}
		NEGF_ASSERT(contact_faces.size()!=0, "There was no contact face.");
		
		// determine for each face the distance to the vertex and pick the minimum
		distance = this->get_distance(vertex, contact_faces[0]);
		for (uint ii=1; ii<contact_faces.size(); ii++)
			distance = min(distance, this->get_distance(vertex, contact_faces[ii]));
		
		return distance;
		break;
	case 2:
		// assemble list of contact edges
		for (uint ii=0; ii<adj_elems.size(); ii++)
		{
			for (uint jj=0; jj<adj_elems[ii]->get_num_edges(); jj++)
			{
				// determine for each edge how many of its vertices are in the contact
				uint num_edge_contactverts = 0;
				if (find(contact_verts.begin(), contact_verts.end(), adj_elems[ii]->get_edge(jj)->get_lower_vertex()) != contact_verts.end())
						num_edge_contactverts++;
				if (find(contact_verts.begin(), contact_verts.end(), adj_elems[ii]->get_edge(jj)->get_upper_vertex()) != contact_verts.end())
						num_edge_contactverts++;
				// if there are both verts in the contact, the edge is a contact edge
				if (num_edge_contactverts ==2) {
					NEGF_ASSERT(find(contact_edges.begin(),contact_edges.end(),adj_elems[ii]->get_edge(jj))==contact_edges.end(),
						"a contact edge belonged to more than one element (how is this possible?).");
					contact_edges.push_back(adj_elems[ii]->get_edge(jj));
				}
			}
		}
		NEGF_ASSERT(contact_edges.size()!=0, "There was no contact edge.");
		
		// determine for each edge the distance to the vertex and pick the minimum
		distance = this->get_distance(vertex, contact_edges[0]);
		for (uint ii=1; ii<contact_edges.size(); ii++)
			distance = min(distance, this->get_distance(vertex, contact_edges[ii]));
		
		return distance;
		break;
	case 1:
		// in 1D a contact consists of one point.
		// NOT IN NEGF!!!
		//NEGF_ASSERT(adj_elems.size()==1 && contact_verts.size()==1, "a 1D contact did not match our expectations.");
		distance = this->get_distance(vertex, contact_verts[0]);
		for (uint ii=1; ii<contact_verts.size(); ii++)
			distance = min(distance, this->get_distance(vertex, contact_verts[ii]));
		return distance;
		break;
	default:
		NEGF_EXCEPTION("Faulty dimensionality.");
		return 0.0;
	}
	NEGF_EXCEPTION("This position is never reached."); return 0.0;
);}


void Geometry::face_plane(uint elem_idx, uint elem_face_idx, double fpt[3], double fdir1[3], double fdir2[3]) const throw (Exception *)
{STACK_TRACE(
    switch(this->get_element(elem_idx)->get_type())
	{
	case element_type::cuboid:
	    if(elem_face_idx >= 6 || elem_face_idx<0)	// numbers in accordance with datexelements.C
			NEGF_EXCEPTION("Invalid face index.");
	    break;
	case element_type::prism:
		if(elem_face_idx >= 5 || elem_face_idx<0)
			NEGF_EXCEPTION("Invalid face index.");
        break;
	case element_type::pyramid:
		if(elem_face_idx >= 5 || elem_face_idx<0)
			NEGF_EXCEPTION("Invalid face index.");
        break;
	case element_type::tetrahedron:
		if(elem_face_idx >= 4 || elem_face_idx<0)
			NEGF_EXCEPTION("Invalid face index.");
        break;
	case element_type::tetrabrick:
		if(elem_face_idx >= 7 || elem_face_idx<0)
			NEGF_EXCEPTION("Invalid face index.");
        break;
	default:
		NEGF_EXCEPTION("Unknown element type.");
	}

	Face * face = this->get_element(elem_idx)->get_face(elem_face_idx);
	for(uint ii=0; ii<3; ii++){
		fpt[ii]   = face->get_vertex(0)->get_coordinate(ii);
		fdir1[ii] = face->get_vertex(1)->get_coordinate(ii) - face->get_vertex(0)->get_coordinate(ii);
		fdir2[ii] = face->get_vertex(2)->get_coordinate(ii) - face->get_vertex(0)->get_coordinate(ii);
    }
);}


/** Given a vertex, find the element containing that vertex.
 *  Note that the routine returns the first element it finds which 
 *  contains the vertex up to some numerical error. If the vertex is on an edge/face
 *  or if it is itself a grid vertex, the result is somewhat random. */
Element * Geometry::get_element_containing(Vertex * vertex) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->get_dimension()==vertex->get_dimension(), "Dimension mismatch between grid and vertex.");
	for (uint ii = 0; ii < this->elements.size(); ii++) {
		if (elements[ii]->contains(vertex)) {
			return elements[ii];
		}
	}
	
	// determine minimum distance to some vertex
	double mindist = 1e100;
	Vertex * minvertex = 0;
	for (uint ii = 0; ii < this->get_num_vertices(); ii++) {
		double dist = this->get_distance(this->get_vertex(ii), vertex);
		if (dist < mindist) {
			mindist = dist;
			minvertex = this->get_vertex(ii);
		}
	}
	NEGF_ASSERT(minvertex!=0 && mindist < 1e100, "somethings wrong.");
	
	// determine minimum distance to faces (3D) / edges (2D) around that vertex
	const vector<Element *> & elems_near_minvertex = this->get_elems_near(minvertex);
	double minfacedist = 1e100; Face * minface = 0;
	double minedgedist = 1e100; Edge * minedge = 0;
	for (uint ii = 0; ii < elems_near_minvertex.size(); ii++) {
		switch (this->get_dimension()) {
		case 3:
			for (uint jj = 0; jj < elems_near_minvertex[ii]->get_num_faces(); jj++) {
				double dist = this->get_distance(vertex, elems_near_minvertex[ii]->get_face(jj));
				if (dist < minfacedist) {
					minfacedist = dist;
					minface = elems_near_minvertex[ii]->get_face(jj);
				}
			}
			break;
		case 2:
			for (uint jj = 0; jj < elems_near_minvertex[ii]->get_num_edges(); jj++) {
				double dist = this->get_distance(vertex, elems_near_minvertex[ii]->get_edge(jj));
				if (dist < minedgedist) {
					minedgedist = dist;
					minedge = elems_near_minvertex[ii]->get_edge(jj);
				}
			}
			break;
		case 1:
		default:
			// do nothing;
			break;
		}
	}
	if (this->get_dimension()==3) {
		cout << "face with minimum distance has the following vertices:" << endl;
		for (uint jj=0; jj < minface->get_num_vertices(); jj++) {
			cout << "    vertex " << jj << "(";
			for (int kk=0; kk < minface->get_vertex(jj)->get_dimension()-1; kk++) {
				cout <<  minface->get_vertex(jj)->get_coordinate(kk) << ", ";
			}
			cout << minface->get_vertex(jj)->get_coordinate(minface->get_vertex(jj)->get_dimension()-1) << ")" << endl;
		}
	}
		
	// fallback: if vertex is really close (<1pm) to some element, take that one
	
	Element * close_elem = 0;
	switch(this->get_dimension())
	{
	case 3:
		if (minfacedist < constants::convert_from_SI(units::length, 1e-12)) {
//			if (this->get_elems_near(minface).size()==1) {
//				logmsg->emit(LOG_WARN,"\nget_element_containing failed for a 3D vertex but it was really close to some element so we took that one.");
//				close_elem = this->get_elems_near(minface)[0];
//			} else {
//				logmsg->emit(LOG_WARN,"\na vertex appeared to be really close to an element but the closest face belongs to more than 1 element.");
//		
			logmsg->emit(LOG_WARN,"\nget_element_containing failed for a 3D vertex but it was really close to some element.");
			close_elem = this->get_elems_near(minface)[0];		
		}
		break;
	case 2:
		if (minedgedist < constants::convert_from_SI(units::length, 1e-12)) {
//			if (this->get_elems_near(minedge).size()==1) {
//				logmsg->emit(LOG_WARN,"\nget_element_containing failed for a 2D vertex but it was really close to some element so we took that one.");
//				close_elem = this->get_elems_near(minedge)[0];
//			} else {
//				logmsg->emit(LOG_WARN,"\na vertex appeared to be really close to an element but the closest edge belongs to more than 1 element.");
//			}
			logmsg->emit(LOG_WARN,"\nget_element_containing failed for a 2D vertex but it was really close to some element.");
			close_elem = this->get_elems_near(minedge)[0];		
		}
		break;
	default:
		break;
	}
	if (close_elem != 0) {
		//return close_elem;
		
		vector<double> coeffs;
		close_elem->get_linear_combination_coefficients(vertex, coeffs);
		NEGF_ASSERT(coeffs.size()==close_elem->get_num_vertices(), "something is weird.");
		for (uint ii = 0; ii < coeffs.size(); ii++) {
			Vertex * v = close_elem->get_vertex(ii);
			cout << "vertex " << ii << "(";
			for (int jj = 0; jj < v->get_dimension()-1; jj++) {
				cout << v->get_coordinate(jj) << ",";
			}
			cout << v->get_coordinate(v->get_dimension()-1) << ")";
			cout << " has linear combination coefficient " << coeffs[ii] << endl;
		}
	}
	
	
	// screen output and error
	switch(this->get_dimension())
	{
	case 3:
		logmsg->emit(LOG_ERROR,"\nget_element_containing failed for 3D vertex (%.12e, %.12e, %.12e)",
				vertex->get_coordinate(0), vertex->get_coordinate(1), vertex->get_coordinate(2));
		logmsg->emit(LOG_ERROR,"    vertex is closest (%e) to grid vertex %d (%.12e, %.12e, %.12e)",
						mindist, minvertex->get_index_global(), minvertex->get_coordinate(0), minvertex->get_coordinate(1),
						minvertex->get_coordinate(2));
		logmsg->emit(LOG_ERROR,"    vertex has a minimum distance of %e (or %.3e[m]) to faces near that vertex.", 
						minfacedist, minfacedist / constants::convert_from_SI(units::length, 1.0));
		break;
	case 2:
		logmsg->emit(LOG_ERROR,"\nget_element_containing failed for 2D vertex (%.12e, %.12e)",
				vertex->get_coordinate(0), vertex->get_coordinate(1));
		logmsg->emit(LOG_ERROR,"    vertex is closest (%e) to grid vertex %d (%.12e, %.12e)",
						mindist, minvertex->get_index_global(), minvertex->get_coordinate(0), minvertex->get_coordinate(1));
		logmsg->emit(LOG_ERROR,"    vertex has a minimum distance of %e (or %.3e[m]) to edges near that vertex.", 
						minedgedist, minedgedist / constants::convert_from_SI(units::length, 1.0));
		break;
	case 1:
		logmsg->emit(LOG_ERROR,"\nget_element_containing failed for 1D vertex (%.12e)",vertex->get_coordinate(0));
		logmsg->emit(LOG_ERROR,"    vertex is closest (%e) to grid vertex %d (%.12e)",
						mindist, minvertex->get_index_global(), minvertex->get_coordinate(0));
		break;
	default:
		NEGF_EXCEPTION("Odd dimensionality.");
	}
	NEGF_EXCEPTION("No element was found containing the given vertex.");
	return NULL;
);}


void Geometry::merge_regions(Region * region1, Region * region2, const string & merged_regionname) throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(region1!=NULL && region2!=NULL, "null pointer encountered.");
	NEGF_ASSERT(region1!=region2, "trying to merge identical regions.");
	
	// search for indices in Geometry's region list
	uint ridx1 = 888888;
	uint ridx2 = 888888;
	for (uint ii = 0; ii < this->nregions; ii++)
	{
		if (this->regions[ii]==region1)
			ridx1 = ii;
		if (this->regions[ii]==region2)
			ridx2 = ii;
	}
	NEGF_ASSERT(ridx1!=888888 && ridx2 != 888888, "one of the regions was not found in the geometry.");
	
	// assign new name to region1
	region1->set_name(merged_regionname);
	
	// assign region1 to all elements that have region2 stored
	for (uint ii = 0; ii < this->elements.size(); ii++) {
		if (this->elements[ii]->get_region()==region2){
			this->elements[ii]->set_region(region1);
		}
	}
	
	// delete region2
	NEGF_ASSERT(regions[ridx2]->get_index()==ridx2, "index inconsistency 1.");
	delete this->regions[ridx2]; this->regions[ridx2] = 0;
	
	// move the last region to the place in the list where region2 has been, including assignment of new index
	if (ridx2!=regions.size()-1) {
		NEGF_ASSERT(regions[regions.size()-1]->get_index()==regions.size()-1, "index inconsistency 2.");
		regions[ridx2] = regions[regions.size()-1];
		regions[ridx2]->set_index(ridx2);
	}
	
	// reduce nregions and the region vector
	this->regions.pop_back();
	this->nregions--;
	
	// user must call prepare() again
	this->prepared = false;
);}


void Geometry::round_vertex_coordinates(uint precision) const throw (Exception *)
{STACK_TRACE(
	double prec = negf_math::pow(10.0,-(double)precision);
	logmsg->emit(LOG_INFO,"|-----------------------------------------------------------------------------------|");
	logmsg->emit(LOG_INFO,"|             rounding vertex coordinates to precision %d or %.1e[m]              |",precision,
						prec / constants::convert_from_SI(units::length, 1.0));
	logmsg->emit(LOG_INFO,"|                      DO NOT CALL THIS DURING A SIMULATION!!!                      |");
	logmsg->emit(LOG_INFO,"|-----------------------------------------------------------------------------------|");
	
	logmsg->init_progress_bar(LOG_INFO_L1,"    looping over vertices", this->get_num_vertices());
	uint count = 0;
	uint next = 0;
	for (uint ii = 0; ii < this->get_num_vertices(); ii++)
	{
		Vertex * v = this->get_vertex(ii);
		NEGF_ASSERT(v!=NULL, "null vertex pointer encountered.");
		for (uint dd = 0; dd < v->get_dimension(); dd++) {
			double scaled = v->get_coordinate(dd)/prec;
			double scaled_floor = floor(scaled);	// floor(11.8) = 11,   floor(-12.2) = -13
			double rounded = (scaled-scaled_floor < 0.5) ? scaled_floor : scaled_floor+1.0; 
			// rounded(11.2)=11, rounded(11.8)=12, rounded(-11.2)=-11, rounded(-11.8)=-12
			// don't wanna work with longs because of size limitations...
			v->set_coordinate(dd, rounded*prec);
		}
		for (uint jj = 0; jj < ii; jj++) {
			NEGF_FASSERT(this->get_distance(v,this->get_vertex(jj)) > 1e-12, "vertices %d and %d have become the same.",ii,jj);
		}
		if(next == count++)  next = logmsg->set_progress_bar(count, this->get_num_vertices());
	}
	logmsg->end_progress_bar();
);}

