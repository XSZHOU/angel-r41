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
#ifndef _GEOMETRY_H_NEGF
#define _GEOMETRY_H_NEGF

#include "all.h"

#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "Element.h"
#include "Contact.h"
#include "Region.h"

using namespace std;

namespace negf {

	/** Top level container for everything related to geometrical issues.
	 *  A Geometry is made out of Regions, Contacts, Elements, Faces, Edges and Vertices. */
	class Geometry {
		
	public:
			
		Geometry(uint nvertices_, uint nedges_, uint nfaces_, uint nelements_) throw (Exception *);
		~Geometry();
		
		// ---------------------
		// setup functions
		// ---------------------

		void 		set_dimension(uint dim) 			{ this->dimension = dim; }				 //!< set dimension
		void		set_num_dfise_elems(uint num_elems) { this->num_dfise_elems = num_elems; } 	 //!< ONLY for mapping to DF-ISE!
		void		set_num_dfise_regions(uint num_regs) { this->num_dfise_regions = num_regs; } //!< ONLY for mapping to DF-ISE!
	
		uint 		add_vertex(Vertex *);	//!< does NOT alter nvertices.
		uint 		add_edge(Edge *);		//!< does NOT alter nedges.
		uint 		add_face(Face *);		//!< does NOT alter nfaces.
		uint 		add_element(Element *);	//!< does NOT alter nelems.
		uint 		add_region(Region *);	//!< does NOT alter nregions.
		uint		add_contact(Contact *);	//!< does NOT alter ncontacts.
		
		bool 		verify() const throw (Exception *);
		void 		prepare() throw (Exception *);
		
		/** Merge region2 into region1. the properties of the new region are inherited from region1. The boolean "ready" is set to false in the end. */
		void		merge_regions(Region * region1, Region * region2, const string & merged_regionname) throw (Exception *);
		
		void 		round_vertex_coordinates(uint precision) const throw (Exception *);
		
		// -------------------------
		// access functions
		// -------------------------
		
		uint 		get_dimension() 		   		const { return this->dimension; }	//!< get dimension
		uint		get_max_connectivity() 	   		const { return this->max_connectivity; } //!< get the maximum number of edges connecting to a Vertex
		
		uint 		get_num_internal_vertices() 	const { return this->num_internal_vertices; }	//!< get number of internal vertices (see "location" property)
		uint 		get_global_vertex_index(uint internal_index) const;
		uint 		get_internal_vertex_index(uint global_index) const;
		
		Vertex * 	get_vertex (uint vert_idx) 		const throw (Exception *);		//!< throws Exception when vert_idx is out of range or when NULL is encountered.
		Edge *		get_edge   (uint edge_idx) 		const throw (Exception *);		//!< throws Exception when edge_idx is out of range or when NULL is encountered.
		Face *      get_face   (uint face_idx) 		const throw (Exception *);		//!< throws Exception when face_idx is out of range or when NULL is encountered.
		Element *   get_element(uint elem_idx) 		const throw (Exception *);		//!< throws Exception when elem_idx is out of range or when NULL is encountered.
		Region *    get_region (uint region_idx)  	const throw (Exception *);		//!< throws Exception when region_idx is out of range or when NULL is encountered.
		Region *    get_region (string name) 	  	const throw (Exception *);		//!< throws Exception when name is not found or when NULL is encountered.
		Contact *	get_contact(uint contact_idx) 	const throw (Exception *);		//!< throws Exception when contact_idx is out of range or when NULL is encountered.
		Contact *	get_contact(string name) 	  	const throw (Exception *);		//!< throws Exception when name is not found or when NULL is encountered.
		
		uint 		get_num_vertices() 				const { return this->nvertices; }		//!< returns number of vertices, regardless whether they are known
		uint 		get_num_edges() 				const { return this->nedges; }			//!< returns number of edges, regardless whether they are known
		uint		get_num_elements() 				const { return this->nelems; }			//!< returns number of elements, regardless whether they are known
		uint		get_num_dfise_elems() 			const { return this->num_dfise_elems; } //!< ONLY for mapping to DF-ISE!
		uint 		get_num_faces() 				const { return this->nfaces; }			//!< returns number of faces, regardless whether they are known
		uint 		get_num_interfaces() 			const { return this->ninterfaces; }		//!< returns number of interfaces (not used)
		uint		get_num_regions() 				const { return this->nregions; }		//!< returns number of regions
		uint		get_num_dfise_regions() 		const { return this->num_dfise_regions; } //!< ONLY for mapping to DF-ISE!
		uint		get_num_contacts() 				const { return this->ncontacts; }		//!< returns number of contacts
		
		bool		is_prepared() 					const { return this->prepared; }		//!< checks whether prepare() was called before.
		Edge *		find_edge(Vertex * v1, Vertex * v2) const throw (Exception *);			//!< tries to find an Edge connecting the given Vertices, throws an error otherwise.
		double		get_area_from_edges(Edge * edge1, Edge * edge2) const throw (Exception *);
		double		get_distance(const Vertex * const v1,     const Vertex * const v2)       const throw (Exception *);//!< get distance between two Vertices
		double		get_distance(const Vertex * const vertex, const Contact * const contact) const throw (Exception *);
		double		get_distance(const Vertex * const vertex, const Face    * const face)    const throw (Exception *);
		double		get_distance(const Vertex * const vertex, const Edge    * const edge)    const throw (Exception *);
		void		face_plane(uint elem_idx, uint elem_face_idx, double fpt[3], double fdir1[3], double fdir2[3]) const throw (Exception *);
		Element *	get_element_containing(Vertex * vertex) const throw (Exception *);

		
		// "secondary" functions which invert the hierarchy geometry-region-element-face-edge-vertex
		uint 			 	     get_num_edges_near(const Vertex * vertex) const { return edges_near_vertex[vertex->get_index_global()].size(); }	//!< get number of edges connected to a Vertex
		const vector<Edge*> &    get_edges_near    (const Vertex * vertex) const { return edges_near_vertex[vertex->get_index_global()]; }			//!< get edges connected to a Vertex
		uint 			 	     get_num_faces_near(const Vertex * vertex) const { return faces_near_vertex[vertex->get_index_global()].size(); }	//!< get number of faces having a certain Vertex as corner
		const vector<Face*> &    get_faces_near    (const Vertex * vertex) const { return faces_near_vertex[vertex->get_index_global()]; }			//!< get faces having a certain Vertex as corner
		uint 		 		     get_num_elems_near(const Edge * edge)     const { return elems_near_edge[edge->get_index_global()].size(); }		//!< get number of elements around an Edge
		const vector<Element*> & get_elems_near    (const Edge * edge)     const { return elems_near_edge[edge->get_index_global()]; }				//!< get all elements around an Edge
		uint 		 		     get_num_elems_near(const Face * face)     const { return elems_near_face[face->get_index_global()].size(); }		//!< get number of elements having a certain Face as face
		const vector<Element*> & get_elems_near    (const Face * face)     const { return elems_near_face[face->get_index_global()]; }				//!< get all elements having a certain Face as face
		uint 			 	     get_num_elems_near(const Vertex * vertex) const { return elems_near_vertex[vertex->get_index_global()].size(); }	//!< get number of elements having a certain vertex as corner
		const vector<Element*> & get_elems_near    (const Vertex * vertex) const { return elems_near_vertex[vertex->get_index_global()]; }			//!< get elements having a certain vertex as corner
		uint 			 	     get_num_regions_near(const Vertex * vertex) const { return regions_near_vertex[vertex->get_index_global()].size(); } 	//!< get number of Regions around a Vertex
		const vector<Region*> &  get_regions_near    (const Vertex * vertex) const { return regions_near_vertex[vertex->get_index_global()]; }		//!< get Regions around a Vertex
		
	protected:	
	
		bool prepared;
	
		vector<Vertex*>   vertices;			//!< stores Vertices
		vector<Edge*> 	  edges;			//!< stores Edges
		vector<Face*> 	  faces;			//!< stores Faces
		vector<Element*>  elements;			//!< stores Elements
		vector<Region*>   regions;			//!< stores Regions
		vector<Contact *> contacts;			//!< stores Contacts
				 
		uint nvertices;				//!< number of Vertices
		uint nedges;				//!< number of Edges
		uint nfaces;				//!< number of Faces
		uint nelems;			//!< number of Elements
		uint ninterfaces;		//!< number of Interfaces (not used)
		uint nregions;			//!< number of Regions
		uint ncontacts;			//!< number of Contacts

		uint dimension;			//!< 1, 2 or 3
		uint max_connectivity;	//!< maximum Vertex connectivity
		uint num_dfise_elems; //!< ONLY for mapping to DF-ISE!
		uint num_dfise_regions; //!< ONLY for mapping to DF-ISE!
		uint num_internal_vertices;
		
		vector<uint> internal_to_global_vertex_indices;
		
		vector< vector<Edge *> >    edges_near_vertex;	 //!< vector of Edge-vectors
		vector< vector<Face *> > 	faces_near_vertex;	 //!< vector of Face-vectors
		vector< vector<Element *> > elems_near_edge; 	 //!< vector of Element-vectors
		vector< vector<Element *> > elems_near_face; 	 //!< vector of Element-vectors
		vector< vector<Element *> > elems_near_vertex; 	 //!< vector of Element-vectors
		vector< vector<Region *> >  regions_near_vertex; //!< vector of Region-vectors
		
	};

} // end of namespace

#endif /*GEOMETRY_H_NEGF*/

