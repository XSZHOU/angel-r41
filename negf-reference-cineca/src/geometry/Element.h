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
#ifndef ELEMENT_H_NEGF
#define ELEMENT_H_NEGF

#include "all.h"

#include "Vertex.h"
#include "Edge.h"
#include "Face.h"

using namespace std;

namespace negf {
	
	class Region;

	// note: distinguish between the vertex index of this element (0,1,2,3 for example)
	//       and the global vertex index (4259, 876, 320, 23091 for example)
	//       and the global element index

	/** Encapsulates everything associated with an element, most importantly its verices, edges and faces. <BR>
	 *  Note: there is no possibility to overwrite a vertex, edge or face once it has been assigned. */
	class Element 
	{

	public:

		Element(uint index_global_, element_type::ElementType element_type_);	//!< number of vertices, edges and faces is determined by the Element type
		~Element();																//!< deallocates inv_jacobian

		// ----------------------------------------------
		// setup functions
		// note: there is no possibility to overwrite a vertex, edge or face once it has been assigned
		// ----------------------------------------------
		void      			add_vertex(Vertex * new_vertex);	//!< automatically sets the local index. throws an exception if the list of vertices is already full.
		void	  			add_edge  (Edge   * new_edge);		//!< automatically sets the local index. throws an exception if the list of vertices is already full.
		void	  			add_face  (Face   * new_face);		//!< automatically sets the local index. throws an exception if the list of vertices is already full.
		void      			set_region(Region * region);		//!< only one Region per element. Gives no warnings when overwriting exisitng Region.
		void	  			set_index_external(uint index_external_) { index_external = index_external_; }	//!< relevant only for DF-ISE stuff
		void	  			set_index_global  (uint index_global_)   { index_global   = index_global_; } 	//!< don't even think about changing this during the simulation

		void      			prepare();		//!< computes the volume, inv_jacobian and the diameter of the element
		bool 	  			verify() const;

		// ----------------------------------------------
		// access functions		
		// ----------------------------------------------
		inline usint 		get_num_vertices() 		const { return num_verts; }	//!< get number of corners constituting the element
		inline usint 		get_num_edges()    		const { return num_edges; }	//!< get number of edges of the element
		inline usint 		get_num_faces()    		const { return num_faces; }	//!< get number of element faces (0 except in 3D)
		double				get_volume()       		const;						//!< get volume. prepare() needs to be called before.
		double				get_diameter()			const;						//!< get diameter. prepare() needs to be called before.
		element_type::ElementType get_type()       	const { return type; }		//!< get type (interval, triangle, tetrahedron, ...)
		Region* 			get_region()       		const;						//!< get the Region that the element lies in (kinda inverts the hierarchy, but very useful)
		void				get_vertex_indices(uint index_array[]) const;		//!< get an array with the global indices of the Element's vertices
		void				get_edge_indices(uint index_array[]) const;			//!< get an array with the global indices of the Element's edges
		void				get_face_indices(uint index_array[]) const;			//!< get an array with the global indices of the Element's faces
		Vertex* 			get_vertex(uint lid) 	const;						//!< get the Vertex with local index lid (=0, 1, maybe 2, maybe 3)
		Edge*   			get_edge(uint lid) 		const;						//!< get the Edge with local index lid (=0, 1, ...)
		Face*   			get_face(uint lid) 		const;						//!< get the Face with local index lid (=0, 1, ...)
		const vector<Vertex *>&	get_vertices() 		const { return vertices; }	//!< get a vector contaning the Element's Vertices
		const vector<Edge *>& 	get_edges()	   		const { return edges; }		//!< get a vector contaning the Element's Edges
		const vector<Face *>& 	get_faces()      	const { return faces; }		//!< get a vector contaning the Element's Faces
		uint				get_local_index(const Vertex * vertex) const;
		uint				get_local_index(const Edge * edge) const;
		uint				get_index_global() 		const { return index_global; }	 //!< get the global index of the Element (very frequently used)
		uint				get_index_external() 	const { return index_external; } //!< get the external index of the Element (used only for input/output stuff)
		const double *		get_invjacobian() 		const;							//!< get an array containing coefficients of the inverse Jacobian
		void				get_linear_combination_coefficients(Vertex * vertex, vector<double> & coeffs) const;

		// ----------------------------------------
		// info functions
		// ----------------------------------------
		bool				is_ready() 				const { return this->ready; }	//!< checks whether prepare() was called
		bool				contains(Vertex * vertex) const;						//!< checks whether the given Vertex is located inside the Element

		// functions used solely for Box Method
		bool				is_semiconductor() 		const { return true; } 		//!< solely used for Box Method. always true
		
	protected:
		void      			compute_volume();			//!< compute Element volume
		void	  			compute_invjacobian();
		void	  			compute_diameter();			//!< compute Element diameter, i.e. maximum distance of 2 points inside Element
		
		
		usint				num_verts;					//!< number of Vertices, set in constructor
		usint				num_edges;					//!< number of Edges, set in constructor
		usint 				num_faces;					//!< number of Faces, set in constructor

		vector<Vertex*>	 	vertices;					//!< contains the Vertices, created by consecutive add_vertex() calls
		vector<Edge*>		edges;						//!< contains the Edges, created by consecutive add_edge() calls
		vector<Face*>		faces;						//!< contains the Faces, created by consecutive add_face() calls
		Region *			region;						//!< set by set_region(); default NULL

		element_type::ElementType type;					//!< interval, triangle, tetrahedron, ...

		uint 				index_global;				//!< "normal" index
		uint				index_external; 			//!< for mapping element indices here (-> e.g.DF-ISE)

		bool         		ready; 						//!< element is ready (prepared)
		double       		element_volume;				//!< volume of the element
		double				diameter;					//!< element diameter (maximum distance of 2 vertices)

		double *			inv_jacobian; 				//!< inverse of Jacobian of affine transformation reference element-->real element
	};

}

#endif /*ELEMENT_H_NEGF*/
