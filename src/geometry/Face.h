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
#ifndef FACE_H_NEGF
#define FACE_H_NEGF

#include "all.h"

#include "Vertex.h"
#include "Edge.h"

using namespace std;

namespace negf
{

	/** Encapsulates everything associated with a face, most importantly the corner Vertices and the Edges between them. <BR>
	 *  Note: There is no possibility to overwrite a vertex, edge or face once it has been assigned. */
	class Face
	{

	public:
		Face(uint index_global_, usint nvertex_, usint nedge_);	//!< number of edges and vertices must be consistent.
		~Face() {}

		// ----------------------------------------------
		// setup functions	
		// note: there is no possibility to overwrite a vertex, edge or face once it has been assigned	
		// ----------------------------------------------
		void      			add_vertex(Vertex * new_vertex);	//!< automatically sets the local index. throws an exception if the list of vertices is already full.
		void	  			add_edge  (Edge   * new_edge);		//!< automatically sets the local index. throws an exception if the list of edges is already full.
		void				set_index_external(uint index_external_) { this->index_external = index_external_; }	//!< relevant for DF-ISE stuff only.
		void	  			set_index_global(uint index_global_) { index_global = index_global_; } //!< don't even think about changing this during the simulation
		bool     			verify() const;						//!< perform some routine checks that everything is correct
		void				prepare();							//!< at the moment only computes the area

		// ----------------------------------------------
		// access functions		
		// ----------------------------------------------
		inline usint 		get_num_vertices() 		const { return nvertex; }		//!< get number of Vertices.
		inline usint 		get_num_edges()    		const { return nedge; }			//!< get number of Edges.
		face_type::FaceType get_type()         		const { return type; }			//!< probably triangle or rectangle.
		double				get_area()				const;							//!< get the area. can only be called after prepare().
		bool				is_ready()				const { return this->ready; }	//!< checks whether prepare() was already called.

		void				get_vertex_indices(long index_array[]) const;			//!< yields an array containing the global indices of the Face's Vertices.
		void				get_edge_indices(long index_array[]) const;				//!< yields an array containing the global indices of the Face's Edges.
		Vertex*   			get_vertex(usint lid) 	const;							//!< get Vertex with local index lid (=0,1,2 or maybe 3).
		Edge*   			get_edge(usint lid) 	const;							//!< get Edge with local index lid (=0,1,2 or maybe 3).
		vector<Vertex *>	get_vertices()     		const { return vertices; }		//!< get the Vertex list
		vector<Edge *>		get_edges()		   		const { return edges; }			//!< get the Edge list
		uint				get_index_global() 		const { return index_global; }	//!< very frequently used function returning the global index
		uint				get_index_external() 	const { return index_external; }//!< very seldomly used function returning the external index

	protected:
		usint				nvertex;		//!< number of Vertices.
		usint				nedge;			//!< number of Edges.
	
		vector<Vertex *>	vertices;		//!< list of Vertices.
		vector<Edge *>		edges;			//!< list of Edges.

		face_type::FaceType	type;			//!< triangle, rectangle etc.

		uint 				index_global;	//!< "normal" index
		uint				index_external;	//!< for correspondence to e.g. DF-ISE

		double				area;			//!< Face area, computed by prepare()

		bool				ready;			//!< false by default, set to true by prepare()
	};

}

#endif /*FACE_H_NEGF*/
