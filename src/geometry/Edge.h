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
#ifndef _EDGE_H_NEGF
#define _EDGE_H_NEGF

#include "all.h"

#include "Vertex.h"

using namespace std;

namespace negf {
	
	/** Encapsulates everything relevant for an edge, most importantly the two Vertices which act as end points ("lower" and "upper" vertex). <BR>
	 *  There are two indices: The global index, which you want to use, and the external index that was read in from the file. */
	class Edge {
	
	public:
		// ---------------------------
		// setup functions
		// ---------------------------
		Edge(uint index_global_, Vertex * vertex1_, Vertex * vertex_2);	//!< end points must be known at the time of creation
		~Edge() {}

		void		set_voronoi_face(double area)	{ voronoi_face = area; } //!< set the area of the dual Voronoi face

		void		set_index_external(uint index_external_) { this->index_external = index_external_; } //!< not very relevant
		void	  	set_index_global(uint index_global_) { index_global = index_global_; } //!< don't even think about changing this during the simulation

		bool 		verify() const;				//!< perform some routine checks that everything is correct

		// ---------------------------
		// access functions
		// ----------------------------
		uint 		get_lower_vertex_index() const;	//!< get the global index of the first end point
		uint 		get_upper_vertex_index() const;	//!< get the global index of the second end point
		Vertex *	get_lower_vertex()       const;	//!< get the first end point
		Vertex * 	get_upper_vertex()       const; //!< get the second end point
		Vertex * 	get_other_vertex(Vertex * v) const; //!< returns the vertex which is not v
		uint		get_index_global()       const { return index_global; }		//!< get the global index - used very frequently
		uint		get_index_external()     const { return index_external; }	//!< get the external index - not very useful
		double		get_length()             const { return length; }			//!< get the distance between the two end points
		double		get_voronoi_face()		 const { return voronoi_face; }		//!< get the area of the dual Voronoi face

	protected:
		Vertex * 	vertex1;		//!< first end point ("lower" vertex)
		Vertex * 	vertex2;		//!< second end point ("upper" vertex)

		uint 		index_global;	//!< "standard" index
		uint		index_external;	//!< used for correspondence e.g. to DF-ISE
		double		length;			//!< distance between the two end points
		double		voronoi_face;	//!< area of the voronoi face associated to this edge
	};
	
}

#endif /*_EDGE_H_NEGF*/
