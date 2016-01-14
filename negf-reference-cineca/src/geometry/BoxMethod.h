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
#ifndef BOXMETHOD_H_NEGF
#define BOXMETHOD_H_NEGF

#include "all.h"

#include "Geometry.h"

using namespace std;
namespace negf { class BoxMethodTest; }

namespace negf {

	/** Class for computing and storing information related to the box method (finite volume discretization). <BR>
	 *  measure is a double pointer containing the Voronoi volumes around each vertex, partitioned into the surrounding elements. <BR>
	 *  coefficient is a double pointer containing the Voronoi cell surface associated with an edge, partitioned into the surrounding elements. <BR>
	 *  Right now only the (trivial) 1D case is implemented. */
	class BoxMethod
	{

	public:
		BoxMethod(const Geometry * const geom);
		 virtual ~BoxMethod();

		const double * const * 	get_measure() const 		{ return measure; }		//!< get partitioned Voronoi volumes
		const double * const * 	get_coefficient() const 	{ return coefficient; } 	//!< get partitioned Voronoi surfaces
		const double * 			get_node_measure()    const { return node_measure; }	//!< get total Voronoi volume around each Vertex
		double 					get_measure(uint elem_idx, uint local_vert_idx) const;
		double 					get_coefficient(uint elem_idx, uint local_edge_idx) const;	
		double 					get_node_measure(uint vert_idx) const;	
					
		/** Helper function providing for a given vertex the fraction by which he is
		 *  surrounded by a given region. We find this by looking at the surface of the 
		 *  voronoi box around the vertex which located is within the region, divided by 
		 *  the entire voronoi surface */
		double get_fraction_of_surroundment(const Vertex * vertex, const Region * region) const;
		 
		 /** the same thing, just for an element */
		double get_fraction_of_surroundment(const Vertex * vertex, const Element * element) const;
		
	protected: 
		
		bool compute_1d();				//!< for 1D grids
		void compute_node_measures();	//!< aggregate individual partitions around a vertex to obtain measure of entire voronoi box arounhd a vertex
	
		const Geometry * const 	geom;	//!< the Geometry object under investigation
	
		/** measure[Number of Elements][Number of vertices per element]
		  * gives the measure (length, area, or volume) of the Voronoi box in each element at each vertex  */
		double** measure;
	
		/** coefficient[Number of Elements][Number of Edges per element]
		  * gives the coefficient of the box in each element at each edge is the coefficient for the matrix evaluation
		  * (but only for this element and only geometrical stuff !)
		  */
		double** coefficient;
	
		double* node_measure; //!< An array of length num_vertices which stores the measure of the voronoi box associated with each vertex
	};
	
}	// end of namespace 

#endif /*BOXMETHOD_H_NEGF*/
