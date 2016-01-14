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
#ifndef _VERTEX_H_NEGF
#define _VERTEX_H_NEGF

#include "all.h"

//#include "Contact.h"

using namespace std;

namespace negf {

	class Contact;

	/** Vertex class encapsulating everything relevant for vertices (points). <BR>
	 *  Note that there are three kinds of indices: The global index, which is the "normal" index and very frequently used,
	 *  the internal index, which only appeals to internal vertices (see the "location" property), and the external
	 *  index which stores the index read from the file. I think global and external are the same. */
	class Vertex {
	
	public:
		// --------------------------
		// setup functions
		// --------------------------
		Vertex(uint index_global_, double x, double y, double z);	//!< create 3D vertex
		Vertex(uint index_global_, double x, double y);				//!< create 2D vertex
		Vertex(uint index_global_, double x);						//!< create 1D vertex
		~Vertex() {}
		
		void	set_index_global(int index_global_) 	 { this->index_global   = index_global_;   } //!< do not change during simulation!!
		void	set_index_internal(int index_internal_)  { this->index_internal = index_internal_; } //!< do not change during simulation!!
		void	set_index_external(uint index_external_) { this->index_external = index_external_; } //!< do not change during simulation!!
		void	set_contact_properties(Contact * contact_, uint idx); //!< assign a contact and the number within the list of vertices of that contact
		void	set_heterointerface(bool yes_or_no)      { this->heterointerface = yes_or_no; }		//!< bool determining whether vertex is at a heterointerface


		void set_coordinate(int xyz, double value);		//!< do not even think of using the following method except for debug purposes!!
		
		// --------------------------
		// access functions
		// --------------------------
		int			get_index_internal() 		const { return this->index_internal; }	 //!< get internal index (seldomly used)
		uint		get_index_external()		const { return this->index_external; }	 //!< get external index (almost never used)
		uint		get_index_global()  		const { return this->index_global; }	 //!< get global index (this is the normal index access routine one wants to use)
		usint		get_dimension() 			const { return this->dimension; }		 //!< 1D-3D
		double		get_coordinate(usint ii) 	const; 									 //!< ii=0,...
		
		// the property at which contact the vertex sits inverts the hierarchy
		// we implemented it anyway because it is often used and speeds up things
		bool		is_at_contact()				const { if (contact==0) return false; else return true; }//!< checks whether vertex is at a contact
		Contact *	get_contact()				const { return contact; }								 //!< returns NULL if the vertex is not at a contact
		uint		get_contact_vertex_index()	const { return contact_vertex_index; }					 //!< returns 0 if the vertex is not at a contact
		
		bool		is_at_heterointerface()		const { return heterointerface; } 		//!< checks whether vertex is at a heterointerface
				
		// additional functions
		double		get_distance_to(Vertex * other_vertex) const;//!< compute distance to another vertex

	protected:

		usint 			dimension;			//!< dimensionality (1/2/3)
		vector<double> 	coord;				//!< stores coordinates (size 1/2/3)
		int 			index_internal;		//!< index that is 0 by default and -1 for "external" vertices (DF-ISE location 'e' or material Gas)
		uint 			index_external;		//!< external index used for correspondence (e.g. DF-ISE index)
		uint 			index_global;		//!< global index of vertex in gridfile
		Contact * 		contact;			//!< pointer to contact, if the vertex belongs to a contact (otherwise null)
		uint 			contact_vertex_index;//!< index of the vertex within the contact list of vertices
		bool 			heterointerface;	//!< boolean that marks if the vertex is located at a heterointerface, false by default
		bool 			at_quantized_region;//!< boolean if vertex is near a quantized region, false by default
		int 			minority_drain;		//!< boolean if vertex is "at the end of the grid" --> will have special recombination when lowdimensional
	};
	
}

#endif /*_VERTEX_H_NEGF*/
