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
#ifndef REGION_H_NEGF
#define REGION_H_NEGF

#include "all.h"

#include "Element.h"
#include "PropertyContainer.h"

using namespace std;

namespace negf {

	/** Encapsulates everything associated with a Region, i.e. name, material (name and property container), elements etc. */
	class Region {
	public:
	
		Region(const char* name);		//!< default material name is "unknown_material". index is set later.
		Region(const string& name);		//!< default material name is "unknown_material". index is set later.
		~Region();

		// access functions
		const string& 						get_name() const;									//!< get the Region (not material!) name
		uint								get_index() const;									//!< get the Region index
		const string&						get_material_name() const;							//!< get the material name
		double								get_material_molefraction() const;					//!< get the material molefraction (throws an error when not assigned)
		bool								has_molefraction() const;							//!< checks whether a molefraction was assigned
		/*const*/ PropertyContainer<double>* 	get_material() const;							//!< get the container storing all properties of the material (throws error when it was not assigned)
		bool 								is_oxide() const;									//!< check (by material name) whether region is not a semiconductor)
		const vector<Element *> & 			get_elements() const { return this->elements; }		//!< get a list of all the Region's Elements
		uint 								get_num_dfise_vertices() const { return num_dfise_vertices; }	//!< get the number of DF-ISE vertices in this Region (important for I/O)
		const vector<uint> & 				get_dfise_vertices() const { return dfise_vertices; }	//!< get a list with DF-ISE vertex numbers (important for I/O)

		// setup functions
		void	set_index(uint index_)		 			 { this->index = index_; this->index_ready = true; }	//!< set the Region index
		void	set_material_name(const char* name);											//!< set the material name
		void	set_material_molefraction(double molefraction);									//!< set the material mole fraction
		void	set_material(PropertyContainer<double>* mat);									//!< set the material itself (container w/ all properties)
		void	set_name(const string & name_) { this->name = name_; }							//!< change the Region name
		void 	add_element(Element * elem);													//!< add an Element to the Region's list of Elements
		void 	set_dfise_region_vertex_numbers(const vector<uint> & vertex_indices);			//!< set the list of DF-ISE vertex numbers
	
		// the material group is of course a porperty of the material and hence should be included
		// in the property container, but its functionality is not yet implemented so we do it
		// here for the sake of simplicity!
		int 	get_materialgroup() const { return 1; }										//!< has no functionality at the moment
				
	protected:
		string						name;					//!< region name
		uint						index;					//!< region index
		string						material_name;			//!< name of the region's material
		double						material_molefraction;  //!< mole fraction of the material (if applicable)
		PropertyContainer<double>*	mat;					//!< material parameters
		bool						index_ready;			//!< flag if region index has been set
		
		vector<Element *> 			elements;				//!< elements within that region
		
		uint 						num_dfise_vertices; 	//!< DF-ISE: df_geom->get_RnumVerts(ridx)
		vector<uint>				dfise_vertices; 		//!< DF-ISE: df_geom->get_Rverts(ridx,num_verts)
	};

}

#endif /*REGION_H_NEGF*/
