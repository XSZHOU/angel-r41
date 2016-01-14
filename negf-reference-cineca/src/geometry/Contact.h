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
#ifndef _CONTACT_H_NEGF
#define _CONTACT_H_NEGF

#include "all.h"

#include "Vertex.h"
#include "Region.h"

using namespace std;


namespace negf {

class Vertex;

	/** A class to store everything assigned to contact (its vertices, boundary conditions, ...)  <BR>
	 *  class Vertex has a contact pointer which is set to NULL if the 
	 *     vertex is not at a contact or to any created contact  <BR>
	 *  basic properties are: name, adjacent material, contained vertices
	 *
	 *  vector boundary_conditions stores the type of the physical quantity whose BC are stored <BR>
	 *  vector bc_ready makes sure these conditions have been set before they are accessed <BR>
	 *  vector bc_values stores the values associated with the boundary condition  <BR>
	 * 
	 *  Every BC has the possibility to be "locked" to some other contact, meaning that
	 *  its value is taken from the other contact (think e.g. of potential at lowdimensional contacts).
	 *  vector bc_locked stores this information, and there are storage places for master- and slave contacts
	 */
	class Contact
	{

	public:

		Contact(string name_);		//!< index will be set later on
		~Contact() {}

		// ------------------------------------
		// access functions
		// ------------------------------------
		
		// basic contact properties
		string					 get_name() 				 const { return name; }
		uint					 get_index() 				 const;
		Region *				 get_adjacent_region() 	 	 const { return adjacent_region; }
		double					 get_area()					 const;
		Vertex *				 get_contact_vertex(uint ii) const;
		const vector<Vertex *> & get_contact_vertices() 	 const { return this->contact_vertices; }
		size_t		 			 get_num_contact_vertices()  const { return this->contact_vertices.size(); }
		bool					 is_vertex_here(Vertex * v)  const;
		
		// master and slave contacts
		Contact *   			 get_master_contact() 		 const;
		Contact *   			 get_slave_contact(uint ii)  const;
		uint					 get_num_slave_contacts() 	 const { return slave_contacts.size(); }
		bool					 is_lowdimensional() 		 const;
		
		// boundary conditions
		bndconds::BndCond 		get_bndcond (quantities::PhysicalQuantity type) const;
		double					get_bc_value(quantities::PhysicalQuantity type) const;
		double					get_bc_value(quantities::PhysicalQuantity type, uint idx) const;
		bool					bc_is_ready (quantities::PhysicalQuantity type) const;

		
		// additional stuff
		uint		get_num_bound_states() 		const { return num_bound_states; }
		
		// ------------------------------------
		// setup functions
		// ------------------------------------
		
		// basic contact properties
		void		set_index(uint index_)		  { index = index_; index_ready = true; } //!< NEVER do this during a simulation!
		void		set_adjacent_region(Region * region_);	//!< Set the name of the material to which the contact is adjacent
		void		add_vertex(Vertex * vertex);			//!< Add a vertex to the contact's list of vertices
		void		set_area(double area_);
		
		// assigning master and slave contacts, locking quantities
		void		assign_master_contact(Contact * master_contact_);	//!< also adds "this" to mastercontact list of slave contacts
		void		add_slave_contact(Contact * contact_); //!< Assign a "master contact" to this contact (makes sense for lowdimensional contacts). Also adds "this" to mastercontact list of slave contacts.
		void		remove_slave_contact(Contact * contact_);
		void		unlock();		//!< Unlock contact: Take away master contact and unlock any physical quantities
		void		lock_quantity(quantities::PhysicalQuantity type);
		
		// boundary conditions
		void		set_bndcond(quantities::PhysicalQuantity type_, bndconds::BndCond cond); //!< Set the boundary condition of a certain physical quantity (e.g. potential->Dirichlet, Neumann etc.)
		void		set_bc_num_values(quantities::PhysicalQuantity type, uint num_values);
		void		set_bnd_value(quantities::PhysicalQuantity type, double value);
		void		set_bnd_value(quantities::PhysicalQuantity type, uint idx, double value);
		
		// additional stuff
		void		set_bound_state_energies(const vector<double> & energies);
		void		clear_bound_state_energies();
		

	protected:
		// basic contact properties
		string				name;					//!< contact name
		uint				index;					//!< contact number
		bool 				index_ready;			//!< true if the contact number has been set
		Region *			adjacent_region;		//!< especially the material of that region is needed!
		double				area;					//!< area (3D) / length (2D) / 1 (1D) of the contact
		bool				area_ready;				//!< true if set_area() was called before
		vector<Vertex *> 	contact_vertices;		//!< list of Vertices of the Contact

		// master and slave contacts
		Contact *			master_contact;			//!< if this contact is 1D/2D: the 3D contact which is physically the same
		vector<Contact *> 	slave_contacts;			//!< if this contact is 3D: stores lowD-contacts which are physically the same

		// boundary conditions
		uint				num_boundary_conditions;	 //!< number of assigned boundary conditions
		vector<quantities::PhysicalQuantity> bc_types;   //!< stores for which physical quantities BCs are available
		vector<bndconds::BndCond> 	boundary_conditions; //!< stores the BC (Dirichlet, Neumann, Ohm, ...)
		vector<bool>				bc_locked;			 //!< if this contact is 1D/2D: stores if a BC (e.g. potential) is locked to the BC of the mastercontact
		vector<bool>				bc_ready;			 //!< true if the BC has been set
		vector<uint>				bc_num_values;		 //!< how many values per BC
		vector< vector<double> > 	bc_values;			 //!< stores the values (Dirichlet->value of potential etc.)
		
		// additional stuff
		vector<double>		bound_state_energies;		 //!< stores energy levels of bound states for lowD-contacts
		uint				num_bound_states;			 //!< number of bound states
		
		// helper functions
		uint				get_type_idx(quantities::PhysicalQuantity type) const;	//!< Helper function
	};

} // end of namespace

#endif /*_CONTACT_H_NEGF*/
