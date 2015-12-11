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
#include "Contact.h"
using namespace negf;


Contact::Contact(string name_):
	name(name_),
	index(0),
	index_ready(false),
	area(0.0),
	area_ready(false)
{STACK_TRACE(
	
	this->bc_types.clear();
	// ---------------------------------------------------------
	// set here which boundary conditions are known
	// ---------------------------------------------------------
	bc_types.push_back(quantities::potential);			// BC 0 - for potential
	bc_types.push_back(quantities::electron_density);	// BC 1 - for electron density
	bc_types.push_back(quantities::hole_density);		// BC 2 - for hole density
	bc_types.push_back(quantities::fermilevel);			// BC 3 - for fermilevel
	
	this->num_boundary_conditions = bc_types.size();
	
	// ------------------------------------------------------------
	// initialization of variables which are not yet known
	// ------------------------------------------------------------

	// essential contact variables
	this->adjacent_region = NULL;
	this->contact_vertices.clear();
	
	// master and slave contacts
	this->master_contact = NULL;
	this->slave_contacts.clear();
	
	// boundary conditions
	this->boundary_conditions.resize(num_boundary_conditions, bndconds::BC_unknown);
	this->bc_locked.resize(num_boundary_conditions, false);
	this->bc_ready.resize(num_boundary_conditions, false);
	this->bc_num_values.resize(num_boundary_conditions, 0);
	vector<double> empty_vector;
	this->bc_values.resize(num_boundary_conditions, empty_vector);

	// additional stuff
	this->num_bound_states = 0;
	this->bound_state_energies.clear();
);}


/** Set the name of the material to which the contact is adjacent */
void Contact::set_adjacent_region(Region * region_)
{STACK_TRACE(
	NEGF_ASSERT(region_!=NULL, "tried to set materialname without a name.");
	if (this->adjacent_region != NULL && this->adjacent_region != region_) {
		logmsg->emit(LOG_WARN, "Warning: overwrote the contact's adjacent region.");
	}
	this->adjacent_region = region_;
);}


/** Add a vertex to the contact's list of vertices */
void Contact::add_vertex(Vertex * vertex)
{STACK_TRACE(
	NEGF_ASSERT(vertex!=0, "tried to add null pointer to cantact vertex list.");
	for (uint ii = 0; ii < contact_vertices.size(); ii++) {
		NEGF_FASSERT(vertex != contact_vertices[ii], "tried to add vertex twice (contact: %s, vertex: %d (%e)).",
				this->get_name().c_str(), vertex->get_index_global(), vertex->get_coordinate(0));
	}
	contact_vertices.push_back(vertex);
	vertex->set_contact_properties(this, contact_vertices.size()-1);
	this->area_ready = false;
);}


void Contact::set_area(double area_)
{STACK_TRACE(
	NEGF_ASSERT(area_>0.0, "invalid area.");
	this->area = area_;
	this->area_ready = true;	
);}


/** Assign a "master contact" to this contact (makes sense for lowdimensional contacts) */
void Contact::assign_master_contact(Contact * master_contact_)
{STACK_TRACE(
	NEGF_ASSERT(master_contact_!=0, "tried to lock null pointer to contact.");
	NEGF_ASSERT(this->master_contact==NULL, "mastercontact was already assigned.");
	this->master_contact = master_contact_;
	this->master_contact->add_slave_contact(this);
);}


void Contact::add_slave_contact(Contact * contact_)
{STACK_TRACE(
	if (find(slave_contacts.begin(), slave_contacts.end(), contact_) != slave_contacts.end())
		NEGF_EXCEPTION("The contact is already in the list of slave contacts.");

	slave_contacts.push_back(contact_);

	if (contact_->get_master_contact()!=this)
			NEGF_EXCEPTION("The contact seems to belong to a different master contact.");
);}


void Contact::remove_slave_contact(Contact * contact_)
{STACK_TRACE(
	vector<Contact *>::iterator ii = find(slave_contacts.begin(), slave_contacts.end(), contact_);
	if (ii==slave_contacts.end())
		NEGF_EXCEPTION("Slave contact was not found in the list of contacts.");
	slave_contacts.erase(ii);
);}


/** Unlock contact: Take away master contact and unlock any physical quantities */
void Contact::unlock()
{STACK_TRACE(
	NEGF_ASSERT(this->master_contact != NULL, "There was no mastercontact to be deassigned.");
	this->master_contact->remove_slave_contact(this);
	this->master_contact = NULL;
	for (uint ii=0; ii < num_boundary_conditions; ii++)
	{
		if (bc_locked[ii]==true) {
			bc_locked[ii] = false;
			bc_ready[ii]  = false;
		}
	}
);}


/** Lock the boundary condition of a certain quantity to some other contact ("master_contact"), 
 *  meaning that any request of that boundary conditions is redirected to the other contact. */
void Contact::lock_quantity(quantities::PhysicalQuantity type)
{STACK_TRACE(
	NEGF_ASSERT(master_contact!=NULL, "master contact was not yet assigned.");
	uint type_idx = this->get_type_idx(type);
	this->bc_locked[type_idx] = true;
	this->bc_ready[type_idx] = true;
);}


/** Set the boundary condition of a certain physical quantity (e.g. potential->Dirichlet, Neumann etc.) */
void Contact::set_bndcond(quantities::PhysicalQuantity type, bndconds::BndCond cond)
{STACK_TRACE(
	uint type_idx = this->get_type_idx(type);
	
	NEGF_ASSERT(bc_locked[type_idx]==false, "trying to assing a locked boundary condition.");
	NEGF_ASSERT(boundary_conditions[type_idx]==bndconds::BC_unknown
				|| boundary_conditions[type_idx]==cond, 
				"a different boundary condition was already set for this type.");
				
	this->boundary_conditions[type_idx] = cond;
	if (cond==bndconds::BC_Neumann)	// little HACK: in the Neumann case, no values need to be assigned to the BC.
		bc_ready[type_idx] = true;
);}


/** Sets how many values correspond to a certain boundary condition
 *  Needs to be called AFTER set_bndcond */
void Contact::set_bc_num_values(quantities::PhysicalQuantity type, uint num_values)
{STACK_TRACE(
	uint type_idx = this->get_type_idx(type);
	
	NEGF_ASSERT(bc_locked[type_idx]==false, "trying to assing a locked boundary condition.");
	NEGF_ASSERT(bc_num_values[type_idx]==0 || bc_num_values[type_idx]==num_values, "bc_num_values has already been set for this type.");
	NEGF_ASSERT(boundary_conditions[type_idx]!=bndconds::BC_unknown, "Set type of boundary condition first.");
	NEGF_ASSERT(boundary_conditions[type_idx]!=bndconds::BC_Neumann, "Values do not make sense for Neumann BC.");
	
	this->bc_num_values[type_idx] = num_values;
	this->bc_values[type_idx].clear();
	this->bc_values[type_idx].resize(num_values, -888888.888888);
);}


/** Assign a value to the boundary condition (e.g. potential, Dirichlet --> phi=1.23V) 
 *  Needs to be called AFTER set_bndcond and set_bc_num_values */ 
void Contact::set_bnd_value(quantities::PhysicalQuantity type, double value)
{STACK_TRACE(
	uint type_idx = this->get_type_idx(type);
	
	NEGF_ASSERT(bc_num_values[type_idx]==1,	"when calling set_bnd_value without a value index, the BC must have exactly 1 value.");
	this->set_bnd_value(type, 0, value);
);}


/** Assign values to the boundary condition (e.g. potential, Dirichlet --> phi=1.23V) 
 *  Needs to be called AFTER set_bndcond and set_bc_num_values */ 
void Contact::set_bnd_value(quantities::PhysicalQuantity type, uint idx, double value)
{STACK_TRACE(
	uint type_idx = this->get_type_idx(type);
	
	NEGF_ASSERT(bc_values[type_idx].size()==bc_num_values[type_idx], "Inconsistency.");
	NEGF_FASSERT(bc_num_values[type_idx]>idx, "invalid index %d.",idx);
	NEGF_ASSERT(this->bc_locked[type_idx]==false, "trying to assing a value to a locked boundary condition.");
	NEGF_ASSERT(boundary_conditions[type_idx]!=bndconds::BC_Neumann, "Values do not make sense for Neumann BC.");

	this->bc_values[type_idx][idx] = value;
	
	// check if every value has been set.
	bool bc_is_readyyy = true;
	for (uint ii = 0; ii < bc_num_values[type_idx]; ii++) {
		if (bc_values[type_idx][ii]==-888888.888888) {
			bc_is_readyyy = false;
			break;
		}
	}
	bc_ready[type_idx] = bc_is_readyyy;
);}





uint Contact::get_index() const
{STACK_TRACE(
	if (!index_ready) {
		NEGF_EXCEPTION("Index has not been set yet."); 
		return 0;
	} else {
		return index;
	}
);}


double Contact::get_area() const
{STACK_TRACE(
	if (!area_ready) {
		NEGF_EXCEPTION("Area was not calculated yet. If you have added vertices recently, the area must be recalculated."); 
		return 0.0;
	} else {
		return this->area;
	}
);}


Vertex * Contact::get_contact_vertex(uint ii) const
{STACK_TRACE(
	NEGF_ASSERT(ii < this->get_num_contact_vertices(), "Invalid index.");
	return contact_vertices[ii];
);}


/** Check if a certain vertex is within the contact */
bool Contact::is_vertex_here(Vertex * vertex) const
{STACK_TRACE(
	for (uint ii = 0; ii < contact_vertices.size(); ii++)
		if (contact_vertices[ii] == vertex)
			return true;
	return false;
);}


/** Determine whether the boundary conditions of the contact are locked to some other contact */
Contact * Contact::get_master_contact() const
{STACK_TRACE(
	if (this->master_contact==NULL) {
		NEGF_EXCEPTION("No master contact was assigned so far.");
		return NULL;
	} else {
		return this->master_contact;
	}
);}


Contact * Contact::get_slave_contact(uint ii) const
{STACK_TRACE(
	NEGF_ASSERT(ii<slave_contacts.size(), "index out of range.")
	return slave_contacts[ii];
);}


bool Contact::is_lowdimensional() const
{STACK_TRACE(
	if (this->master_contact==NULL) {
		for (uint ii = 0; ii < num_boundary_conditions; ii++)
			NEGF_ASSERT(bc_locked[ii]==false, "Inconsistency between locked and mastercontact.");
		return false;
	} else {
		return true;
	}
);}


bndconds::BndCond Contact::get_bndcond(quantities::PhysicalQuantity type) const
{STACK_TRACE(
	uint type_idx = this->get_type_idx(type);
	if (this->bc_locked[type_idx]) {
		return this->get_master_contact()->get_bndcond(type);
	} else {
		NEGF_FASSERT(bc_ready[type_idx]==true, "Boundary condition has not been set in contact %d(%s).",
				this->get_index(), this->get_name().c_str());
		return this->boundary_conditions[type_idx];
	}
);}


double Contact::get_bc_value(quantities::PhysicalQuantity type) const
{STACK_TRACE(
	uint type_idx = this->get_type_idx(type);

	NEGF_ASSERT(this->bc_ready[type_idx]==true, "Set up the boundary condition first.");
	if (this->bc_locked[type_idx]) {
		return this->get_master_contact()->get_bc_value(type);
	} else {
		NEGF_ASSERT(bc_num_values[type_idx]==1, "when calling this function, the bnd.cond. must have exactly 1 value.");	
		return this->get_bc_value(type, 0);
	}
);}


double Contact::get_bc_value(quantities::PhysicalQuantity type, uint idx) const
{STACK_TRACE(
	uint type_idx = this->get_type_idx(type);

	if (this->bc_locked[type_idx]) {
		return this->get_master_contact()->get_bc_value(type, idx);
	} else {
		NEGF_FASSERT(bc_ready[type_idx]==true, "Boundary condition has not been set in contact %d(%s).",
				this->get_index(), this->get_name().c_str());
		NEGF_ASSERT(bc_num_values[type_idx]==bc_values[type_idx].size(), "Inconsistency.");
		NEGF_ASSERT(bc_num_values[type_idx]>idx, "Invalid index.");
		NEGF_ASSERT(boundary_conditions[type_idx]!=bndconds::BC_unknown, "Unknown BC.");
		NEGF_ASSERT(boundary_conditions[type_idx]!=bndconds::BC_Neumann, "Value does not make sense for Neumann BC.");
		
		return bc_values[type_idx][idx];
	}
);}


/** Checks if a BC has been properly set up */
bool Contact::bc_is_ready(quantities::PhysicalQuantity type) const
{STACK_TRACE(
	uint type_idx = this->get_type_idx(type);
	
	if (bc_locked[type_idx]) {
		return this->get_master_contact()->bc_is_ready(type);
	} else {
		return this->bc_ready[type_idx];
	}
);}


/** Helper function */
uint Contact::get_type_idx(quantities::PhysicalQuantity type) const
{STACK_TRACE(
	for (uint ii = 0; ii < this->num_boundary_conditions; ii++) {
		if (this->bc_types[ii]==type) 
			return ii;
	}
	NEGF_EXCEPTION("Type not found in the list of possible types for boundary conditions.");
)};



void Contact::set_bound_state_energies(const vector<double> & energies)
{STACK_TRACE(
	this->bound_state_energies = energies;
	this->num_bound_states = energies.size();
);}


void Contact::clear_bound_state_energies()
{STACK_TRACE(
	this->bound_state_energies.clear();
	this->num_bound_states = 0;
);}

