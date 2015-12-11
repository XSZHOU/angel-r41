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
#include "PoissonBM.h"
using namespace negf;


PoissonBM::PoissonBM(const Geometry * const grid_, const BoxMethod * const box_method_):
	Poisson(grid_),
	box_method(box_method_)
{STACK_TRACE(
	NEGF_ASSERT(grid_!=0 && box_method_!=0, "tried to create Poisson equation with null pointer(s).");

	// compute the scaling factor
	this->scaling_factor = 0.0;
	double epsilon_guess = constants::convert_from_SI(units::dielectric, 10.0 * constants::SIeps0);
	const double * const * coefficient = box_method->get_coefficient();
	uint counter = 0;
	for (uint elem_idx = 0; elem_idx < grid->get_num_elements(); elem_idx++) {
		for (uint ii = 0; ii < grid->get_element(elem_idx)->get_num_edges(); ii++) {
			this->scaling_factor += epsilon_guess * coefficient[elem_idx][ii];
			counter++;
		}
	}
	this->scaling_factor /= counter;
);}


double PoissonBM::get_scaling_factor() const 
{STACK_TRACE(
	if (scaling_factor==0)
		NEGF_EXCEPTION("The scaling factor was not set up yet.");
	return scaling_factor;
);}


/** Return the Newton function (whose roots are searched)
 * @param line the variable. line+offset = line number in Jacobian */
double PoissonBM::get_newton_function(uint line) const
{STACK_TRACE(
	const Vertex * 		   vertex       = grid->get_vertex(line);
	const double * const * coefficient  = box_method->get_coefficient();
	const double *         node_measure = box_method->get_node_measure();
	const double ec = constants::convert_from_SI(units::charge, constants::SIec);

	NEGF_ASSERT(vertex!=0 && epsilon!=0 && edensity!=0 && hdensity!=0 && doping!=0
				&& coefficient!=0 && node_measure!=0 && scaling_factor!=0, "Some quantity is missing.");

	// treatment of Dirichlet contact vertices (Poisson eqn is replaced with explicit value)
	if (vertex->is_at_contact() && vertex->get_contact()->get_bndcond(quantities::potential)==bndconds::BC_Dirichlet) 
	{
		return scaling_factor * (
				  this->get_value(line) 
				- vertex->get_contact()->get_bc_value(quantities::potential, vertex->get_contact_vertex_index()));
		// vertex->get_contact_vertex_index() returns the index that the vertex has in the list of the contact's vertices
	}
	
	// NEGF: we also have to do something in case of contacts and Neumann conditions
	if (vertex->is_at_contact() && vertex->get_contact()->get_bndcond(quantities::potential)==bndconds::BC_Neumann) 
	{
		const vector<Edge *> & near_edges = grid->get_edges_near(vertex);
		bool a_device_vertex_is_connected = false;
		Vertex * device_interior_vertex = 0;
		for (uint ii=0; ii < near_edges.size(); ii++) 
		{
			Vertex * v = (near_edges[ii]->get_lower_vertex()==vertex) 
				? near_edges[ii]->get_upper_vertex() : near_edges[ii]->get_lower_vertex();
			if (!v->is_at_contact()) {
				a_device_vertex_is_connected = true;
				device_interior_vertex = v;
				break;
			}
		}
		
		if (!a_device_vertex_is_connected) 
		{
			// make sure that a "contact-interior" vertex gets the same value as a 
			// vertex directly adjacent to the actual device
			const vector<Vertex *> & contact_verts = vertex->get_contact()->get_contact_vertices();
			Vertex * v_adjacent_to_device = 0;
			for (uint ii=0; ii < contact_verts.size(); ii++) {
				const vector<Edge *> & edges_near_ii = grid->get_edges_near(contact_verts[ii]);
				for (uint jj=0; jj < edges_near_ii.size(); jj++) {
					Vertex * v = (edges_near_ii[jj]->get_lower_vertex()==contact_verts[ii]) 
						? edges_near_ii[jj]->get_upper_vertex() : edges_near_ii[jj]->get_lower_vertex();
					if (!v->is_at_contact()) {
						v_adjacent_to_device = contact_verts[ii];
						break;
					}
				}
				if (v_adjacent_to_device!=0) break;
			}
			NEGF_ASSERT(v_adjacent_to_device!=0, "no vertex was found which is adjacent to the device interior !?");
			
			// also in case of nonzero neumann field, we absolutely need flat potenetial WITHIN the contacts
			// because potential is added to hamiltonian which enters boundary self-energies, and potential MUST be flat there!
			return scaling_factor * (this->get_value(line) - this->get_value(v_adjacent_to_device->get_index_global()));
		} else {
			// vertex which is at contact but directly connected to the device interior: 
			// continue with normal treatment, normal formula
			
			// NO! same value as device-internal vertex!
			// desired value = value at device interface - distance * E-field
			if (fabs(this->neumann_field) > 1e-10) NEGF_ASSERT(grid->get_dimension()==1, "Neumann field only supported for d=1!");
			double dist = vertex->get_coordinate(0) - device_interior_vertex->get_coordinate(0);
			return scaling_factor * (this->get_value(line) - this->get_value(device_interior_vertex->get_index_global()) - dist * this->neumann_field);
		}
	}

	double result = 0.0;

	vector<Edge *> edges = grid->get_edges_near(vertex);
	for (uint ii = 0; ii < edges.size(); ii++)
	{
		uint other_vertex = (edges[ii]->get_lower_vertex()==vertex)
							? edges[ii]->get_upper_vertex()->get_index_global()
							: edges[ii]->get_lower_vertex()->get_index_global();
		vector<Element *> elems_near_edge = grid->get_elems_near(edges[ii]);
		for (uint jj = 0; jj < elems_near_edge.size(); jj++)
		{
			uint elem_idx = elems_near_edge[jj]->get_index_global();
			uint local_edge_idx = elems_near_edge[jj]->get_local_index(edges[ii]);

			result += coefficient[elem_idx][local_edge_idx] 
					* epsilon->get_value(elem_idx) 
					* ( this->get_value(line) - this->get_value(other_vertex) );
		}
	}

	result -= node_measure[line] * ec * (hdensity->get_value(line) - edensity->get_value(line) + doping->get_value(line));

	for (uint nu = 4; nu < dependencies.size(); ++nu)
	{
		int charge_sign;
		if (dependencies[nu]->get_type() == quantities::electron_density)
			charge_sign = -1;
		else if (dependencies[nu]->get_type() == quantities::hole_density)
			charge_sign = +1;
		else
			NEGF_EXCEPTION("Poisson dependency not supported.");
		
		if (vertex->is_at_contact() && vertex->get_contact()->get_bndcond(quantities::potential)==bndconds::BC_Neumann) {
			// treatment of Neumann contact vertices: no space charge at contact!
		} else {
			result -= node_measure[line] * ec * charge_sign * dependencies[nu]->get_value(line);
		}
	}
	
	// check: total number of particles from dependencies
	if (line==this->number_of_variables-1)
	{
		for (uint nu = 4; nu < dependencies.size(); ++nu)
		{
			double N = 0.0;
			for (uint ii=0; ii < this->number_of_variables; ii++)
				N += node_measure[ii] * dependencies[nu]->get_value(ii);
			logmsg->emit(LOG_INFO_L2, "%s: %7.3e carriers in total (timestamp %d)", 
					dependencies[nu]->get_name().c_str(), N, dependencies[nu]->get_timestamp());
		}
	}

	return result;
);}


/** Get the partial derivative of the Newton function w.r.t. a certain dependency
 * @param line which variable
 * @param eqn the dependency w.r.t which the derivative is computed
 * @param nonzeros number storing the number of array entries used
 * @param indices array storing the indices of the dependency equation where the derivative is nonzero
 * @param coeff array storing the corresponding derivatives */
void PoissonBM::get_newton_direct_derivatives(uint line, const Equation * eqn, 
								uint & nonzeros, uint indices[], double coeffs[]) const
{STACK_TRACE(
	NEGF_FASSERT(line<grid->get_num_vertices(), "invalid vertex index (%d). this->get_num_variables()=%d",line,this->get_num_variables());
	const Vertex *         vertex   	= grid->get_vertex(line);
	const double * const * coefficient  = box_method->get_coefficient();
	const double *         node_measure = box_method->get_node_measure();
	const double ec = constants::convert_from_SI(units::charge, constants::SIec);

	NEGF_ASSERT(vertex!=0 && epsilon!=0 && edensity!=0 && hdensity!=0 && doping!=0 
				&& coefficient!=0 && node_measure!=0 && scaling_factor!=0, "Some quantity is missing.");

	NEGF_ASSERT(eqn!=NULL, "the dependency requested was NULL!");

	nonzeros = 0;

	// treatment of Dirichlet contact vertices (Poisson eqn is replaced with explicit value)
	if (vertex->is_at_contact() && vertex->get_contact()->get_bndcond(quantities::potential)==bndconds::BC_Dirichlet)
	{
		if (eqn==this) {
			nonzeros   = 1;
			indices[0] = line;
			if (coeffs!=NULL)
				coeffs[0]  = scaling_factor;
			return;
		}
		for (uint nu = 0; nu < dependencies.size(); nu++) {
			if (eqn==dependencies[nu]) {
				nonzeros = 0;
				return;
			}
		}
		NEGF_EXCEPTION("Unknown equation.");
		return;
	}
	
	// NEGF: we also have to do something in case of contacts and Neumann conditions
	if (vertex->is_at_contact() && vertex->get_contact()->get_bndcond(quantities::potential)==bndconds::BC_Neumann) 
	{
		const vector<Edge *> & near_edges = grid->get_edges_near(vertex);
		bool a_device_vertex_is_connected = false;
		Vertex * device_interior_vertex = 0;
		for (uint ii=0; ii < near_edges.size(); ii++) 
		{
			Vertex * v = (near_edges[ii]->get_lower_vertex()==vertex) 
				? near_edges[ii]->get_upper_vertex() : near_edges[ii]->get_lower_vertex();
			if (!v->is_at_contact()) {
				a_device_vertex_is_connected = true;
				device_interior_vertex = v;
				break;
			}
		}
		
		if (!a_device_vertex_is_connected) 
		{
			// make sure that a "contact-interior" vertex gets the same value as a 
			// vertex directly adjacent to the actual device
			const vector<Vertex *> & contact_verts = vertex->get_contact()->get_contact_vertices();
			Vertex * v_adjacent_to_device = 0;
			for (uint ii=0; ii < contact_verts.size(); ii++) 
			{
				const vector<Edge *> & edges_near_ii = grid->get_edges_near(contact_verts[ii]);
				for (uint jj=0; jj < edges_near_ii.size(); jj++) {
					Vertex * v = (edges_near_ii[jj]->get_lower_vertex()==contact_verts[ii]) 
						? edges_near_ii[jj]->get_upper_vertex() : edges_near_ii[jj]->get_lower_vertex();
					if (!v->is_at_contact()) {
						v_adjacent_to_device = contact_verts[ii];
						break;
					}
				}
				if (v_adjacent_to_device!=0) break;
			}
			NEGF_ASSERT(v_adjacent_to_device!=0, "no vertex was found which is adjacent to the device interior !?");
			
			if (eqn==this) {
				nonzeros   = 2;
				indices[0] = line;
				indices[1] = v_adjacent_to_device->get_index_global();
				if (coeffs!=NULL) {
					coeffs[0]  = scaling_factor;
					coeffs[1]  = -scaling_factor;
				}
				return;
			} else {
				nonzeros = 0;
				return;
			}
		} else {
			// vertex which is at contact but directly connected to the device interior: 
			// continue with normal treatment
			// NO!!! same value as device-interior!
			if (eqn==this) {
				nonzeros   = 2;
				indices[0] = line;
				indices[1] = device_interior_vertex->get_index_global();
				if (coeffs!=NULL) {
					coeffs[0]  = scaling_factor;
					coeffs[1]  = -scaling_factor;
				}
				return;
			} else {
				nonzeros = 0;
				return;
			}
		}
	}
	
	
	if (eqn==this)  // dependence of Poisson equation on potential
	{
		indices[0] = line;
		if (coeffs!=NULL)
			coeffs[0] = 0;
		nonzeros++;

		vector<Edge *> edges = grid->get_edges_near(vertex);
		for (uint ii=0; ii < grid->get_num_edges_near(vertex); ii++)
		{
			int other_vert_idx = -1;
			if (edges[ii]->get_lower_vertex()==vertex)
				other_vert_idx =  edges[ii]->get_upper_vertex()->get_index_global();
			if (edges[ii]->get_upper_vertex()==vertex)
				other_vert_idx =  edges[ii]->get_lower_vertex()->get_index_global();
			NEGF_ASSERT(other_vert_idx != -1, "Faulty list get_edges_near(vertex)");

			indices[nonzeros] = other_vert_idx;
			if (coeffs!=NULL) {
			coeffs[nonzeros]  = 0;
			vector<Element *> elems_near_edge = grid->get_elems_near(edges[ii]);
			for (uint jj = 0; jj < grid->get_num_elems_near(edges[ii]); jj++)
			{
				uint local_edge_idx = elems_near_edge[jj]->get_local_index(edges[ii]);

				uint elem_idx  = elems_near_edge[jj]->get_index_global(); // epsilon is defined on elements
				
				coeffs[0]         += epsilon->get_value(elem_idx) * coefficient[elem_idx][local_edge_idx];
				coeffs[nonzeros]  -= epsilon->get_value(elem_idx) * coefficient[elem_idx][local_edge_idx];
			}
			}	// coeffs!=NULL
			nonzeros++;
		}
		return;
	}
	if(eqn==epsilon) {
		vector<Element *> elems = grid->get_elems_near(vertex);
		for (uint ii=0; ii < grid->get_num_elems_near(vertex); ii++)
		{
			uint elem_idx     = elems[ii]->get_index_global();
			indices[nonzeros] = elem_idx;
			if (coeffs!=NULL) {
			coeffs[nonzeros]  = 0;
			uint irrelevant_edges = 0; // just for debugging
			for (uint jj = 0; jj < elems[ii]->get_num_edges(); jj++)
			{
				int other_vert_idx = -1;
				if (elems[ii]->get_edge(jj)->get_lower_vertex()==vertex)
					other_vert_idx =  elems[ii]->get_edge(jj)->get_upper_vertex()->get_index_global();
				if (elems[ii]->get_edge(jj)->get_upper_vertex()==vertex)
					other_vert_idx =  elems[ii]->get_edge(jj)->get_lower_vertex()->get_index_global();
				if (other_vert_idx != -1) // other_idx = -1 is easily possible
					coeffs[nonzeros] += coefficient[elem_idx][jj] * (this->get_value(line) - this->get_value(other_vert_idx));
				else
					irrelevant_edges++;
			}
			switch (elems[ii]->get_type())  // just to make sure...
			{
			case element_type::tetrahedron:
				NEGF_ASSERT(irrelevant_edges==3, "number of irrelevant edges is not correct."); break;
			case element_type::triangle:
				NEGF_ASSERT(irrelevant_edges==1, "number of irrelevant edges is not correct."); break;
			case element_type::rectangle:
				NEGF_ASSERT(irrelevant_edges==2, "number of irrelevant edges is not correct."); break;
			case element_type::cuboid:
				NEGF_ASSERT(irrelevant_edges==9, "number of irrelevant edges is not correct."); break;
			default:
				break; // no checking
			}
			} // coeffs!=NULL
			nonzeros++;
		}
		return;
	}
	if(eqn==edensity) {
		nonzeros   = 1;
		indices[0] = line;
		if (coeffs!=NULL) {
			if (vertex->is_at_contact() && vertex->get_contact()->get_bndcond(quantities::potential)==bndconds::BC_Neumann) {
				// treatment of Neumann contact vertices: no space charge at contact!
				//coeffs[0]  = 0.0;
				coeffs[0]  = node_measure[line] * ec;
			} else {
				coeffs[0]  = node_measure[line] * ec;
			}
		}
		return;
	}
	if(eqn==hdensity) {
		nonzeros   = 1;
		indices[0] = line;
		if (coeffs!=NULL) {
			if (vertex->is_at_contact() && vertex->get_contact()->get_bndcond(quantities::potential)==bndconds::BC_Neumann) {
				// treatment of Neumann contact vertices: no space charge at contact!
				//coeffs[0]  = 0.0;
				coeffs[0]  = -1 * node_measure[line] * ec;
			} else {
				coeffs[0]  = -1 * node_measure[line] * ec;
			}
		}
		return;
	}
	if(eqn==doping) {
		nonzeros   = 1;
		indices[0] = line;
		if (coeffs!=NULL) {
			if (vertex->is_at_contact() && vertex->get_contact()->get_bndcond(quantities::potential)==bndconds::BC_Neumann) {
				// treatment of Neumann contact vertices: no space charge at contact!
				//coeffs[0]  = 0.0;
				coeffs[0]  = -1 * node_measure[line] * ec;
			} else {
				coeffs[0]  = -1 * node_measure[line] * ec;
			}
		}
		return;
	}
	for (uint nu = 4; nu <  dependencies.size(); nu++) // interpolated charge from other grids	for (nu = 4; nu < dependencies.size(); ++nu)
	{
		if (eqn==dependencies[nu])
		{
			int charge_sign;
			if (dependencies[nu]->get_type() == quantities::electron_density)
				charge_sign = -1;
			else if (dependencies[nu]->get_type() == quantities::hole_density)
				charge_sign = +1;
			else
				NEGF_EXCEPTION("Poisson dependency not supported.");
			nonzeros   = 1;
			indices[0] = line;
			if (coeffs!=NULL){
				if (vertex->is_at_contact() && vertex->get_contact()->get_bndcond(quantities::potential)==bndconds::BC_Neumann) {
					// treatment of Neumann contact vertices: no space charge at contact!
					coeffs[0]  = 0.0;
					//coeffs[0]  = -1 * node_measure[line] * charge_sign * ec;
				} else {
					coeffs[0]  = -1 * node_measure[line] * charge_sign * ec;
				}
			}
			return;
		}
	}

	NEGF_EXCEPTION("Poisson equation does not seem to depend on this equation.");
	return;
);}
