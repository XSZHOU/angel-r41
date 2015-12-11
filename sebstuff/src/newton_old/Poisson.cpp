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
#include "Poisson.h"
using namespace negf;


Poisson::Poisson(const Geometry * const grid_):
	grid(grid_),
	epsilon(0),
	edensity(0),
	hdensity(0),
	doping(0)
{STACK_TRACE(
	NEGF_ASSERT(grid!=0, "tried to create Poisson equation with null pointer(s).");
	this->number_of_variables = grid->get_num_vertices();
	this->its_type = quantities::potential;

	this->neumann_field = 0.0;
	
	dependencies.clear();
	dependencies.push_back(0);	// mandatory dependence: dependencies[0] = epsilon
	dependencies.push_back(0);	// mandatory dependence: dependencies[1] = electron density
	dependencies.push_back(0);	// mandatory dependence: dependencies[2] = hole density
	dependencies.push_back(0);	// mandatory dependence: dependencies[3] = net ionized dopant density

	current_variable_values.resize(this->number_of_variables, 0.0);
);}


void Poisson::add_source(Equation * density_eqn) 
{STACK_TRACE(
	NEGF_ASSERT(density_eqn!=0, "Tried to allocate null pointer.");
	NEGF_ASSERT(density_eqn->get_type() == quantities::electron_density 
			 || density_eqn->get_type() == quantities::hole_density, "Wrong equation.");
	NEGF_ASSERT(density_eqn->get_num_variables() == this->get_num_variables(), "Density eqn has wrong # variables.");
	dependencies.push_back(density_eqn); 
);}


void Poisson::allocate_epsilon(Equation * epsilon_eqn) 
{STACK_TRACE(
	NEGF_ASSERT(epsilon_eqn!=0, "Tried to allocate null pointer.");
	NEGF_ASSERT(epsilon_eqn->get_type() == quantities::dielectric_coeff, "Wrong equation.");
	NEGF_ASSERT(epsilon_eqn->get_num_variables() == this->grid->get_num_elements(), "Epsilon is on different grid.");
	dependencies[0] = epsilon_eqn;
	epsilon			= epsilon_eqn;
);}


/** Allocate the net ionized dopants.
 *  ATTENTION: POSITIVE VALUES MEAN N-DOPING, SINCE THEY ARE POSITIVELY CHARGED IONIZED DOPANTS */
void Poisson::allocate_doping(Equation * doping_eqn) 
{STACK_TRACE(
	NEGF_ASSERT(doping_eqn!=0, "Tried to allocate null pointer.");
	NEGF_ASSERT(doping_eqn->get_type() == quantities::hole_density, "Wrong equation.");
	NEGF_ASSERT(doping_eqn->get_num_variables() == this->get_num_variables(), 
					"Doping is on different grid.");
	dependencies[3] = doping_eqn;
	doping			= doping_eqn;
);}


void Poisson::allocate_owncharge(Equation * edensity_eqn, Equation * hdensity_eqn) 
{STACK_TRACE(
	NEGF_ASSERT(edensity_eqn!=0 && hdensity_eqn!=0, "Tried to allocate null pointer.");
	NEGF_ASSERT(   edensity_eqn->get_type() == quantities::electron_density
				&& hdensity_eqn->get_type() == quantities::hole_density, "Wrong equation.");
	NEGF_ASSERT(edensity_eqn->get_num_variables() == this->get_num_variables(), 
					"el. density is on different grid.");
	NEGF_ASSERT(hdensity_eqn->get_num_variables() == this->get_num_variables(), 
					"hole density is on different grid.");
	dependencies[1] = edensity_eqn;
	dependencies[2] = hdensity_eqn; 
	edensity		= edensity_eqn;
	hdensity		= hdensity_eqn;
);}


/** Compute initial guess of electrostatic potential. <BR>
 *  This is done by a linear interpolation of the Dirichlet values at the contacts. <BR>
 *  If one of the contacts is not Dirichlet but Neumann or so, the initial guess is zero.
 *  @param tstamp The timestamp the equation gets after setting the initial guess */
void Poisson::initial_guess(uint tstamp)
{STACK_TRACE(
	// check that boundary conditions for potential are set
	for (uint cc = 0; cc < this->grid->get_num_contacts(); cc++)
		NEGF_ASSERT( this->grid->get_contact(cc)->bc_is_ready(quantities::potential),
					"Please set the boundary condition or all the contacts first." );
	
	vector<double> potential_init; 
	potential_init.resize(this->grid->get_num_vertices(), 0.0);
	
	// if not all contacts have Dirichlet BC, intitialize with zero
	bool calc_potential = true;
	for (uint cc = 0; cc < this->grid->get_num_contacts(); cc++)
		if (this->grid->get_contact(cc)->get_bndcond(quantities::potential) != bndconds::BC_Dirichlet)
			calc_potential = false;
	if (!calc_potential) {
		logmsg->emit(LOG_WARN, "Initial guess for potential is =0 because not all contacts ahve Dirichlet BC.");
		this->set_values(potential_init, tstamp);
		return;
	}
	
	vector<double> distance_to_contacts; 
	vector<uint> nearest_vertex; 		// stores the nearest contact vertex for all contacts (index of contact list of vertices)
	distance_to_contacts.resize(this->grid->get_num_contacts(), 0.0);
	nearest_vertex.resize(this->grid->get_num_contacts(), 0); 
	for (uint ii = 0; ii < this->grid->get_num_vertices(); ii++)
	{		
		// get distance of vertex to contacts
		for (uint cc = 0; cc < this->grid->get_num_contacts(); cc++) {
			const vector<Vertex *> & contact_vertices = grid->get_contact(cc)->get_contact_vertices();
			NEGF_ASSERT(contact_vertices.size()>0, "no contact vertices ?!");
			
			distance_to_contacts[cc] = grid->get_distance(grid->get_vertex(ii), contact_vertices[0]);
			nearest_vertex[cc] = 0;
			for (uint jj=1; jj < contact_vertices.size(); jj++) {
				double distance = grid->get_distance(grid->get_vertex(ii), contact_vertices[jj]);
				if (distance < distance_to_contacts[cc] ) {
					distance_to_contacts[cc] = distance;
					nearest_vertex[cc] = jj;
				}
			}
		}

		// old version was to interpolate linearly via 
		// phi(x) = \frac{1}{\sum_{i,j}} \sum_{i,j} \frac{d_i}{d_i+d_j} \phi_j + \frac{d_j}{d_i+d_j} \phi_i
		// where i-j is the contact pair i-j and d_i is the distance of vertex x to contact i
		// not implemented anymore
		 
		// -------------------------------------
		// case where we have 0, 1 or 2 contacts
		// -------------------------------------
		if (grid->get_num_contacts()==0)
			NEGF_EXCEPTION("Must have at least one contact to make an initial guess for the potential.");
		if (grid->get_num_contacts()==1) {
			potential_init[ii] = grid->get_contact(0)->get_bc_value(quantities::potential, 0);
			continue;
		}		
		if (grid->get_num_contacts()==2) {
			double sum = distance_to_contacts[0]+distance_to_contacts[1];
			potential_init[ii] =  distance_to_contacts[0]/sum * grid->get_contact(1)->get_bc_value(quantities::potential, 0)
								+ distance_to_contacts[1]/sum * grid->get_contact(0)->get_bc_value(quantities::potential, 0);
			continue;
		}
		
		
		// -------------------------------------
		// case where we have 3 or more contacts
		// -------------------------------------
		Element * helper;
		uint num_nearest_contacts;
		switch (grid->get_dimension())
		{
		case 1:
			NEGF_EXCEPTION("It is hard to imagine a 1D grid with more than 2 contacts.");
			break;
		case 2: /** GET 3 NEAREST CONTACTS **/
			helper = new Element(0, element_type::triangle);
			num_nearest_contacts = 3;
			break;
		case 3: /** GET 4 NEAREST CONTACTS IF POSSIBLE **/
			// if we only have 3 contacts: the fourth point will be VERY far away
			helper = new Element(0, element_type::tetrahedron);
			num_nearest_contacts = 4;
			break;
		default:
			NEGF_EXCEPTION("Strange dimensionality.");
		}
		
		vector<uint> nearest_contacts;
		for (uint jj = 0; jj < grid->get_num_contacts(); jj++) {
			uint counter = 0;
			for (uint kk = 0; kk <grid->get_num_contacts(); kk++) {
				if (distance_to_contacts[kk] < distance_to_contacts[jj])
					counter++;
			}
			if (counter<num_nearest_contacts) { 
				// contact jj is amongst the num_nearest_contacts closest contacts
				nearest_contacts.push_back(grid->get_contact(jj)->get_index());
			}
		}
		
		/** assign vertices to the helper element
		 *  for dim=3, num_contacts=3 we use the HACK that the fourth contact is far, far away */
		uint nvertices;
		if (grid->get_dimension()==3 && grid->get_num_contacts()==3)
		{
			NEGF_ASSERT(nearest_contacts.size()==3, "nearest_contacts.size()!=3");
			nvertices = 3;
			const double very_large_distance = constants::convert_from_SI(units::length, 100.0); // 100 meter
			Vertex * helpvert = new Vertex(0, very_large_distance, very_large_distance, very_large_distance);
			helper->add_vertex(helpvert);
		} else {
			NEGF_ASSERT(nearest_contacts.size()==num_nearest_contacts, "nearest_contacts.size()!=num_nearest_contacts");
			nvertices = num_nearest_contacts;
		}
		for (uint jj = 0; jj < nvertices; jj++) {
			helper->add_vertex(grid ->get_contact(nearest_contacts[jj])
							 			->get_contact_vertex(nearest_vertex[nearest_contacts[jj]]) );
		}					   		
		helper->prepare();
		
		vector<double> interpol_coeffs;
		helper->get_linear_combination_coefficients(grid->get_vertex(ii), interpol_coeffs);
		if (grid->get_dimension()==3 && grid->get_num_contacts()==3)
		{
			NEGF_ASSERT(interpol_coeffs.size()==4 && nearest_contacts.size()==3, "interpol_coeffs.size()!=4");
			// interpol_coeffs[3] will not be used
		} else {		
			NEGF_ASSERT(interpol_coeffs.size()==nearest_contacts.size(), "interpol_coeffs.size()!=nearest_contacts.size()");
		}
		
		// CALCULATE INTERPOLATED POTENTIAL
		double check = 0.0;
		for (uint jj = 0; jj < nearest_contacts.size(); jj++) {
			potential_init[ii] += interpol_coeffs[jj] * grid->get_contact(nearest_contacts[jj])->get_bc_value(quantities::potential, 0);
			check += interpol_coeffs[jj];
		}
		
		// some checking
		if (fabs(check - 1.0) > 1e-5) {
			logmsg->emit(LOG_INFO,"vertex %d: check failed (%13.10e instead of 1.0)", ii, check);
			for (uint jj = 0; jj < interpol_coeffs.size(); jj++)
				logmsg->emit(LOG_INFO,"   interpol_coeffs[%d]=%e",jj,interpol_coeffs[jj]);
			if (grid->get_dimension()==3 && grid->get_num_contacts()==3) {
				NEGF_EXCEPTION("It is likely that the helper vertex was put too close to the structure.");
			} else {
				NEGF_EXCEPTION("Check your interpol_coeffs.");
			}
		}
		if (grid->get_vertex(ii)->is_at_contact())
		{
			Contact * contact = grid->get_vertex(ii)->get_contact();
			if (fabs(potential_init[ii] - contact->get_bc_value(quantities::potential, 0)) > 1e-6) {
				logmsg->emit(LOG_INFO,"vertex %d at contact %d should have value %e.",
							ii, contact->get_index(), contact->get_bc_value(quantities::potential, 0) );
				for (uint jj = 0; jj < nearest_contacts.size(); jj++) {
					uint cidx = nearest_contacts[jj];
					logmsg->emit(LOG_INFO,"   nearest contact %d is number %d, value %e, coeff %e; nearest_vertex_global_idx=%d.",
							jj, cidx, 
							grid->get_contact(cidx)->get_bc_value(quantities::potential, 0),
							interpol_coeffs[jj], 
							grid->get_contact(cidx)->get_contact_vertex(nearest_vertex[cidx])->get_index_global());
				}
				NEGF_EXCEPTION("For contact vertices the guess should be exactly the contact potential.");
			}
		}
	}
	this->set_values(potential_init, tstamp);
);}


void Poisson::set_neumann_electric_field(double field)
{STACK_TRACE(
	this->neumann_field = field;
);}

