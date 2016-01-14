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
#include "PoissonFEM.h"

#ifndef NOTDKP

using namespace negf;

#include "tdkp/common/Logger.h"
#include "tdkp/probdefs/PoissonEquationAQUA.h"
#define VEPPO_ASSERT(cond, msg)	   if(!(cond)) { ostringstream _sout; _sout << msg; if(_sout.str().size() == 0) { _sout << #cond; } NEGF_EXCEPTION(_sout.str().c_str()); }


PoissonFEM::PoissonFEM(
	const Geometry * const grid_,
	const MaterialDatabase& material_database, 
	const char* gridfile, 
	//const char* polarization_charge_file
	const vector<double> & static_rhs_
) :	Poisson(grid_),
	piezo_decreaser(1.0)
{//STACK_TRACE(
	/** the following variables have already been set in the Poisson() constructor:
	 *  - number_of_variables = grid->get_num_vertices
	 *  - its_type = quantities::potential
	 *  - grid
	 */
	NEGF_ASSERT(grid!=0, "null pointer encountered.");
	
	try {
	if (mpi->get_rank()==constants::mpi_master_rank) {
		tdkp::Logger::get_instance()->set_level(tdkp::LOG_INFO);
	} else {
		tdkp::Logger::get_instance()->set_level(tdkp::LOG_ERROR); 
	}
 	tdkp::PoissonEquationAQUA tdkp_poisson_wrapper(
 		gridfile, /*polarization_charge_file*/"",
 		constants::convert_from_SI(units::dielectric, constants::SIeps0), // this gives eps0 in negf units
 		constants::convert_from_SI(units::length, 1.0e-6)         // and this should give length in negf units, grid is in micrometers
 	);
 
 	for(unsigned int ii = 0; ii < material_database.get_num_materials(); ii++) {
 		tdkp_poisson_wrapper.set_permittivites(
 			material_database.get_material(ii)->get_name(),
 			material_database.get_material(ii)->get("static_dielectric_constant")
 		);
 	} 	
 	tdkp_poisson_wrapper.assemble();
  	tdkp_poisson_wrapper.get_stiffness_matrix(stiffness_matrix);
  	tdkp_poisson_wrapper.get_mass_matrix(mass_matrix);
  	tdkp_poisson_wrapper.get_icol(icol);
  	tdkp_poisson_wrapper.get_prow(prow);
  	//tdkp_poisson_wrapper.get_static_rhs(static_rhs);
  	
  	// assign polarization
  	if (static_rhs_.size()==grid->get_num_vertices()) {
  		this->static_rhs = static_rhs_;
  	} else {
  		NEGF_ASSERT(static_rhs_.size()==0, "size of handed over static_rhs should be either 0 or num_vertices.");
  		this->static_rhs.assign(grid->get_num_vertices(), 0.0);
  	}
  	
  	ostringstream sout;
  	sout << "prow.size() - 1 (" << prow.size()-1 << ") == grid->get_num_vertices() (" << grid->get_num_vertices() << ")";
  	NEGF_ASSERT(prow.size()             == grid->get_num_vertices()+1, sout.str().c_str());
	NEGF_ASSERT(prow[grid->get_num_vertices()] == (signed)stiffness_matrix.size(), "prow[grid_->get_num_vertices()] == stiffness_matrix.size()");
	NEGF_ASSERT(stiffness_matrix.size() == mass_matrix.size(), "stiffness_matrix.size() == mass_matrix.size()");
	NEGF_ASSERT(icol.size()             == stiffness_matrix.size(), "icol.size() == stiffness_matrix.size()");
	NEGF_ASSERT(static_rhs.size()       == (unsigned int)grid->get_num_vertices(), "static_rhs.size() == grid->get_num_vertices()");  
	
	/*logmsg->emit(LOG_INFO,"TDKP Prow: ");
	for (uint ii=0; ii<prow.size(); ii++) {
		logmsg->emit_noendl(LOG_INFO,"%d  ",prow[ii]);
	}
	logmsg->emit(LOG_INFO,"");
	
	logmsg->emit(LOG_INFO,"TDKP Icol: ");
	for (uint ii=0; ii<icol.size(); ii++) {
		logmsg->emit_noendl(LOG_INFO,"%d  ",icol[ii]);
	}
	logmsg->emit(LOG_INFO,"");
	
	logmsg->emit(LOG_INFO,"TDKP Mass: ");
	for (uint ii=0; ii<mass_matrix.size(); ii++) {
		logmsg->emit_noendl(LOG_INFO,"%.3f  ",mass_matrix[ii]);
	}
	logmsg->emit(LOG_INFO,"");
	
	logmsg->emit(LOG_INFO,"TDKP Stiffness: ");
	for (uint ii=0; ii<stiffness_matrix.size(); ii++) {
		logmsg->emit_noendl(LOG_INFO,"%.3f  ",stiffness_matrix[ii]);
	}
	logmsg->emit(LOG_INFO,"");*/
	
	
	bool debug = false;
	if (debug && grid->get_dimension()==2) {
		// get vector containing int Ni rho for static rho (piezo/polarization), all shape funcrtions Ni 
		vector<double> static_charges;
		tdkp_poisson_wrapper.get_static_rhs(static_charges);
		NEGF_ASSERT(static_charges.size()==grid->get_num_vertices(), "expected num_vertices entries in static charge vector.");
		// static_charges[i] stores int Ni(x) rho(x) dx. for quasi-1D structures, rho(x,y)=rho*delta(y); hence static_charges[i] = rho * int_xi^xj Ni(x,0) dx
		// furthermore, on the edge i-j the linear shape function is Ni(x,0)=(xj-x)/(xj-xi). Therefore:
		// static_charges[i] = rho * ( xj - 1/(xj-xi) * 0.5(xj^2-xi^2) ) = rho * (xj - 0.5(xj-xi)(xj+xi)/(xj-xi)) = rho * (xj - 0.5(xj+xi)) 
		//                   = rho * 0.5(xj-xi).
		
		// we expect charges at material boundaries
		for (uint ii=0; ii<grid->get_num_vertices(); ii++) {
			Vertex * v = grid->get_vertex(ii);
			const vector<Region *> & adj_regions = grid->get_regions_near(v);
			if (adj_regions.size()<2) continue;
			
			const vector<Edge *> & adj_edges = grid->get_edges_near(v);
			Edge   * boundary_edge = 0;
			Vertex * v2 = 0;
			for (uint jj=0; jj<adj_edges.size(); jj++) {
				Vertex * vert = (adj_edges[jj]->get_lower_vertex()==v) ? adj_edges[jj]->get_upper_vertex() : adj_edges[jj]->get_lower_vertex();
				const vector<Region *> & adj_regions2 = grid->get_regions_near(vert);
				if (adj_regions2.size()>=2) {
					NEGF_ASSERT(boundary_edge==0, ">1 boundary edges???");
					boundary_edge = adj_edges[jj];
					v2 = vert;
				}
			}
			double rho1 = static_charges[v->get_index_global()];
			double rho2 = static_charges[v2->get_index_global()];
			double x1 = v->get_coordinate(0);
			double x2 = v2->get_coordinate(0);
			if (x2 < x1) continue; // double counting 
			
			logmsg->emit(LOG_INFO,"v(%e,%e), v2(%e,%e): rho1=%10.3e, rho2=%10.3e, rho1/(0.5(x2-x1))=%10.3e, rho2/(0.5(x1-x2))=%10.3e.",
					x1, v->get_coordinate(1), x2, v2->get_coordinate(1), rho1, rho2, rho1/(0.5*(x2-x1)), rho2/(0.5*(x1-x2)));
		}
		logmsg->emit(LOG_INFO, "rho/(0.5(xi-xj) should have units (%.1gec) / (%.1gnm)^2.",
					1.0/constants::convert_from_SI(units::charge, constants::SIec), 1.0/constants::convert_from_SI(units::length, 1e-9));
	}
	
	
	} catch (std::string s) { NEGF_EXCEPTION(s.c_str()); }
//);  	
}

void PoissonFEM::set_piezo_decreaser(double new_decreaser) 
{STACK_TRACE(
	NEGF_ASSERT(new_decreaser >= 0.0 && new_decreaser < 10.0, "invalid decreaser value");
	this->piezo_decreaser = new_decreaser;
);}
	

/** Return the Newton function (whose roots are searched) evaluated for variable ii
 * @param ii the global vertex index 
*/
double PoissonFEM::get_newton_function(uint line) const
{//STACK_TRACE(
	const Vertex * vertex = grid->get_vertex(line);
	const double ec       = constants::convert_from_SI(units::charge, constants::SIec);
	double result = 0.0;

	NEGF_ASSERT(vertex!=0 && epsilon!=0 && edensity!=0 && hdensity!=0 && doping!=0, "Some quantity is missing.");

	// -------------------------------------------------------------------------
	// newton function is: - nabla^2 u - rho = 0; and u - u_bc at the boundary
	// -------------------------------------------------------------------------
	
	// -------------------------------------------------
	// first, check if vertex is at contact
	// -------------------------------------------------
	if (vertex->is_at_contact() && 
	    vertex->get_contact()->get_bndcond(quantities::potential) == bndconds::BC_Dirichlet) {
	    // vertex->get_contact_vertex_index() returns the index of the vertex in the list of the contact's vertices
		return this->get_value(line) - vertex->get_contact()->get_bc_value(quantities::potential, vertex->get_contact_vertex_index());		
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
			return (this->get_value(line) - this->get_value(v_adjacent_to_device->get_index_global()));
		} else {
			// vertex which is at contact but directly connected to the device interior: 
			// continue with normal treatment, normal formula
			
			// NO! same value as device-internal vertex!
			// desired value = value at device interface - distance * E-field
			if (fabs(this->neumann_field) > 1e-10) NEGF_ASSERT(grid->get_dimension()==1, "Neumann field only supported for d=1!");
			double dist = vertex->get_coordinate(0) - device_interior_vertex->get_coordinate(0);
			return (this->get_value(line) - this->get_value(device_interior_vertex->get_index_global()) - dist * this->neumann_field);
		}
	}
	
	
	
	// -------------------------------------------------
	// so, its not a contact, but a standard equation,	
	// -------------------------------------------------
	for(int kk = prow[line]; kk < prow[line + 1]; kk++) {
		VEPPO_ASSERT((signed)stiffness_matrix.size() > kk, "");
		VEPPO_ASSERT((signed)mass_matrix.size() > kk, "");
		// -------------------------------------------------
		// we assemble the function starting with stiffness matrix contribution
		// -------------------------------------------------
		result += stiffness_matrix[kk] // epsilon grad Ni * grad Nkk 
		        * this->get_value(icol[kk]);
		// -------------------------------------------------        
		// now we subtract the mass matrix contributions
		// containing the nodal charge field	
		// -------------------------------------------------
		result -= mass_matrix[kk] * ec * (
		       		hdensity->get_value(icol[kk])  
			      - edensity->get_value(icol[kk])
				  + doping->get_value(icol[kk])  
			    );
		// -------------------------------------------------
		// and we also subtract the mass matrix contributiosn
		// caused by charge from quantized charge density 
		// -------------------------------------------------
		for(unsigned int nu = 4; nu < dependencies.size(); nu++) {
			int charge_sign;
			if(dependencies[nu]->get_type() == quantities::electron_density) {
				charge_sign = -1;
			} else if(dependencies[nu]->get_type() == quantities::hole_density) {
				charge_sign = +1;
			} else {
				NEGF_EXCEPTION("Poisson dependency not supported.");
			}
			result -= mass_matrix[kk] * ec * charge_sign * dependencies[nu]->get_value(icol[kk]);
		}
	}
	
	// -------------------------------------------------
	// finally we have to subtract the static rhs 
	// -------------------------------------------------
	if (fabs(piezo_decreaser-1.0)>1e-10 && line==30) cout << "piezo_decreaser=" << this->piezo_decreaser << "    ";
	result -= ec * this->piezo_decreaser * static_rhs[line]; 
	
	// output debug info if sheet charge exceeds 0.01 C/m2
	/*#pragma omp critical
	{
	if (fabs(static_rhs[line]) > 0.01 * constants::convert_from_SI(units::charge, 1.0) / constants::convert_from_SI(units::area, 1.0)) {
		cout.precision(4);
		double tmp = 0.0;
		for(int kk = prow[line]; kk < prow[line + 1]; kk++) {
			tmp -= mass_matrix[kk] * ec * (
		       		hdensity->get_value(icol[kk])  
			      - edensity->get_value(icol[kk])
				  + doping->get_value(icol[kk])  
			    );
			for(unsigned int nu = 4; nu < dependencies.size(); nu++) {
				int charge_sign;
				if(dependencies[nu]->get_type() == quantities::electron_density) {
					charge_sign = -1;
				} else if(dependencies[nu]->get_type() == quantities::hole_density) {
					charge_sign = +1;
				} else {
					NEGF_EXCEPTION("Poisson dependency not supported.");
				}
				tmp -= mass_matrix[kk] * ec * charge_sign * dependencies[nu]->get_value(icol[kk]);
			}
		}
		cout << "vertex " << line << " (y=" << vertex->get_coordinate(1) << "): charge_dens=" << 	tmp 
			 << ", static_rhs=" << -ec*hack_factor*static_rhs[line] << endl;  
	}
	}*/
	
	return result;
//);
}


/** @param ii which variable
 * @param eqn the dependency w.r.t which the derivative is computed
 * @param nonzeros number storing the number of array entries used
 * @param indices array storing the indices of the dependency equation where the derivative is nonzero
 * @param coeffs array storing the corresponding derivatives. coeffs==NULL --> sparsity pattern only
*/
void PoissonFEM::get_newton_direct_derivatives(uint line, const Equation * eqn, 
								uint & nonzeros, uint indices[], double coeffs[]) const
{
//STACK_TRACE(
	const Vertex * vertex = grid->get_vertex(line);
	const double ec = constants::convert_from_SI(units::charge, constants::SIec);

	NEGF_ASSERT(vertex!=0 && epsilon!=0 && edensity!=0 && hdensity!=0 && doping!=0, "Some quantity is missing.");
	NEGF_ASSERT(eqn!=NULL, "the dependency requested was NULL!");

	nonzeros = 0;
	
	// -------------------------------------------------------------------------------------
	// treatment of Dirichlet contact vertices (Poisson eqn is replaced with explicit value)
	// -------------------------------------------------------------------------------------
	if(vertex->is_at_contact() && vertex->get_contact()->get_bndcond(quantities::potential)==bndconds::BC_Dirichlet) {		
		if(eqn == this) {
			nonzeros   = 1;
			indices[0] = line;
			// return if this is only a structure inquiry
			if(coeffs == NULL) {
				return;
			}
			// equation is u - u_bc
			coeffs[0]  = 1.0;
			return;
		}
		// -------------------------------------		
		// don't depend on anything else
		// -------------------------------------
		for(unsigned int nu = 0; nu < dependencies.size(); nu++) {
			if(eqn == dependencies[nu]) {
				nonzeros = 0;
				return;
			}
		}
		NEGF_EXCEPTION("Unknown equation.");
		return;
	}
	
	// -------------------------------------------------------------------------------------
	// NEGF: we also have to do something in case of contacts and Neumann conditions
	// -------------------------------------------------------------------------------------
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
					coeffs[0]  = 1.0;
					coeffs[1]  = -1.0;
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
					coeffs[0]  = 1.0;
					coeffs[1]  = -1.0;
				}
				return;
			} else {
				nonzeros = 0;
				return;
			}
		}
	}
		
	
	
	// ----------------------------------------------
	// the pattern does not change for any equation
	// ----------------------------------------------
	nonzeros = prow[line + 1] - prow[line];
	for(unsigned int kk = 0; kk < nonzeros; kk++) {		
		indices[kk] = icol[prow[line] + kk];
	}
	// ----------------------------------------------
	// structure queries are fine with that
	// ----------------------------------------------
	if(coeffs == 0) {
		return;	
	}
	
	// ----------------------------------------------
	// direct dependence of Poisson equation on potential
	// given by int eps grad Ni grad Nj
	// ----------------------------------------------
	if(eqn == this) { 
		nonzeros = prow[line + 1] - prow[line];
		// copy stiffness matrix entries
		for(unsigned int kk = 0; kk < nonzeros; kk++) {		
			coeffs[kk]  = stiffness_matrix[prow[line] + kk];
		}
		return;
	}
	// ----------------------------------------------
	// that's for the kassel guy ... i assume epsilon to be constant
	// ----------------------------------------------
	if(eqn == epsilon) {
		nonzeros = 0;			
		return;
	}
	// ----------------------------------------------
	// for the rest: i just need to know the vorzeichen of the charge
	// ----------------------------------------------
	int charge_sign = 0;
	if(eqn == edensity) {
		charge_sign = - 1;				
	} else if(eqn==hdensity || eqn == doping) {
		charge_sign = + 1;	
	} else {
		// ----------------------------------------------
		// huston, check quantized spreads for their sign
		// ----------------------------------------------
		for(unsigned int nu = 4; nu < dependencies.size(); nu++) {	
			if(dependencies[nu] == eqn) {
				if(dependencies[nu]->get_type() == quantities::electron_density) {
					charge_sign = - 1;
				} else if(dependencies[nu]->get_type() == quantities::hole_density) {
					charge_sign = + 1;
				} else {
					NEGF_EXCEPTION("Poisson dependency not supported.");
				}
			}
		}
	}
	if(charge_sign == 0) {
		NEGF_EXCEPTION("Poisson equation does not seem to depend on this equation.");
		return;
	}
	// --------------------------------------
	// so, now its simply int ec Ni Nj
	// --------------------------------------
	for(unsigned int kk = 0; kk < nonzeros; kk++) {		
		coeffs[kk] = (-1.0) * ec * charge_sign * mass_matrix[prow[line] + kk];
	}
	return;
	
//);
}

#endif // ndef NOTDKP


