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
#include "InterfaceEffMass.h"
using namespace negf;

InterfaceEffMass::InterfaceEffMass(const Geometry * xspace_,
						 		const MaterialDatabase * db,
						 		const double temperature) throw (Exception *):
	xspace(xspace_)
{STACK_TRACE(
	NEGF_ASSERT(xspace!=NULL, "null pointer encountered.");
	NEGF_ASSERT(xspace->get_dimension()==1, "class is only designed for 1D.");
	uint Nvert = xspace->get_num_vertices();
	
	const double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
	
	// construct the FEM-discretized basis function overlap matrix and hamiltonian
	this->overlap = DGEMatrix(Nvert,Nvert);
	this->core_hamiltonian = GEMatrix(Nvert,Nvert);
	for (uint xx=1; xx<=Nvert; xx++) 
	{
		Vertex * vx = xspace->get_vertex(xx-1);
		const vector<Element *> elems_near_x = xspace->get_elems_near(vx);
		NEGF_ASSERT(elems_near_x.size()==1 || elems_near_x.size()==2, "unexpected number of edges near x.");
		for (uint ii = 0; ii < elems_near_x.size(); ii++) 
		{
			NEGF_ASSERT(elems_near_x[ii]->get_num_edges()==1, "expected exactly 1 edge.");
			Edge * edge = elems_near_x[ii]->get_edge(0);
			Vertex * vy = (edge->get_lower_vertex()==vx) ? edge->get_upper_vertex() : edge->get_lower_vertex();
			uint yy = vy->get_index_global() + 1;
			double lxy =  edge->get_length();
			
			const PropertyContainer<double> * mat = elems_near_x[ii]->get_region()->get_material();
			double m   = constants::convert_from_SI(units::mass, mat->get("electron_effective_mass") * constants::SIm0);
			double fac = hbar*hbar / (2.0 * m);
			double Ec  = TdkpInfoDesk::get_cbedge(mat, temperature, db);
			
			// CONVERSION TO TDKP UNITS: Ec is energy, must be eV; fac is  lxy is length, must be nm
			Ec  = 1.0/constants::SIec * Ec/constants::convert_from_SI(units::energy, 1.0);
			lxy = 1e9                 * lxy/constants::convert_from_SI(units::length, 1.0);
			// fac is energy*length^2, must be eV*nm^2
			fac = 1.0/constants::SIec*1e18 * fac/(constants::convert_from_SI(units::energy, 1.0)*constants::convert_from_SI(units::area, 1.0));
			
			this->overlap(xx,xx) += lxy / 3.0;
			this->overlap(xx,yy) = lxy / 6.0;
			
			// band edge - element-wise constant
			this->core_hamiltonian(xx,xx) += Ec * lxy / 3.0;
			this->core_hamiltonian(xx,yy) += Ec * lxy / 6.0;
			
			// kinetic part of Hamiltonian
			this->core_hamiltonian(xx,xx) += fac / lxy;
			this->core_hamiltonian(xx,yy) += -fac / lxy;
		}
	}
		
	// checks
	Matd tmp2(this->overlap);
	Matd tmp1(Nvert,Nvert); trans(tmp2, tmp1); // tmp1 = overlap^T
	tmp1 -= tmp2;
	NEGF_ASSERT(negf_math::matrix_norm(tmp1) < 1e-14, "M should be symmetric.");
	Matc tmp4(this->core_hamiltonian);
	Matc tmp3(Nvert,Nvert);	trans(tmp4, tmp3); // tmp3 = core_hamiltonian^T
	tmp3 -= tmp4;
	NEGF_ASSERT(negf_math::matrix_norm(tmp3) < 1e-14, "H should be symmetric.");
	
	this->potential_augmented_hamiltonian = this->core_hamiltonian;
	this->potential_and_k_augmented_hamiltonian = this->core_hamiltonian;
	
);}


// add FEM-discretized phi!
void InterfaceEffMass::set_potential(const vector<double>& node_potential)
{STACK_TRACE(
	uint Nvert = xspace->get_num_vertices();
	NEGF_ASSERT(node_potential.size()==Nvert, "inconsistent potential vector.");
	this->potential_augmented_hamiltonian = this->core_hamiltonian;
	for (uint xx=1; xx<=Nvert; xx++) 
	{
		Vertex * vx = xspace->get_vertex(xx-1);
		const vector<Edge *> edges_near_x = xspace->get_edges_near(vx);
		NEGF_ASSERT(edges_near_x.size()==1 || edges_near_x.size()==2, "unexpected number of edges near x.");
		for (uint ii = 0; ii < edges_near_x.size(); ii++) 
		{
			Vertex * vy = (edges_near_x[ii]->get_lower_vertex()==vx) ? edges_near_x[ii]->get_upper_vertex() : edges_near_x[ii]->get_lower_vertex();
			uint yy = vy->get_index_global() + 1;
			double lxy =  edges_near_x[ii]->get_length();
			
			// V(x) = sum_i Vi Ni(x)
			//this->potential_augmented_hamiltonian(xx,xx) += node_potential[vx->get_index_global()] * lxy;
			//this->potential_augmented_hamiltonian(xx,yy) += node_potential[vx->get_index_global()] * lxy / 6.0;
			
			// element-wise constant V(x)
			double V_el = 0.5 * (node_potential[vx->get_index_global()] + node_potential[vy->get_index_global()]);
			
			// CONVERSION TO TDKP UNITS: V_el is energy, must be eV; lxy is length, must be nm
			V_el = 1.0/constants::SIec * V_el/constants::convert_from_SI(units::energy, 1.0);
			lxy = 1e9                  * lxy/constants::convert_from_SI(units::length, 1.0);
			
			this->potential_augmented_hamiltonian(xx,xx) += V_el * lxy / 3.0;
			this->potential_augmented_hamiltonian(xx,yy) += V_el * lxy / 6.0;
		}
	}
);}


// H(k) = H(0) + hbar^2 k^2 / (2m) * M
void InterfaceEffMass::assemble_hamiltonian(const double& k_transversal_nm)
{STACK_TRACE(
	double k_transversal = constants::convert_from_SI(units::density_1d, 1e9*k_transversal_nm);
	uint Nvert = xspace->get_num_vertices();
	const double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
	
	this->potential_and_k_augmented_hamiltonian = this->potential_augmented_hamiltonian;
	for (uint xx=1; xx<=Nvert; xx++) 
	{
		Vertex * vx = xspace->get_vertex(xx-1);
		const vector<Element *> elems_near_x = xspace->get_elems_near(vx);
		NEGF_ASSERT(elems_near_x.size()==1 || elems_near_x.size()==2, "unexpected number of elements near x.");
		for (uint ii = 0; ii < elems_near_x.size(); ii++) 
		{
			NEGF_ASSERT(elems_near_x[ii]->get_num_edges()==1, "expected exactly 1 edge.");
			Edge * edge = elems_near_x[ii]->get_edge(0);
			double lxy =  edge->get_length();
			Vertex * vy = (edge->get_lower_vertex()==vx) ? edge->get_upper_vertex() : edge->get_lower_vertex();
			uint yy = vy->get_index_global() + 1;
			
			const PropertyContainer<double> * mat = elems_near_x[ii]->get_region()->get_material();
			double m   = constants::convert_from_SI(units::mass, mat->get("electron_effective_mass") * constants::SIm0);
			double fac = hbar*hbar * k_transversal*k_transversal / (2.0 * m);
			
			// CONVERSION TO TDKP UNITS: fac is energy, must be eV; lxy is length, must be nm
			fac = 1.0/constants::SIec * fac/constants::convert_from_SI(units::energy, 1.0);
			lxy = 1e9                 * lxy/constants::convert_from_SI(units::length, 1.0);
						
			 // fac is element-wise constant
			this->potential_and_k_augmented_hamiltonian(xx,xx) += fac * lxy / 3.0;
			this->potential_and_k_augmented_hamiltonian(xx,yy) += fac * lxy / 6.0;
		}
	}
);}


void InterfaceEffMass::set_strain(StrainPolarization * strainpol)
{STACK_TRACE(
    NEGF_EXCEPTION("InterfaceEffMass abort: NYI.");
);}

