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
#ifndef CONTACTFERMILEVEL_H_
#define CONTACTFERMILEVEL_H_

#include "all.h"

#include "Geometry.h"
#include "Kspace.h"
#include "Options.h"
#include "Hamiltonian.h"
#include "Energies.h"
#include "Overlap.h"
#include "SEContacts.h"

#include "MaterialDatabase.h"

namespace negf {
	
	/** Computes and stores the Fermilevel at contact 0.
	    The Fermilevel EF is found numerically by populating the DOS of the contact with the Fermi function 1/(1+exp((E-EF)/kT)). 
		However, the contact DOS is constructed from the contact Hamiltonian assuming zero electrostatic potential!
		So the implicit assumpion here is that contact 0 has el.stat. potential 0. */
	class ContactFermilevel
	{
	public:
	
		ContactFermilevel(const Energies * energies_, Hamiltonian * ham_, const Overlap * ov_, SEContacts * se_contacts_) throw (Exception *);
		~ContactFermilevel() {}
		
		// obtain the fermilevel at contact 0 (where phi=0, calculated from the DOS)
		// note: doping>0 means n-doping 
		double get_contact_0_fermilevel() const { return contact_0_fermilevel; }
		void compute_contact_0_fermilevel(double doping) throw (Exception *);
		void calculate_dos_contact0();
		void calculate_contact0_dos_from_device_GR();
		
		const Geometry 						* get_xspace() 		const { return this->xspace; }
		const Kspace 						* get_kspace() 		const { return this->kspace; }
		const Options 						* get_options() 	const { return this->options; }	
		const Energies 						* get_energies() 	const { return this->energies; }	
		const Hamiltonian					* get_hamiltonian() const { return this->ham; }		
		
	protected:
		
		// helpers
		void calculate_contact_density_from_fermilevel(double EF, bool calculate_derivative, vector<double> & n_p) const;
		
		
		Hamiltonian * 					ham;
		const Geometry * 				xspace;
		const Kspace   * 				kspace;
		const Options  * 				options;
		const Energies * 				energies;
		const Overlap * 				ov;
		SEContacts * 					se_contacts;
		const uint 						Nx;
		const uint 						NxNn;
		
		double 							contact_0_fermilevel;
		
		vector<double> 					contact_0_dos_n;
		vector<double> 					contact_0_dos_p;
		
	};
	
} // end namespace negf

#endif /*CONTACTFERMILEVEL_H_*/
