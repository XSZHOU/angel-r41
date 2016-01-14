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
#ifndef SEBUETTIKER_H_
#define SEBUETTIKER_H_

#include "all.h"

#include "PropertyContainer.h"
#include "Geometry.h"
#include "Kspace.h"
#include "Energies.h"
#include "Options.h"
#include "NEGFObject.h"
#include "Hamiltonian.h"
#include "Overlap.h"
#include "SelfEnergy.h"
#include "GreenFunctions.h"

namespace negf {
	
	/** Self-energy associated to Buettiker probes */
	class SEBuettiker : public SelfEnergy {
	public:
		
		SEBuettiker(const Hamiltonian * ham_,
					const Overlap * ov_,
					const Geometry * xspace_, 
					const Kspace * kspace_, 
					const Energies * energies_, 
					const Options * options_,
					const GreenFunctions * gf_,
					const SelfEnergy * se_contact_,
					const double & energy_parameter_);
		
		~SEBuettiker() {}
		
		void 	set_energy_parameter(const double & new_parameter);
		double  get_energy_parameter() const { return this->energy_parameter; }
		
		// inherited from SelfEnergy class
		void calculate();
						
	protected:
		
		// -----------------------
		// helper functions
		// -----------------------
		void 				calculate_retarded();		// only needs to be executed once
		void 				calculate_lesser_greater();	// for FIXED fermilevels
		
		void 				interpolate_fermilevels();  // linear interpolation between contacts
		void 				compute_fermilevels();		// Newton iteration s.th. current is conserved
		bool 				newton_step();
		void				compute_newton_function(Vecd & result) const;
		void 				compute_newton_derivative(Matd & result) const;
				
		void 				compute_buettiker_current();
		
		// -----------------------
		// class variables
		// -----------------------
		const Hamiltonian * ham;
		const Overlap * 	ov;
		const GreenFunctions * gf;
		const SelfEnergy *  se_contact;	// contact self-energy
		double 				energy_parameter;
		const double 		kT; 		// temperature, in energy units 
		bool 				i_am_master;		// true for the master thread
		const uint 			how_many_neighbours; // Hamiltonian / overlap is assumed to be zero for |ii-jj| > how_many_neighbours
		vector<Contact *> 	fixed_fermilevels;	// ... for sites connected to contacts
		bool 				first_calculation;
		
		vector<double> 		site_fermilevels;		
	};
	
	
} // end of namespace

#endif /*SEBUETTIKER_H_*/
