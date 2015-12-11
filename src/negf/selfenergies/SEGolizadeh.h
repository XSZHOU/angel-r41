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
#ifndef SEGOLIZADEH_H_
#define SEGOLIZADEH_H_

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
	
	/** Self-energy associated to Golizadeh probes:
	 *  SigmaR(i,j) = D(i,j) GR(i,j), SigmaL(i,j) = D(i,j) GL(i,j)
	 *  ... where D(i,j) = delta_ij d_m for momentum relaxation or
	 *            D(i,j) = d_p          for phase relaxation */
	class SEGolizadeh : public SelfEnergy {
	public:
		
		SEGolizadeh(const Hamiltonian * ham_,
					const Overlap * ov_,
					const Geometry * xspace_, 
					const Kspace * kspace_, 
					const Energies * energies_, 
					const Options * options_,
					const GreenFunctions * gf_,
					const double & energy_parameter_,
					const bool momentum_relaxation_);
		
		~SEGolizadeh() {}
		
		void 	set_energy_parameter(const double & new_parameter);
		double  get_energy_parameter() const { return this->energy_parameter; }
		
		// inherited from SelfEnergy class
		void calculate();
						
	protected:
		
		// -----------------------
		// helper functions
		// -----------------------
		void 				calculate_retarded();	
		void 				calculate_lesser_greater();	
						
		// -----------------------
		// class variables
		// -----------------------
		const Hamiltonian * ham;
		const Overlap * 	ov;
		const GreenFunctions * gf;
		double 				energy_parameter;	
		const bool 			momentum_relaxation; // if true, D(i,j) = delta_ij d_m
												 // if false, D(i,j) = d_p
		const bool security_checking;
	};
	
	
} // end of namespace

#endif /*SEGOLIZADEH_H_*/
