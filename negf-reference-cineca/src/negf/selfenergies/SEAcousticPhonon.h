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
#ifndef SEACOUSTICPHONON_H_
#define SEACOUSTICPHONON_H_

#include "all.h"

#include "PropertyContainer.h"
#include "MaterialDatabase.h"
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
	
	/** Acoustic phonon self-energy for planar structures with radial-approximation quantities (G(k)=G(|k|))
	 *  Exactly formula (A13) in Lake, Klimeck, Bowen, Jovanovic, J.Appl.Phys. 81, 7845 (1997) */
	class SEAcousticPhonon : public SelfEnergy {
	public:
		
		SEAcousticPhonon(const Overlap * ov_,
					const Geometry * xspace_, 
					const Kspace * kspace_, 
					const Energies * energies_, 
					const Options * options_,
					const GreenFunctions * gf_,
					const MaterialDatabase * db);
		
		~SEAcousticPhonon() {}
		
		// inherited from SelfEnergy class
		void calculate();
		
		void set_scaling(double new_scaling);
						
	protected:
		
		// helper function
		uint find_delta_index(uint xx, uint yy) const;
		
		// -----------------------
		// class variables
		// -----------------------
		const Overlap        * ov;
		const GreenFunctions * gf;
		
		vector<double> prefactor;	// stores for every band D^2*kT/(rho*c^2*a) with the corresponding deformation potential
		double scaling;				// multiplier that gives possibility to scale scattering down
		
		const bool new_version;
		uint Nl;
		vector<double> delta_ij;
	};
	
	
} // end of namespace

#endif /*SEACOUSTICPHONON_H_*/
