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
#ifndef TRANSMISSION_H_
#define TRANSMISSION_H_

#include "all.h"

#include "Geometry.h"
#include "Kspace.h"
#include "Energies.h"
#include "NEGFObject.h"
#include "SelfEnergy.h"
#include "GreenFunctions.h"

namespace negf {
	
	/** Transmission spectrum T(E) = Tr(Gam2*GR*Gam1*GA)
	 *  Gam = i(SR-SA) = i(SG-SL)
	 *  Gam2 is the SE of the left contact, Gam1 the SE of the right contact 
	 *  note: only k=0 should be taken since the lateral (k-)energy just shifts everything up and down! */
	class Transmission {
	public:
		
		Transmission(const Geometry * xspace_, 
					const Kspace * kspace_, 
					const Energies * energies_, 
					const GreenFunctions * gf_,
					const SelfEnergy * se_contact_) throw (Exception *);
		
		~Transmission() {}
		
		void calculate() throw (Exception *);
		
		void write_to_file(const char * filename) throw (Exception *);
								
	protected:
		
		const Geometry       *  xspace;
		const Kspace         *  kspace;
		const Energies       *  energies;
		const GreenFunctions *  gf;
		const SelfEnergy     *  se_contact;
		uint 					Nx;
		uint 					NxNn;
		uint 					Nk;
		uint 					NE;
		uint 					myNE;
		
		// emission spectrum (only master thread)
		vector<double> transmission_spectrum;		// all k's
		vector<double> transmission_spectrum_k0;	// k=0 only
	};
	
} // end of namespace

#endif /*TRANSMISSION_H_*/
