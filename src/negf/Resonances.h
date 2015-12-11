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
#ifndef RESONANCES_H_
#define RESONANCES_H_

#include "all.h"

#include "InnerLoop.h"

#ifndef NOTDKP
#include "TdkpInfoDesk.h"
#include "tdkp/interface/Interface.h"
#include "tdkp/interface/InterfaceNEGFWell.h"
#endif

namespace negf {
	
	
	/** Determine bound-state resonances of a well in the middle of the structure
	 *  First, simple strategy: a single resonance at the maximum of LDOS(k=0)
	 *  Second strategy: use PML's and TDKP to find quasi-bound energies and widths 
	 *  The PML functionality is only provided when compiling with TDKP. */
	class Resonances {
	public:
	
		Resonances(InnerLoop * inner_) throw (Exception *);
		~Resonances() {}
		
		void determine_new_energy_grid(bool divergence_avoiding_only, const vector<double> & fermilevels_) throw (Exception *);  // for energy grid
		
	protected:

		void 			determine_divergences();
		void 			determine_ldos_resonances();
		void 			determine_all_ldos_resonances();
		void 			determine_pml_resonances();
		
		double 			monotonic_function(const double & E); 			 //!< for energy grid; only used by master thread
		double 			monotonic_function_derivative(const double & E); //!< for energy grid; only used by master thread

		vector<double>  fermilevels;    //!< there will be grid refinement around these energies
        vector<double>  divergences;    //!< only used by master thread. stores energies which are to be avoided by the energy grid.
		vector<double>  resonances;		//!< only used by master thread. stores energies around which a refinement is made
		vector<double>  broadenings;	//!< only used by master thread. stores the sharpness of the refinement function for a resonance

		// piecewise-linear monotonic function stuff
		const bool 		linear;			//!< determines whether refinement around Fermilevels is linear or of (1-f)-type
		const double 	lin_ref_num;	//!< number determining the width of refinement around FLs
		const double 	lin_delta;		//!< unitless number to determine sharpness of transition linear - flat
		double 			ca;
		double 			cb;
		double 			cc;
		double 			cd;
		double 			ce;
		double 			cf;
		
		// LDOS resonance finder (for all k) neglects resonances that are outside [EF_min-n*kT, EF_max+n*kT]
		static const double reson_k_neglect_kT = 20.0; // was 10
		
		
		Energies 		* energies;
		const Options	* options;
		const Geometry	* xspace;
		const Kspace 	* kspace;
		GreenFunctions  * gf;
		SelfEnergies 	* se;
		Hamiltonian 	* ham;
		PostProcessing 	* pp;
		
#ifndef NOTDKP
		TdkpInfoDesk 			  * infodesk;
		tdkp::InterfaceWellRadial * pml;
#endif 
		Geometry * pmlgrid;
	};
	
} // end namespace negf


#endif /*RESONANCES_H_*/
