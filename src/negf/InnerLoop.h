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
#ifndef INNERLOOP_H_
#define INNERLOOP_H_

#include "all.h"

#include "Hamiltonian.h"
#include "ContactFermilevel.h"
#include "Overlap.h"
#include "NEGFObject.h"
#include "GreenFunctions.h"
#include "SelfEnergies.h"
#include "SEContacts.h"
#include "PostProcessing.h"
#include "MaterialDatabase.h"

namespace negf {
	
	
	/** Reach self-consistency for a fixed electrostatic field 
	 *  The calculation of the GF is performed explicitly in this class, whereas
	 *  The SE calculation is outsourced to SelfEnergies */
	class InnerLoop {
	public:
	
		InnerLoop(Hamiltonian * ham_, Overlap * ov_, GreenFunctions * gf_, 
				SelfEnergies * se_, PostProcessing * pp_, MaterialDatabase * material_) throw (Exception *);
		~InnerLoop() {}
		
		// steering
		void 			set_max_inner_iterations(uint new_number);
	
		// the argument is used for adjusting the max. number of inner iterations and 
		// not terminating the program in the first outer iteration
		// return value marks whether convergeenc was reached
		bool perform(uint outer_iteration, const double err_crit) throw (Exception *);
		
		void initial_guess() throw(Exception);	// computes retarded and advanced GF
		
		Hamiltonian	   * get_hamiltonian()     const { return ham; }
		Overlap		   * get_overlap()         const { return ov; }
		GreenFunctions * get_green_functions() const { return gf; }
		SelfEnergies   * get_self_energies()   const { return se; }
		PostProcessing * get_post_processing() const { return pp; }
		ContactFermilevel * get_contact_0_fermilevel() const { return contact_0_fermi; }
		
		const Geometry * get_xspace() 		   const { return xspace; }
		const Kspace   * get_kspace() 		   const { return kspace; }
		const Options  * get_options() 		   const { return options; }
		Energies * 		 get_energies() 	   const { return energies; }
	
	protected:
		
		void 			iterate();
		bool 			converged(const double err_crit);
	
		void 			perform_frey_iteration();
		void 			broaden_contact_states();
				
		uint 			max_inner_iterations;
		Matd	 		old_spectral_edensity;
		Matd	 		old_spectral_hdensity;
		Matd	 		old_spectral_ecurrent;
		Matd	 		old_spectral_hcurrent;
		vector<double>	old_edensity;
		vector<double> 	old_hdensity;
		vector<double> 	old_ecurrent;
		vector<double> 	old_hcurrent;
	
		Hamiltonian	   * ham;
		ContactFermilevel * contact_0_fermi;
		Overlap	       * ov;
		GreenFunctions * gf;
		SelfEnergies   * se;
		PostProcessing * pp;
		MaterialDatabase * material;
		
		const Geometry * xspace;
		const Kspace   * kspace;
		const Options  * options;
		Energies * 		 energies;
		
		// objects for Frey method of finding self-energies in contact
		Geometry       * left_freyspace;
		Hamiltonian    * left_freyham;
		Overlap        * left_frey_overlap;
		GreenFunctions * left_frey_gf;
		SelfEnergies   * left_frey_energies;
		Geometry       * right_freyspace;
		Hamiltonian    * right_freyham;
		Overlap        * right_frey_overlap;
		GreenFunctions * right_frey_gf;
		SelfEnergies   * right_frey_energies;
		vector< vector< vector< cplx > > >  left_total_frey_broadening; // only used by master process
		vector< vector< vector< cplx > > > right_total_frey_broadening; // only used by master process
	};
	
} // end namespace negf
	
#endif /*INNERLOOP_H_*/
