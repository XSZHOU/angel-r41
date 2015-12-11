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
#ifndef PPEMISSION_H_
#define PPEMISSION_H_

#include "all.h"

#include "Geometry.h"
#include "Kspace.h"
#include "Energies.h"
#include "Options.h"
#include "NEGFObject.h"
#include "Hamiltonian.h"
#include "Overlap.h"
#include "SelfEnergy.h"
#include "GreenFunctions.h"

#include "SEPhotonSpontaneous.h"

namespace negf {
	
	/** Luminescence spectrum from spontaneous photon emission
	 *  Please refer to the report of S.Steiger for the meaning of the quantities */
	class PPEmission {
	public:
		
		PPEmission(const Overlap * ov_,
					const Geometry * xspace_, 
					const Kspace * kspace_, 
					const Energies * energies_, 
					const Options * options_,
					const Hamiltonian * ham_,
					const GreenFunctions * gf_,
					const SEPhotonSpontaneous * spont_);
		
		~PPEmission() {}
		
		void calculate();
		
		void write_recombination_to_file(const char * filename);
		void write_spectrum_to_file(const char * filename);
								
	protected:
		
		// helper functions
		void determine_mpi_stuff();
		void communicate_As();
		void compute_QQ_Jphot_spectrum();
		void determine_needed_processes(
					vector< vector<int> > & processes_needed,
					vector<int> & pp_E0_idx,
					vector<int> & pp_E1_idx,
					vector<int> & pp_E0_plus_Emin_idx,
					vector<int> & pp_E1_plus_Emax_idx) const;		
		
		// -----------------------
		// class variables
		// -----------------------
		const Overlap        * ov;
		const Geometry       * xspace;
		const Kspace         * kspace;
		const Energies       * energies;
		const Options  		 * options;
		const Hamiltonian    * ham;
		const GreenFunctions * gf;
		const SEPhotonSpontaneous * spont;
		uint 			Nx;
		uint 			Nn;
		uint 			Nk;
		uint 			NE;
		uint 			myNE;
		
		// discretization of spectrum
		vector<double> 	hw_grid;
		uint 			Nhw;
		double 			hw_min;
		double 			hw_max;
		
		// helper quantities
		vector< vector< vector< cplx > > >  PP;    // PP[x][E][E']
		vector< vector< vector< cplx > > >  QQ;    // QQ[x][E][hw]
		vector< vector< vector< cplx > > >  QQ_total; // only used by master process
		vector< vector< vector< cplx > > >  PP_VB;    // PP[x][E][E']
		vector< vector< vector< cplx > > >  QQ_VB;    // QQ[x][E][hw]
		vector< vector< vector< cplx > > >  QQ_total_VB; // only used by master process
		
		// photocurrent (only master thread)
		GEMatrix 		jphot; // jphot[x][hw]
		GEMatrix 		jphot_VB;
		
		// emission spectrum (only master thread)
		DGEMatrix 		spectrum; // spectrum[a][hw]
		DGEMatrix 		spectrum_VB;
		double 			area;
		vector<double>  output_spectrum;
		vector<double>  output_spectrum_VB;
		
		// MPI related
		uint 			E0_idx;	// boundary indices of own energies
		uint 			E1_idx;	
		uint 			E0_plus_Emin_idx;	// needed for determining which energies of other processes are needed
		uint 			E1_plus_Emax_idx;
		
	};
	
	
} // end of namespace

#endif /*PPEMISSION_H_*/
