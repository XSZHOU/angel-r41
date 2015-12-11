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
#ifndef LUMINESCENCE_H_
#define LUMINESCENCE_H_

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

	enum photon_current_mode {
		lake,
		galerpin
	};
	
	/** Luminescence spectrum from spontaneous photon emission
	 *  Please refer to the report of S.Steiger for the meaning of the quantities */
	class Luminescence {
	public:
		
		Luminescence(const Overlap * ov_,
					const Geometry * xspace_, 
					const Kspace * kspace_, 
					const Energies * energies_, 
					const Options * options_,
					const Hamiltonian * ham_,
					const GreenFunctions * gf_,
					const SEPhotonSpontaneous * spont_,
					photon_current_mode mode_) throw (Exception *);
		
		~Luminescence() {}
		
		void calculate() throw (Exception *);
		
		void write_recombination_to_file(const char * filename) throw (Exception *);
		void write_spectrum_to_file(const char * filename) throw (Exception *);
		
		void snapshot(double voltage) throw (Exception *);
		void write_power_to_file(const char * filename) throw (Exception *);
								
	protected:
		
		// helper functions
		void 	determine_mpi_stuff();
		void 	communicate_As();
		void 	compute_PP(const uint ee, const vector<SEMat> & ALGmat, const vector<double> & ALGnorm, 
							const vector< vector<Matc> > & HminusEM, bool below);
		void 	compute_QQ_Jphot_spectrum();
		
		void 	determine_hw_grid();
		double 	monotonic_function(const double & hw); 			  // for energy grid; only used by master thread
		double 	get_monotonic_function_contribution(const double & hw, const double & Emin, const double & Emax);
		double 	monotonic_function_derivative(const double & hw); // for energy grid; only used by master thread
		double  get_monotonic_function_contribution_derivative(const double & hw, const double & Emin, const double & Emax);
		
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
		uint 			NxNn;
		uint 			Nk;
		uint 			NE;
		uint 			myNE;
		
		photon_current_mode mode;
		
		// discretization of spectrum
		vector<double> 	hw_grid;
		uint 			Nhw;
		double 			hw_min;
		double 			hw_max;
		double 			Egap_min;
		
		// helper quantities
		vector< vector< vector< cplx > > >  PPplus;     // [x][E][E'], E+hwmax >= E' >= E+hwmin
		vector< vector< vector< cplx > > >  PPminus;    // [x][E][E'], E-hwmax <= E' <= E-hwmin
		vector< vector< vector< cplx > > >  QQ;    		// QQ[x][E][hw]
		vector< vector< vector< cplx > > >  QQ_total; 	// only used by master process
		vector< vector< vector< cplx > > >  PPplus_VB;  
		vector< vector< vector< cplx > > >  PPminus_VB;  
		vector< vector< vector< cplx > > >  QQ_VB;   
		vector< vector< vector< cplx > > >  QQ_total_VB; 
				
		// emission spectrum (only master thread)
		double 			area;
		Matd	 		spectrum; 			// 1st column --> position, 2nd column --> spectrum of emission at that place
		Matd	 		spectrum_VB; 
		Matd            output_spectrum;	// 1st column --> emission energy, 2nd column --> spectrum
		Matd            output_spectrum_VB;
		double  		output_power;
		double  		output_power_VB;
		
		// snapshot arrays, store the voltage and the computed power for that voltage
		Matd			power_snapshots;
		
		// photocurrent (only master thread)
		Matc	 		jphot; 				// 1st column --> position, 2nd column --> radiative recombination at that place
		Matc	 		jphot_VB;
		
		// MPI related
		uint 			E0_idx;				// boundary indices of own energies
		uint 			E1_idx;	
		uint 			E0_plus_Emin_idx;	// needed for determining which energies of other processes are needed
		uint 			E1_plus_Emax_idx;
		uint 			E0_minus_Emax_idx;
		uint 			E1_minus_Emin_idx;
		
	};
	
	
} // end of namespace

#endif /*LUMINESCENCE_H_*/
