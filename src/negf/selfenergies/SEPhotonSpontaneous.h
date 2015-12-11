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
#ifndef SEPHOTONSPONTANEOUS_H_
#define SEPHOTONSPONTANEOUS_H_

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

#include "TdkpInfoDesk.h"	 // for cbedge...
#include "CSRMatrix.h"

namespace negf {
	
	/** Photon self-energy from spontaneous emission */
	class SEPhotonSpontaneous : public SelfEnergy {
	public:
		
		SEPhotonSpontaneous(const Overlap * ov_,
					const Geometry * xspace_, 
					const Kspace * kspace_, 
					const Energies * energies_, 
					const Options * options_,
					const GreenFunctions * gf_,
					const MaterialDatabase * db) throw (Exception *);
		
		~SEPhotonSpontaneous();
		
		// inherited from SelfEnergy class
		void calculate() throw (Exception *);
		
		void set_scaling(double new_scaling) throw (Exception *);
		
		double get_min_emission_energy() const { return this->hwmin; }
		double get_max_emission_energy() const { return this->hwmax; }
		// note that there is no fixed photon energy grid in this class
		// hw = E1-E2 is flexible for available energies hwmin<=E1-E2<=hwmax
		
		// will be used for computation of luminescence spectrum
		// stored only for own energies!
		const vector<unsigned long>   & get_AL_real_comp_size()  const { return this->AL_real_comp_size; }
		const vector<unsigned long>   & get_AL_imag_comp_size()  const { return this->AL_imag_comp_size; }
		const vector<unsigned char *> & get_AL_real_compressed() const { return this->AL_real_compressed; }
		const vector<unsigned char *> & get_AL_imag_compressed() const { return this->AL_imag_compressed; }
		const vector<unsigned long>   & get_AG_real_comp_size()  const { return this->AG_real_comp_size; }
		const vector<unsigned long>   & get_AG_imag_comp_size()  const { return this->AG_imag_comp_size; }
		const vector<unsigned char *> & get_AG_real_compressed() const { return this->AG_real_compressed; }
		const vector<unsigned char *> & get_AG_imag_compressed() const { return this->AG_imag_compressed; }
		const vector<double>		  & get_AL_norm(uint ee2) const;
		const vector<double>		  & get_AG_norm(uint ee2) const;
		
		
		// also used in class Luminescence
		void determine_needed_processes(
					vector< vector<int> > & processes_needed,
					vector<int> & pp_E0_minus_hwmax_idx,
					vector<int> & pp_E1_minus_hwmin_idx,
					vector<int> & pp_E0_idx,
					vector<int> & pp_E1_idx,
					vector<int> & pp_E0_plus_hwmin_idx,
					vector<int> & pp_E1_plus_hwmax_idx) const;	
		
		void output_debug_info();	
	protected:
		
		// helper functions
		void determine_mpi_stuff();
		void communicate_As();
		void calculate_lesser_greater();	
		void calculate_retarded();
		void get_CB_VB_norms(const Matc & A, vector<double> & result); // for debug only
		
		// -----------------------
		// class variables
		// -----------------------
		const Overlap        * ov;
		const GreenFunctions * gf;
		
		double hwmin;	// maximum emission energy
		double hwmax;	// minimum emission energy
				
		vector< vector<SEMat> > AL;	// helper array; AL[i,j] stores GL at some other processes energy, k-point j
		vector< vector<SEMat> > AG;
		vector< vector<SEMat> > AR;
		vector< vector<double> > AL_norm;
		vector< vector<double> > AG_norm;
		vector< vector<double> > AR_norm;
		
		vector< vector< CSRMatrix<double> *> > Mt;	// M_abcd stores M_ik * M_lj * \int d\phi \int d\theta ...
		
		double prefactor;	// stores for every band D^2*kT/(rho*c^2*a) with the corresponding deformation potential
		double scaling;		// multiplier that gives possibility to scale scattering down
		
		const bool complicated_retarded;
		
		const bool security_checking;
		
		// MPI related
		uint E0_idx;	// boundary indices of own energies
		uint E1_idx;	
		uint E0_minus_hwmax_idx;	// needed for determining which energies of other processes are needed
		uint E1_minus_hwmin_idx;
		uint E0_plus_hwmin_idx;
		uint E1_plus_hwmax_idx;
		uint nE_below;
		uint nE_self;
		uint nE_above;

		// temporary array used during decompression
		// always have the same size --> need to be allocated only once
		unsigned char * compressed_array_real;	
		unsigned char * compressed_array_imag;
		unsigned char * decompress_array_real;	
		unsigned char * decompress_array_imag;
		
		// zipped versions of own AL, AG
		// only for own energies!
		vector<double *> AL_real_data;	// AL_real_data[ee2] will store copies of the complex AL[ee2] (for all k) split into real and imag parts
		vector<double *> AG_real_data;
		vector<double *> AR_real_data;
		vector<double *> AL_imag_data;
		vector<double *> AG_imag_data;
		vector<double *> AR_imag_data;
		vector<unsigned char *> AL_real_char;		// will store pointer to same memory as AL_...._data, but different type
		vector<unsigned char *> AG_real_char;
		vector<unsigned char *> AR_real_char;
		vector<unsigned char *> AL_imag_char;
		vector<unsigned char *> AG_imag_char;
		vector<unsigned char *> AR_imag_char;
		vector<unsigned long> num_chars;			// will store how many chars the compressed array of complex antihermitian matrices corresponds to
		vector<unsigned long> num_chars_AR;			// will store how many chars the compressed array of complex matrices corresponds to
		vector<unsigned long> AL_real_comp_size;	// will store the size of the zipped data
		vector<unsigned long> AG_real_comp_size;
		vector<unsigned long> AR_real_comp_size;
		vector<unsigned long> AL_imag_comp_size;
		vector<unsigned long> AG_imag_comp_size;
		vector<unsigned long> AR_imag_comp_size;
		vector<unsigned char *> AL_real_compressed;	// will store pointers to the zipped data
		vector<unsigned char *> AG_real_compressed;
		vector<unsigned char *> AR_real_compressed;
		vector<unsigned char *> AL_imag_compressed;
		vector<unsigned char *> AG_imag_compressed;
		vector<unsigned char *> AR_imag_compressed;

	};
	
	
} // end of namespace

#endif /*SEPHOTONSPONTANEOUS_H_*/
