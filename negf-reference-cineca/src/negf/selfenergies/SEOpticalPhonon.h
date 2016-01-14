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
#ifndef SEOPTICALPHONON_H_
#define SEOPTICALPHONON_H_

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
	
	/** Polar Optical Phonon self-energy for planar structures with radial-approximation quantities (G(k)=G(|k|)) */
	class SEOpticalPhonon : public SelfEnergy {
	public:
		
		SEOpticalPhonon(const Overlap * ov_,
					const Geometry * xspace_, 
					const Kspace * kspace_, 
					const Energies * energies_, 
					const Options * options_,
					const GreenFunctions * gf_,
					const MaterialDatabase * db);
		
		~SEOpticalPhonon() {}
		
		// inherited from SelfEnergy class
		void calculate();
		
		void determine_mpi_stuff(); // needs to be called whenever energy grid is changed - called automatically in calculate()
		
		void set_scaling(double new_scaling);

		void output_debug_info();

		void test_operation() throw (Exception *);
		
	protected:
		
		// -----------------------
		// helper functions
		// -----------------------
		void determine_interpolations(); // sets up downcoupling, upcoupling
		
		void calculate_lesser_greater();
		void communicate_As();
		
		void calculate_retarded();
		
		void calculate_F();
		void get_CB_VB_norms(const Matc & A, vector<double> & result);

		uint find_F_index(uint xx, uint yy) const;
		
		// -----------------------
		// class variables
		// -----------------------
		const Overlap        * ov;
		const GreenFunctions * gf;
		
		double lattice_constant;
		double hw;					// optical phonon energy
		double prefactor;			// Froehlich scattering strength prefactor
		double Nphonon;				// Bose-Einstein distribution with hw
		
		double scaling;				// multiplier that gives possibility to scale scattering down
		
		const bool complicated_retarded;
		const bool nemo_retarded;
		
		const bool security_checking;

		bool diagonals_only;		// set in constructor. true-> SE is diagonal IN BANDS (NOT in space!)
		
		const double q0;
		
		vector< vector< vector<double> > > F; // form factor
		double Fmax;
		const double F_neglect;		// F < F_neglect * Fmax --> contribution is not added
		uint Nl;
		
		// linear interpolation of E+-hw
		// Sigma(Ei) = A sum_j downcoupling(i,j) G(Ej) + B sum_k upcoupling(i,k) G(Ek), 
		// Ej~Ei-hw, Ek~Ei+hw
		// A is N+1 for SL and N for SG
		// B is N+1 for SG and N for SL
		Matd downcoupling;		
		Matd upcoupling;
		
		// MPI related
		uint E0_idx;
		uint Emax_idx;
		uint E0_minus_hw_idx;
		uint Emax_minus_hw_idx;
		uint E0_plus_hw_idx;
		uint Emax_plus_hw_idx;
		uint nE_below;
		uint nE_self;
		uint nE_above;
		vector< vector<SEMat> > AL;	// helper array; AL[i,j] stores AL at energy i, k-point j
		vector< vector<SEMat> > AG;
		vector< vector<SEMat> > AR;
		vector< uint > 			control;
		uint shift;
		
	};
	
	
} // end of namespace

#endif /*SEOPTICALPHONON_H_*/
