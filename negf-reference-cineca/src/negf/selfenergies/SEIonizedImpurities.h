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
#ifndef SEIONIZEDIMPURITIES_H_
#define SEIONIZEDIMPURITIES_H_

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
//#include "Equation.h"

namespace negf {
	
	/** Ionized impurity self-energy for planar structures with radial-approximation quantities (G(k)=G(|k|)) */
	class SEIonizedImpurities : public SelfEnergy {
	public:
		
		SEIonizedImpurities(const Overlap * ov_,
					const Geometry * xspace_, 
					const Kspace * kspace_, 
					const Energies * energies_, 
					const Options * options_,
					const GreenFunctions * gf_) throw (Exception *);
		
		~SEIonizedImpurities() {}
		
		// needs to be called before anything
		/*void set_up(const Equation * doping,
					const Equation * dielectric,
					const Equation * Nc,
					const Equation * Nv) throw (Exception *);*/
        void set_up(const vector<double> & doping,
                    const vector<double> & dielectric,
                    const vector<double> & Nc,
                    const vector<double> & Nv) throw (Exception *);
		
		// inherited from SelfEnergy class
		void calculate();
		void set_scaling(double new_scaling);
				
	protected:
		
		// -----------------------
		// helper functions
		// -----------------------
		
		void calculate_lesser_greater();
		void calculate_retarded();
		
		void calculate_F(const vector<double> & doping, const vector<double> & dielectric, const vector<double> & nu);
		uint find_F_index(uint xx, uint yy) const;
		
		// -----------------------
		// class variables
		// -----------------------
		const Overlap        * ov;
		const GreenFunctions * gf;
		
		bool prepared;
		
		double scaling;				// multiplier that gives possibility to scale scattering down
		
		// F-related
		uint Nl;
		vector< vector< vector<double> > > F; // form factor
		double Fmax;
		const double F_neglect;		// F < F_neglect * Fmax --> contribution is not added	
	};
	
	
} // end of namespace

#endif /*SEIONIZEDIMPURITIES_H_*/
