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
#ifndef STRAINPOLARIZATION_H_NEGF
#define STRAINPOLARIZATION_H_NEGF

#include "all.h"

#include "Geometry.h"
#include "MaterialDatabase.h"


using namespace std;

namespace negf {

	/** Sets up planar strain and polarization.
	 *  Basic assumption: growth is ALWAYS along [001]
	 *  The material at the contacts is taken to be the intrinsic material (must be the same for both contacts).
	 *  --> The strain tensor is diagonal and the biaxial strain is given gy exx=eyy=(a-ai)/ai, ezz=-2*C13/C33*exx.
	 *      ai is the lattice constant of the intrinsic material. 
	 *  --> The polarization of a material is given by:
	 *      wurtzite:   Pi = e31*(exx+eyy) + e33*ezz + psp
	 *      zincblende: Pi = 0                               (would be 2*e14*exy, but exy=0)
	 *  --> The sheet charge at material interfaces is given by rho(x) = -(P1-P2)*delta(x)
	 *  --> For the FEM-Poisson equation we need int rho(x) Ni(x) = -(P1-P2).
	 *
	 * The polarization, and also the sheet charges, which are calculated from the MaterialDatabase
	 * may be scaled by a factor "pol_decreaser".
	 *  */
	class StrainPolarization
	{
		public:
			StrainPolarization(const Geometry * grid_, 
								const MaterialDatabase * db_, 
								const double pol_decreaser = 1.0);
			~StrainPolarization() {}

			// access functions
			double get_exx_eyy				(const Region * reg) const { return exx[reg->get_index()]; }
			double get_ezz					(const Region * reg) const { return ezz[reg->get_index()]; }
			double get_polarization			(const Region * reg) const { return pol[reg->get_index()]; }
			double get_sheet_charge			(const Vertex * v)   const { return sheet_charge[v->get_index_global()]; }
			double get_sheet_charge_between (const Region * reg1, const Region * reg2) const;

		protected:
			const Geometry * 		 	grid;
			const MaterialDatabase * 	db;

			vector<double> 				exx;			// for each region
			vector<double> 				ezz;			// for each region
			vector<double> 				pol;			// for each region, =piezo+psp
			vector<double> 				piezo;			// for each region
			vector<double> 				sheet_charge;	// for each vertex
		
			vector<bool> 				is_wurtzite;	// all regionwise
			vector<double> 				alattice;
			vector<double> 				C11;			// zincblende material parameters
			vector<double> 				C12;
			vector<double> 				C13;			// wurtzite material parameters
			vector<double> 				C33;
			vector<double> 				e31;
			vector<double> 				e33;
			vector<double> 				psp;

	};

}	// end namespace gain

#endif /*STRAINPOLARIZATION_H_NEGF*/
