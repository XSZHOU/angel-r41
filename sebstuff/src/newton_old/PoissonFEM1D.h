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
#ifndef POISSONFEM1D_H_NEGF
#define POISSONFEM1D_H_NEGF

#include "all.h"

#include "Geometry.h"
#include "ImplicitEquation.h"
#include "Poisson.h"
#include "MaterialDatabase.h"

using namespace std;

namespace negf {

	/** 1D Poisson Equation using FEM discretization */
	class PoissonFEM1D : public Poisson {
		public:
			PoissonFEM1D(
				const Geometry * const grid_,
				const MaterialDatabase& material_database, 
				const char* gridfile, 
				const vector<double> & static_rhs_			//!< contains int Ni(x) rho(x), rho are the polarization sheet charges
			);
			~PoissonFEM1D() {}	
			
			void set_piezo_decreaser(double new_decreaser_value);
			
		protected:
			/** virtual functions of the ImplicitEquation class */
			double 	get_newton_function(uint line) const;
			void 	get_newton_direct_derivatives(uint ii, const Equation * eqn, uint & nonzeros, 
												  uint indices[], double coeff[]) const;
		private:
			
			vector<double> stiffness_matrix;	//!< epsilon grad Ni * grad Nj in csr
			vector<double> mass_matrix;			//!< Ni * Nj in csr
			vector<int> icol;					//!< CSR pattern for mass and stiffness matrix
			vector<int> prow;					//!< CSR pattern for mass and stiffness matrix
			vector<double> static_rhs;			//!< static rhs contribution (int Ni * rho_static)
			
			double piezo_decreaser;				//!< factor that multiplies the static rhs, to decrease sheet densities
	};

}	// end of namespace

#endif /*POISSONFEM1D_H_NEGF*/
