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
#ifndef POISSON_H_NEGF
#define POISSON_H_NEGF

#include "all.h"

#include "Geometry.h"
#include "BoxMethod.h"
#include "ImplicitEquation.h"

using namespace std;

namespace negf {

	/** Poisson Equation for the electrostatic potential, defined on the grid vertices. <BR>
	 *  mandatory dependencies: dielectric constant, electron density, hole density, dopant density <BR>
	 *  optional dependencies:  sources of electron or hole density (e.g. interpolated low-dimensional density) <BR>
	 *  This is only the base class, independent of the discretization chosen.
	 */
	class Poisson : public ImplicitEquation
	{
		public:
			Poisson(const Geometry * const grid_);
			~Poisson() {}

			void set_neumann_electric_field(double field);  //!< set electric field for Neumann contacts (default: 0) */
			
			const Geometry  * get_grid() const { return grid; }

			// functions to allocate dependencies
			void add_source(Equation * charge_eqn); 		//!< adds a density source. eqn must be on the same grid
			void allocate_epsilon(Equation * epsilon_eqn); 	//!< adds the dielectric constant and computes the scaling factor
			void allocate_doping(Equation * doping_eqn);	//!< adds the net ionized dopant density (same grid!)
			void allocate_owncharge(Equation * edensity_eqn, Equation * hdensity_eqn); //!< adds electron and hole charge
			
			void initial_guess(uint tstamp);				//!< compute initial guess
			
		protected:
			// virtual functions of the ImplicitEquation class
			virtual double 	get_newton_function(uint line) const = 0;
			virtual void 	get_newton_direct_derivatives(uint ii, const Equation * eqn, uint & nonzeros, 
												  uint indices[], double coeff[]) const = 0;
			
			const Geometry  * const grid;	//!< grid on whose vertices the electrostatic potential is defined

			double 			neumann_field;	//!< electric field for Neumann contacts (default: 0)
												  
			// aliases for dependencies
			Equation * epsilon;		//!< dielectric constant
			Equation * edensity;	//!< electron density
			Equation * hdensity;	//!< hole density
			Equation * doping;		//!< doping density (>0 --> n-doping)
	};

}	// end of namespace

#endif /*POISSON_H_NEGF*/
