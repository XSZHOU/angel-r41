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
#ifndef EFFECTIVEDOS_H_NEGF
#define EFFECTIVEDOS_H_NEGF

#include "all.h"

#include "ExplicitEquation.h"

using namespace std;

namespace negf {

	/** Calculates the Effective Density of States from the temperature, mass and dimensionality:
	 *  \f$ N_c(3D) = 2 * pow(m_c*kT/(2*\pi*\hbar^2), 1.5)  \f$ <BR>
	 *  \f$ N_c(2D) = m_c*kT/(2*\pi*\hbar^2)                \f$ <BR>
	 *  \f$ N_c(1D) = sqrt(m_c*kT/(2*\pi*\hbar^2))          \f$ <BR>
	 *  ... and similar for \f$N_v\f$ (with \f$m_v\f$)
	 * */
	class EffectiveDOS : public ExplicitEquation
	{
		public:
			EffectiveDOS(Equation * temperature_, Equation * mass_, uint dos_dim_,
							quantities::PhysicalQuantity electron_or_hole);
			~EffectiveDOS() {};

		protected:
			
			/** virtual functions from the ExplicitEquation class */
			virtual double  compute_value(uint line) const;
			virtual void 	direct_derivatives(uint line, const Equation * eqn, uint & nonzeros, 
										uint indices [], double coeff []) const;

			uint dos_dim; //!< stores dimensionality of the density of states
			
			// aliases for dependencies
			Equation * temperature;
			Equation * mass;

	};

}	// end of namespace

#endif /*EFFECTIVEDOS_H_NEGF*/
