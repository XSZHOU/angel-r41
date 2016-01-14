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
#ifndef SEMICLASSICALDENSITY_H_NEGF
#define SEMICLASSICALDENSITY_H_NEGF

#include "all.h"

#include "ExplicitEquation.h"

using namespace std;

namespace negf {

	/** Semiclassical density which is fixed and cannot be computed itself (use set_values() instead)
	 *  but which possesses semiclassical derivatives w.r.t. potential.
	 * */
	class SemiclassicalDensity : public ExplicitEquation
	{
		public:
			SemiclassicalDensity(Equation * potential_, 
								Equation * bandedge_,
								Equation * effdos_, 
								Equation * temperature_);
			~SemiclassicalDensity() {};

		protected:
			
			/** virtual functions from the ExplicitEquation class */
			virtual double  compute_value(uint line) const;
			virtual void 	direct_derivatives(uint line, const Equation * eqn, uint & nonzeros, 
										uint indices [], double coeffs []) const;
			
			/** aliases for dependencies */
			Equation * potential;
			Equation * bandedge;
			Equation * effdos;
			Equation * temperature;
	};

}	// end of namespace

#endif /*SEMICLASSICALDENSITY_H_NEGF*/
