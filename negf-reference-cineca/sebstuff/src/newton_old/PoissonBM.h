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
#ifndef POISSONBM_H_NEGF
#define POISSONBM_H_NEGF

#include "all.h"

#include "Geometry.h"
#include "BoxMethod.h"
#include "ImplicitEquation.h"
#include "Poisson.h"

using namespace std;

namespace negf {

	/** Poisson Equation using box method (finite volume) discretization */
	class PoissonBM : public Poisson
	{

		public:
			PoissonBM(const Geometry * grid_, const BoxMethod * const box_method_);
			~PoissonBM() {}

			const BoxMethod * get_boxmethod()      const { return box_method; }
			double            get_scaling_factor() const;
			
			
		protected:
			// virtual functions of the ImplicitEquation class
			double 	get_newton_function(uint line) const;
			void 	get_newton_direct_derivatives(uint ii, const Equation * eqn, uint & nonzeros, 
												  uint indices[], double coeff[]) const;
			
			const BoxMethod * box_method;

			double			scaling_factor;	//!< scaling factor for equation for Dirichlet boundary vertices, computed only once
	};

}	// end namespace

#endif /*POISSONBM_H_NEGF*/
