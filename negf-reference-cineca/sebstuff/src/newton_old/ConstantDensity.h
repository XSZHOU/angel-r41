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
#ifndef CONSTANTDENSITY_H_NEGF
#define CONSTANTDENSITY_H_NEGF

#include "all.h"

#include "Geometry.h"
#include "ExplicitEquation.h"

using namespace std;

namespace negf {

	/** Implements a constant density (electrons or holes), defined on the vertices of a grid.
	 *  The Equation timestamp is set to 99999. */
	class ConstantDensity : public ExplicitEquation
	{
	public:

		ConstantDensity(quantities::PhysicalQuantity type_, const Geometry * grid_, double density_);
		virtual ~ConstantDensity() {}

		// additional functions
		const Geometry * get_grid() const { return grid; }


	protected:
	
		// virtual functions of the ExplicitEquation class
		virtual void   direct_derivatives(uint line, const Equation * eqn, uint & nonzeros, 
							uint indices [], double coeff []) const;
		virtual double compute_value(uint line) const
				 { NEGF_EXCEPTION("The class is not designed for this."); }

		const Geometry * grid;		//!< how much density
		double     		density; 	//!< the globally defined, constant density
	};

} // end of namespace

#endif /*CONSTANTDENSITY_H_NEGF*/
