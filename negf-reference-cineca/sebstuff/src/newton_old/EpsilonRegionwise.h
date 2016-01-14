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
#ifndef EPSILONREGIONWISE_H_NEGF
#define EPSILONREGIONWISE_H_NEGF

#include "all.h"

#include "Region.h"
#include "Geometry.h"
#include "MaterialDatabase.h"
#include "ExplicitEquation.h"

using namespace std;

namespace negf {

	/** Implements a regionwise constant dielectric constant. <BR>
	 *  The values are read in from the material database. epsilon is defined on elements */
	class EpsilonRegionwise : public ExplicitEquation
	{
		public:
			EpsilonRegionwise(const Geometry * grid_,
							  const MaterialDatabase * db_);
			~EpsilonRegionwise() {}

			// further functions
			const Geometry * get_grid() const { return this->grid; }
			double get_epsilon(Region * region) const; //!< epsilon is constant within Region because it's the same material everywhere

		protected:
			const Geometry * grid;
			const MaterialDatabase * db;
			map<Region *, double> region_values;

			// virtual functions from the ExplicitEquation class
			virtual void direct_derivatives(uint line, const Equation * eqn, uint & nonzeros, 
										uint indices [], double coeff []) const
				 { NEGF_EXCEPTION("The class is not designed for this."); }
			virtual double compute_value(uint line) const
				 { NEGF_EXCEPTION("The class is not designed for this."); }

	};

}	// end of namespace

#endif /*EPSILONREGIONWISE_H_NEGF*/
