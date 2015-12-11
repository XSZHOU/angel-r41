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
#ifndef ARBITRARYDATA_H_NEGF
#define ARBITRARYDATA_H_NEGF

#include "all.h"

#include "ExplicitEquation.h"

using namespace std;

namespace negf {

	/** Implements a container for arbitrary, constant data which is assigned in the constructor.
	 *  The timestamp of the Equation is set to 99999. */
	class ArbitraryData : public ExplicitEquation
	{
		public:
			ArbitraryData(quantities::PhysicalQuantity type_, const vector<double> & values);
			~ArbitraryData() {}

		protected:

			// virtual functions from the ExplicitEquation class
			virtual void direct_derivatives(uint line, const Equation * eqn, uint & nonzeros, 
										uint indices [], double coeff []) const
				 { NEGF_EXCEPTION("The class is not designed for this."); }
			virtual double compute_value(uint ii) const
				 { NEGF_EXCEPTION("The class is not designed for this."); }

	};

}	// end of namespace

#endif /*ARBITRARYDATA_H_NEGF*/
