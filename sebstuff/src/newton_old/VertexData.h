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
#ifndef VERTEXDATA_H_NEGF
#define VERTEXDATA_H_NEGF


#include "all.h"
#include "Exception.h"
#include "Logger.h"

#include "Geometry.h"
#include "ExplicitEquation.h"


using namespace std;

namespace negf {

	/** Implements a container for data defined on vertices. The values are passed in the constructor. The timestamp is set to 99999. */
	class VertexData : public ExplicitEquation
	{
		public:
			VertexData(const Geometry * const grid_, quantities::PhysicalQuantity type_, vector<double> & values);
			VertexData(const Geometry * const grid_, quantities::PhysicalQuantity type_);
			~VertexData() {}

			// overwritten functions from class Equation
			void 			newton_update(const double * update, uint tstamp) {}

			// further functions
			const Geometry * const get_grid() const { return this->grid; }

		protected:
			const Geometry * const grid;

			// virtual functions from the ExplicitEquation class
			virtual void 	direct_derivatives(uint line, const Equation * eqn, uint & nonzeros, 
										uint indices [], double coeff []) const; 
			virtual double compute_value(uint ii) const
				 { NEGF_EXCEPTION("The class is not designed for this."); }

	};

}	// end namespace 

#endif /*VERTEXDATA_H_NEGF*/
