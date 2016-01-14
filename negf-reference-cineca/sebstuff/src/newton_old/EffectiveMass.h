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
#ifndef EFFECTIVEMASS_H_NEGF
#define EFFECTIVEMASS_H_NEGF

#include "all.h"

#include "Element.h"
#include "Geometry.h"
#include "BoxMethod.h"
#include "MaterialDatabase.h"
#include "ExplicitEquation.h"

using namespace std;

namespace negf {

	/** Stores the (electron or hole) effective mass either on vertices or elements. <BR>
	 *  Naturally the quantity is defined on elements, for vertices interpolation w/ voronoi mesures is done 
	 */
	class EffectiveMass : public ExplicitEquation
	{
		public:
			EffectiveMass(const Geometry * grid_, const BoxMethod * boxmethod_,
					 const MaterialDatabase * db_, 
					 quantities::PhysicalQuantity electrons_or_holes_, string verts_or_elems_);
			~EffectiveMass() {};

			double get_element_value(uint elem_idx) const;	//!< get the mass defined on an element
			
		protected:
			
			// virtual functions from the ExplicitEquation class
			virtual double  compute_value(uint line) const 
				{ NEGF_EXCEPTION("class is not designed for this."); }
			virtual void 	direct_derivatives(uint line, const Equation * eqn, uint & nonzeros, 
										uint indices [], double coeff []) const
				{ NEGF_EXCEPTION("class is not designed for this."); }

			const Geometry         * grid;
			const BoxMethod        * boxmethod;
			const MaterialDatabase * db;
			
			vector<double>				 elem_values;		 //!< element values are stored anyway
			quantities::PhysicalQuantity electrons_or_holes; //!< stores if its the conduction or valence band
			string 						 verts_or_elems;	 //!< stores if data is defined on vertices or elements

	};

}	// end of namespace

#endif /*EFFECTIVEMASS_H_NEGF*/
