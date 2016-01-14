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
#ifndef BANDEDGE_H_NEGF
#define BANDEDGE_H_NEGF

#include "all.h"

#include "Element.h"
#include "Region.h"
#include "Geometry.h"
#include "BoxMethod.h"
#include "MaterialDatabase.h"
#include "PropertyContainer.h"
#include "TernaryPropertyContainer.h"
#include "ExplicitEquation.h"
#include "InputParser.h"

using namespace std;

namespace negf {

	/** Stores the (conduction or valence) band edge either on vertices or elements. <BR>
	 *  Note that this quantity is really the band edge, without potential! <BR>
	 *  Naturally the quantity is defined on elements, for vertices interpolation w/ voronoi mesures is done 
	 *  For vertices at the interface of a quantized and an unquantized region, only the unquantized region value is taken. <BR>
	 *  The bandedges are computed as: Ec = -electron_affinity, Ev = -electron_affinity-bandgap
	 */
	class Bandedge : public ExplicitEquation
	{
		public:
			Bandedge(const Geometry * grid_, const BoxMethod * boxmethod_,
					 const MaterialDatabase * db_, 
					 quantities::PhysicalQuantity electrons_or_holes_, 
					 string verts_or_elems_,
					 double temperature_);
			~Bandedge() {};
			
			double		get_element_value(uint elem_idx) const;	//!< uncorrected; used in Heterointerface class
			
			bool 		is_corrected() 							const { return corrected; }
			Equation * 	get_corrected_edge() 					const;
			void 		assign_corrected(Equation * corrected_eqn_);
			
			/** strain corrections */
			void 		assign_strain_correction(Equation * correction_);	//!< alters current_variable_values
			Equation *  get_strain_correction() 				const { return strain_correction; }
			bool 		is_strained()							const { return strained; } //!< true when assign_strain.. was called before
			
			//!< Computation of CB edge is nontrivial.
			static double get_cbedge(const PropertyContainer<double> * mat, const double & T/*emperature*/, 
									 const MaterialDatabase * const database);
			
		protected:
					
			/** virtual functions from the ExplicitEquation class */
			virtual double  compute_value(uint line) const 
			{NEGF_EXCEPTION("class is not designed for this."); }
			virtual void 	direct_derivatives(uint line, const Equation * eqn, uint & nonzeros, 
										uint indices [], double coeff []) const
				{ NEGF_EXCEPTION("class is not designed for this."); }

			const Geometry         * grid;					//!< the grid on whose vertices the band edge is defined
			const BoxMethod        * boxmethod;				//!< the box method, important for material boundaries
			const MaterialDatabase * db;					//!< edges are calculated w/ these parameters
			
			quantities::PhysicalQuantity electrons_or_holes; //!< stores if its the conduction or valence band
			string 						 verts_or_elems;	 //!< stores if data is defined on vertices or elements
			double 						 temperature;		 //!< (glocally constant) temperature
			
			// correction for quantum regions
			bool 						corrected;
			vector<double> 				elem_values;	  	//!< stores values on elements anyway.
			vector<double> 				corrected_values; 	//!< stores values of the quantum-corrected band edges
			Equation * 					corrected_eqn;
			
			// correction due to strain
			Equation * 					strain_correction;		//!< the Equation storing the strain-corrected values
			bool 						strained;				//!< is there strain present
	};

}	// end of namespace

#endif /*BANDEDGE_H_NEGF*/
