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
#ifndef DOPING_H_NEGF
#define DOPING_H_NEGF

#include "all.h"

#include "Geometry.h"
#include "BoxMethod.h"
#include "ExplicitEquation.h"
#include "InputParser.h"


using namespace std;

namespace negf {

	/** Reads in the Doping Concentration 
	 * @param grid       the grid on which the doping is defined
	 * @param bm         the box method for that grid 
	 * @param filename   the filename from which the fields are read
	 * @param DopingUnit doping in .dat-file is assumed to be in cm-3 (0) or m-3 (1)
	 * @param DopingSign doping is multiplied by +1 or -1
	 * @param UseGrainBoundaryDoping sets if the fields BoronGrainBoundaryConcentration etc are read in
	 * @param dim        determines conversion factor (uints::density_3d, units::density_2d etc) */
	class Doping : public ExplicitEquation
	{
		public:

#ifndef NODFISE
            //! OLD - read in from DF-ISE .dat-file
			Doping( const  Geometry * grid_,
					const  BoxMethod * bm_, 
					const  string filename, 
					uint   options_DopingUnit,
					int    options_DopingSign,
					bool   options_UseGrainBoundaryDoping,
					uint   options_dim);
#endif

			// NEW - information is in .cmd-file
            Doping( const  Geometry * grid_,
                    const  BoxMethod * bm_,
                    const map< string, PropertyContainer<double> * > * cmdfile);

			~Doping() {}

			// further functions
			const Geometry * get_grid() const { return this->grid; }

		protected:
			const Geometry * grid;
			const BoxMethod * bm;

			// virtual functions from the ExplicitEquation class
			virtual void 	direct_derivatives(uint line, const Equation * eqn, uint & nonzeros, 
										uint indices [], double coeff []) const
				 { NEGF_EXCEPTION("The class is not designed for this."); }
			virtual double compute_value(uint ii) const
				 { NEGF_EXCEPTION("The class is not designed for this."); }

	};

}	// end namespace gain

#endif /*DOPING_H_NEGF*/
