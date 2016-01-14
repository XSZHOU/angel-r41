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
#ifndef KSPACE_H_
#define KSPACE_H_

#include "all.h"

#include "Options.h"
#include "DomainPoint.h"
#include "DomainMaster.h"
#include "Domain.h"
#include "Geometry.h"
#include "MaterialDatabase.h"

namespace negf {
	
	class Kspace {
		
	public:
		
		// sets up k-space grid, dimensionality depends on x-grid (D(k) = 3-D(x))
		// remember: for a given energy, all k-points are assigned to the same MPI process
		Kspace(const Geometry * xspace, const Options * options_, const MaterialDatabase * db) throw (Exception *);
		~Kspace() {}
	
		// ------------------
		// access functions
		// ------------------
		uint 				get_number_of_points() 	const { return Nk; }
		
		/** please note: for a 2D k-space with axial approximation (only |k| is considered),
		    point.get_weight() will include the factor 2*pi! */
		const DomainPoint & get_point(uint idx) 	const { return *(k_grid[idx]); } 
		
		uint 				get_dimension() 		const { return dimension; } 
		double 				get_kmin() 				const { return kmin; } 
		double 				get_kmax() 				const { return kmax; } 
		
		enum spacing {
		    equal,
		    square_root
		};
		
		enum integration_rule {
		    trapez,
		    romberg_simpson,
		    three_point
		};

	protected:
	
		void test_accuracy(const double & m);
		
		// old
		// DomainMaster k_grid;
		
		// new
		vector<DomainPoint *> k_grid;	// stores a DomainVertex and the weights 
		uint 				dimension;
		uint 				Nk;
		double 				kmin;
		double 				kmax;
		const Options *     options;
		
		spacing          my_spacing;
		integration_rule my_rule;

	};
	
} // end of namespace


#endif /*KSPACE_H_*/
