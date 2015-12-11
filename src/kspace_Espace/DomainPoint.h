/*
Copyright (c) 2009 Sebastian Steiger, Ratko Veprek, Integrated Systems Laboratory, ETH Zurich.
Comments, suggestions, criticism or bug reports are welcome: steiger@purdue.edu. 

This file is part of ANGEL, a simulator for LEDs based on the NEGF formalism.
The software is distributed under the Lesser GNU General Public License (LGPL).
ANGEL is free software: you can redistribute it and/or modify it under the terms 
of the Lesser GNU General Public License v3 or later. ANGEL is distributed
without any warranty; without even the implied warranty of merchantability or 
fitness for a particular purpose. See also <http://www.gnu.org/licenses/>.
*/
#ifndef _DOMAINPOINT_H
#define _DOMAINPOINT_H

#include "all.h"

using namespace std;

namespace negf {

	class DomainVertex {
	
	public:
	
		// setup functions
		DomainVertex(int index_, double x, double y, double z);	// 3D
		DomainVertex(int index_, double x, double y);			// 2D
		DomainVertex(int index_, double x);						// 1D
		~DomainVertex() {}
		void set_index(uint index_);
		void set_coordinate(int xyz, double value); // do not even think of using the following method except for debug purposes!!!!!!!!!
		
		// access functions
		int			get_index()  				const;
		usint		get_dimension() 			const 	{ return dimension; }
		double		get_coordinate(usint ii) 	const; // ii=0,...
				
		// additional functions
		double		get_distance_to(DomainVertex * other_vertex) const;

	private:
		usint dimension;		// dimensionality
		vector<double> coord;	// stores coordinates
		int index;				// index
	};
	
	
	/** basic building unit: a vertex with some weight */
	class DomainPoint {
		
	public:
	
		DomainPoint(const double& x);	
		DomainPoint(const double& x, const double& y);
		DomainPoint(const double& x, const double& y, const double& z);
		// standard copy constructor is o.k.!
		
		const double 	get_coord(uint idx) const { return vertex.get_coordinate(idx); }
		const double 	get_weight() 		const { return weight; }
		double 			get_coord_abs() 	const;
		int  			get_index() 		const { return vertex.get_index(); }
		uint  			get_dimension() 	const { return vertex.get_dimension(); }
		
		void 			set_weight(const double& weight_) { weight = weight_; }
		void 			set_index(uint index_)    { vertex.set_index(index_); }
		
	protected:
	
		DomainVertex  vertex;	// has index and coordinates
		double 	 	  weight;
	};
	
	
	ostream& operator<<(ostream& out, const DomainPoint& point);

} // namespace negf

#endif
