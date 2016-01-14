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
#ifndef DOMAIN_H_
#define DOMAIN_H_

#include "all.h"

#include "DomainMaster.h"

namespace negf {


	/** TO CREATE THESE DOMAINS, CALL THE FOLLOWING FUNCTIONS: */
	
	/** create 1D regular grid */
	void create_1D_domain_trapezoidal(DomainMaster& domain, const double& x0, const double& x1, const uint num);
	/** create radial along [1 0] */
	void create_2D_domain_radial(DomainMaster& domain, const double& r0, const double& r1, const uint num);
	/** create radial along [1 0 0] */
	void create_3D_domain_radial(DomainMaster& domain, const double direction[3], const double& r0, const double& r1, const uint num);
	
	
	
	
	/** The following is able to interpolate between 2 grids (domains), e.g. if the target is a finer grid as the source
	 *  which values are taken from the old gridpoints for a new gridpoint is stored in contribution */
	
	class DomainPointContribution {
	public:	
		DomainPointContribution() : point_idx(0), contribution(0.0) {}
		DomainPointContribution(unsigned int point_idx_, const double& contribution_)
		  : point_idx(point_idx_), contribution(contribution_) {}
		
		/** point idx in from Domain Master */
		unsigned int point_idx;
		/* usually inside ]0,1] (except for extrapolations) */
		double contribution;		
	};
	
	class DomainMap {
	public:
		typedef list<DomainPointContribution>::const_iterator PointMapIterator;
		
		DomainMap(const DomainMaster& source_, const DomainMaster& target_);
		
		PointMapIterator begin(unsigned int target_point_idx) {
			NEGF_ASSERT(point_map.size() > target_point_idx, "point_map.size() > target_point_idx");
			return point_map[target_point_idx].begin();	
		}
	
		PointMapIterator end(unsigned int target_point_idx) {
			NEGF_ASSERT(point_map.size() > target_point_idx, "point_map.size() > target_point_idx");
			return point_map[target_point_idx].end();
		}
	
		const DomainMaster& get_source_domain() { return source; }
		const DomainMaster& get_target_domain() { return target; }	
				
	private:
		const DomainMaster& source;
		const DomainMaster& target;
		vector<list<DomainPointContribution> > point_map;
	};
	

} // end of namespace

#endif /*DOMAIN_H_*/
