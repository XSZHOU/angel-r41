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
#include "Domain.h"
using namespace negf;

// =========================================================
// Implementations of domain constructors
// =========================================================

/** build trapezoidal 1D domain */
void negf::create_1D_domain_trapezoidal(DomainMaster& domain, const double& x0, const double& x1, const unsigned int num) 
{STACK_TRACE(
	NEGF_ASSERT(num > 1, "num > 1");
	domain.clean(); 
	
	const double dx  = (x1 - x0) / static_cast<double>(num - 1);
	const double wdx = negf_math::abs(dx);
	
	for(unsigned int ii = 0; ii < num; ii++) {
		double weight = (ii > 0 && ii < num - 1) ? wdx : wdx / 2.0;
		domain.add_node(new DomainNodeSingularPoint(false, weight, x0 + dx * ii)); 
	}
	domain.update(); 		
);}

/** build radial 2D domain in [10] direction */
void negf::create_2D_domain_radial(DomainMaster& domain, const double& r0, const double& r1, const unsigned int num) 
{STACK_TRACE(
	NEGF_ASSERT(r1 > r0, "r1 > r0");

	domain.clean();

	// ------------------------------------------
	// o.k. we build it that way: 
	// ------------------------------------------
	// num_k_values - 2 is the number of normal radial intervals
	// interval_length = (r1 - r0) / (num k values - 1)	
	// intervals:
	// first  [r0, r0+dr/2]
	// then   [r0 + dr/2 + ii * dr, +dr]
	// last   [r1 - dr/2, r1]
	// ------------------------------------------- 
			
	const double dr  = (r1 - r0) / static_cast<double>(num - 1);
	
	// 2pi ri dr
	double weight = constants::pi * (((r0 + dr / 2.0) * (r0 + dr / 2.0)) - (r0 * r0));
	
	// add first node at r0 as singular point node
	domain.add_node(new DomainNodeSingularPoint(true, weight, r0, 0.0));
	// add subsequent nodes
	for(unsigned int ii = 0; ii < num - 2; ii++) {
		domain.add_node(new DomainNodeRadialPlane(1.0, 0.0, r0 + dr/2.0 + dr*ii, r0 + dr/2.0 + dr*ii + dr)); 		
	}
	// add last node
	weight = constants::pi * ((r1 * r1) - ((r1 - dr / 2.0) * (r1 - dr / 2.0)));
	domain.add_node(new DomainNodeSingularPoint(true, weight, r1, 0));
		  
	domain.update();
);}

void negf::create_3D_domain_radial(DomainMaster& domain, const double direction[3], const double& r0, const double& r1, const unsigned int num) 
{STACK_TRACE(
	logmsg->emit(LOG_WARN,"WARNING: 3D RADIAL DOMAIN DOES NOT WORK: WEIGHTS ARE NOT SET TO POINTS!");
	NEGF_ASSERT(r1 >= r0, "r1 >= r0"); 

	// normalize direction vector
	double mydir[3] = { direction[0], direction[1], direction[2] };
	negf_math::vector_normalize_3d(mydir);
	
	domain.clean();

	if(num == 1) { 
		domain.add_node(new DomainNodeSingularPoint(true, 1.0e-100, r0 * mydir[0], r0 * mydir[1], r0 * mydir[2]));
	} else {	
		const double dr  = (r1 - r0) / static_cast<double>(num - 1);
		
		double rr = r0 + dr / 2.0;
		double rl = r0;
		double weight = 4.0 / 3.0 * constants::pi * ((rr * rr) - (rl * rl));
		
		// add first node at r0 as singular point node
		domain.add_node(new DomainNodeSingularPoint(true, weight,  r0 * mydir[0], r0 * mydir[1], r0 * mydir[2]));
		// add subsequent nodes		
		for(unsigned int ii = 0; ii < num - 2; ii++) {
			rl = r0 + dr / 2 + dr * ii;
			rr = r0 + dr / 2 + dr * (ii + 1);
			weight = 4.0 / 3.0 * constants::pi * ((rr * rr) - (rl * rl));	
			double rii = (rr + rl) / 2.0;
			domain.add_node(new DomainNodeSingularPoint(true, weight, rii * mydir[0], rii * mydir[1], rii * mydir[2]));					
		}
		// add last node
		rl = r1 - dr / 2;
		rr = r1;
		weight = 4.0 / 3.0 * constants::pi * ((rr * rr) - (rl * rl));
		domain.add_node(new DomainNodeSingularPoint(true, weight, rr * mydir[0], rr * mydir[1], rr * mydir[2]));
	}
	domain.update();
);}



// =====================================
// domain map implementatiion
// =====================================

DomainMap::DomainMap(const DomainMaster& source_, const DomainMaster& target_) 
: source(source_), 
  target(target_),
  point_map(target.get_number_of_points())
{STACK_TRACE(
	// -----------------------------------------------
	// check that both domains have the same dimension
	// -----------------------------------------------
	NEGF_ASSERT(source.get_dimension() == target.get_dimension(), "source.get_dimension() == target.get_dimension()");
	NEGF_ASSERT(source.get_number_of_points() > 1, "source.get_number_of_points() > 1");  	
	NEGF_ASSERT(target.get_number_of_points() > 1, "target.get_number_of_points() > 1");	

	// -----------------------------------------------
	// choose strategy
	// -----------------------------------------------
	if(source.radial() && target.radial()) {
		vector<double> source_place(source.get_number_of_points());
		vector<double> target_place(target.get_number_of_points());
		// --------------------------------------------
		// check if domain points are ordered properly
		// --------------------------------------------
		unsigned int idx = 0;
		DomainMaster::point_const_iterator it = source.begin();
		source_place[idx++] = (*it)->get_coord_abs();		
		for(it++; it != source.end(); it++, idx++) {
			source_place[idx] = (*it)->get_coord_abs();			
			NEGF_ASSERT(source_place[idx - 1] < source_place[idx], "source_place[idx - 1] < source_place[idx]"); 			    
		}
		idx = 0;
		it = target.begin();
		target_place[idx++] = (*it)->get_coord_abs();		
		for(it++; it != target.end(); it++, idx++) {
			target_place[idx] = (*it)->get_coord_abs();			
			NEGF_ASSERT(target_place[idx - 1] < target_place[idx], "target_place[idx - 1] < target_place[idx]"); 			    
		}
		// --------------------------------------------
		// loop over target values and build map
		// --------------------------------------------
		// extrapolation 
		// f(x1 + dx) = f(N1) + (f(N2) - f(N1)) / (N2 - N1) dx
		//            = f(N1) * (1 - dx/(N2 - N1)) + f(N2) * (dx / (N2 - N1))
		idx = 0;
		double dx; 
		double dN2N1 = source_place[1] - source_place[0];
		while(target_place[idx] < source_place[0]) {
			dx = target_place[idx] - source_place[0];
			point_map[idx].push_back(
				DomainPointContribution(0, 1.0 - dx / dN2N1)
			);
			point_map[idx].push_back(
				DomainPointContribution(1, dx / dN2N1)
			);
			idx++;			
		}
		// ---------------------------------------------
		// interpolate in between
		// ---------------------------------------------
		unsigned int source_frame_idx = 0;
		double fraction;
		while(idx < target_place.size() && target_place[idx] <= source_place.back()) {
			
  
			
			// shift source frame if necessary
			while(source_frame_idx < source_place.size() - 1 &&
			      target_place[idx] >= source_place[source_frame_idx + 1]) {
				source_frame_idx++;
			}
			/*
			cout << "target: " << idx << "(" << target_place.size() << ") = " << target_place[idx] << " source: " 
			     << source_frame_idx << "( " << source_place.size() << ") = " << source_place[source_frame_idx] << " next = " << source_place[source_frame_idx + 1] << "\n";*/  
			// catch exact match
			if(source_place[source_frame_idx] == target_place[idx]) {
				point_map[idx].push_back(
					DomainPointContribution(source_frame_idx, 1.0)
				);
			} else {
				NEGF_ASSERT(source_frame_idx < source_place.size() - 1, "source_frame_idx < source_place.size() - 1"); 
				// in between
				fraction = (target_place[idx] - source_place[source_frame_idx]) 
				         / (source_place[source_frame_idx + 1] - source_place[source_frame_idx]); 
				point_map[idx].push_back(
					DomainPointContribution(source_frame_idx, 1.0 - fraction)
				);				
				point_map[idx].push_back(
					DomainPointContribution(source_frame_idx + 1, fraction)
				);	
			}
			idx++;			       			
		}				
		// same with end points
		// f(xn + dx) = f(xn) + (f(xn) - f(xn-1)) / dxnxn-1 * dx
		//            = f(xn) * (1 + dx / dxnxn-1) - f(xn-1) * dx/dxnxn-1
		idx = target_place.size() - 1;
		dN2N1 = source_place[source_place.size() - 1] - source_place[source_place.size() - 2];
		while(target_place[idx] > source_place.back()) {
			dx = target_place[idx] - source_place.back();
			point_map[idx].push_back(
				DomainPointContribution(source_place.size() - 1, 1.0 + dx / dN2N1)
			);
			point_map[idx].push_back(
				DomainPointContribution(source_place.size() - 2, - dx / dN2N1)
			);
			idx--;							
		}
	} else {
		NEGF_EXCEPTION("sorry, mapping is currently only implemented for radial domains");	
	}
);}

