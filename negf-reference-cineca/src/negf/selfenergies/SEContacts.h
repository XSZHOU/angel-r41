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
#ifndef SECONTACTS_H_
#define SECONTACTS_H_

#include "all.h"

#include "PropertyContainer.h"
#include "Geometry.h"
#include "Kspace.h"
#include "Energies.h"
#include "Options.h"
#include "NEGFObject.h"
#include "Hamiltonian.h"
#include "Overlap.h"
#include "SelfEnergy.h"

namespace negf {
	
	/** Contact self-energy */
	class SEContacts : public SelfEnergy {
		
	friend class ContactFermilevel;
	
	public:
		
		SEContacts(const Hamiltonian * ham_,
					const Overlap * ov_,
					const Geometry * xspace_, 
					const Kspace * kspace_, 
					const Energies * energies_, 
					const Options * options_);
		
		~SEContacts() {}
		
		// inherited from SelfEnergy class
		void calculate();
		
		bool was_calculated() const { return this->calculated; }
		const vector<Vertex *> & get_interface_vertices(uint cc) const { return this->interface_vertices[cc]; }	// needed in OuterLoop.cpp
		const vector<Vertex *> & get_second_row_vertices(uint cc) const { return this->second_row_vertices[cc]; }	// needed in OuterLoop.cpp
		double 					 get_dx(uint cc) const { return this->dx[cc]; }
		
		void 					 assign_broadening(vector< vector< vector< vector<cplx> > > > & frey_broadening_);

		void 					 assign_bandedges(const vector<double> & cb_, const vector<double> & vb_);
		
		double 					 get_fermi(uint contact_idx, uint global_Eidx) const;	// also used in InnerLoop (frey broadening)
		
	protected:
		
		/** helper functions */
		void calculate_retarded();
		void calculate_easy_retarded();
		void calculate_lesser_greater();
		void get_internal_idx_of_vertices_near_contact(uint contact_idx, vector<uint> & verts) const;
		double get_F0(uint contact_idx, uint global_Eidx) const;
				
		const Hamiltonian * ham;
		const Overlap * ov;
		double temperature;
		
		vector< vector<Vertex *> > interface_vertices;	// gives for each contact the vertices at the interface to the device
		vector< vector<Vertex *> > second_row_vertices;	// gives for each interface_vertex the connected vertex lying perpendicular inside the CONTACT
		vector< vector<Vertex *> > device_vertices;		// gives for each interface_vertex the connected vertex lying perpendicular inside the DEVICE
		vector<double> dx;								// gives for each contact the distance between an interface vertex and its second-row vertex
		
		bool calculated;
		
		const bool security_checking;
		
		// stores for each own energy ee2, k-point kk, band index nn and contact cc (in this order) the complex optical potential
		vector< vector< vector< vector<cplx> > > > frey_broadening;

		// band edges to cut off some injecting states, if .cmd-file option "InjectingStatesCutoff" is set
		vector<double> cb;
		vector<double> vb;
	};
	
	
} // end of namespace

#endif /*SECONTACTS_H_*/
