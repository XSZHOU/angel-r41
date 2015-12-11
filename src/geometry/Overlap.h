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
#ifndef OVERLAP_H_NEGF
#define OVERLAP_H_NEGF

#include "all.h"

#include "Geometry.h"
#include "Options.h"
#include "Hamiltonian.h"

namespace negf {
	
	/** gives FEM overlap matrices, inflated by the ndegrees of freedom, either for all vertices or the internal vertices only */
	class Overlap
	{
	public:
	
		Overlap(Hamiltonian * ham_, const Geometry * xspace_) throw (Exception *);
		~Overlap() {}
	
		// obtain the (Nvert*Nn)^2 overlap matrix at ALL x-points (with global vertex indices)
		const OVMat & get_overlap() const;
		
		// obtain the (Nx*Nn)^2 overlap matrix between INTERNAL x-points (with internal vertex indices)
		const OVMat & get_internal_overlap() const;
		
		// obtain the Nvert*Nvert overlap matrix at ALL x-points but NOT INFLATED
		const Matd &  get_vertex_overlap() const { return this->overlap; }
		
		// obtain the Nx*Nx overlap matrix at ALL x-points but NOT INFLATED
		const Matd &  get_vertex_internal_overlap() const { return this->small_overlap; }
					
	protected:
		
		Hamiltonian * 	 ham;
		const Geometry * xspace;
		const uint 		 Nx;			// internal vertices where NEGF problem is solved
		const uint 		 NxNn;
		const uint 		 Nvert; 		// all vertices of the geometry, including those in the contact
		
		Matc 	global_overlap;			// size (Nvert*Nn)^2, Nvert=grid->get_num_vertices()
		Matc 	internal_overlap;		// size (Nx*Nn)^2, Nx=grid->get_num_internal_vertices()
		BMatc   global_overlap_banded;
		BMatc   internal_overlap_banded;
		Matd 	overlap;
		Matd 	small_overlap;
		
	};
	
} // end namespace negf

#endif /*OVERLAP_H_NEGF*/
