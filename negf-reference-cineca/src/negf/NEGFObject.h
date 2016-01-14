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
#ifndef NEGFOBJECT_H_
#define NEGFOBJECT_H_

#include "all.h"

#include "Geometry.h"
#include "Kspace.h"
#include "Energies.h"
#include "Options.h"


namespace negf {
	
	/** base class for NEGF objects which have arguments x,x';k,E
	 * it should be noted that actually it is x,n,x',n';k,E, but the band index n is lumped into x
	 * so the object is a matrix for each k,E
	 * numbering: the matrix corresponding to k(i), E(j) is stored in real_space_matrices[i*(myNE-1)+j];
	 * 
	 * note on number of x-points: NEGF is only solved on the vertices NOT being contact vertices!
	 * */
	class BandedNEGFObject {
	public:
	
		BandedNEGFObject(const Geometry * xspace_, const Kspace * kspace_, 
				   const Energies * energies_, const Options * opts_, uint offdiags);
		~BandedNEGFObject();
		
		// get a locally stored matrix
		BMatc & get(uint kidx, uint global_Eidx);
		
		// get all locally stored matrices to some energy index (every k)
		vector<BMatc> & get(uint global_Eidx);
		
	protected:
		
		vector< vector<BMatc> > matrices;
		
		const Geometry * xspace;
		const Kspace   * kspace;
		const Energies * energies;
		const Options  * opts;
		
		const uint Nx;  // number of NEGF x-points
		const uint NxNn;// Nx*Nn
		const uint Nk;	// number of k points
		const uint NE;	// number of energy points in total
		uint myNE;		// number of energy points of the calling process, may change during the simulation
		
	};	
		
	/** same, but with full matrices */
	class FullNEGFObject {
	public:
		FullNEGFObject(const Geometry * xspace_, const Kspace * kspace_, 
				   const Energies * energies_, const Options * opts_);
		~FullNEGFObject();
		
		Matc 		 & get(uint kidx, uint global_Eidx);
		vector<Matc> & get(uint global_Eidx);
		
	protected:
		
		vector< vector<Matc> > matrices;
		
		const Geometry * xspace;
		const Kspace   * kspace;
		const Energies * energies;
		const Options  * opts;
		
		const uint Nx;  // number of NEGF x-points
		const uint NxNn;// Nx*Nn
		const uint Nk;	// number of k points
		const uint NE;	// number of energy points in total
		uint myNE;		// number of energy points of the calling process, may change during the simulation
	};	
	

} // end namespace negf


#endif /*NEGFOBJECT_H_*/
