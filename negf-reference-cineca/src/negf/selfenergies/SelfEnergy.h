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
#ifndef SELFENERGY_H_
#define SELFENERGY_H_

#include "all.h"

#include "Geometry.h"
#include "Kspace.h"
#include "Energies.h"
#include "Options.h"
#include "NEGFObject.h"

namespace negf {

#ifdef USE_BANDED
	typedef BandedNEGFObject SL_object_type;
#else 
	typedef FullNEGFObject SL_object_type;
#endif

	/** Base class for a Self-Energy (i.e., the retarded, lesser and greater matrices) 
	 * template parameter: matrix type (full of banded <-- variable number of offdiagonals) */
	class SelfEnergy {
	public:
		/** TO BE CALLED BY CONSTRUCTORS OF DERIVED CLASSES BEFORE DOING ANYTHING ELSE */
		SelfEnergy(const Geometry * xspace_, const Kspace * kspace_, const Energies * energies_, const Options * options_, uint num_offdiags);
		
		virtual ~SelfEnergy();
		
		/** TO BE IMPLEMENTED BY DERIVED CLASSES */
		virtual void calculate() = 0;
		
		/** functions to access the whole self-energy object */
		SL_object_type * get_retarded() const { return this->SigmaR; }
		SL_object_type * get_lesser() 	const { return this->SigmaL; }
		SL_object_type * get_greater() 	const { return this->SigmaG; }
		
		/** functions to access locally stored matrices */
		SEMat & 			get_retarded(uint kidx, uint global_Eidx) const { return this->SigmaR->get(kidx, global_Eidx); }
		SEMat & 			get_lesser  (uint kidx, uint global_Eidx) const { return this->SigmaL->get(kidx, global_Eidx); }
		SEMat & 			get_greater (uint kidx, uint global_Eidx) const { return this->SigmaG->get(kidx, global_Eidx); }
	
		/** other trivial access functions */
		const Geometry 	* 	get_xspace() 	const { return this->xspace; }
		const Kspace 	* 	get_kspace() 	const { return this->kspace; }
		const Options 	* 	get_options() 	const { return this->options; }
		const Energies 	* 	get_energies()  const { return this->energies; }
				
	protected:
				
		const Options  * options;	// stores what scattering is used, etc.
		const Geometry * xspace;	// real-space grid
		const Kspace   * kspace;	// k-space grid (including integration weights for each point)
		const Energies * energies;	// 1D energy space grid (including integration weights for each point)
		
		const uint Nx;		// number of x-points
		const uint NxNn;	// product Nx*Nn
		const uint Nk;		// number of k-points
		const uint NE;		// number of energy points in total
		uint       myNE;	// number of energy points belonging to the current MPI thread, may change during the simulation
		const uint Nvert;	// total number of vertices, including contact vertices
	
		SL_object_type * SigmaR;	// retarded component
		SL_object_type * SigmaL;	// lesser component
		SL_object_type * SigmaG;	// greater component
	};


} // end of namespace

#endif /*SELFENERGY_H_*/
