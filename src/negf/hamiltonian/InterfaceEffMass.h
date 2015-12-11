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
#ifndef INTERFACEEFFMASS_H_
#define INTERFACEEFFMASS_H_

#include "all.h"

#include "Geometry.h"
#include "PropertyContainer.h"
#include "MaterialDatabase.h"
#include "StrainPolarization.h"

#ifndef NOTDKP
	#include "tdkp/interface/Interface.h" // must have TDKP in #include-path
	#include "tdkp/interface/InterfaceNEGFWell.h" 
#else
	#include "InterfaceNEGFWell.h" 
#endif
#include "TdkpInfoDesk.h"
/* InterfaceNEGFWell.h defines the following typedefs:
typedef GeMatrix<FullStorage<complex<double>, ColMajor> >   Hamiltonian;
typedef GeMatrix<FullStorage<double, ColMajor> >            Overlap;
*/

namespace negf {
	
	/** Effective Mass Hamiltonian and FEM stiffness matrix... */
#ifndef NOTDKP
	class InterfaceEffMass : public tdkp::InterfaceNEGFWell
#else
	class InterfaceEffMass : public InterfaceNEGFWell
#endif
	{
	public:
	
		InterfaceEffMass(const Geometry * xspace_,
						 const MaterialDatabase * db,
						 const double temperature) throw (Exception *);
		~InterfaceEffMass() {}
		
		/** inherited functions */
		void set_potential(const vector<double>& node_potential);	// energetic potential, not electrostatic potential!
		void assemble_hamiltonian(const double& k_transversal_nm);
        void set_strain(StrainPolarization * strainpol_);
		
		void get_hamiltonian(GEMatrix & target_matrix) { target_matrix = potential_and_k_augmented_hamiltonian; }
		void get_overlap(DGEMatrix & target_matrix) { target_matrix = overlap; }
		
	protected:
				
		const Geometry * xspace;
		
		DGEMatrix overlap;
		GEMatrix core_hamiltonian;
		GEMatrix potential_augmented_hamiltonian;
		GEMatrix potential_and_k_augmented_hamiltonian;

        StrainPolarization * strainpol;
		
	};
	
} // end namespace negf

#endif /*INTERFACEEFFMASS_H_*/
