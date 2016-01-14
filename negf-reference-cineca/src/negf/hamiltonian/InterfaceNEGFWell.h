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
#ifndef INTERFACENEGFWELL_H_NEGF
#define INTERFACENEGFWELL_H_NEGF

// ==============================================================================
// This is a replacement for tdkp::InterfaceNEGFWell in case TDKP is not present
// ==============================================================================

#include "all.h"
#include "StrainPolarization.h"
using namespace negf;

namespace negf {

	//typedef flens::GeMatrix<flens::FullStorage<complex<double>, flens::ColMajor> >   Hamiltonian;
	//typedef flens::GeMatrix<flens::FullStorage<double, flens::ColMajor> >            Overlap;
	
	/** abstract negf hamiltonian builder interface */
	class InterfaceNEGFWell {
	public:
	
		//static  InterfaceNEGFWell* factory(const string& id_string, const InterfaceConfiguration& config, InformationDesk& information_desk, const char* gridfile);
		virtual ~InterfaceNEGFWell() {}	
		virtual void set_potential(const vector<double>& node_potential)     = 0;
		virtual void assemble_hamiltonian(const double& k_transversal)       = 0;	
		virtual void get_hamiltonian(/*Hamiltonian*/GEMatrix& target_matrix) = 0;
		virtual void get_overlap(/*Overlap*/DGEMatrix& overlap)              = 0;

        //! assign strain corrections. the Hamiltonian needs to re-implement this function.
        virtual void set_strain(StrainPolarization * strainpol) { NEGF_EXCEPTION("strain corrections not implemented."); }
		    		
	private:

};

}

#endif /*INTERFACENEGFWELL_H_NEGF*/
