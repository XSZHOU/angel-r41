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
#ifndef HAMILTONIAN_H_
#define HAMILTONIAN_H_

#include "all.h"

#include "Geometry.h"
#include "Kspace.h"
#include "Options.h"
#include "PropertyContainer.h"
#include "TernaryPropertyContainer.h"
#include "MaterialDatabase.h"
#include "DomainMaster.h"
#include "DomainPoint.h"
#include "StrainPolarization.h"

#ifndef NOTDKP
#include "tdkp/interface/Interface.h" // must have TDKP in #include-path
#include "tdkp/interface/InterfaceNEGFWell.h" 
#endif
#include "InterfaceEffMass.h"
#include "InterfaceEffMassOrtho.h"
#include "InterfaceEffMassOrtho2Band.h"
#include "TdkpInfoDesk.h"

namespace negf {
	
	/** gives KP (or whatever)  Hamiltonian for a given TRANSVERSAL (not in transport direction) k-point
	 * (without electrostatic potential, i.e. only the kinetic part) */
	class Hamiltonian
	{
	public:
	
		Hamiltonian(const Geometry * xspace_, const Kspace * kspace_, const Options * options_,
						 const MaterialDatabase * db_, const char * xgridfilename_) throw (Exception *);
		~Hamiltonian() {}

        // -----------------------------
        // setup
        // -----------------------------

		//! set the new electrostatic potential
		void set_electrostatic_potential(const vector<double> & elstat_potential) throw (Exception *);

		void set_strain(StrainPolarization * strainpol_);

		// -----------------------------
		// access
		// -----------------------------
	
		// obtain the Hamiltonian at ALL x-points (not only the internal) for a given transversal k-point 
		void get(const DomainPoint & kk, Matc & result) const throw (Exception *);
		void get(double kk, Matc & result) const { NEGF_EXCEPTION("do not use!"); }
		void get(int kk, Matc & result) const { NEGF_EXCEPTION("do not use!"); }
		
		// obtain the Hamiltonian at the internal x-points only for a given transversal k-point 
		void get_internal(const DomainPoint & kk, Matc & result) const throw (Exception *);
		void get_internal(double kk, Matc & result) const { NEGF_EXCEPTION("do not use!"); }
		void get_internal(int kk, Matc & result) const { NEGF_EXCEPTION("do not use!"); }
		
		// obtain the overlap matrix of the basis functions - will be called by Overlap class - reinitializes result matrix
		void get_overlap(Matd & result) const throw (Exception *);
		
		// other trivial access functions
		const vector<double> & 				  get_electrostatic_potential() const { return this->electrostatic_potential; }
		const Geometry 						* get_xspace() 		const { return this->xspace; }
		const Kspace 						* get_kspace() 		const { return this->kspace; }
		const Options 						* get_options() 	const { return this->options; }
		const MaterialDatabase 				* get_material_db() const { return this->db; }
#ifndef NOTDKP
		const tdkp::InterfaceConfiguration  & get_tdkp_config() const { return this->config; }
		TdkpInfoDesk 						* get_tdkp_infodesk() const { return this->infodesk; }
		tdkp::InterfaceNEGFWell 			* get_tdkp_interface() const { return this->interface; }
#endif		
	
	protected:
		
		// helpers
#ifndef NOTDKP
		void configure_interface(tdkp::InterfaceConfiguration & conf, bool is_contact);
#endif		

		vector<double> 					electrostatic_potential;
		
		const Geometry * 				xspace;
		const Kspace   * 				kspace;
		const Options  * 				options;
		const MaterialDatabase * 		db;
		string 							xgridfilename;
		uint 							Nx;
		uint 							NxNn;
		
		TdkpInfoDesk * 					infodesk;
		
#ifndef NOTDKP
		tdkp::InterfaceConfiguration 	config;
		tdkp::InterfaceNEGFWell *		interface;
		tdkp::InterfaceConfiguration	contact_0_config;
		tdkp::InterfaceBulkRadial * 	contact_0_bandstructure;
#else
		InterfaceNEGFWell *				interface;
#endif
		vector<double> 					contact_0_dos_n;
		vector<double> 					contact_0_dos_p;
		
		Matc 							work_matrix;	// used in get_internal		
	};
	
} // end namespace negf

#endif /*HAMILTONIAN_H_*/
