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
#ifndef POSTPROCESSING_H_
#define POSTPROCESSING_H_

#include "all.h"

#include "Geometry.h"
#include "Kspace.h"
#include "Energies.h"
#include "Options.h"

#include "SelfEnergies.h"
#include "GreenFunctions.h"
#include "Current.h"
#include "Transmission.h"
#include "QuasiFermilevel.h"

#include "TdkpInfoDesk.h" // for get_spin_degeneracy

#include "SEPhotonSpontaneous.h"
#include "Luminescence.h"

namespace negf {
	
	/** wrapper for quantities that can be obtained from the Green functions and self energies */
	class PostProcessing
	{
	public:
	
		PostProcessing(const Hamiltonian * ham_, 
						const Overlap * ov_,
						const GreenFunctions * gf_, 
						const SelfEnergies * se_, 
						const Options * opts_, 
						const Geometry * xspace_,
						const Kspace * kspace_, 
						const Energies * energies_) throw (Exception *);
		~PostProcessing();
	
		/** compute LDOS(x,E) = -1/pi*sum_{n,k} Im GR(x,n,x',n';k,E) */
		void 				compute_local_dos() throw (Exception *);
		
		/** compute n(x,E) = sum_{nc,k} -i*GL(x,nc,x,nc;k,E) or 
		 *          p(x,E) = sum_{nv,k} +i*GG(x,nv,x,nv;k,E)       */
		void 				compute_spectral_edensity() throw (Exception *);
		void 				compute_spectral_hdensity() throw (Exception *);
		
		/** compute n(x) = sum_E dE/2pi n(x,E) or 
		 *          p(x) = sum_E dE/2pi p(x,E)
		 *  the result will be stored in the master process only! */
		void 				compute_edensity() throw (Exception *);
		void 				compute_hdensity() throw (Exception *);
		
		void 				compute_spectral_current() const throw (Exception *);
		void 				compute_current() const throw (Exception *);
		
		void 				compute_transmission() const throw (Exception *);
		
        void                compute_quasi_fermilevel() const throw (Exception *);

		void 				compute_scattering_current(SelfEnergyType type) const throw (Exception *); 
		
		void 				check_conservation(SelfEnergyType type) const throw (Exception *); 
		
		void 				output_selfenergies(double Etarget, const char * filename);
		
		/** get the computed results */
		const Matd & 	get_local_dos() 		const { return LDOS; }
		const Matd & 	get_local_dos_k0() 	 	const { return LDOS_k0; }
		const Matd & 	get_local_dos_VB() 	 	const { return LDOS_VB; }
		const Matd & 	get_spectral_edensity() const { return spectral_edensity; }
		const Matd & 	get_spectral_hdensity() const { return spectral_hdensity; }
		const Matd & 	get_spectral_ecurrent() const { return this->ecurrent->get_spectral_current(); }
		const Matd & 	get_spectral_hcurrent() const { return this->hcurrent->get_spectral_current(); }
		const Matd & 	get_entire_spectral_edensity()  const throw (Exception *);
		const Matd & 	get_entire_spectral_hdensity()  const throw (Exception *);
		const Matd & 	get_entire_spectral_ecurrent()  const throw (Exception *) { return this->ecurrent->get_entire_spectral_current(); }
		const Matd & 	get_entire_spectral_hcurrent()  const throw (Exception *) { return this->hcurrent->get_entire_spectral_current(); }
		const Matd & 	get_entire_spectral_ecurrent2() const throw (Exception *) { return this->ecurrent->get_entire_spectral_current2(); }
		const Matd & 	get_entire_spectral_hcurrent2() const throw (Exception *) { return this->hcurrent->get_entire_spectral_current2(); }
		const Matd & 	get_entire_scattering_ecurrent(SelfEnergyType type) const throw (Exception *) { return this->ecurrent->get_entire_scattering_current(type); }
		const Matd & 	get_entire_scattering_hcurrent(SelfEnergyType type) const throw (Exception *) { return this->hcurrent->get_entire_scattering_current(type); }
		const Matd & 	get_entire_local_dos()  const throw (Exception *);
		const Matd & 	get_entire_local_dos_k0()  const throw (Exception *);
		const Matd & 	get_entire_local_dos_VB()  const throw (Exception *);
		const vector<double> & get_edensity()    const throw (Exception *);
		const vector<double> & get_hdensity()    const throw (Exception *);
		const vector<double> & get_ecurrent()    const { return this->ecurrent->get_current(); }
		const vector<double> & get_hcurrent()    const { return this->hcurrent->get_current(); }
		ContactCurrent *  	get_contact_current() const throw (Exception *);

		Transmission *      get_transmission() const throw (Exception *);

        QuasiFermilevel *   get_quasi_fermilevel() const throw (Exception *);
		
		/** other trivial access functions */
		const Geometry 	* 	get_xspace() 	const { return this->xspace; }
		const Kspace 	* 	get_kspace() 	const { return this->kspace; }
		const Options 	* 	get_options() 	const { return this->opts; }
		const Energies 	* 	get_energies()  const { return this->energies; }
		
		/** luminescence-related */
		void 				set_up_luminescence(SEPhotonSpontaneous * spont_) throw (Exception *);
		void 				compute_luminescence() const throw (Exception *);
		Luminescence * 		get_luminescence() const throw (Exception *);	// Lake-style
		Luminescence * 		get_luminescence2() const throw (Exception *);	// Galerpin-style
		
	protected:
		
		/* helpers */
		void compute_spectral_xdensity(const vector<uint> & bands, bool e_or_h);
		void compute_xdensity(vector<double> & density, const Matd & spectral_density);
	
		const Geometry	 	 * xspace;
		const Kspace 		 * kspace;
		const Energies 		 * energies;
		const Options 		 * opts;
		const Hamiltonian    * ham;
		const Overlap  		 * ov;
		const GreenFunctions * gf;
		const SelfEnergies 	 * se;
		const uint 			   Nx;
		const uint 			   NxNn;
		
		Matd  LDOS;					// size myNE*Nx
		Matd  LDOS_k0;				// size myNE*Nx
		Matd  LDOS_VB;				// size myNE*Nx
		Matd  spectral_edensity;	// size myNE*Nx
		Matd  spectral_hdensity;	// size myNE*Nx
		Matd  entire_spectral_edensity;	// size NE*Nx (root only)
		Matd  entire_spectral_hdensity;	// size NE*Nx (root only)
		Matd  entire_LDOS;			// size NE*Nx (root only)
		Matd  entire_LDOS_k0;		// size NE*Nx (root only)
		Matd  entire_LDOS_VB;		// size NE*Nx (root only)
		
		Current * ecurrent;
		Current * hcurrent;
		ContactCurrent * contact_current;
		
		Transmission * transmission;
		
		QuasiFermilevel * qfl;

		Luminescence * luminescence;	// Lake-style
		Luminescence * luminescence2;	// Galerpin-style
		
		vector<double> edensity;		// will be filled only in the master process!
		vector<double> hdensity;		// will be filled only in the master process!
		
		bool 	  i_am_master;
		
	};
	
} // end namespace negf

#endif /*POSTPROCESSING_H_*/
