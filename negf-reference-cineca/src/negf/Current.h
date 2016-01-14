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
#ifndef CURRENT_H_
#define CURRENT_H_

#include "all.h"

#include "Hamiltonian.h"
#include "Overlap.h"
#include "NEGFObject.h"
#include "GreenFunctions.h"
#include "SelfEnergies.h"
#include "TdkpInfoDesk.h" // for get_spin_degeneracy


namespace negf {
	
	
	class ContactCurrent {
		public:
			ContactCurrent(const Geometry * xspace_, const Energies * energies_);
			~ContactCurrent() {}
		
			void set_entire_spectral_ecurrent(const Matd * entire_spectral_ecurrent_);
			void set_entire_spectral_hcurrent(const Matd * entire_spectral_hcurrent_);
			void set_eJleft_eJright(const double * eJleft_, const double * eJright_) { this->eJleft = eJleft_; this->eJright = eJright_; }
			void set_hJleft_hJright(const double * hJleft_, const double * hJright_) { this->hJleft = hJleft_; this->hJright = hJright_; }
			
			void snapshot(double voltage); 		// associates the current variable value with "voltage" and stores it
			void save_to_file(const char * filename); // writes all stored voltages and current into the given filename

	        // from Equation
			uint get_num_variables() const { return this->number_of_variables; }
	        uint get_timestamp() const { return timestamp; }
	        void set_timestamp(uint new_timestamp) { this->timestamp = new_timestamp; }
	        void compute_values(uint new_timestamp);

		protected:
			virtual double  compute_value(uint line) const;
			
			vector<double> 			 voltages;
			vector<double> 			 Jlefts;
			vector<double> 			 Jrights;
			const double * 			 eJleft;
			const double *		 	 eJright;
			const double * 			 hJleft;
			const double *		 	 hJright;
			vector< vector<double> > contact_currents;
			
			const Geometry * xspace;
			const Energies * energies;
			const Matd * entire_spectral_ecurrent;
			const Matd * entire_spectral_hcurrent;

	        // from Equation
	        uint timestamp;
	        uint number_of_variables;
	        vector<double> current_variable_values;
	};
	
	/** Spectral current and current at contacts */
	class Current {
	public:
	
		Current(const Hamiltonian * ham_, const Overlap * ov_, const GreenFunctions * gf_,
				const SelfEnergies * se_, quantities::PhysicalQuantity e_or_h_) throw (Exception *);
		~Current() {}
			
		/** compute 
		             J(x,E) = 1/hbar * Tr(H*GL - GL*H) 
		*/
		void 				compute_spectral_current() throw (Exception *);
		
		/** compute 
		            dJ/dx(x,E) = 1/hbar *       Tr_perp(SL*GG - GL*SG)_ii 
		    or
			        dJ/dx(x,E) = 1/hbar * 0.5 * Tr_perp(SL*GG - GL*SG + GG*SL - SG*GL)_ii
		*/
		void 				compute_spectral_current2() throw (Exception *);
		
		/** compute 
		           dJ_{scatt}/dx(x,E) = 1/hbar * Tr_perp(SLs*GG - GL*SGs)_ii 
		    where SLGs is only the scattering self-energy */
		void 				compute_scattering_current(SelfEnergyType scattering_type) throw (Exception *);
				
		/** integrate J(x,E) over energy to obtain J(x), calculated from HG-GH */
		void 				compute_current() throw (Exception *);
		
		/** obtain J(x,E), calculated from HG-GH, for the own energy points. size: myNE x Nx-1 */
		const Matd & 		get_spectral_current() const { return this->spectral_current; }
		
		/** obtain J(x,E), calculated from HG-GH, for all energy points. size: NE x Nx-1. master process only! */
		const Matd & 		get_entire_spectral_current() const throw (Exception *);
		
		/** obtain J(x), calculated from HG-GH. size: Nx-1. master process only! */
		const vector<double> & get_current() const throw (Exception *);
		
		/** obtain J(x,E), calculated from SLGG-GGSL, for all energy points. size: NE x Nx. master process only! */
		const Matd & 		get_entire_spectral_current2() const throw (Exception *);
		
		const double &		get_Jleft () const { return this->Jleft; }
		const double &		get_Jright() const { return this->Jright; }
		
		const Matd & 		get_entire_scattering_current(SelfEnergyType scattering_type) const throw (Exception *);
	
	protected:
		
		quantities::PhysicalQuantity e_or_h;			// electron or hole current?
		
		// current computed from H*GL - GL*H
		Matd  		spectral_current;			// size myNE*(Nx-1)
		Matd  		entire_spectral_current;	// size NE*(Nx-1) (root only)
		vector<double>  current;				// size Nx-1
		
		// current computed from SL*GG - GL*SG
		Matd  		spectral_current2;			// size myNE*(Nx-1)
		Matd  		entire_spectral_current2;	// size NE*(Nx-1) (root only)
	
		const Hamiltonian	 * ham;
		const Overlap	     * ov;
		const GreenFunctions * gf;
		const SelfEnergies   * se;
		
		const Geometry * xspace;
		const Kspace   * kspace;
		const Options  * options;
		const Energies * energies;
		
		const uint Nx;
		const uint NxNn;
		
		bool i_am_master;
		
		//ContactCurrent * contact_current;
		
		// matrices to store individual scattering currents (master process only)
		Matd  Jcoh;			// Coherent current (contact SE only)
		double Jleft;		// current at left contact, calculated from coherent current
		double Jright;		// current at left contact, calculated from coherent current
		Matd  Jbuettiker;	// Scattering current from Buettiker probes
		Matd  Jgolim;		// Golizadeh momentum relaxation
		Matd  Jgolip;		// Golizadeh phase relaxation
		Matd  Jpop;			// Polar optical phonon scattering
		Matd  Jac;			// Acoustic phonon scattering
		Matd  Jspont;		// Spontaneous photon emission
		Matd  Jionimp;		// Ionized impurities
		
		cplx J_FisherLee;
	};
	
} // end namespace negf
	
#endif /*CURRENT_H_*/
