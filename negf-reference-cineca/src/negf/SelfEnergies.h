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
#ifndef SELFENERGIES_H_
#define SELFENERGIES_H_

#include "all.h"

#include "Geometry.h"
#include "Kspace.h"
#include "Energies.h"
#include "Options.h"
#include "NEGFObject.h"
#include "SelfEnergy.h"
#include "Hamiltonian.h"
#include "Overlap.h"
#include "MaterialDatabase.h"
#include "GreenFunctions.h"

#include "SEContacts.h"
#include "SEBuettiker.h"
#include "SEGolizadeh.h"
#include "SEOpticalPhonon.h"
#include "SEAcousticPhonon.h"
#include "SEPhotonSpontaneous.h"
#include "SEIonizedImpurities.h"

namespace negf {
		
	enum SelfEnergyType { 
		SEtype_contact, 
		SEtype_buettiker,
		SEtype_golizadeh_momentum,
		SEtype_golizadeh_phase,
		SEtype_optical_phonon,
		SEtype_acoustic_phonon,
		SEtype_spont_photon,
		SEtype_ion_imp,
		SEtype_all
	}; 	
	
	/** Stores and computes the TOTAL retarded, lesser and greater electronic self-energies */
	class SelfEnergies : public SelfEnergy
	{
	public:
	
		SelfEnergies(const Hamiltonian * ham,
					const Overlap * ov,
					const Geometry * xspace_, 
					const Kspace * kspace_, 
					const Energies * energies_, 
					const Options * options_,
					const GreenFunctions * gf_,
					const MaterialDatabase * db_) throw (Exception *);
		~SelfEnergies() {}
		
		void 					calculate() throw (Exception *);
		
		bool 					is_ballistic() const;
		bool 					has_self_energy(SelfEnergyType type) const;
		SelfEnergy * 			get_selfenergy(SelfEnergyType type) const;
		
		void 					initial_guess() throw (Exception *); // calculates contact self-energy
		
		// get access to the individual self-energies (throws exception if not existing)
		SEContacts          *	get_contact_selfenergy() const { return this->contact_selfenergy; }
		SEBuettiker         * 	get_buettiker_selfenergy() const throw (Exception *);
		SEGolizadeh         * 	get_golizadeh_m_selfenergy() const throw (Exception *);
		SEGolizadeh         * 	get_golizadeh_p_selfenergy() const throw (Exception *);
		SEOpticalPhonon     * 	get_optical_phonon_selfenergy() const throw (Exception *);
		SEAcousticPhonon    * 	get_acoustic_phonon_selfenergy() const throw (Exception *);
		SEPhotonSpontaneous * 	get_spontaneous_photon_selfenergy() const throw (Exception *);
		SEIonizedImpurities * 	get_ionized_impurities_selfenergy() const throw (Exception *);
	
	protected:
		
		void add_self_energy(SelfEnergy * selfenergy);
		
		const GreenFunctions * gf;
		const MaterialDatabase * db;
		
		// INDIVIDUAL self-energies
		vector<SelfEnergy *> selfenergies;
		
		// separate pointer for every self-energy
		SEContacts  		* contact_selfenergy;
		SEBuettiker 		* buettiker_selfenergy;
		SEGolizadeh 		* golizadeh_m_selfenergy;
		SEGolizadeh 		* golizadeh_p_selfenergy;
		SEOpticalPhonon 	* optical_phonon_selfenergy;
		SEAcousticPhonon	* acoustic_phonon_selfenergy;
		SEPhotonSpontaneous * spont_photon_selfenergy;
		SEIonizedImpurities * ion_imp_selfenergy;
		 
	};
	
	
} // end of namespace

#endif /*SELFENERGIES_H_*/
