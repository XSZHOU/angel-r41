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
#include "SelfEnergies.h"
using namespace negf;


SelfEnergies::SelfEnergies(const Hamiltonian * ham,
						   const Overlap * ov,
						   const Geometry * xspace_, 
						   const Kspace * kspace_, 
						   const Energies * energies_,
						   const Options * options_,
						   const GreenFunctions * gf_,
						   const MaterialDatabase * db_) throw (Exception *):
	SelfEnergy(xspace_,kspace_,energies_,options_, constants::odSE),
	gf(gf_),
	db(db_)
{STACK_TRACE(
	logmsg->emit(LOG_INFO,"finished allocating memory for total self energy.");
	NEGF_ASSERT(gf!=NULL && db!=NULL, "null pointer encountered.");
	
	// initialize pointers to zero
	this->contact_selfenergy         = NULL;
	this->buettiker_selfenergy       = NULL;
	this->golizadeh_m_selfenergy     = NULL;
	this->golizadeh_p_selfenergy     = NULL;
	this->optical_phonon_selfenergy  = NULL;
	this->acoustic_phonon_selfenergy = NULL;
	this->spont_photon_selfenergy    = NULL;
	this->ion_imp_selfenergy         = NULL;
	this->selfenergies.clear();
	
	// the contact self-energy (for all contacts) is ALWAYS included!
	logmsg->emit_header("setting up contact self-energy");
	this->contact_selfenergy = new SEContacts(ham, ov, xspace, kspace, energies, options);
	this->selfenergies.push_back(contact_selfenergy);
		
	// check whether Buettiker probes are turned on
	if (options->exists("Buettiker") && options->get("Buettiker")==1) {
		logmsg->emit_header("setting up Buettiker self-energy");
		NEGF_ASSERT(options->exists("BuettikerParameter"), "Must have \"BuettikerParameter\"!");
		const double energy_parameter = constants::convert_from_SI(units::energy, constants::SIec * 
							options->get("BuettikerParameter"));
		SEBuettiker * buettiker = new SEBuettiker(ham, ov, xspace, kspace, energies, options,
											gf, this->contact_selfenergy, energy_parameter);
		this->selfenergies.push_back(buettiker);
		this->buettiker_selfenergy = buettiker;
	}
	
	// check for Golizadeh momentum (and phase) relaxation
	if (options->exists("GolizadehMomentumRelaxation") && options->get("GolizadehMomentumRelaxation")==1) {
		logmsg->emit_header("setting up Golizadeh momentum relaxation self-energy");
		NEGF_ASSERT(options->exists("GolizadehMomentumParameter"), "Must have \"GolizadehMomentumParameter\"!");
		const double energy_parameter = negf_math::pow(constants::convert_from_SI(units::energy, constants::SIec), 2.0) * 
							options->get("GolizadehMomentumParameter");
		SEGolizadeh * golizadeh_m = new SEGolizadeh(ham, ov, xspace, kspace, energies, options, gf, energy_parameter, true);
		this->selfenergies.push_back(golizadeh_m);
		this->golizadeh_m_selfenergy = golizadeh_m;
	}
	
	// check for Golizadeh phase relaxation (dephasing)
	if (options->exists("GolizadehDephasing") && options->get("GolizadehDephasing")==1) {
		logmsg->emit_header("setting up Golizadeh dephasing self-energy");
		NEGF_ASSERT(options->exists("GolizadehDephasingParameter"), "Must have \"GolizadehDephasingParameter\"!");
		const double energy_parameter = negf_math::pow(constants::convert_from_SI(units::energy, constants::SIec), 2.0) * 
							options->get("GolizadehDephasingParameter");
		SEGolizadeh * golizadeh_p = new SEGolizadeh(ham, ov, xspace, kspace, energies, options, gf, energy_parameter, false);
		this->selfenergies.push_back(golizadeh_p);
		this->golizadeh_p_selfenergy = golizadeh_p;
	}
	
	// check for optical phonon self-energy
	if (options->exists("OpticalPhonons") && options->get("OpticalPhonons")==1) {
		logmsg->emit_header("setting up optical phonon self-energy");
		SEOpticalPhonon * opt = new SEOpticalPhonon(ov, xspace, kspace, energies, options, gf, db);
		this->selfenergies.push_back(opt);
		this->optical_phonon_selfenergy = opt;
		
		//opt->test_operation();
	}
	
	// check for acoustic phonon self-energy
	if (options->exists("AcousticPhonons") && options->get("AcousticPhonons")==1) {
		logmsg->emit_header("setting up acoustic phonon self-energy");
		SEAcousticPhonon * ac = new SEAcousticPhonon(ov, xspace, kspace, energies, options, gf, db);
		this->selfenergies.push_back(ac);
		this->acoustic_phonon_selfenergy = ac;
	}
	
	// check for photon self-energy (spontaneous emission only)
	if (options->exists("SpontaneousPhotons") && options->get("SpontaneousPhotons")==1) {
		logmsg->emit_header("setting up spontaneous photon self-energy");
		SEPhotonSpontaneous * spont = new SEPhotonSpontaneous(ov, xspace, kspace, energies, options, gf, db);
		this->selfenergies.push_back(spont);
		this->spont_photon_selfenergy = spont;
	}
	
	// check for ionized impurities self-energy
	if (options->exists("IonizedImpurities") && options->get("IonizedImpurities")==1) {
		logmsg->emit_header("setting up ionized impurities self-energy");
		SEIonizedImpurities * ionimp = new SEIonizedImpurities(ov, xspace, kspace, energies, options, gf);
		this->selfenergies.push_back(ionimp);
		this->ion_imp_selfenergy = ionimp;
	}
);}

						   
bool SelfEnergies::is_ballistic() const
{STACK_TRACE(
	if (this->selfenergies.size()==1) { // contact self-energy only
		return true;
	} else {
		return false;
	}
);}
						   
bool SelfEnergies::has_self_energy(SelfEnergyType type) const
{STACK_TRACE(
	switch(type) {
	case SEtype_contact:   			return true; break;
	case SEtype_buettiker: 			return (this->buettiker_selfenergy==NULL) ? false : true; break;
	case SEtype_golizadeh_momentum: return (this->golizadeh_m_selfenergy==NULL) ? false : true; break;
	case SEtype_golizadeh_phase: 	return (this->golizadeh_p_selfenergy==NULL) ? false : true; break;
	case SEtype_optical_phonon: 	return (this->optical_phonon_selfenergy==NULL) ? false : true; break;
	case SEtype_acoustic_phonon: 	return (this->acoustic_phonon_selfenergy==NULL) ? false : true; break;
	case SEtype_spont_photon: 		return (this->spont_photon_selfenergy==NULL) ? false : true; break;
	case SEtype_ion_imp:	 		return (this->ion_imp_selfenergy==NULL) ? false : true; break;
	default: NEGF_EXCEPTION("Self-energy not treated."); return false;
	}
);}
						   
SelfEnergy * SelfEnergies::get_selfenergy(SelfEnergyType type) const
{STACK_TRACE(
	NEGF_ASSERT(this->has_self_energy(type), "Self energy does not exist.");
	switch(type) {
	case SEtype_contact:   			return this->contact_selfenergy; break;
	case SEtype_buettiker: 			return this->buettiker_selfenergy; break;
	case SEtype_golizadeh_momentum: return this->golizadeh_m_selfenergy; break;
	case SEtype_golizadeh_phase: 	return this->golizadeh_p_selfenergy; break;
	case SEtype_optical_phonon: 	return this->optical_phonon_selfenergy; break;
	case SEtype_acoustic_phonon: 	return this->acoustic_phonon_selfenergy; break;
	case SEtype_spont_photon: 		return this->spont_photon_selfenergy; break;
	case SEtype_ion_imp: 			return this->ion_imp_selfenergy; break;
	case SEtype_all: 				NEGF_EXCEPTION("Cannot access total self-energy."); break;
	default: NEGF_EXCEPTION("Self-energy not treated."); return false;
	}
);}

SEBuettiker * SelfEnergies::get_buettiker_selfenergy() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->has_self_energy(SEtype_buettiker), "Buettiker was not turned on.");
	return this->buettiker_selfenergy;
);}

SEGolizadeh * SelfEnergies::get_golizadeh_m_selfenergy() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->has_self_energy(SEtype_golizadeh_momentum), "Golizadeh momentum relaxation was not turned on.");
	return this->golizadeh_m_selfenergy;
);}

SEGolizadeh * SelfEnergies::get_golizadeh_p_selfenergy() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->has_self_energy(SEtype_golizadeh_phase), "Golizadeh dephasing was not turned on.");
	return this->golizadeh_p_selfenergy;
);}

SEOpticalPhonon * SelfEnergies::get_optical_phonon_selfenergy() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->has_self_energy(SEtype_optical_phonon), "Optical phonons were not turned on.");
	return this->optical_phonon_selfenergy;
);}

SEAcousticPhonon * SelfEnergies::get_acoustic_phonon_selfenergy() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->has_self_energy(SEtype_acoustic_phonon), "Acoustic phonons were not turned on.");
	return this->acoustic_phonon_selfenergy;
);}

SEPhotonSpontaneous * SelfEnergies::get_spontaneous_photon_selfenergy() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->has_self_energy(SEtype_spont_photon), "Sponanteous photons were not turned on.");
	return this->spont_photon_selfenergy;
);}

SEIonizedImpurities * SelfEnergies::get_ionized_impurities_selfenergy() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->has_self_energy(SEtype_ion_imp), "Ionized impurities were not turned on.");
	return this->ion_imp_selfenergy;
);}

void SelfEnergies::calculate() throw (Exception *)
{STACK_TRACE(
	logmsg->emit_header("calculating self energies");
	
	// *** DAMPING FACTOR: Sigma = (1-alpha)*Sigma_new + alpha*Sigma_old ***
	bool damping = false;
	double alpha = 0.0;
	if (options->exists("SelfEnergyUnderrelaxation")) {
		alpha = options->get("SelfEnergyUnderrelaxation");
		NEGF_ASSERT(alpha>=0.0 && alpha<=1.0, "expected self energy underrelaxation to be between 0 and 1.");
		if (alpha>0.0) {
			damping = true;
		}
	}
	
	// calculate individual self-energies
	for (uint ii = 0; ii < selfenergies.size(); ii++) {
		
		if (damping) {
			for (uint ee2 = 0; ee2 < energies->get_my_number_of_points(); ee2++) {	
				uint ee = energies->get_global_index(ee2);
				for (uint kk = 0; kk < Nk; kk++) {	
					SEMat & SR = this->SigmaR->get(kk, ee);
					SEMat & SL = this->SigmaL->get(kk, ee);
					SEMat & SG = this->SigmaG->get(kk, ee);
					if (ii==0) {// boundary self-energy
						SR = SEMat_create(NxNn);
						SL = SEMat_create(NxNn);
						SG = SEMat_create(NxNn);
					} else {
						add(selfenergies[ii]->get_retarded(kk, ee), alpha, SR);
						add(selfenergies[ii]->get_lesser  (kk, ee), alpha, SL);
						add(selfenergies[ii]->get_greater (kk, ee), alpha, SG);
					}
				}
			}
		}
			
		double t1 = MPI_Wtime();
		selfenergies[ii]->calculate();
		double t2 = MPI_Wtime();
		logmsg->emit(LOG_INFO,"Needed %.2e[s] for the entire self-energy.", t2-t1);
		
		if (damping) {
			for (uint ee2 = 0; ee2 < energies->get_my_number_of_points(); ee2++) {	
				uint ee = energies->get_global_index(ee2);
				for (uint kk = 0; kk < Nk; kk++) {	
					SEMat & SR = this->SigmaR->get(kk, ee);
					SEMat & SL = this->SigmaL->get(kk, ee);
					SEMat & SG = this->SigmaG->get(kk, ee);
					if (ii==0) {	// boundary self-energy
						SR += selfenergies[ii]->get_retarded(kk, ee);
						SL += selfenergies[ii]->get_lesser  (kk, ee);
						SG += selfenergies[ii]->get_greater (kk, ee);
					} else {
						add(selfenergies[ii]->get_retarded(kk, ee), 1.0-alpha, SR);
						add(selfenergies[ii]->get_lesser  (kk, ee), 1.0-alpha, SL);
						add(selfenergies[ii]->get_greater (kk, ee), 1.0-alpha, SG);
					}
				}
			}
		}
	}
	mpi->synchronize_processes();
	
	if (damping) {
		logmsg->emit(LOG_INFO,"Underrelaxed scattering self energies by S=%g*S_new+%g*S_old",1-alpha,alpha);
	} else {
		logmsg->emit_small_header("adding up self energies");
		
		// add everything to the total self-energy (for the energies of the calling process)
		for (uint ee2 = 0; ee2 < energies->get_my_number_of_points(); ee2++)
		{	
			uint ee = energies->get_global_index(ee2);
			logmsg->emit_noendl_all(LOG_INFO_L2,"  p%d:E=%d...",mpi->get_rank(),ee);
			for (uint kk = 0; kk < Nk; kk++)
			{	
				SEMat & SR = this->SigmaR->get(kk, ee);
				SEMat & SL = this->SigmaL->get(kk, ee);
				SEMat & SG = this->SigmaG->get(kk, ee);
				
				// re-initialize to zero
				SR = SEMat_create(NxNn);
				SL = SEMat_create(NxNn);
				SG = SEMat_create(NxNn);
				// add
				for (uint ii = 0; ii < selfenergies.size(); ii++) {
					SR += selfenergies[ii]->get_retarded(kk, ee);
					SL += selfenergies[ii]->get_lesser(kk, ee);
					SG += selfenergies[ii]->get_greater(kk, ee);
				}
			}
		}
	}
	mpi->synchronize_processes();
);}


void SelfEnergies::add_self_energy(SelfEnergy * selfenergy)
{STACK_TRACE(
	// security checks
	NEGF_ASSERT(selfenergy!=NULL && selfenergy->get_xspace()==this->xspace && selfenergy->get_kspace()==this->kspace 
			&& selfenergy->get_energies()==this->energies && selfenergy->get_options()==this->options,
			"self energy cannot be added.");
	for (uint ii=0; ii < this->selfenergies.size(); ii++) {
		NEGF_ASSERT(selfenergies[ii]!=selfenergy, "self energy was already added.");
	}
	this->selfenergies.push_back(selfenergy);
);}


void SelfEnergies::initial_guess() throw (Exception *)
{STACK_TRACE(
	logmsg->emit_header("Initial guess for self-energies");
	
	this->contact_selfenergy->calculate();
	
	// add contact SE to the total SE (for the energies of the calling process)
	logmsg->emit(LOG_INFO,"Assigning contact self-energy to total self-energy...");
	for (uint ee2 = 0; ee2 < energies->get_my_number_of_points(); ee2++)
	{	
		uint ee = energies->get_global_index(ee2);
		for (uint kk = 0; kk < Nk; kk++)
		{	
			this->SigmaR->get(kk, ee) = this->contact_selfenergy->get_retarded(kk, ee);
			this->SigmaL->get(kk, ee) = this->contact_selfenergy->get_lesser(kk, ee);
			this->SigmaG->get(kk, ee) = this->contact_selfenergy->get_greater(kk, ee);
		}
	}
);}
