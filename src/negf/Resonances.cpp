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
#include "Resonances.h"
using namespace negf;


Resonances::Resonances(InnerLoop * inner) throw (Exception *):
	linear(false),	// switches between linear and Fermi-like refinement around Fermilevels. IWCE example only converges w/ Fermi
	lin_ref_num(2.0),
	lin_delta(0.1),
	energies(inner->get_energies()),
	options (inner->get_options()),
	xspace  (inner->get_xspace()),
	kspace  (inner->get_kspace()),
	gf      (inner->get_green_functions()),
	se      (inner->get_self_energies()),
	ham     (inner->get_hamiltonian()),
	pp      (inner->get_post_processing())
{STACK_TRACE(
	NEGF_ASSERT(energies!=NULL && options!=NULL && xspace!=NULL && kspace!=NULL
			 && gf!=NULL && se!=NULL && ham!=NULL && pp!=NULL, "null pointer encountered.");
	this->fermilevels.resize(2, 0.0);
	this->pmlgrid = NULL;
	
	// set up some polynomial coefficients (see mathematica file for the calculation of these)
	const double lin_delta2 = lin_delta*lin_delta; 
	const double lin_delta3 = lin_delta*lin_delta*lin_delta; 
	const double lin_delta4 = lin_delta*lin_delta*lin_delta*lin_delta; 
	const double lin_ref_num2 = lin_ref_num*lin_ref_num; 
	const double lin_ref_num3 = lin_ref_num*lin_ref_num*lin_ref_num; 
	const double lin_ref_num4 = lin_ref_num*lin_ref_num*lin_ref_num*lin_ref_num; 
	/* this->ca = 0.0; 
	this->cb = 0.0; 
	this->cc = 0.0;
	this->cd = 1.0 / (8.0 * lin_ref_num2 * lin_delta);
	this->ce = (1.0 + lin_delta) / (4.0 * lin_ref_num * lin_delta);
	this->cf = (1.0 + 2.0*lin_delta + lin_delta2) / (8.0 * lin_delta);	*/
	this->ca = 0.0;
	this->cb =                                                     -1.0 / (32.0 * lin_ref_num4 * lin_delta3);
	this->cc =                                                     -1.0 / ( 8.0 * lin_ref_num3 * lin_delta3);
	this->cd =                                   3.0*(lin_delta2 - 1.0) / (16.0 * lin_ref_num2 * lin_delta3);
	this->ce =                  (2.0*lin_delta3 + 3.0*lin_delta2 - 1.0) / ( 8.0 * lin_ref_num  * lin_delta3);
	this->cf = (3.0*lin_delta4 + 8.0*lin_delta3 + 6.0*lin_delta2 - 1.0) / (32.0 * lin_delta3);
	
	if (options->exists("PMLResonances") && fabs(options->get("PMLResonances")-1.0)<1e-10) {
		NEGF_ASSERT(fabs(options->get("kp_method")-0.0)>1e-10 && fabs(options->get("kp_method")-1.0)>1e-10
					&& fabs(options->get("kp_method")-3.0)>1e-10, "The KP model must be one with a TDKP Hamiltonian (not 0,1 or 3)");
#ifndef NOTDKP
		const tdkp::InterfaceConfiguration  & config = ham->get_tdkp_config();
		this->infodesk = ham->get_tdkp_infodesk();
		string pmlgridfile = fnames->get_filename();
		pmlgridfile.append("_pml");
		
		InputParser parser;
		this->pmlgrid = parser.read_dfise_grd(pmlgridfile.c_str());
		pmlgridfile.append(".grd");
		char logf[1000];
		sprintf(logf,"%sthread%d.PMLlog", fnames->get_outfiledirectory().c_str(),mpi->get_rank());
		string logfilename(logf);
		this->pml = tdkp::InterfaceFactory::create_well_radial_pml(logfilename, config, *infodesk, pmlgridfile.c_str());
#endif
	}
);}


void Resonances::determine_new_energy_grid(bool divergence_avoiding_only, const vector<double> & fermilevels_) throw (Exception *)
{STACK_TRACE(
	logmsg->emit_header("Determining new energy grid");
	this->fermilevels = fermilevels_;
	
	uint nE = this->energies->get_number_of_points();
	double Emin = energies->get_energy_from_global_idx(0);
	double Emax = energies->get_energy_from_global_idx(nE-1);
	double kT = constants::convert_from_SI(units::energy, constants::SIkb * this->options->get("temperature"));
	
	vector<double> new_energy_grid;
	new_energy_grid.resize(nE, 0.0);
	
	vector<bool> new_ignore_array;
	new_ignore_array.resize(nE, false);
	
	bool midres_k0   = (options->exists("Resonances")     && fabs(options->get("Resonances")    -1.0)<1e-10) ? true : false;
	bool midres_kall = (options->exists("ResonancesAllK") && fabs(options->get("ResonancesAllK")-1.0)<1e-10) ? true : false;
	bool pmlres      = (options->exists("PMLResonances")  && fabs(options->get("PMLResonances") -1.0)<1e-10) ? true : false;
	
	// the new energy grid is determined solely by the master process and then communicated via MPI
	// otherwise we could run into problems e.g. when the Newton procedure stops very close to the 
	// convergence criterion, which is also influenced by numerical noise...
	if (mpi->get_rank()==constants::mpi_master_rank)
	{
		new_energy_grid[0]    = Emin;
		new_energy_grid[nE-1] = Emax;
		
		// ----------------------------------------------------------------------------------------------------
		// if the option is turned on, find RTD resonance by looking at LDOS in the middle of the structure
		// ----------------------------------------------------------------------------------------------------
		this->resonances.clear();
		this->broadenings.clear();
		if (midres_k0 && !divergence_avoiding_only) 
		{
			// determine resonances using LDOS@k=0 in the middle of the structure. 
			// will assign resonances to vectors "resonances", "broadenings"
			this->determine_ldos_resonances();
		}
		if (midres_kall && !divergence_avoiding_only)
		{
			// determine resonances using LDOS at all k-vectors in the middle of the structure. 
			// will assign resonances to vectors "resonances", "broadenings"
			this->determine_all_ldos_resonances();
		}
		if (pmlres && !divergence_avoiding_only) 
		{
			// determine resonances using PMLs. will assign resonances to vectors "resonances", "broadenings"
			this->determine_pml_resonances();
		}
		
		// find 1D-LDOS-divergences in the lead Hamiltonians which must be avoided --> array "divergences"
		if ((!options->exists("IncoherentContacts")) || 
			( options->exists("IncoherentContacts") && options->exists("IncoherentContacts")==0) ||
			( options->exists("IncoherentContacts") && options->get("IncoherentContactBroadening")==0.0) ) {
			this->determine_divergences();
		}
		
		double mmin = this->monotonic_function(Emin);
		double mmax = this->monotonic_function(Emax);
		//uint num_contributions = 3 + resonances.size();
		
		// energy must not be in interval [divergence, divergence+Ecrit]
		const double Ecrit  = constants::convert_from_SI(units::energy, constants::SIec * 
								constants::divergence_avoid); 
		// energy will be replaced by divergence-dE_div
		const double dE_div = constants::convert_from_SI(units::energy, constants::SIec * 
								constants::divergence_distance);
			
		for (uint ii=1; ii < nE-1; ii++) 
		{
			// initial guess: as if monotonic_function was linear
			double E = Emin + ii * (Emax-Emin)/(nE-1);
			
			LoggerLevel level = LOG_INFO_L3;
			if (divergence_avoiding_only) {
				E = energies->get_energy_from_global_idx(ii);
			} else {
				logmsg->emit(level,"ii=%d:",ii);
				// we search for energy E such that monotonic_function(E) = mmin + ii * (mmax-mmin)/(nE-1)
				// we do this by Newton iteration
				const uint max_iter = 10000;
				double dE_max = constants::convert_from_SI(units::energy, constants::SIec * 0.001); 
				const double desired_result = mmin + ii * (mmax-mmin)/(nE-1);
				if (broadenings.size()>0) dE_max = min(dE_max, broadenings[0]); // in case of resonances need smaller dE_max!
				
				double dE = 1e10;
				uint iter = 0;
				while (fabs(dE) > constants::convert_from_SI(units::energy, constants::SIec * 1e-6) && iter < max_iter) {
					double    F = this->monotonic_function(E) - desired_result;
					double dFdE = this->monotonic_function_derivative(E);
					dE = -F/dFdE;
					// limit update
					if (fabs(dE) > dE_max) {
						logmsg->emit_noendl(level,"limiting dE=+-kT instead of %e...",dE);
						dE = negf_math::sign(dE) * dE_max;
					}
					
					E += dE;
					logmsg->emit/*_noendl*/(level,"now E=%.6e, dE=%.2e, F=%.3e, dFdE=%.3e    ",E,dE,F,dFdE);
					iter++;
				}
				NEGF_FASSERT(iter<max_iter-1, "Newton for new energy grid did not converge for energy point %d: E=%e, dE=%e, desired_result=%e, function=%e, dFdE=%e",
						ii, E, dE, desired_result, this->monotonic_function(E), this->monotonic_function_derivative(E));
			}
			
			for (uint dd=0; dd < this->divergences.size(); dd++) {
				//if (E>divergences[dd] && fabs(E-divergences[dd])<Ecrit) {	// only E>divergences[dd] is a problem!
				if (fabs(E-divergences[dd])<Ecrit) {	
				// E<divergences[dd] might also be a problem because resonance might be a little broadened.
				// AND IN THE VALENCE BAND, E<divergences[dd] IS THE PROBLEM!!!
					double new_value = divergences[dd] - dE_div;
					if (E<1.0) {
						new_value = divergences[dd] + dE_div;
					}
					if (ii==0 || (ii>0 && new_energy_grid[ii-1] < new_value)) {
						logmsg->emit(LOG_INFO,"Replaced energy %.8g with %.8g due to divergence at %.8g", E, divergences[dd]-dE_div, divergences[dd]);
						E = new_value;
					} else {
						logmsg->emit(LOG_WARN,"Warning: energy %.8g is close to divergence %.8g and cannot be replaced!", E, divergences[dd]);
						//if (ii>0) logmsg->emit(LOG_WARN,"         will ignore this point in the energy integration. next-lower point: %.8g", new_energy_grid[ii-1]);
						//new_ignore_array[ii] = true;
						if (ii>0) {
							E = new_energy_grid[ii-1] + constants::convert_from_SI(units::energy, constants::SIec * 1e-9);
							logmsg->emit(LOG_WARN,"         Setting it to %.8g.", E);
						}
					}
				}
			}
			new_energy_grid[ii] = E;
			logmsg->emit(level,"");
		}
						
		// send result to all other threads
		for (int proc=0; proc<mpi->get_num_procs(); proc++) {
			if (proc==constants::mpi_master_rank) continue;
			int tag  = 654;
			int dest = proc;
			mpi->send(new_energy_grid, dest, tag);
			mpi->send(new_ignore_array, dest, tag);
		}
	} else {
		// determine_all_ldos_resonances() relies on all processes
		if (midres_kall && !divergence_avoiding_only) {
			this->determine_all_ldos_resonances();
		}
		
		// receive resonances, broadenings from master thread
		int tag  = 654;
		int source = constants::mpi_master_rank;
		mpi->recv(new_energy_grid, nE, source, tag);
		mpi->recv(new_ignore_array, nE, source, tag);
	}
	mpi->synchronize_processes();		
	
	energies->set_ignore_array(new_ignore_array);	// needs to be called BEFORE assign_new_energy_grid because of the wights determination in there
	energies->assign_new_energy_grid(new_energy_grid);
);}


double Resonances::monotonic_function(const double & E)
{STACK_TRACE(
	uint nE = this->energies->get_number_of_points();
	double Emin = energies->get_energy_from_global_idx(0);
	double Emax = energies->get_energy_from_global_idx(nE-1);
	double kT = constants::convert_from_SI(units::energy, constants::SIkb * this->options->get("temperature"));

	kT = constants::kT_broad_hack * kT;
	
	/*if (mpi->get_rank()==0) {
		for (uint ii=0; ii<1000; ii++) {
			double nu = -1.2*lin_ref_num + double(ii)/1000.0 * 2.4*lin_ref_num;
			double tmp = 888.888;
			int count = -1;
			if (nu <= -(1+lin_delta)*lin_ref_num) { tmp = 0.0; count=0; }
			if (nu >=  (1+lin_delta)*lin_ref_num) { tmp = 1.0; count=1; }
			if (nu >= -(1-lin_delta)*lin_ref_num && nu <= (1-lin_delta)*lin_ref_num) { tmp = (nu+lin_ref_num) / (2.0*lin_ref_num); count=2; }
			if (nu >  -(1+lin_delta)*lin_ref_num && nu < -(1-lin_delta)*lin_ref_num) { tmp = ca*nu*nu*nu*nu*nu + cb*nu*nu*nu*nu + cc*nu*nu*nu + cd*nu*nu + ce*nu +       cf; count=4; }
			if (nu >   (1-lin_delta)*lin_ref_num && nu <  (1+lin_delta)*lin_ref_num) { tmp = ca*nu*nu*nu*nu*nu - cb*nu*nu*nu*nu + cc*nu*nu*nu - cd*nu*nu + ce*nu + (1.0-cf); count=5; }
			cout << tmp << "(" << count << ")   ";
		}
		
		cout << endl;
		for (uint ii=0; ii<1000; ii++) {
			double nu = -1.2*lin_ref_num + double(ii)/1000.0 * 2.4*lin_ref_num;
			double tmp = 888.888;
			int i = -1;
			if (nu <= -(1+lin_delta)*lin_ref_num) { tmp = 0.0; i=0; }
			if (nu >=  (1+lin_delta)*lin_ref_num) { tmp = 0.0; i=1; }
			if (nu >= -(1-lin_delta)*lin_ref_num && nu <= (1-lin_delta)*lin_ref_num) { tmp = 1.0 / (2.0*lin_ref_num); i=2; }
			if (nu >  -(1+lin_delta)*lin_ref_num && nu < -(1-lin_delta)*lin_ref_num) { tmp = ca*5.0*nu*nu*nu*nu + cb*4.0*nu*nu*nu + cc*3.0*nu*nu + cd*2.0*nu + ce; i=3; }
			if (nu >   (1-lin_delta)*lin_ref_num && nu <  (1+lin_delta)*lin_ref_num) { tmp = ca*5.0*nu*nu*nu*nu - cb*4.0*nu*nu*nu + cc*3.0*nu*nu - cd*2.0*nu + ce; i=4; }
			cout << tmp << "(" << i << ")   ";
		}
		
	}
	mpi->synchronize_processes();
	NEGF_EXCEPTION("Stop.");*/
	
	// ---------------------------------------------------------------
	// first contribution to monotonic function: constant slope 0...1!
	// ---------------------------------------------------------------
	double result = E / (Emax - Emin);
	
	//return result; // uncomment this line for a uniform energy grid
	
	// ---------------------------------------------------------------
	// second contribution to monotonic function: 1-f_L
	// ---------------------------------------------------------------
	double x2 = constants::egrid_weight_fL;
	NEGF_ASSERT(this->xspace->get_num_contacts()==2, "2 contacts expected!");
	double nu = (E - this->fermilevels[0]) / kT;
	if (linear) {
		double tmp = 888.888;
		if (nu <= -(1+lin_delta)*lin_ref_num) { tmp = 0.0; }
		if (nu >=  (1+lin_delta)*lin_ref_num) { tmp = 1.0; }
		if (nu >= -(1-lin_delta)*lin_ref_num && nu <= (1-lin_delta)*lin_ref_num) { tmp = (nu+lin_ref_num) / (2.0*lin_ref_num); }
		if (nu >  -(1+lin_delta)*lin_ref_num && nu < -(1-lin_delta)*lin_ref_num) { tmp = ca*nu*nu*nu*nu*nu + cb*nu*nu*nu*nu + cc*nu*nu*nu + cd*nu*nu + ce*nu +       cf; }
		if (nu >   (1-lin_delta)*lin_ref_num && nu <  (1+lin_delta)*lin_ref_num) { tmp = ca*nu*nu*nu*nu*nu - cb*nu*nu*nu*nu + cc*nu*nu*nu - cd*nu*nu + ce*nu + (1.0-cf); }
		NEGF_ASSERT(tmp!=888.888, "a case was not considered.");
		result += x2*tmp;
	} else {
		result += x2 * (1.0 - 1.0/(1.0 + negf_math::exp(nu)));
	}
		
	// ---------------------------------------------------------------
	// third contribution to monotonic function: 1-f_R
	// ---------------------------------------------------------------
	double x3 = constants::egrid_weight_fR;
	nu = (E - this->fermilevels[1]) / kT;
	if (linear) {
		double tmp = 888.888;
		if (nu <= -(1+lin_delta)*lin_ref_num) { tmp = 0.0; }
		if (nu >=  (1+lin_delta)*lin_ref_num) { tmp = 1.0; }
		if (nu >= -(1-lin_delta)*lin_ref_num && nu <= (1-lin_delta)*lin_ref_num) { tmp = (nu+lin_ref_num) / (2.0*lin_ref_num); }
		if (nu >  -(1+lin_delta)*lin_ref_num && nu < -(1-lin_delta)*lin_ref_num) { tmp = ca*nu*nu*nu*nu*nu + cb*nu*nu*nu*nu + cc*nu*nu*nu + cd*nu*nu + ce*nu +       cf; }
		if (nu >   (1-lin_delta)*lin_ref_num && nu <  (1+lin_delta)*lin_ref_num) { tmp = ca*nu*nu*nu*nu*nu - cb*nu*nu*nu*nu + cc*nu*nu*nu - cd*nu*nu + ce*nu + (1.0-cf); }
		NEGF_ASSERT(tmp!=888.888, "a case was not considered.");
		result += x3*tmp;
	} else {
		result += x3 * (1.0 - 1.0/(1.0 + negf_math::exp(nu)));		
	}
	
	// -----------------------------------------------------------------------
	// fourth contribution to monotonic function: density-of-states-related
	// not yet implemented!
	// -----------------------------------------------------------------------
	
	// -----------------------------------------------------------------------
	// fifth contribution to monotonic function: resonances of Hamiltonian
	// -----------------------------------------------------------------------
	for (uint ii=0; ii < this->resonances.size(); ii++) {
		const double & Eres = this->resonances[ii];
		const double & Gam = this->broadenings[ii];
		
		// Integral of 1/(pi*Gam) Gam^2/((E' - E0)^2 + Gam^2) from -infty to E
		// gives 1/(2pi) * [ 2arctan((E-E0)/Gam) + i*log(-i/Gam)-log(i/Gam) ] = 1/(2pi) * [ -2arctan((4-E)/Gam) + pi ]
		double x4 = constants::egrid_weight_res / resonances.size();
		result += x4 * (1.0/constants::pi * negf_math::atan((E-Eres)/Gam) + 1.0/2.0);
	}
	
	return result;
);}


double Resonances::monotonic_function_derivative(const double & E)
{STACK_TRACE(
	uint nE = this->energies->get_number_of_points();
	double Emin = energies->get_energy_from_global_idx(0);
	double Emax = energies->get_energy_from_global_idx(nE-1);
	double kT = constants::convert_from_SI(units::energy, constants::SIkb * this->options->get("temperature"));
	
	kT = constants::kT_broad_hack * kT;
	
	// ---------------------------------------------------------------
	// first contribution to monotonic function: constant slope 0...1!
	// ---------------------------------------------------------------
	double result = 1.0 / (Emax - Emin);
	
	//return result; // uncomment this line for a uniform energy grid
	
	// ---------------------------------------------------------------
	// second contribution to monotonic function: 1-f_L
	// ---------------------------------------------------------------
	double x2 = constants::egrid_weight_fL;
	NEGF_ASSERT(this->xspace->get_num_contacts()==2, "2 contacts expected!");
	double nu = (E - this->fermilevels[0]) / kT;
	if (linear) {
		double tmp = 888.888;
		int i = -1;
		if (nu <= -(1+lin_delta)*lin_ref_num) { tmp = 0.0; i=0; }
		if (nu >=  (1+lin_delta)*lin_ref_num) { tmp = 0.0; i=1; }
		if (nu >= -(1-lin_delta)*lin_ref_num && nu <= (1-lin_delta)*lin_ref_num) { tmp = 1.0 / (2.0*lin_ref_num); i=2; }
		if (nu >  -(1+lin_delta)*lin_ref_num && nu < -(1-lin_delta)*lin_ref_num) { tmp = ca*5.0*nu*nu*nu*nu + cb*4.0*nu*nu*nu + cc*3.0*nu*nu + cd*2.0*nu + ce; i=3; }
		if (nu >   (1-lin_delta)*lin_ref_num && nu <  (1+lin_delta)*lin_ref_num) { tmp = ca*5.0*nu*nu*nu*nu - cb*4.0*nu*nu*nu + cc*3.0*nu*nu - cd*2.0*nu + ce; i=4; }
		NEGF_ASSERT(fabs(tmp-888.888)>1e-13, "a case was not considered.");
		logmsg->emit_noendl(LOG_INFO_L3,"f_L contrib (%d): %.3e   ", i, x2*tmp);
		result += x2*tmp * 1.0/kT;
	} else {
		result += x2 * 1.0/(2.0 + negf_math::exp(nu) + negf_math::exp(-nu)) * 1.0/kT;
	}
	
	// ---------------------------------------------------------------
	// third contribution to monotonic function: 1-f_R
	// ---------------------------------------------------------------
	double x3 = constants::egrid_weight_fR;
	nu = (E - this->fermilevels[1]) / kT;
	if (linear) {
		double tmp = 888.888;
		int i=-1;
		if (nu <= -(1+lin_delta)*lin_ref_num) { tmp = 0.0; i=0; }
		if (nu >=  (1+lin_delta)*lin_ref_num) { tmp = 0.0; i=1; }
		if (nu >= -(1-lin_delta)*lin_ref_num && nu <= (1-lin_delta)*lin_ref_num) { tmp = 1.0 / (2.0*lin_ref_num); i=2; }
		if (nu >  -(1+lin_delta)*lin_ref_num && nu < -(1-lin_delta)*lin_ref_num) { tmp = ca*5.0*nu*nu*nu*nu + cb*4.0*nu*nu*nu + cc*3.0*nu*nu + cd*2.0*nu + ce; i=3; }
		if (nu >   (1-lin_delta)*lin_ref_num && nu <  (1+lin_delta)*lin_ref_num) { tmp = ca*5.0*nu*nu*nu*nu - cb*4.0*nu*nu*nu + cc*3.0*nu*nu - cd*2.0*nu + ce; i=4; }
		NEGF_ASSERT(fabs(tmp-888.888)>1e-13, "a case was not considered.");
		logmsg->emit_noendl(LOG_INFO_L3,"f_R contrib (%d): %.3e   ", i, x3*tmp);
		result += x3*tmp * 1.0/kT;
	} else {
		result += x3 * 1.0/(2.0 + negf_math::exp(nu) + negf_math::exp(-nu)) * 1.0/kT;
	}
	
	// -----------------------------------------------------------------------
	// fourth contribution to monotonic function: density-of-states-related
	// not yet implemented!
	// -----------------------------------------------------------------------
	
	// -----------------------------------------------------------------------
	// fifth contribution to monotonic function: resonances of Hamiltonian
	// -----------------------------------------------------------------------
	for (uint ii=0; ii < this->resonances.size(); ii++) {
		const double & Eres = this->resonances[ii];
		const double & Gam = this->broadenings[ii];
		
		// derivative (wrt E) of 1/(2pi) * [ 2arctan((E-Eres)/Gam) + pi ]
		// --> 1/pi * arctan'((E-Eres)/Gam)*(1/Gam) = 1/(pi*Gam)*1.0/(1.0+((E-Eres)/Gam)^2)
		double x4 = constants::egrid_weight_res / resonances.size();
		double tmp = (E-Eres)/Gam;
		result += x4 * 1.0/(Gam*constants::pi) * 1.0/(1.0+tmp*tmp);
	}
		
	return result;
);}


/** determine energetic position of divergences due to singularity of 1D-DOS for a given k-vector
    part of this code is copied form SEContacts.cpp
    called by master thread only */
void Resonances::determine_divergences()
{STACK_TRACE(
	this->divergences.clear();

	uint Nk    = this->kspace->get_number_of_points();
	uint Nvert = this->xspace->get_num_vertices();
	const OVMat & M = this->gf->get_overlap()->get_overlap(); // defined on all indices (not just internal), (Nvert*Nn)^2
	NEGF_ASSERT(M.num_rows()==Nvert*Nn && M.num_cols()==Nvert*Nn, "inconsistent overlap matrix size.");
	for (uint kk = 0; kk < Nk; kk++)
	{	
		// get Hamiltonian including energy from transversal k-vector and electrostatic potential!
		Matc H(Nvert*Nn,Nvert*Nn);
		this->ham->get(kspace->get_point(kk), H);	
		NEGF_ASSERT(H.num_rows()==Nvert*Nn && H.num_cols()==Nvert*Nn, "inconsistent Hamiltonian size.");
		
		for (uint cc=0; cc < this->xspace->get_num_contacts(); cc++) 
		{
			// ---------------------------------------------------------------------------------------------------------
			// construct the matrix Hc of size 2NcNn*2NcNn, Nc = number of vertices at the interface contact-device
			// Hc(      1:NcNn,      1:NcNn) will be denoted H00
			// Hc(NcNn+1:2NcNn,      1:NcNn) will be denoted H10
			// Hc(      1:NcNn,NcNn+1:2NcNn) will be denoted H01
			// Hc(NcNn+1:2NcNn,NcNn+1:2NcNn) will be denoted H11
			// ---------------------------------------------------------------------------------------------------------
			//uint cidx = 0; <ss 7.7.09> that's a bug, isn't it? 
			uint cidx = cc;
			uint Nc = this->se->get_contact_selfenergy()->get_interface_vertices(cidx).size();
			Matc Hc(2*Nc*Nn, 2*Nc*Nn);
			Matc Mc(2*Nc*Nn, 2*Nc*Nn);
			for (uint ii=1; ii<=Nc; ii++) 
			{
				uint gii1 = this->se->get_contact_selfenergy()->get_interface_vertices(cidx)[ii-1]->get_index_global() + 1;
				uint gii2 = this->se->get_contact_selfenergy()->get_second_row_vertices(cidx)[ii-1]->get_index_global() + 1;
				
				// ordering is always (xx-1)*Nn+nn in Hc, Mc!
				for (uint mm=1; mm<=Nn; mm++) {
					for (uint nn=1; nn<=Nn; nn++) {
						Hc(      (ii-1)*Nn+mm,       (ii-1)*Nn+nn) = H(get_mat_idx(gii1,mm,Nvert), get_mat_idx(gii1,nn,Nvert));
						Hc(      (ii-1)*Nn+mm, Nc*Nn+(ii-1)*Nn+nn) = H(get_mat_idx(gii1,mm,Nvert), get_mat_idx(gii2,nn,Nvert));
						Hc(Nc*Nn+(ii-1)*Nn+mm,       (ii-1)*Nn+nn) = H(get_mat_idx(gii2,mm,Nvert), get_mat_idx(gii1,nn,Nvert));
						Hc(Nc*Nn+(ii-1)*Nn+mm, Nc*Nn+(ii-1)*Nn+nn) = H(get_mat_idx(gii2,mm,Nvert), get_mat_idx(gii2,nn,Nvert));
						
						Mc(      (ii-1)*Nn+mm,       (ii-1)*Nn+nn) = M(get_mat_idx(gii1,mm,Nvert), get_mat_idx(gii1,nn,Nvert));
						Mc(      (ii-1)*Nn+mm, Nc*Nn+(ii-1)*Nn+nn) = M(get_mat_idx(gii1,mm,Nvert), get_mat_idx(gii2,nn,Nvert));
						Mc(Nc*Nn+(ii-1)*Nn+mm,       (ii-1)*Nn+nn) = M(get_mat_idx(gii2,mm,Nvert), get_mat_idx(gii1,nn,Nvert));
						Mc(Nc*Nn+(ii-1)*Nn+mm, Nc*Nn+(ii-1)*Nn+nn) = M(get_mat_idx(gii2,mm,Nvert), get_mat_idx(gii2,nn,Nvert));
					}
				}
			}			
			Matc Hc_00(Nc*Nn,Nc*Nn); Hc.get_submatrix(      1,   Nc*Nn,       1,   Nc*Nn, Hc_00);	
			Matc Hc_01(Nc*Nn,Nc*Nn); Hc.get_submatrix(      1,   Nc*Nn, Nc*Nn+1, 2*Nc*Nn, Hc_01); // UPPER diagonal
			Matc Hc_10(Nc*Nn,Nc*Nn); Hc.get_submatrix(Nc*Nn+1, 2*Nc*Nn,       1,   Nc*Nn, Hc_10); // LOWER diagonal
			Matc Hc_11(Nc*Nn,Nc*Nn); Hc.get_submatrix(Nc*Nn+1, 2*Nc*Nn, Nc*Nn+1, 2*Nc*Nn, Hc_11);
			
			Matc Mc_00(Nc*Nn,Nc*Nn); Mc.get_submatrix(      1,   Nc*Nn,       1,   Nc*Nn, Mc_00);	
			Matc Mc_01(Nc*Nn,Nc*Nn); Mc.get_submatrix(      1,   Nc*Nn, Nc*Nn+1, 2*Nc*Nn, Mc_01);
			Matc Mc_10(Nc*Nn,Nc*Nn); Mc.get_submatrix(Nc*Nn+1, 2*Nc*Nn,       1,   Nc*Nn, Mc_10);
						
			// ----------------------------------------------
			// construct Htot = Hc_00 + Hc_01 + Hc_10
			//           Mtot = Mc_00 + Mc_01 + Mc_10
			// ----------------------------------------------
			Matc Htot(Nc*Nn, Nc*Nn);
			Htot += Hc_00;
			Htot += Hc_01;
			Htot += Hc_10;
			Matc Mtot(Nc*Nn, Nc*Nn);
			Mtot += Mc_00;
			Mtot += Mc_01;
			Mtot += Mc_10;
			
			// ----------------------------------------------
			// invert Mtot
			// ----------------------------------------------
			invert(Mtot);
			
			// ----------------------------------------------
			// construct tmp = Mtot^-1 * Htot
			// ----------------------------------------------
			Matc tmp(Nc*Nn, Nc*Nn);
			mult(Mtot, Htot, tmp);
			
			// ----------------------------------------------
			// solve eigenvalue problem tmp*psi=lambda*psi
			// ----------------------------------------------
			uint n = tmp.num_rows();
			cplx AA[n*n];
			for (uint ii=1; ii <= n; ii++) {
				for (uint jj=1; jj <= n; jj++) {
					AA[(jj-1)*n+(ii-1)] = tmp(ii,jj);
				}
			}
			cplx lambda[n];			// will store eigenvalues
			cplx vr[n*n];
			negf_math::zgeev(n, AA, lambda, vr); 
			
			// -----------------------------------------------------
			// add eigenvalues to the "list of forbidden energies"
			// -----------------------------------------------------
			for (uint ii=0; ii<n; ii++) {
				NEGF_FASSERT(fabs(lambda[ii].imag()) < constants::imag_err, "imaginary eigenvalue detected: (%e,%e)", lambda[ii].real(), lambda[ii].imag());
				bool divergence_already_present = false;
				for (int jj=divergences.size()-1; jj>=0; jj--) {
					if (fabs(divergences[jj]-lambda[ii].real()) < constants::convert_from_SI(units::energy, constants::SIec * 1e-9)) {
						divergence_already_present = true;
						break;
					}
				}
				if (!divergence_already_present) {
					this->divergences.push_back(lambda[ii].real());
				}
				logmsg->emit(LOG_INFO_L2,"found divergence E=%.8g (kk=%d,cc=%d)",lambda[ii].real(),kk,cc);
			}
		}
	}
	logmsg->emit(LOG_INFO_L3,"Divergences: ");
	std::sort(this->divergences.begin(), this->divergences.end()); 
	for (uint ii=0; ii < this->divergences.size(); ii++) {
		logmsg->emit_noendl(LOG_INFO_L3,"E=%.8g   ", this->divergences[ii]);
	}
	logmsg->emit(LOG_INFO_L3,"");
);}



/**  determine resonances from shape of LDOS@k=0 in the middle of the structure
     at the moment only a single resonance is sought
     called by master thread only */
void Resonances::determine_ldos_resonances()
{STACK_TRACE(
	uint NE = energies->get_number_of_points();
	// get LDOS at k=0 in the middle of the structure
	const Matd & LDOS_k0 = this->pp->get_entire_local_dos_k0();
	uint xmid = LDOS_k0.num_cols() / 2;	// integer division
	
	// find maximum of LDOS
	int Eidx_with_max_LDOS = 0;
	double max_LDOS = -1.0;
	for (uint ii=1; ii<=LDOS_k0.num_rows(); ii++) {
		if (LDOS_k0(ii, xmid) > max_LDOS) {
			max_LDOS = LDOS_k0(ii, xmid);
			Eidx_with_max_LDOS = ii;
		}
	}
	
	// debug info
	for (uint xx=2; xx<=this->xspace->get_num_internal_vertices()-2; xx++) {
		// find maximum at this place
		uint Eidx_with_max_LDOS_xx = 0;
		double max_LDOS_xx = -1.0;
		for (uint ii=1; ii<=NE; ii++) {
			if (LDOS_k0(ii,xx) > max_LDOS_xx) {
				max_LDOS_xx = LDOS_k0(ii,xx);
				Eidx_with_max_LDOS_xx = ii;
			}
		}
		logmsg->emit(LOG_INFO,"x=%d: max_LDOS=%e (E=%g)",xx,max_LDOS_xx, energies->get_energy_from_global_idx(Eidx_with_max_LDOS_xx));
	}

	// assign that energy point as resonance
	double reson = energies->get_energy_from_global_idx(Eidx_with_max_LDOS);
	double broad = constants::convert_from_SI(units::energy, 1e-3 * constants::SIec); // hard coded width at the moment
	this->resonances.push_back(reson);
	this->broadenings.push_back(broad); 
	logmsg->emit(LOG_INFO,"Found resonance at E=%.6g. Will apply broadening of %.3g.",reson, broad);
);}


/** same as determine_ldos_resonances(), but now for all k-vectors, not only k=0
 *  must be called by all MPI processes since LDOS at an arbitrary k is not saved anywhere else and
 *  must be computed */
void Resonances::determine_all_ldos_resonances()
{STACK_TRACE(
	logmsg->emit(LOG_INFO,"Computing resonances of LDOS for all k.");
	uint NE   = energies->get_number_of_points();
	uint myNE = energies->get_my_number_of_points();
	uint Nk   = kspace->get_number_of_points();
	uint Nx = xspace->get_num_internal_vertices();
	uint xmid = xspace->get_num_internal_vertices() / 2 + 1;	// integer division, +1->FLENS index
	double spin = get_spin_degeneracy(options->get("kp_method"), quantities::electron_density); // assumption: all bands have same degeneracy..
		                 
	// compute LDOS in the middle of the structure -> for own energies --> for all k-vectors
	logmsg->emit(LOG_INFO_L2,"Computing local part.");
	Matd midLDOSlocal(myNE,Nk);
	for (uint ee2=0; ee2<myNE; ee2++) {
		uint ee = energies->get_global_index(ee2);
		for (uint kk=0; kk<Nk; kk++) {
			const Matc & GR = gf->get_retarded(kk, ee);
			
			// determine maximum of LDOS(k,E) around (+- 10 vertices) xmid
			midLDOSlocal(ee2+1,kk+1) = 0.0;
			for (int xx=negf_math::max(5, int(xmid)-10); xx<=min(int(Nx)-5, (int)xmid+10); xx++) {
				double tmp = 0.0;
				for (uint nn = 1; nn <= Nn; nn++) {
					// LDOS = i/2pi * (GR-GA)(x,x)
					tmp += (-1.0/constants::pi) * spin * GR((xx-1)*Nn+nn,(xx-1)*Nn+nn).imag();
				}
				NEGF_ASSERT(tmp >= -1e-5, "negative LDOS?");
				midLDOSlocal(ee2+1,kk+1) = negf_math::max(midLDOSlocal(ee2+1,kk+1), tmp);
			}
			
		}
	}
	mpi->synchronize_processes();
	
	// communicate to master thread
	if (mpi->get_rank()==constants::mpi_master_rank) 
	{
		Matd midLDOS(NE,Nk);
		
		logmsg->emit(LOG_INFO_L2,"Computing own energy points (master)");
		// own energy points
		for (uint ee2=0; ee2<myNE; ee2++) {	
			uint ee = energies->get_global_index(ee2);
			for (uint kk=1; kk<=Nk; kk++) {
				midLDOS(ee+1,kk) = midLDOSlocal(ee2+1, kk);
			}
		}
		
		// energy points of other processors
		logmsg->emit(LOG_INFO_L2,"Computing other energy points (master)");
		for (int ppp=0; ppp<mpi->get_num_procs(); ppp++) {
			if (ppp==constants::mpi_master_rank) continue;
			
			// receive from other process
			uint num_energy_points = energies->get_number_of_points(ppp);
			Matd tmp_dos(num_energy_points, Nk);
			int tag = ppp;
			mpi->recv(tmp_dos, ppp, tag);
			
			// add to midLDOS
			uint start_idx = energies->get_start_global_idx(ppp) + 1;
			uint stop_idx = energies->get_stop_global_idx(ppp) + 1;
			for (uint ee=start_idx; ee<=stop_idx; ee++) {
				for (uint kk=1; kk<=Nk; kk++) {
					midLDOS(ee,kk) = tmp_dos(ee-start_idx+1, kk);
				}
			}
		}
		
		// compute!
		logmsg->emit(LOG_INFO_L2,"Computing total LDOS(E,k,x=xmid) (master)");
		double max_LDOS_k0 = -1.0;
		for (uint kk=1; kk<=Nk; kk++)
		{	
			// find maximum of LDOS for that k-point
			uint Eidx_with_max_LDOS = 0;
			double max_LDOS = -1.0;
			for (uint ee=1; ee<=midLDOS.num_rows(); ee++) {
				if (midLDOS(ee,kk) > max_LDOS) {
					max_LDOS = midLDOS(ee,kk);
					Eidx_with_max_LDOS = ee-1;			// ii is 1-based, Eidx_with_max_LDOS is 0-based
				}
			}
			if (kk==1) max_LDOS_k0 = max_LDOS;
			NEGF_ASSERT(max_LDOS_k0>1e-10, "max_LDOS_k0 was smaller than 1e-10.");
			
			if (fabs(max_LDOS/max_LDOS_k0) < 0.0001) continue; // can be the case for high k-vectors --> either completely zero or unbound --> not resonant

			double reson = energies->get_energy_from_global_idx(Eidx_with_max_LDOS);
			double broad = constants::convert_from_SI(units::energy, 5e-3 * constants::SIec); // hard coded width at the moment
			
			// determine maximal, minimal fermilevels
			// skip resonance if it is too far off and will not contribute substantially to current
			NEGF_ASSERT(this->fermilevels.size()==2, "expected 2 contacts/fermilevels");
			double EF_min = min(this->fermilevels[0], this->fermilevels[1]);
			double EF_max = max(this->fermilevels[0], this->fermilevels[1]);
			double kT = constants::convert_from_SI(units::energy, constants::SIkb * this->options->get("temperature"));
			if (reson > EF_max + reson_k_neglect_kT*kT || reson < EF_min - reson_k_neglect_kT*kT) {
				continue;
			}

			// skip resonance if it is already in the list (from another kk)
			bool skip_resonance = false;
			uint other_resonance = 0;
			for (uint ii=0; ii<resonances.size(); ii++) {
				if (fabs(resonances[ii]-reson) < constants::convert_from_SI(units::energy, 0.0001 * constants::SIec)) {
					other_resonance = ii;
					skip_resonance = true;
					break;
				}
			}
			if (skip_resonance) {
				logmsg->emit(LOG_INFO,"  kk=%d: Skipping resonance at E=%.6g (idx=%d) because it is close to E=%.6g.",
					kk-1, reson, Eidx_with_max_LDOS, resonances[other_resonance]);
				continue;
			}
			
			// assign that energy point as resonance
			this->resonances.push_back(reson);
			this->broadenings.push_back(broad);
			logmsg->emit(LOG_INFO,"  kk=%d: Found resonance at E=%.6g (ee=%d). Assume broadening of %.3g.", kk-1, reson, Eidx_with_max_LDOS, broad);
		}
	} else {
		// send to master process
		int dest = constants::mpi_master_rank;
		int tag = mpi->get_rank();
        mpi->send(midLDOSlocal, dest, tag);
	}
	mpi->synchronize_processes();
);}


/** 
class Interface
	void set_minimal_edges(const double& cb_edge, const double& vb_edge);
	void set_maximal_edges(const double& max_cb, const double& min_vb);        
	const double& get_maximal_cb_edge() ;              
	const double& get_maximal_vb_edge() ;
	bool maximal_edges_set() ;
	void set_potential(const vector<double>& potential);
	
	void calculate() ;
	bool is_calculated() ;

	// returns the number of valid cb bands in the sense that they were not affected by the barrier sorting. 
	// the barrier sorting does only affect kp models that treat electron and holes in the same hamiltonian.
	uint get_number_of_valid_cb_subbands() ;
	uint get_number_of_valid_vb_subbands() ;
	
	// returns the number of bound bands
	// if the band barriers are set, the number of bound subbands are given by the number of bands between minimal cb edge and the cb _barrier
	uint get_number_of_bound_cb_subbands() ;
	uint get_number_of_bound_vb_subbands() ;

	// returns the number of subbands that should be calculated;   bound_subbands <= valid_subbands <= max_num_subbands
	uint get_max_number_of_cb_subbands() ;
	uint get_max_number_of_vb_subbands() ;

    void write_bandstructure(const char* filename) ;
	const double& get_cb_base_energy(unsigned int subband) ;
	const double& get_vb_base_energy(unsigned int subband) ;
	virtual	const vector<double>& get_cb_base_probability(unsigned int subband) ;
	const vector<double>& get_vb_base_probability(unsigned int subband) ;

	void calculate_bandedges(vector<double>& node_cb_edge, vector<double>& node_vb_edge);	

class InterfaceSingleDispersion : public Interface 
	bool cb_is_parabolic() ;
	bool vb_is_parabolic() ;

	int get_cb_num_k_values() ;
	int get_vb_num_k_values() ;
	int get_max_num_k_values() ;
	double get_k_value(unsigned int kk) ;
		
	const double& get_cb_energy(unsigned int band, unsigned int kk) ;
	const double& get_vb_energy(unsigned int band, unsigned int kk) ;
		
	// returns the imaginary part of the energy for a given subband at given kk index (only nonzero for PML calculations) 
	double get_cb_imag_energy(unsigned int band, unsigned int kk) ;
	double get_vb_imag_energy(unsigned int band, unsigned int kk) ;    
	
	// evaluates the bandstructure at the given k values 
	void evaluate_cb_energy(unsigned int band, const vector<double>& k_values, vector<double>& energies) ;
	void evaluate_vb_energy(unsigned int band, const vector<double>& k_values, vector<double>& energies) ;        
	// finds cb inverse a given k value 
	void get_cb_k_of_energy(unsigned int band, const double& energy, vector<double>& k_values) ;
	void get_vb_k_of_energy(unsigned int band, const double& energy, vector<double>& k_values) ;
		
	const vector<double>& get_cb_probability(unsigned int band, unsigned int kk) ;
	const vector<double>& get_vb_probability(unsigned int band, unsigned int kk) ;
	void write_probability(bool cb_band, const char* filename, unsigned int subband, unsigned int kk) ;
	void write_bandstructure(const char* filename) ;
    
this->pml is InterfaceWellRadial which is derived from InterfaceSingleDispersion
*/

/** called by master thread only 
 *  look for "uint nk = ..." definition to check whether all k-vectors are taken or only k=0 */
#ifdef NOTDKP
void Resonances::determine_pml_resonances()
{STACK_TRACE(
	NEGF_EXCEPTION("NEGF was compiled with -DNOTDKP flag --> no PML resonance functionality.");
);}
#else
void Resonances::determine_pml_resonances()
{
	//try {
	// -----------------------------------------------------------
	// Determine energetic potential
	// -----------------------------------------------------------
    const double ec = constants::convert_from_SI(units::charge, constants::SIec);
	vector<double> potential = ham->get_electrostatic_potential(); 
	for (uint ii=0; ii<potential.size(); ii++) {
		potential[ii] = -ec*potential[ii];
	}
	
	// -----------------------------------------------------------
	// Interpolate energetic potential onto PML grid
	// -----------------------------------------------------------
	vector<double> pml_potential; pml_potential.resize(this->pmlgrid->get_num_vertices(), 0.0);
	for (uint ii=0; ii<pml_potential.size(); ii++) {
		double x = this->pmlgrid->get_vertex(ii)->get_coordinate(0);
		
		// left end
		if (x <= xspace->get_vertex(0)->get_coordinate(0)) {
			pml_potential[ii] = potential[0];
			continue;
		}
		// right end
		uint Nvert = xspace->get_num_vertices();
		if (x >= xspace->get_vertex(Nvert-1)->get_coordinate(0)) {
			pml_potential[ii] = potential[Nvert-1];
			continue;
		}
		
		// in between
		uint lower_idx = 0;
		while (lower_idx < Nvert-1 && xspace->get_vertex(lower_idx+1)->get_coordinate(0) < x) {
			lower_idx++;
		}
		NEGF_ASSERT(lower_idx<Nvert-1, "expected lower_idx<Nvert-1");		
		uint upper_idx = lower_idx+1;
		double xlow = xspace->get_vertex(lower_idx)->get_coordinate(0);
		double xupp = xspace->get_vertex(upper_idx)->get_coordinate(0);
		double frac = (x-xlow) / (xupp-xlow); // frac=0 --> x=xlow
		pml_potential[ii] = (1.0-frac)*potential[lower_idx] + frac*potential[upper_idx];
	}
	/*cout << "pml_potential = ";
    for (uint ii=0; ii<pml_potential.size(); ii++) {
    	cout << pml_potential[ii] << "  ";
    }
    cout << endl;*/
    
	// -----------------------------------------------------------
	// Determine band edges (including potential)
	// -----------------------------------------------------------
    vector<double> cbedge; cbedge.resize(xspace->get_num_vertices(), 0.0);
    vector<double> vbedge; vbedge.resize(xspace->get_num_vertices(), 0.0);
    const double T = options->get("temperature");
    const MaterialDatabase * db = ham->get_material_db();
    for (uint ii=0; ii<xspace->get_num_vertices(); ii++) {
    	Vertex * v = xspace->get_vertex(ii);
    	for (uint jj=0; jj<xspace->get_num_regions_near(v); jj++) {
    		const PropertyContainer<double> * mat = xspace->get_regions_near(v)[jj]->get_material();
    		double cb = this->infodesk->get_cbedge(mat, T, db);
    		double vb = this->infodesk->get_property(mat->get_name(), "valence_band_edge");
    		cbedge[ii] += cb;
    		vbedge[ii] += vb;
    	}
    	cbedge[ii] /= xspace->get_num_regions_near(v);
    	vbedge[ii] /= xspace->get_num_regions_near(v);
    	cbedge[ii] += potential[ii];
    	vbedge[ii] += potential[ii];
    }
    
    // -----------------------------------------------------------
    // Determine vertices where bound states are supposed to sit
    // -----------------------------------------------------------
    vector<Vertex *> bound_vertices;
    for (uint ii=0; ii<xspace->get_num_vertices(); ii++) {
    	Vertex * v = xspace->get_vertex(ii);
    	const vector<Region *> regs_near_v = xspace->get_regions_near(v);
    	for (uint jj=0; jj<regs_near_v.size(); jj++) {
	    	const string & regname = regs_near_v[jj]->get_name();
	    	if (regname.length()>6 && regname.substr(regname.length()-6)=="_bound") {
		    	bound_vertices.push_back(v);
		    	break;
	    	}
    	}
    }
    NEGF_ASSERT(bound_vertices.size()>0, "No vertices w/ adjacent regions ending \"_bound\" found!");
    
	// -----------------------------------------------------------
	// Determine bound energy intervals
	// -----------------------------------------------------------
    double Ecmin = 1e100;
    double Ecmax = -1e100;
    double Evmin = 1e100;
    double Evmax = -1e100;
    for (uint ii=0; ii<bound_vertices.size(); ii++) {
    	double Ec = cbedge[bound_vertices[ii]->get_index_global()];
    	double Ev = vbedge[bound_vertices[ii]->get_index_global()];
    	if (Ec<Ecmin) Ecmin = Ec;
    	if (Ec>Ecmax) Ecmax = Ec;
    	if (Ev<Evmin) Evmin = Ev;
    	if (Ev>Evmax) Evmax = Ev;
    }
	NEGF_ASSERT(Ecmin!=1e100 && Ecmax!=-1e100 && Evmin!=1e100 && Evmax!=-1e100, "could not determine bound energy intervals.");
    logmsg->emit(LOG_INFO,"Bound energy intervals: CB=[%e,%e], VB=[%e,%e]",Ecmin,Ecmax,Evmin,Evmax);
	//} catch (Exception * e) { NEGF_EXCEPTION(e->get_reason().c_str()); }
    
    logmsg->emit(LOG_INFO,"Please ignore stinky solution warnings if there are."); 
    logmsg->emit(LOG_INFO,"The bound energy interval was calculated near the bound regions only whereas ");
    logmsg->emit(LOG_INFO,"solutions located in the PML will appear where the electrostatic might be lower.");
    try {
    this->pml->set_potential(pml_potential);
	this->pml->set_minimal_edges(Ecmin, Evmax);
    this->pml->set_maximal_edges(Ecmax, Evmin);  
        
	// -----------------------------------------------------------
	// Calculate imaginary energies using PML boundaries
	// -----------------------------------------------------------
    this->pml->calculate();
    } catch (std::string e) { NEGF_EXCEPTION(e.c_str()); }
    
	// -----------------------------------------------------------
	// Set up vector of resonance energies of bound states
	// Do some screen output
	// -----------------------------------------------------------
	vector<double> & res   = this->resonances;
	vector<double> & broad = this->broadenings;
	const double state_broadening = constants::convert_from_SI(units::energy, 1e-3 * constants::SIec);
	// note vectors are NOT cleared, i.e. resonances are appended to existing resonances
			
    uint num_cbbands = this->pml->get_number_of_bound_cb_subbands();
    uint num_vbbands = this->pml->get_number_of_bound_vb_subbands();
    uint nk_valid = min(this->pml->get_cb_num_k_values(), this->pml->get_vb_num_k_values());
    uint nk = /*1*/nk_valid;
    logmsg->emit(LOG_INFO,"Searching %d k-values",nk);
    for (uint kk=0; kk<nk; kk++)
    {
	    logmsg->emit(LOG_INFO,"PML resonances @k=%d:",kk);
	    for (uint ii=0; ii<num_cbbands; ii++) 
		{
	    	// determine vertex w/ maximum probability
	    	Vertex * v = 0; 
	    	const vector<double> & prob = this->pml->get_cb_probability(ii,kk);
	    	double maxprob = -1e100;
	    	for (uint jj=0; jj<prob.size(); jj++) {
	    		if (prob[jj]>maxprob) {
	    			maxprob = prob[jj];
	    			v = this->pmlgrid->get_vertex(jj);
	    		}
	    	}
	    	NEGF_ASSERT(v!=NULL, "v=0");
	    	// determine if that vertex is near a bound region
	    	bool v_within_bound_region = false;
	    	const vector<Region *> regs_near_v = this->pmlgrid->get_regions_near(v);
	    	for (uint jj=0; jj<regs_near_v.size(); jj++) {
		    	const string & regname = regs_near_v[jj]->get_name();
		    	if (regname.length()>6 && regname.substr(regname.length()-6)=="_bound") {
			    	v_within_bound_region = true;
			    	break;
		    	}
	    	}
	    	
			// screen and file output
			// add to resonance vector if bound
			LoggerLevel level;
			string bub;
	    	char filename[1000];
	    	sprintf(filename,"pml_cb%d.dat",ii);
			if (v_within_bound_region) {
				res.push_back(this->pml->get_cb_energy(ii, kk));
				//broad.push_back(this->pml->get_cb_imag_energy(ii, kk));
				broad.push_back(state_broadening);
	
				level = LOG_INFO;
				bub = "BOUND";
	    		this->pml->write_probability(true, filename, ii, kk);
			} else {
				level = LOG_INFO_L2;
				bub = "unbound";
			}
	    	logmsg->emit(level,"CB-%d  (%s)  E=(%e,%e)  xmax=%e (%s)", ii, filename, this->pml->get_cb_energy(ii, kk), 
	    			this->pml->get_cb_imag_energy(ii, kk), v->get_coordinate(0), bub.c_str());
	    }
	    for (uint ii=0; ii<num_vbbands; ii++) 
		{
	    	// determine vertex w/ maximum probability
	    	Vertex * v = 0; 
	    	const vector<double> & prob = this->pml->get_vb_probability(ii,kk);
	    	double maxprob = -1e100;
	    	for (uint jj=0; jj<prob.size(); jj++) {
	    		if (prob[jj]>maxprob) {
	    			maxprob = prob[jj];
	    			v = this->pmlgrid->get_vertex(jj);
	    		}
	    	}
	    	NEGF_ASSERT(v!=NULL, "v=0");
	    	// determine if that vertex is near a bound region
	    	bool v_within_bound_region = false;
	    	const vector<Region *> regs_near_v = this->pmlgrid->get_regions_near(v);
	    	for (uint jj=0; jj<regs_near_v.size(); jj++) {
		    	const string & regname = regs_near_v[jj]->get_name();
		    	if (regname.length()>6 && regname.substr(regname.length()-6)=="_bound") {
			    	v_within_bound_region = true;
			    	break;
		    	}
	    	}
	    
			// screen and file output
			// add to resonance vector if bound
			LoggerLevel level;
			string bub;
	    	char filename[1000];
	    	sprintf(filename,"pml_vb%d.dat",ii);
			if (v_within_bound_region) {
				res.push_back(this->pml->get_vb_energy(ii, kk));
				//broad.push_back(this->pml->get_vb_imag_energy(ii, kk));
				broad.push_back(state_broadening);
	
				level = LOG_INFO;
				bub = "BOUND";
	    		this->pml->write_probability(false, filename, ii, kk);
			} else {
				level = LOG_INFO_L2;
				bub = "unbound";
			}
	    	logmsg->emit(level,"VB-%d    E=(%e,%e)  xmax=%e (%s)", ii,this->pml->get_vb_energy(ii, kk), 
	    			this->pml->get_vb_imag_energy(ii, kk), v->get_coordinate(0), bub.c_str());
	    }
    }
}//);}
#endif // def NOTDKP
