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
#include "Luminescence.h"
using namespace negf;

// see all.h whether DIAGONAL_ONLY is #defined

Luminescence::Luminescence(const Overlap * ov_,
						const Geometry * xspace_, 
						const Kspace * kspace_, 
						const Energies * energies_, 
						const Options * options_,
						const Hamiltonian * ham_,
						const GreenFunctions * gf_,
						const SEPhotonSpontaneous * spont_,
						photon_current_mode mode_) throw (Exception *):
	ov(ov_),
	xspace(xspace_),
	kspace(kspace_),
	energies(energies_),
	options(options_),
	ham(ham_),
	gf(gf_),
	spont(spont_),
	mode(mode_),
	area(constants::convert_from_SI(units::area, 1e-4)), // 1cm2
	output_power(0.0),
	output_power_VB(0.0),
	jphot(Matc(1,1)),
	jphot_VB(Matc(1,1))
{STACK_TRACE(
	NEGF_ASSERT(ov!=NULL && xspace!=NULL && kspace!=NULL && energies!=NULL && options!=NULL && ham!=NULL && spont!=NULL, "null pointer.");
	this->Nx = xspace->get_num_internal_vertices();
	this->NxNn = Nx*Nn;
	this->Nk = kspace->get_number_of_points();
	this->NE = energies->get_number_of_points();
	this->myNE = energies->get_my_number_of_points();
	this->Nhw = 100;
	
	// ----------------------------------------------------------
	// prepare Egap_min, to be used for hw_grid discretization
	// ----------------------------------------------------------
	this->Egap_min = 1e100;
	for (uint ii=0; ii<xspace->get_num_regions(); ii++) 
	{
		const PropertyContainer<double> * mat = xspace->get_region(ii)->get_material();		
		double Ec = constants::convert_from_SI(units::energy, constants::SIec * 
						TdkpInfoDesk::get_cbedge(mat, options->get("temperature"), ham->get_material_db()));
		double Ev = constants::convert_from_SI(units::energy, constants::SIec * mat->get("valence_band_edge"));
		double Egap = Ec - Ev;
		
		this->Egap_min = min(this->Egap_min, Egap);
	}
		
	// -----------------------------------------------
	// initially regular discretization of spectrum
	// -----------------------------------------------
	this->hw_min = spont->get_min_emission_energy();
	this->hw_max = spont->get_max_emission_energy();
	double dhw = (hw_max-hw_min)/(Nhw-1);
	this->hw_grid.resize(Nhw, 0.0);
	for (uint ii=0; ii<Nhw; ii++) {
		this->hw_grid[ii] = hw_min + ii*dhw;
	}
	
	// --------------------------------
	// allocate memory
	// --------------------------------
	vector<cplx>           tmp1; tmp1.resize(Nhw, 0.0);
	vector< vector<cplx> > tmp2; tmp2.resize(myNE, tmp1);
	this->QQ   .resize(Nx, tmp2); 
	this->QQ_VB.resize(Nx, tmp2); 
	this->spectrum    = Matd(Nhw, Nx); 
	this->spectrum_VB = Matd(Nhw, Nx); 
	
	this->power_snapshots = Matd(1,1);
);}


void Luminescence::determine_mpi_stuff()
{STACK_TRACE(
	logmsg->emit_small_header("Determining MPI communication stuff for photons");
	const uint nE = energies->get_number_of_points();
	
	// MPI: set up the array of matrices GL, GG which is needed to calculate the own self-energies
	
	this->E0_idx = energies->get_start_global_idx(mpi->get_rank());
	double E0    = energies->get_energy_from_global_idx(E0_idx);
	this->E1_idx = energies->get_stop_global_idx(mpi->get_rank());
	double E1    = energies->get_energy_from_global_idx(E1_idx);
	
	// -----------------------------------------------------------------------------
	// compute the following quantities:
	// 1. E0_minus_Emax_idx - the energy just below E0-Emax, at least zero
	// 2. E1_minus_Emin_idx - the energy just above E1-Emin, at most E0
	// 3. E0_plus_Emin_idx  - the energy just below E0+Emin, at least E1
	// 4. E1_plus_Emax_idx  - the energy just above E1+Emax, at most the highest energy
	// 5. nE_below = E1_minus_Emin_idx-E0_minus_Emax_idx+1
	// 6. nE_above = E1_plus_Emax_idx-E0_plus_Emin_idx+1
	// -----------------------------------------------------------------------------
	
	// determine lowest needed global energy index BELOW interval
	this->E0_minus_Emax_idx = E0_idx;
	while (E0_minus_Emax_idx > 0 
		   && energies->get_energy_from_global_idx(E0_minus_Emax_idx) > E0-this->hw_max) {
		E0_minus_Emax_idx--;
	}
	// determine highest needed global energy index BELOW interval
	this->E1_minus_Emin_idx = E0_idx;
	while (E1_minus_Emin_idx > 0 
		   && energies->get_energy_from_global_idx(E1_minus_Emin_idx) > E1-this->hw_min) {
		E1_minus_Emin_idx--;
	}
	if (E1_minus_Emin_idx<E0_idx) E1_minus_Emin_idx++;
	// determine lowest needed global energy index ABOVE interval
	this->E0_plus_Emin_idx = E1_idx;
	while (E0_plus_Emin_idx < nE-1 
		   && energies->get_energy_from_global_idx(E0_plus_Emin_idx) < E0+this->hw_min) {
		E0_plus_Emin_idx++;
	}
	if (E0_plus_Emin_idx > E1_idx) E0_plus_Emin_idx--;
	// determine highest needed global energy index ABOVE interval
	this->E1_plus_Emax_idx = E1_idx;
	while (E1_plus_Emax_idx < nE-1
		   && energies->get_energy_from_global_idx(E1_plus_Emax_idx) < E1+this->hw_max) {
		E1_plus_Emax_idx++;
	}
	logmsg->emit_all(LOG_INFO_L3,"p%d (energy indices %d...%d = %.3geV...%.3geV) needs indices %d(%.3geV)...%d(%.3geV) and %d(%.3geV)...%d(%.3geV)",
			mpi->get_rank(), E0_idx, E1_idx, E0, E1,
			E0_minus_Emax_idx,energies->get_energy_from_global_idx(E0_minus_Emax_idx),
			E1_minus_Emin_idx,energies->get_energy_from_global_idx(E1_minus_Emin_idx),
			E0_plus_Emin_idx, energies->get_energy_from_global_idx(E0_plus_Emin_idx),
			E1_plus_Emax_idx, energies->get_energy_from_global_idx(E1_plus_Emax_idx));
	
	uint nE_below = E1_minus_Emin_idx-E0_minus_Emax_idx+1;
	uint nE_above = E1_plus_Emax_idx-E0_plus_Emin_idx+1;
	logmsg->emit_all(LOG_INFO_L2,"p%d: nE_above=%d, nE_below=%d",mpi->get_rank(), nE_above, nE_below);
	
	// OBSOLETE: we want to get rid of the points E0_idx and Emax_idx in the auxiliary arrays because we don't need them	
	
	// -----------------
	// allocate space
	// -----------------
	vector<cplx> tmp1; tmp1.resize(nE_above, 0.0);
	vector< vector<cplx> > tmp2; tmp2.resize(myNE, tmp1);
	this->PPplus   .clear(); this->PPplus   .resize(Nx, tmp2);
	this->PPplus_VB.clear(); this->PPplus_VB.resize(Nx, tmp2);
	
	vector<cplx> tmp3; tmp3.resize(nE_below, 0.0);
	vector< vector<cplx> > tmp4; tmp4.resize(myNE, tmp3);
	this->PPminus   .clear(); this->PPminus   .resize(Nx, tmp4);
	this->PPminus_VB.clear(); this->PPminus_VB.resize(Nx, tmp4);
	
	
	mpi->synchronize_processes();
);}


void Luminescence::determine_hw_grid()
{STACK_TRACE(
	this->hw_grid.resize(this->Nhw, 0.0);
	
	// equal discretization
	if (false) {
		double dhw = (this->hw_max - this->hw_min) / (this->Nhw - 1);
		for (uint ww=0; ww<this->Nhw; ww++) {
			this->hw_grid[ww] = hw_min + ww*dhw;
		}
		return;
	}
	
	// general discretization using the same monotonic function trick as in class Resonances
	if (mpi->get_rank()==constants::mpi_master_rank)
	{
		for (uint ww=0; ww<Nhw; ww++) 
		{
			// initial guess
			double hw =  hw_min + ww* (hw_max-hw_min) / (Nhw-1);
			
			// things needed for Newton iteration
			const uint   max_iter = 1000;
			const double kT = constants::convert_from_SI(units::energy, constants::SIkb * options->get("temperature"));
			const double dhw_max = kT; 
			const double mmin = this->monotonic_function(hw_min);
			const double mmax = this->monotonic_function(hw_max);
			const double desired_result = mmin + ww * (mmax-mmin)/(Nhw-1);
			
			// perform Newton iteration
			double dhw = 1e10;
			uint iter = 0;
			while (fabs(dhw) > constants::convert_from_SI(units::energy, constants::SIec * 1e-6) && iter < max_iter) {
				double    F = this->monotonic_function(hw) - desired_result;
				double dFdhw = this->monotonic_function_derivative(hw);
				dhw = -F/dFdhw;
				// limit update
				if (fabs(dhw) > dhw_max) { dhw = negf_math::sign(dhw) * dhw_max; }
				
				hw += dhw;
				iter++;
			}
			NEGF_FASSERT(iter<max_iter-1, "Newton for new hw-grid did not converge for energy point %d: hw=%e, dhw=%e, desired_result=%e, function=%e, dFdhw=%e",
						ww, hw, dhw, desired_result, this->monotonic_function(hw), this->monotonic_function_derivative(hw));
			
			// assign value
			this->hw_grid[ww] = hw;
		}
		
		// send result to all other threads
		for (int proc=0; proc<mpi->get_num_procs(); proc++) {
			if (proc==constants::mpi_master_rank) continue;
			int tag  = 654;
			int dest = proc;
			mpi->send(hw_grid, dest, tag);
		}
	} else {
		// receive new grid from master thread
		int tag  = 654;
		int source = constants::mpi_master_rank;
		mpi->recv(hw_grid, Nhw, source, tag);
	}
);}


double Luminescence::monotonic_function(const double & hw)
{STACK_TRACE(	
	// ---------------------------------------------------------------
	// first contribution to monotonic function: constant slope 0...1!
	// ---------------------------------------------------------------
	double result = hw / (hw_max - hw_min);
	
	// -----------------------------------------------------------------------------
	// second contribution to monotonic function: constant slope Egap...Egap+Edelta
	// -----------------------------------------------------------------------------
	const double x2 = 0.5;
	const double Edelta = constants::convert_from_SI(units::energy, constants::SIec * 0.1);
	result += x2 * this->get_monotonic_function_contribution(hw, this->Egap_min, Egap_min+Edelta);
	
	// -----------------------------------------------------------------------------------------
	// third contribution to monotonic function: constant slope Euser-Ewidth/2...Euser+Ewidth/2
	// -----------------------------------------------------------------------------------------
	if (options->exists("LuminescenceRefinementEnergy") && options->exists("LuminescenceRefinementWidth")) {
		const double x3 = 1.0;
		const double Euser  = constants::convert_from_SI(units::energy, options->get("LuminescenceRefinementEnergy") * constants::SIec);
		const double Ewidth = constants::convert_from_SI(units::energy, options->get("LuminescenceRefinementWidth")  * constants::SIec);
		result += x3 * this->get_monotonic_function_contribution(hw, Euser-Ewidth/2.0, Euser+Ewidth/2.0);
	}
		
	return result;
);}

double Luminescence::get_monotonic_function_contribution(const double & hw, const double & Emin, const double & Emax)
{STACK_TRACE(
	if (hw > Emin && hw < Emax) {
		return (hw-Emin) / (Emax-Emin);
	} else if (hw >= Emax) {
		return 1.0;
	} else {
		return 0.0;
	}
);}


double Luminescence::monotonic_function_derivative(const double & hw)
{STACK_TRACE(	
	// ---------------------------------------------------------------
	// first contribution to monotonic function: constant slope 0...1!
	// ---------------------------------------------------------------
	double result = 1.0 / (hw_max - hw_min);
	
	// --------------------------------------------------------------------------------
	// second contribution to monotonic function: constant slope Egap...Egap+Edelta
	// --------------------------------------------------------------------------------
	const double x2 = 0.5;
	const double Edelta = constants::convert_from_SI(units::energy, constants::SIec * 0.1);
	result += x2 * this->get_monotonic_function_contribution_derivative(hw, this->Egap_min, Egap_min+Edelta);
	
	// -----------------------------------------------------------------------------------------
	// third contribution to monotonic function: constant slope Euser-Ewidth/2...Euser+Ewidth/2
	// -----------------------------------------------------------------------------------------
	if (options->exists("LuminescenceRefinementEnergy") && options->exists("LuminescenceRefinementWidth")) {
		const double x3 = 1.0;
		const double Euser  = constants::convert_from_SI(units::energy, options->get("LuminescenceRefinementEnergy") * constants::SIec);
		const double Ewidth = constants::convert_from_SI(units::energy, options->get("LuminescenceRefinementWidth")  * constants::SIec);
		result += x3 * this->get_monotonic_function_contribution_derivative(hw, Euser-Ewidth/2.0, Euser+Ewidth/2.0);
	}
	return result;
);}


double Luminescence::get_monotonic_function_contribution_derivative(const double & hw, const double & Emin, const double & Emax)
{STACK_TRACE(
	if (hw > Emin && hw < Emax) {
		return 1.0 / (Emax-Emin);
	} else {
		return 0.0;
	}
);}

void Luminescence::calculate() throw (Exception *)
{STACK_TRACE(
	switch (this->mode) {
	case galerpin:  logmsg->emit_header("Computing luminescence, MODE: GALERPIN"); break;
	case lake: 		logmsg->emit_header("Computing luminescence, MODE: LAKE"); break;
	}
	double t0 = MPI_Wtime();
	
	this->determine_mpi_stuff();
	
	// ---------------------------------------------------------------
	// determine new discretization hw_grid (hwmin, hmax are fixed)
	// ---------------------------------------------------------------
	this->determine_hw_grid();
	
	// --------------------------------------------------
	// initialize the matrices PPminus and PPplus to zero
	// --------------------------------------------------
	NEGF_ASSERT(PPplus.size()==Nx,"PPplus.size()==Nx");
	for (uint xx=0; xx<Nx; xx++) {
		NEGF_ASSERT(PPplus[xx].size()==myNE, "PPplus[xx].size()==myNE");
		for (uint ee2=0; ee2<myNE; ee2++) {
			NEGF_ASSERT(this->E1_plus_Emax_idx-this->E0_plus_Emin_idx+1 == PPplus[xx][ee2].size(), "PPplus[xx][ee2].size()==nE_above");
			for (uint ee=this->E0_plus_Emin_idx; ee<=this->E1_plus_Emax_idx; ee++) {
				PPplus   [xx][ee2][ee-this->E0_plus_Emin_idx] = 0.0;
				PPplus_VB[xx][ee2][ee-this->E0_plus_Emin_idx] = 0.0;
			}
		}
	}
	NEGF_ASSERT(PPminus.size()==Nx,"PPminus.size()==Nx");
	for (uint xx=0; xx<Nx; xx++) {
		NEGF_ASSERT(PPminus[xx].size()==myNE, "PPminus[xx].size()==myNE");
		for (uint ee2=0; ee2<myNE; ee2++) {
			NEGF_ASSERT(this->E1_minus_Emin_idx-this->E0_minus_Emax_idx+1 == PPminus[xx][ee2].size(), "PPminus[xx][ee2].size()==nE_below");
			for (uint ee=this->E0_minus_Emax_idx; ee<=this->E1_minus_Emin_idx; ee++) {
				PPminus   [xx][ee2][ee-this->E0_minus_Emax_idx] = 0.0;
				PPminus_VB[xx][ee2][ee-this->E0_minus_Emax_idx] = 0.0;
			}
		}
	}
	
	this->communicate_As();
	
	this->compute_QQ_Jphot_spectrum();
	
	double t1 = MPI_Wtime();
	logmsg->emit(LOG_INFO,"Total time needed for the computation of the luminescence spectrum: %5.2g[s]",t1-t0);
);}


void Luminescence::communicate_As()
{STACK_TRACE(
	            
	// ========================================================================
	// get AL from processes storing E0_plus_Emin_idx...E1_plus_Emax_idx
	// and compute the helper quantity "PP" (for an explanation-->report of S.Steiger)
	// PP = array of (Nx-1)*myNE*(E1_plus_Emax_idx-E0_plus_Emin_idx+1) complex numbers
	// ========================================================================
	int root = constants::mpi_master_rank;
	int my_rank = mpi->get_rank();
	
	// --------------------------------------------------------------------------------------------------------
	// set up an array which will store for every process all other processes it relies on (EXcluding itself)
	// same processes as in SEPhotonSpontaneous!
	// --------------------------------------------------------------------------------------------------------
	logmsg->emit(LOG_INFO,"Setting up processes_needed...");
	vector< vector<int> > processes_needed; 
	vector<int> pp_E0_plus_Emin_idx; 	
	vector<int> pp_E1_plus_Emax_idx; 
	vector<int> pp_E0_idx;			
	vector<int> pp_E1_idx;			
	vector<int> pp_E0_minus_Emax_idx; 
	vector<int> pp_E1_minus_Emin_idx;
	spont->determine_needed_processes(processes_needed, pp_E0_minus_Emax_idx, pp_E1_minus_Emin_idx, pp_E0_idx, pp_E1_idx, pp_E0_plus_Emin_idx, pp_E1_plus_Emax_idx);
			
	// ----------------------------------------------------------------
	// zipped data arrays are already available in SEPhotonSpontaneous!
	// ----------------------------------------------------------------
	const vector<unsigned long>   & AL_real_comp_size  = spont->get_AL_real_comp_size();
	const vector<unsigned long>   & AL_imag_comp_size  = spont->get_AL_imag_comp_size();
	const vector<unsigned char *> & AL_real_compressed = spont->get_AL_real_compressed();
	const vector<unsigned char *> & AL_imag_compressed = spont->get_AL_imag_compressed();
	const vector<unsigned long>   & AG_real_comp_size  = spont->get_AG_real_comp_size();
	const vector<unsigned long>   & AG_imag_comp_size  = spont->get_AG_imag_comp_size();
	const vector<unsigned char *> & AG_real_compressed = spont->get_AG_real_compressed();
	const vector<unsigned char *> & AG_imag_compressed = spont->get_AG_imag_compressed();
	NEGF_ASSERT(AL_real_compressed.size()==myNE, "AL_real_data was not yet set up");
	for (uint ee2=0; ee2<myNE; ee2++) {
		NEGF_ASSERT(AL_real_compressed[ee2]!=NULL, "AL_real_data was not yet set up!");
	}
	mpi->synchronize_processes();
	
	// --------------------------------------------------------------------------------------
	// determine how much data is sent in total
	// we do this by going through the same algorithm as when sending the data, just w/o sending
	// necessary for the amount of buffer that needs to be allocated for a safe operation
	// --------------------------------------------------------------------------------------
	logmsg->emit(LOG_INFO,"Creating total amount of data to be sent...");
	NEGF_FASSERT(this->E1_idx-this->E0_idx+1==myNE,"E1_idx=%d, E0_idx=%d, myNE=%d",this->E1_idx,this->E0_idx,myNE);
	long buffer_size_needed = 0;
	vector<bool> process_was_computed; 
	process_was_computed.resize(mpi->get_num_procs(), false); 
	while (true)
	{
		vector<bool> receiver; receiver.resize(mpi->get_num_procs(), false);
		vector<bool> sender;   sender.resize(mpi->get_num_procs(), false);
				
		// -------------------------------------
		// set up arrays w/ receiver and sender
		// -------------------------------------
		mpi->determine_senders_receivers(processes_needed, process_was_computed, sender, receiver);
		
		// -------------------------------------------
		// we're done if there are no more receivers
		// -------------------------------------------
		uint num_receivers = 0;
		uint num_senders = 0;
		for (int pp=0; pp < mpi->get_num_procs(); pp++) {
			if (receiver[pp]) num_receivers++;
			if (sender[pp]) num_senders++;
		}
		if (num_receivers == 0) {
			break;
		}
		
		bool i_am_receiver = receiver[my_rank];
		bool i_am_sender   = sender[my_rank];
		
		if (i_am_sender) {
			NEGF_FASSERT(!i_am_receiver, "p%d is both sender and receiver!",my_rank);
		}
		
		if (i_am_receiver) {
			// is not of our concern right now
		} else 
		{
			if (!i_am_sender) { // this is possible!
			} else {
						
			// PART 1: sender sends his AL's to processes storing energies BELOW his energy interval
			for (uint ee=this->E0_idx; ee<=this->E1_idx; ee++) 
			{
				// determine receiver BELOW current process - if there is none, continue to next energy
				for (int pp=0; pp<mpi->get_num_procs(); pp++) 
				{
					bool pp_receives = receiver[pp] && pp_E0_plus_Emin_idx[pp]<=int(ee) && pp_E1_plus_Emax_idx[pp]>=int(ee);
					if (!pp_receives) continue;
					
					// so the current process needs to send AL[ee-this->E0_idx] plus the size information !
					buffer_size_needed += AL_real_comp_size[ee-this->E0_idx] + AL_imag_comp_size[ee-this->E0_idx] + 2*sizeof(unsigned long) + 8*MPI_BSEND_OVERHEAD;
				}	
			}
			
			// PART 2: sender sends his AG's to processes storing energies ABOVE his energy interval
			for (uint ee=this->E0_idx; ee<=this->E1_idx; ee++) 
			{			
				// determine receiver ABOVE current process - if there is none, continue to next energy
				for (int pp=0; pp<mpi->get_num_procs(); pp++) 
				{
					bool pp_receives = receiver[pp] && pp_E0_minus_Emax_idx[pp]<=int(ee) && pp_E1_minus_Emin_idx[pp]>=int(ee);
					if (!pp_receives) continue;
					
					// so the current process needs to send AG[ee-this->E0_idx] plus the size information !
					buffer_size_needed += AG_real_comp_size[ee-this->E0_idx] + AG_imag_comp_size[ee-this->E0_idx] + 2*sizeof(unsigned long) + 8*MPI_BSEND_OVERHEAD;
				}	
			}
			} // if(i_am_sender)
		}
	
		// mark receivers as computed
		// needs to be performed in ALL threads (senders, receivers and those which are neither)
		for (int pp=0; pp < mpi->get_num_procs(); pp++) {
			if (receiver[pp]) {
				process_was_computed[pp] = true;
			}
		}
	} // while(true)
	logmsg->emit_all(LOG_INFO_L3,"p%d will need to send %d chars in total.", mpi->get_rank(), buffer_size_needed);
	mpi->synchronize_processes();
	
	// ------------------------------
	// allocate buffer!
	// ------------------------------
	logmsg->emit(LOG_INFO,"Allocating MPI buffer...");
	unsigned long buffersize_long = buffer_size_needed+1000;
	NEGF_ASSERT(buffersize_long < 2147483648UL, "Buffer size does not fit into an int!!@!");
	int buffersize = int(buffersize_long);
	char * buffer = new char[buffersize];
	int err = MPI_Buffer_attach(buffer, buffersize);
	NEGF_FASSERT(err==0, "MPI_Buffer_attach gave error %d",err);
	mpi->synchronize_processes();
	
	// ------------------------------
	// prepare Hamiltonian
	// ------------------------------
	logmsg->emit(LOG_INFO,"Preparing Hamiltonians...");
	const OVMat & M = ov->get_internal_overlap();
	NEGF_ASSERT(M.num_cols()==NxNn && M.num_rows()==M.num_cols(), "wrong overlap matrix.");
	Matc Hsmall(NxNn,NxNn);
	Matc HminusEM_tmp1(NxNn,NxNn);
	vector< Matc >         HminusEM_tmp2; HminusEM_tmp2.resize(Nk, HminusEM_tmp1);
	vector< vector<Matc> > HminusEM;      HminusEM.resize(myNE, HminusEM_tmp2);
	if (this->mode==lake) // otherwise not needed
	{
		for (uint kk=0; kk<Nk; kk++) 
		{
			ham->get_internal(kspace->get_point(kk), Hsmall);
			NEGF_ASSERT(Hsmall.num_rows()==Nx*Nn, "something went wrong.");
		
			NEGF_ASSERT(xspace->get_dimension()==1, "only 1D is implemented!");
			for (uint ee2 = 0; ee2 < myNE; ee2++) 
			{
				uint ee = energies->get_global_index(ee2);
				const double E = energies->get_energy_from_global_idx(ee);
						
				Matc & HmEM = HminusEM[ee2][kk];	
				mult(M, -E, HmEM);
				HmEM += Hsmall;
			}
		}	
	}
	
	// ------------------------------
	// do MPI communication!
	// ------------------------------
	logmsg->emit(LOG_INFO,"Communicate!");
	process_was_computed.clear();
	process_was_computed.resize(mpi->get_num_procs(), false);
	while (true)
	{
		vector<bool> receiver; receiver.resize(mpi->get_num_procs(), false);
		vector<bool> sender;   sender.resize(mpi->get_num_procs(), false);
		
		// -------------------------------------
		// set up arrays w/ receiver and sender
		// -------------------------------------
		mpi->determine_senders_receivers(processes_needed, process_was_computed, sender, receiver);
		
		// -------------------------------------------
		// we're done if there are no more receivers
		// -------------------------------------------
		uint num_receivers = 0;
		uint num_senders = 0;
		for (int pp=0; pp < mpi->get_num_procs(); pp++) {
			if (receiver[pp]) num_receivers++;
			if (sender[pp]) num_senders++;
		}
		
		// some screen output
		if (my_rank==root) {
			logmsg->emit_noendl_all(LOG_INFO, "This time we have %d receivers: ",num_receivers);
			for (int pp=0; pp < mpi->get_num_procs(); pp++) {
				if (receiver[pp]) logmsg->emit_noendl_all(LOG_INFO, "%d ", pp);
			}
			logmsg->emit_noendl_all(LOG_INFO, " and %d senders: ",num_senders);
			for (int pp=0; pp < mpi->get_num_procs(); pp++) {
				if (sender[pp]) logmsg->emit_noendl_all(LOG_INFO, "%d ", pp);
			}
			logmsg->emit_all(LOG_INFO, "");
		}
		
		if (num_receivers == 0) {
			break;
		}
		//mpi->synchronize_processes();
				
		// ---------------------------------------------------------------------------
		// receiver processes receive all their needed energies from sender processes
		// ---------------------------------------------------------------------------
		bool i_am_receiver = receiver[my_rank];
		bool i_am_sender   = sender[my_rank];
		
		if (i_am_sender) {
			NEGF_FASSERT(!i_am_receiver, "p%d is both sender and receiver!",my_rank);
		}
		
		if (i_am_receiver)
		{
			SEMat tmp = SPhotMat_create(NxNn);
			vector<SEMat> ALmat; ALmat.resize(Nk, tmp);
			vector<SEMat> AGmat; AGmat.resize(Nk, tmp);
			vector<double> ALnorm; ALnorm.resize(Nk, 0.0);
			vector<double> AGnorm; AGnorm.resize(Nk, 0.0);
			
			double mpi_time1 = 0.0;
			double comp_time1 = 0.0;
			double mpi_time2 = 0.0;
			double comp_time2 = 0.0;
			
			// PART 1: receiver receives missing AL's ABOVE his energy interval
			for (uint ee=this->E0_plus_Emin_idx; ee<=this->E1_plus_Emax_idx; ee++) 
			{				
				int sender_id = energies->get_process_computing(ee);
				if (sender_id==my_rank) {
					// this is possible, ee could be E0_idx!
					continue;
				}
				NEGF_FASSERT(sender[sender_id] && !receiver[sender_id],"p%d expected sender %d but something went wrong.", my_rank, sender_id);
				int tag = ee;
								
				double t1 = MPI_Wtime();
				logmsg->emit_all(LOG_INFO_L3,"p%d waits for ee=%d from p%d", my_rank, ee, sender_id);
				mpi->recv(ALnorm, Nk, sender_id, tag);
				mpi->recv_antihermitians(ALmat, sender_id, tag);// get stuff from an energy below --> outscattering --> AL
				logmsg->emit_all(LOG_INFO_L3,"p%d just got ee=%d from p%d", my_rank, ee, sender_id);
				double t2 = MPI_Wtime();
				
				this->compute_PP(ee, ALmat, ALnorm, HminusEM, false);
				
				double t3 = MPI_Wtime();
				mpi_time1 += t2-t1;
				comp_time1 += t3-t2;
			}
			logmsg->emit_all(LOG_ERROR,"receiver p%d needed %.3g[s] for MPI and %.3g[s] for the first part",mpi->get_rank(),mpi_time1,comp_time1);
			
			// PART 2: receiver receives missing AG's BELOW his energy interval
			for (uint ee=this->E0_minus_Emax_idx; ee<=this->E1_minus_Emin_idx; ee++) 
			{				
				int sender_id = energies->get_process_computing(ee);
				if (sender_id==my_rank) {
					// this is possible, ee could be E0_idx!
					continue;
				}
				NEGF_FASSERT(sender[sender_id] && !receiver[sender_id],"p%d expected sender %d but something went wrong.", my_rank, sender_id);
				int tag = ee;
								
				double t1 = MPI_Wtime();
				logmsg->emit_all(LOG_INFO_L3,"p%d waits for ee=%d from p%d", my_rank, ee, sender_id);
				mpi->recv(AGnorm, Nk, sender_id, tag);
				mpi->recv_antihermitians(AGmat, sender_id, tag);// get stuff from an energy below --> outscattering --> AG, GG
				logmsg->emit_all(LOG_INFO_L3,"p%d just got ee=%d from p%d", my_rank, ee, sender_id);
				double t2 = MPI_Wtime();
				
				this->compute_PP(ee, AGmat, AGnorm, HminusEM, true);
				
				double t3 = MPI_Wtime();
				mpi_time2 += t2-t1;
				comp_time2 += t3-t2;
			}
			logmsg->emit_all(LOG_ERROR,"receiver p%d needed %.3g[s] for MPI and %.3g[s] for the second part",mpi->get_rank(),mpi_time2,comp_time2);
			
		} else 
		{
			if (!i_am_sender) { // this is possible!
			} else {
				
			// PART 1: sender sends his AL's to processes storing energies BELOW his energy interval
			for (uint ee=this->E0_idx; ee<=this->E1_idx; ee++) 
			{			
				// determine receiver(s) BELOW current process - if there is none, continue to next energy
				for (int pp=0; pp<mpi->get_num_procs(); pp++) 
				{
					bool pp_receives = receiver[pp] && pp_E0_plus_Emin_idx[pp]<=int(ee) && pp_E1_plus_Emax_idx[pp]>=int(ee);
					if (!pp_receives) continue;
					int receiver_id = pp;
						
					int tag = ee;
					int tag2 = ee+1;
				
					logmsg->emit_all(LOG_INFO_L3,"p%d AL sends ee=%d downwards to p%d", my_rank, ee, receiver_id);
					// BUFFERED SEND! mpi->send_antihermitians and mpi->recv_antihermitians send/receive the following 
					// quantities in order: real_char_size, imag_char_size, real_compressed, imag_compressed
					
					mpi->send(spont->get_AL_norm(ee-this->E0_idx) , receiver_id, tag);
					
					// send array lengths
					NEGF_ASSERT(AL_real_comp_size.size()>ee-this->E0_idx, "AL_real_comp_size has wrong size.");
					int real_char_size = int(AL_real_comp_size[ee-this->E0_idx]);
					int imag_char_size = int(AL_imag_comp_size[ee-this->E0_idx]);
					err = MPI_Bsend(&real_char_size, 1, MPI_INT, receiver_id, tag, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of AL_real_char_size failed: err=%d",err);
					err = MPI_Bsend(&imag_char_size, 1, MPI_INT, receiver_id, tag2, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of AL_imag_char_size failed: err=%d",err);
		
					// send arrays
					err = MPI_Bsend(AL_real_compressed[ee-this->E0_idx], real_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of AL_real_compressed failed: err=%d",err);
					err = MPI_Bsend(AL_imag_compressed[ee-this->E0_idx], imag_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag2, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of AL_imag_compressed failed: err=%d",err);
					
					logmsg->emit_all(LOG_INFO_L3,"p%d AL just sent ee=%d to p%d", my_rank, ee, receiver_id);
				}	
			}
			
			// PART 2: sender sends his AG's to processes storing energies ABOVE his energy interval
			for (uint ee=this->E0_idx; ee<=this->E1_idx; ee++) 
			{			
				// determine receiver ABOVE current process - if there is none, continue to next energy
				for (int pp=0; pp<mpi->get_num_procs(); pp++) 
				{
					bool pp_receives = receiver[pp] && pp_E0_minus_Emax_idx[pp]<=int(ee) && pp_E1_minus_Emin_idx[pp]>=int(ee);
					if (!pp_receives) continue;
					int receiver_id = pp;
					
					int tag = ee;
					int tag2 = ee+1;
				
					logmsg->emit_all(LOG_INFO_L3,"p%d sends ee=%d upwards to p%d", my_rank, ee, receiver_id);
					// BUFFERED SEND!
					// mpi->send_antihermitians and mpi->recv_antihermitians send/receive the following 
					// quantities in order: real_char_size, imag_char_size, real_compressed, imag_compressed
					
					mpi->send(spont->get_AG_norm(ee-this->E0_idx) , receiver_id, tag);
					
					// send array lengths
					NEGF_ASSERT(AG_real_comp_size.size()>ee-this->E0_idx, "AG_real_comp_size has wrong size.");
					int real_char_size = int(AG_real_comp_size[ee-this->E0_idx]);
					int imag_char_size = int(AG_imag_comp_size[ee-this->E0_idx]);
					err = MPI_Bsend(&real_char_size, 1, MPI_INT, receiver_id, tag, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of (2) real_char_size failed: err=%d",err);
					err = MPI_Bsend(&imag_char_size, 1, MPI_INT, receiver_id, tag2, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of (2) imag_char_size failed: err=%d",err);
		
					// send arrays
					err = MPI_Bsend(AG_real_compressed[ee-this->E0_idx], real_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of (2) real_compressed failed: err=%d",err);
					err = MPI_Bsend(AG_imag_compressed[ee-this->E0_idx], imag_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag2, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of (2) imag_compressed failed: err=%d",err);
					
					logmsg->emit_all(LOG_INFO_L3,"p%d just sent ee=%d upwards to p%d", my_rank, ee, receiver_id);
				}	
			}
			} // if(i_am_sender)
		}
	
		// mark receivers as computed
		// needs to be performed in ALL threads (senders, receivers and those which are neither)
		for (int pp=0; pp < mpi->get_num_procs(); pp++) {
			if (receiver[pp]) {
				process_was_computed[pp] = true;
			}
		}
		//mpi->synchronize_processes();
	} // while(true)
			
	// security check
	for (int pp=0; pp < mpi->get_num_procs(); pp++) {
		NEGF_ASSERT(process_was_computed[pp], "a process was not computed.");
	}
	
	// ----------------------------------------------------------------
	// release buffer
	// ----------------------------------------------------------------
	logmsg->emit(LOG_INFO,"Deallocating MPI buffer");
	err = MPI_Buffer_detach(&buffer, &buffersize);
	delete [] buffer;
	NEGF_FASSERT(err==0, "MPI_Buffer_attach gave error %d",err);
	
	// ----------------------------------------------------------------
	// no release of zipped data!
	// ----------------------------------------------------------------
	logmsg->emit_noendl_all(LOG_INFO,"p%d   ",mpi->get_rank());	
	mpi->synchronize_processes();
);}


// below = true --> A=AG. below=false --> A=AL.
void Luminescence::compute_PP(const uint ee, const vector<SEMat> & ALGmat, const vector<double> & ALGnorm, const vector< vector<Matc> > & HminusEM, bool below)
{STACK_TRACE(
	NEGF_ASSERT(ALGmat.size()==Nk && ALGnorm.size()==Nk && HminusEM.size()==energies->get_my_number_of_points(), "invalid vector size encountered.");
		
	// set up CB/VB distinction
	vector<uint> cb_bands; options->get_conduction_degrees_of_freedom(cb_bands);
	vector<uint> vb_bands; options->get_valence_degrees_of_freedom(vb_bands);
	for (uint nn=0; nn<Nn; nn++) {
		bool cb = false;
		bool vb = false;
		for (uint ii=0; ii<cb_bands.size(); ii++) {	if (cb_bands[ii]==nn) {	cb = true; break; }	}
		for (uint ii=0; ii<vb_bands.size(); ii++) {	if (vb_bands[ii]==nn) {	vb = true; break; }	}
		NEGF_ASSERT((cb && !vb) || (!cb && vb), "band must be CB or VB!");
	}
	
	SEMat tmp = SPhotMat_create(NxNn);
	if (this->mode==galerpin)
	{
		// ----------------------------------------------------------------------------
		// For all Ei (own energies) and k:
		// 1. Compute C(Ei,ee,k) = AG(ee,k)*GL(Ei,k)
		// 2. Take trace, multiply with k-space weight, and store in PPminus[x][Ei][ee]
		// PP has units length^{-2}*time^{-1}*energy^{-1} (NO divide by dx)
		// ----------------------------------------------------------------------------
		Matc CCtmp(NxNn,NxNn);
		Matc CCtmp2(NxNn,NxNn);
		Matc Aee_full(NxNn,NxNn);
		for (uint ee2=0; ee2<myNE; ee2++) 
		{
			for (uint kk=0; kk<Nk; kk++) 
			{
				double wk = kspace->get_point(kk).get_weight() / (constants::pi*constants::pi*4.0);
				uint Ei_idx = energies->get_global_index(ee2);
				const GLMat & G = (below) ? gf->get_lesser(kk,Ei_idx) : gf->get_greater(kk,Ei_idx);
				const SEMat & Aee = ALGmat[kk]; // ALG = AG (below) or AL (above)
				
				// skip if matrix norm is too small
				if (ALGnorm[kk] < constants::ALGnorm_neglect) continue;
				
				vector<cplx> AeeG; get_diag(Aee, G, AeeG);
				vector<cplx> GAee; get_diag(G, Aee, GAee);
				
				for (uint xx=0; xx<Nx; xx++) 
				{
					
					// take trace!
					cplx trace_CB = 0.0;
					for (uint nn=0; nn<cb_bands.size(); nn++) {
						uint idx = get_mat_idx(xx+1,cb_bands[nn]+1,Nx) - 1;
						trace_CB += 0.5 * (GAee[idx] + AeeG[idx]);
					}
					cplx trace_VB = 0.0;
					for (uint nn=0; nn<vb_bands.size(); nn++) {
						uint idx = get_mat_idx(xx+1,vb_bands[nn]+1,Nx) - 1;
						trace_VB += 0.5 * (GAee[idx] + AeeG[idx]);
					}
					
					// multiply with k-weight and store
					if (below) {
						NEGF_ASSERT(PPminus.size()==Nx, "PPminus.size()==Nx failed.");
						NEGF_ASSERT(PPminus[xx].size()==myNE, "PPminus[xx].size()==myNE failed.");
						NEGF_ASSERT(ee-this->E0_minus_Emax_idx < PPminus[xx][ee2].size(), "ee-this->E0_minus_Emax_idx<PPminus[xx][ee2].size() failed.");
						this->PPminus   [xx][ee2][ee-this->E0_minus_Emax_idx] += wk * trace_CB; 
						this->PPminus_VB[xx][ee2][ee-this->E0_minus_Emax_idx] += wk * trace_VB;
					} else {
						NEGF_ASSERT(PPplus.size()==Nx, "PPplus.size()==Nx failed.");
						NEGF_ASSERT(PPplus[xx].size()==myNE, "PPplus[xx].size()==myNE failed.");
						NEGF_ASSERT(ee-this->E0_plus_Emin_idx < PPplus[xx][ee2].size(), "ee-this->E0_plus_Emin_idx<PPplus[xx][ee2].size() failed.");
						this->PPplus   [xx][ee2][ee-this->E0_plus_Emin_idx] += wk * trace_CB; 
						this->PPplus_VB[xx][ee2][ee-this->E0_plus_Emin_idx] += wk * trace_VB; 
					}
				}
			}
		}
	} else 
	{ 	
		// ------------------
		// mode=lake
		// ------------------
		
		// ----------------------------------------------------------------------------
		// For all Ei (own energies) and k:
		// 1. Compute C(Ei,ee,k) = GR(Ei,k)*ALG(ee,k)*GA(Ei,k)
		// 2. Compute t(x,x+1)*C(x+1,x) - t(x+1,x)*C(x,x+1) for all x
		// 3. Take trace, multiply with k-space weight, and store
		// ----------------------------------------------------------------------------
		double time1 = 0.0;
		double time2 = 0.0;
		double t1, t2, t3;
		
#ifdef USE_BANDED
		// determine sparsity pattern which is the same for all A-matrices
		vector< vector<int> > pattern(NxNn);
		for (uint k = 0; k < NxNn; k++) {
			for (uint l=1; l<=NxNn; l++) {
				if(ALGmat[0].stores(k+1,l)) {
					pattern[k].push_back(l);	// pattern is 1-based (but the index to access is of course 0-based)					
				}
			}					
		}
#endif
				
		Matc CCtmp(NxNn,NxNn);
		for (uint ee2=0; ee2<myNE; ee2++) 
		{
			for (uint kk=0; kk<Nk; kk++) 
			{
				uint Ei_idx = energies->get_global_index(ee2);
				double wk = kspace->get_point(kk).get_weight() / (constants::pi*constants::pi*4.0);
				const Matc & GR = gf->get_retarded(kk,Ei_idx);
				const Matc & GA = gf->get_advanced(kk,Ei_idx);
				const SEMat & Aee = ALGmat[kk];
				
				// skip if matrix norm is too small
				if (ALGnorm[kk] < 1e-20) continue;
				
				t1 = MPI_Wtime();				
#ifdef USE_BANDED
				// afterwards we need CCtmp(xx,mm, xx+1,nn) and CCtmp(xx+1,mm, xx,nn) only! (1<=xx<=Nx-1)
				uint l, k;
				uint i1, j1, i2, j2;
				uint Nnp1 = Nn + 1;
				cplx tmp1, tmp2;
				for (uint xx=1; xx<=Nx-1; xx++) {
					for (uint mm=1; mm<Nnp1; mm++) {
						i1 = get_mat_idx(xx  ,mm,Nx);
						i2 = get_mat_idx(xx+1,mm,Nx);
						for (uint nn=1; nn<Nnp1; nn++) {
							j1 = get_mat_idx(xx+1,nn,Nx);
							j2 = get_mat_idx(xx  ,nn,Nx);
							
							CCtmp(i1,j1) = 0.0;
							CCtmp(i2,j2) = 0.0;
							for (uint kit=0; kit < NxNn; kit++) {
								k = kit + 1;
								
								tmp1 = 0.0; // will store (A*GA)_{kj1} = \sum_l A_{kl}GA_{lj1}
								tmp2 = 0.0; // will store (A*GA)_{kj2} = \sum_l A_{kl}GA_{lj2}
								for (uint lit = 0; lit < pattern[kit].size(); lit++) {
									l = pattern[kit][lit];
										
									tmp1 += Aee.nocheck(k,l)*GA(l,j1);
									tmp2 += Aee.nocheck(k,l)*GA(l,j2);
								}
								CCtmp(i1,j1) += GR(i1,k)*tmp1;
								CCtmp(i2,j2) += GR(i2,k)*tmp2;
							}
						}
					}
				}
#else
				mult(Aee, GA, tmp);   // tmp = Aee * GA; A = AG (below) or AL(above)
				mult(GR, tmp, CCtmp); // CCtmp = GR * tmp;
#endif
				
				t2 = MPI_Wtime();
				const Matc & HmEM = HminusEM[ee2][kk];

				vector<cplx> PPP; PPP.resize(Nn, 0.0);
				for (uint xx=0; xx<Nx-1; xx++) 
				{								
					// compute!
					for (uint nn=1; nn<=Nn; nn++) {
						uint x1n = get_mat_idx(xx+1,nn,Nx);
						uint x2n = get_mat_idx(xx+2,nn,Nx);
						for (uint ii=1; ii<=Nn; ii++) {
							uint x1i = get_mat_idx(xx+1,ii,Nx);
							uint x2i = get_mat_idx(xx+2,ii,Nx);
							PPP[nn-1] += HmEM(x1n,x2i)*CCtmp(x2i,x1n) - HmEM(x2n,x1i)*CCtmp(x1i,x2n);
						}
					}
						
					// take trace!
					cplx trace_CB = 0.0;
					for (uint nn=0; nn<cb_bands.size(); nn++) {
						trace_CB += PPP[cb_bands[nn]];
					}
					cplx trace_VB = 0.0;
					for (uint nn=0; nn<vb_bands.size(); nn++) {
						trace_VB += PPP[vb_bands[nn]];
					}
					
					// multiply with k-weight and store
					if (below) {
						NEGF_ASSERT(PPminus.size()==Nx-1, "xx<PPminus.size() failed.");
						NEGF_ASSERT(PPminus[xx].size()==myNE, "ee2<PPminus[xx].size() failed.");
						NEGF_ASSERT(ee-this->E0_minus_Emax_idx < PPminus[xx][ee2].size(), "ee-this->E0_minus_Emax_idx<PPminus[xx][ee2].size() failed.");
						this->PPminus   [xx][ee2][ee-this->E0_minus_Emax_idx] += wk * trace_CB;  // PPminus was initialized to 0 at the beginning of this routine
						this->PPminus_VB[xx][ee2][ee-this->E0_minus_Emax_idx] += wk * trace_VB;
					} else {
						NEGF_ASSERT(PPplus.size()==Nx-1, "xx<PPplus.size() failed.");
						NEGF_ASSERT(PPplus[xx].size()==myNE, "ee2<PPplus[xx].size() failed.");
						NEGF_ASSERT(ee-this->E0_plus_Emin_idx < PPplus[xx][ee2].size(), "ee-this->E0_plus_Emin_idx<PPplus[xx][ee2].size() failed.");
						this->PPplus   [xx][ee2][ee-this->E0_plus_Emin_idx] += wk * trace_CB;  // PPplus was initialized to 0 at the beginning of this routine
						this->PPplus_VB[xx][ee2][ee-this->E0_plus_Emin_idx] += wk * trace_VB;
					}
				}
				t3 = MPI_Wtime();
				time1 += t2-t1; 
				time2 += t3-t2; 
			}
		}
		logmsg->emit_all(LOG_INFO_L2,"ee=%d (p%d): time1=%.3g=%2.0g%%, time2=%.3g=%2.0g%%",ee,mpi->get_rank(),
				time1, time1/(time1+time2+1e-8)*100, time2, time2/(time1+time2+1e-8)*100);
	}
);}


void Luminescence::compute_QQ_Jphot_spectrum()
{STACK_TRACE(
	// =====================================================================
	// COMPUTE QQ(xx,E_i,hw_j) from PPplus and PPminus by linear combination
	// =====================================================================
	logmsg->emit(LOG_INFO,"Computing QQ");
	for (uint ee2=0; ee2<this->myNE; ee2++) 
	{
		uint ee = energies->get_global_index(ee2);
		double E = energies->get_energy_from_global_idx(ee);
		for (uint ww=0; ww<this->Nhw; ww++) 
		{
			NEGF_ASSERT(ww<hw_grid.size(), "ww<hw_grid.size() failed!");
			double hw = this->hw_grid[ww];
			double E_plus_hw  = E+hw;
			double E_minus_hw = E-hw;
			
			for (uint xx=0; xx<Nx; xx++) 
			{
				NEGF_ASSERT(QQ.size()==Nx && QQ[xx].size()==myNE && QQ[xx][ee2].size()==Nhw, "QQ had wrong size.");
				this->QQ   [xx][ee2][ww] = 0.0;
				this->QQ_VB[xx][ee2][ww] = 0.0;
								
				// --------------------------------------------------------
				// find energy grid points adjacent to E+hw				
				// --------------------------------------------------------
				bool E_plus_hw_inside_energy_grid = (E_plus_hw <= energies->get_energy_from_global_idx(NE-1));
				uint ee_lower = NE-1;
				while (ee_lower>0 && energies->get_energy_from_global_idx(ee_lower)>E_plus_hw) {
					ee_lower--;
				}
				if (ee_lower==NE-1) E_plus_hw_inside_energy_grid = false;
				if (E_plus_hw_inside_energy_grid) {
					uint ee_upper = ee_lower + 1;
					double Elow = energies->get_energy_from_global_idx(ee_lower);
					double Eupp = energies->get_energy_from_global_idx(ee_upper);
					double frac = (E_plus_hw - Elow) / (Eupp-Elow); // frac=1 --> E_plus_hw=Eupp
					NEGF_FASSERT(frac>=0.0 && frac<=1.0, "frac should be between 0 and 1, encountered %e instead (PPplus).", frac);
					
					// linear combination - frac=0 means E_minus_hw=Elow!
					NEGF_ASSERT(PPplus.size()==Nx && PPplus[xx].size()==myNE, "PPplus had wrong size.");
					NEGF_ASSERT(ee_lower>=this->E0_plus_Emin_idx, "ee_lower<this->E0_minus_Emax_idx");
					NEGF_ASSERT(ee_upper<=this->E1_plus_Emax_idx, "ee_upper>this->E1_minus_Emin_idx");
					NEGF_ASSERT(ee_lower-E0_plus_Emin_idx < PPplus[xx][ee2].size(), "ee_lower-E0_minus_Emax_idx<PPplus[xx][ee2].size() failed.");
					NEGF_ASSERT(ee_upper-E0_plus_Emin_idx < PPplus[xx][ee2].size(), "ee_upper-E0_minus_Emax_idx<PPplus[xx][ee2].size() failed.");
					if (this->mode==galerpin || this->mode==lake) 
					{
						this->QQ   [xx][ee2][ww] += (1.0-frac) * this->PPplus   [xx][ee2][ee_lower-E0_plus_Emin_idx] + frac * this->PPplus   [xx][ee2][ee_upper-E0_plus_Emin_idx];				
						this->QQ_VB[xx][ee2][ww] += (1.0-frac) * this->PPplus_VB[xx][ee2][ee_lower-E0_plus_Emin_idx] + frac * this->PPplus_VB[xx][ee2][ee_upper-E0_plus_Emin_idx];				
					}
				}
				
				// --------------------------------------------------------
				// find energy grid points adjacent to E-hw				
				// --------------------------------------------------------
				bool E_minus_hw_inside_energy_grid = (E_minus_hw >= energies->get_energy_from_global_idx(0));
				uint ee_upper = 0;
				while (ee_upper<NE-1 && energies->get_energy_from_global_idx(ee_upper)<E_minus_hw) {
					ee_upper++;
				}
				if (ee_upper==0) E_minus_hw_inside_energy_grid = false;
				if (E_minus_hw_inside_energy_grid) {
					ee_lower = ee_upper - 1;
					double Elow = energies->get_energy_from_global_idx(ee_lower);
					double Eupp = energies->get_energy_from_global_idx(ee_upper);
					double frac = (E_minus_hw - Elow) / (Eupp-Elow); // frac=1 --> E_minus_hw=Eupp
					NEGF_FASSERT(frac>=0.0 && frac<=1.0, "frac should be between 0 and 1, encountered %e instead (PPminus).",frac);
					
					// linear combination - frac=0 means E_minus_hw=Elow!
					NEGF_ASSERT(PPminus.size()==Nx && PPminus[xx].size()==myNE, "PPminus had wrong size.");
					NEGF_ASSERT(ee_lower>=this->E0_minus_Emax_idx, "ee_lower<this->E0_minus_Emax_idx");
					NEGF_ASSERT(ee_upper<=this->E1_minus_Emin_idx, "ee_upper>this->E1_minus_Emin_idx");
					NEGF_ASSERT(ee_lower-E0_minus_Emax_idx < PPminus[xx][ee2].size(), "ee_lower-E0_minus_Emax_idx<PPminus[xx][ee2].size() failed.");
					NEGF_ASSERT(ee_upper-E0_minus_Emax_idx < PPminus[xx][ee2].size(), "ee_upper-E0_minus_Emax_idx<PPminus[xx][ee2].size() failed.");
					if (this->mode==galerpin) 
					{
						this->QQ   [xx][ee2][ww] -= (1.0-frac) * this->PPminus   [xx][ee2][ee_lower-E0_minus_Emax_idx] + frac * this->PPminus   [xx][ee2][ee_upper-E0_minus_Emax_idx];				
						this->QQ_VB[xx][ee2][ww] -= (1.0-frac) * this->PPminus_VB[xx][ee2][ee_lower-E0_minus_Emax_idx] + frac * this->PPminus_VB[xx][ee2][ee_upper-E0_minus_Emax_idx];				
					} else 
					{ 	// lake: =, NOT -=!
						this->QQ   [xx][ee2][ww] = (1.0-frac) * this->PPminus   [xx][ee2][ee_lower-E0_minus_Emax_idx] + frac * this->PPminus   [xx][ee2][ee_upper-E0_minus_Emax_idx];				
						this->QQ_VB[xx][ee2][ww] = (1.0-frac) * this->PPminus_VB[xx][ee2][ee_lower-E0_minus_Emax_idx] + frac * this->PPminus_VB[xx][ee2][ee_upper-E0_minus_Emax_idx];				
					}
					// note the minus sign!!!!!!
				}
			}
		}
	}
	
	mpi->synchronize_processes();
	
	// =====================================================================
	// Send QQ to master thread
	// =====================================================================
	logmsg->emit(LOG_INFO,"Sending QQ to master thread");
	const double init_num = -8888.0;
	// send QQ-array to master thread
	if (mpi->get_rank()==constants::mpi_master_rank) 
	{
		vector<cplx> tmp1; tmp1.resize(Nhw, init_num);
		vector< vector<cplx> > tmp2; tmp2.resize(NE, tmp1);
		this->QQ_total   .assign(Nx, tmp2); // array of size Nx*NE*Nhw
		this->QQ_total_VB.assign(Nx, tmp2); // array of size Nx*NE*Nhw
		
		// own contribution
		for (uint ee2=0; ee2<myNE; ee2++) {
			uint ee = energies->get_global_index(ee2);
			for (uint xx=0; xx<Nx; xx++) {
				for (uint ww=0; ww<Nhw; ww++) {
					NEGF_ASSERT(QQ.size()==Nx && QQ[xx].size()==myNE && QQ[xx][ee2].size()==Nhw, "wrong size of QQ.");
					NEGF_ASSERT(QQ_total.size()==Nx && QQ_total[xx].size()==NE && QQ_total[xx][ee].size()==Nhw, "wrong size of QQ_total.");
					this->QQ_total   [xx][ee][ww] = this->QQ   [xx][ee2][ww];
					this->QQ_total_VB[xx][ee][ww] = this->QQ_VB[xx][ee2][ww];
				}
			}
		}
		
		// collect the pieces
		for (int pp=0; pp<mpi->get_num_procs(); pp++) 
		{
			if (pp==constants::mpi_master_rank) continue;
			uint ppNE = energies->get_number_of_points(pp);
			uint start_idx = energies->get_start_global_idx(pp);
			
			// receive from other process
			uint length = Nx*ppNE*Nhw;
			vector<double> tmp_real   ; tmp_real   .resize(length, init_num);
			vector<double> tmp_imag   ; tmp_imag   .resize(length, init_num);
			vector<double> tmp_real_VB; tmp_real_VB.resize(length, init_num);
			vector<double> tmp_imag_VB; tmp_imag_VB.resize(length, init_num);
			int tag = pp;
			mpi->recv(tmp_real   , length, pp, tag);
			mpi->recv(tmp_imag   , length, pp, tag);
			mpi->recv(tmp_real_VB, length, pp, tag);
			mpi->recv(tmp_imag_VB, length, pp, tag);
			
			// add to total matrix
			for (uint ee2=0; ee2<ppNE; ee2++) {
				uint ee = start_idx + ee2;
				for (uint xx=0; xx<Nx; xx++) {
					for (uint ww=0; ww<Nhw; ww++) {
						uint idx = xx*Nhw*ppNE + ee2*Nhw + ww;
						NEGF_ASSERT(idx<tmp_real.size() && idx<tmp_imag.size(), "tmp_real and/or tmp_imag has wrong size.");
						NEGF_ASSERT(QQ_total.size()==Nx && QQ_total[xx].size()==NE && QQ_total[xx][ee].size()==Nhw, "wrong size of QQ_total.");
						NEGF_ASSERT(abs(QQ_total[xx][ee][ww]-init_num) < 1e-14, "expected free place QQ_total[xx][ee][ww]");
						NEGF_ASSERT(fabs(tmp_real[idx]-init_num) > 1e-12 && fabs(tmp_imag[idx]-init_num) > 1e-12, "expected filled tmp_real, tmp_imag");
						QQ_total   [xx][ee][ww] = tmp_real   [idx] + constants::imag_unit * tmp_imag   [idx];
						QQ_total_VB[xx][ee][ww] = tmp_real_VB[idx] + constants::imag_unit * tmp_imag_VB[idx];
					}
				}
			}
		}
	} else {		
		// prepare array to send
		vector<double> tmp_real   ; tmp_real   .resize(Nx*myNE*Nhw, init_num);
		vector<double> tmp_imag   ; tmp_imag   .resize(Nx*myNE*Nhw, init_num);
		vector<double> tmp_real_VB; tmp_real_VB.resize(Nx*myNE*Nhw, init_num);
		vector<double> tmp_imag_VB; tmp_imag_VB.resize(Nx*myNE*Nhw, init_num);
		for (uint ee2=0; ee2<myNE; ee2++) {
			for (uint xx=0; xx<Nx; xx++) {
				for (uint ww=0; ww<Nhw; ww++) {
					uint idx = xx*Nhw*myNE + ee2*Nhw + ww;
					NEGF_ASSERT(idx<tmp_real.size() && idx<tmp_imag.size(), "tmp_real and/or tmp_imag has wrong size.");
					NEGF_ASSERT(QQ.size()==Nx && QQ[xx].size()==myNE && QQ[xx][ee2].size()==Nhw, "wrong size of QQ.");
					NEGF_ASSERT(fabs(tmp_real[idx]-init_num) < 1e-12 && fabs(tmp_imag[idx]-init_num) < 1e-12, "expected free place in tmp_real and tmp_imag.");
					tmp_real   [idx] = this->QQ   [xx][ee2][ww].real();
					tmp_imag   [idx] = this->QQ   [xx][ee2][ww].imag();
					tmp_real_VB[idx] = this->QQ_VB[xx][ee2][ww].real();
					tmp_imag_VB[idx] = this->QQ_VB[xx][ee2][ww].imag();
				}
			}
		}
		
		// send to master process
		int dest = constants::mpi_master_rank;
		int tag = mpi->get_rank();
		mpi->send(tmp_real   , dest, tag);
		mpi->send(tmp_imag   , dest, tag);
		mpi->send(tmp_real_VB, dest, tag);
		mpi->send(tmp_imag_VB, dest, tag);
	}
 
	// =====================================================================
	// COMPUTE Jphot, spectrum! (result is in master process only)
	// =====================================================================
	// master thread performs energy integration
	logmsg->emit(LOG_INFO,"Computing Jphot, spectrum");
	if (mpi->get_rank()==constants::mpi_master_rank) 
	{
		const double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
		
		// spin!
		// only add to jphot, spectrum and power are computed from this quantity
		double spin    = get_spin_degeneracy(options->get("kp_method"), quantities::electron_density); 
		double spin_VB = get_spin_degeneracy(options->get("kp_method"), quantities::hole_density); 
		NEGF_ASSERT(fabs(spin-spin_VB) < 1e-14,
				"A band structure which has different spin degeneracy in the CB and VB makes no sense in my opinion.");
		//spin = 1; // only for the moment...
		
		this->jphot       		 = Matc(Nhw, Nx);	// mode=lake: xx=Nx is not used
		this->jphot_VB    		 = Matc(Nhw, Nx);	
		this->spectrum    		 = Matd(Nhw, Nx);		
		this->spectrum_VB		 = Matd(Nhw, Nx);
		this->output_spectrum    = Matd(2, Nhw);		
		this->output_spectrum_VB = Matd(2, Nhw);
		this->output_power       = 0.0;
		this->output_power_VB    = 0.0;
		for (uint ww=0; ww<Nhw; ww++) 
		{
			double hw = this->hw_grid[ww];
			this->output_spectrum   (1,ww+1) = hw;
			this->output_spectrum_VB(1,ww+1) = hw;
			this->output_spectrum   (2,ww+1) = 0.0;
			this->output_spectrum_VB(2,ww+1) = 0.0;
			
			if (mode==galerpin) 
			{
				for (uint xx=0; xx<Nx; xx++) {
					for (uint ee=0; ee<NE; ee++) {
						double dE = energies->get_weight_from_global_idx(ee);
						this->jphot   (ww+1,xx+1) += spin * dE/(4.0*constants::pi*constants::pi*hbar) * hw *        this->QQ_total   [xx][ee][ww];
						this->jphot_VB(ww+1,xx+1) += spin * dE/(4.0*constants::pi*constants::pi*hbar) * hw * (-1) * this->QQ_total_VB[xx][ee][ww];
					}
				}
				
				for (uint xx=0; xx<Nx; xx++) {
					cplx Jphot    = this->jphot   (ww+1,xx+1);
					cplx Jphot_VB = this->jphot_VB(ww+1,xx+1);
					NEGF_FASSERT(fabs(Jphot   .imag()) < 100*constants::imag_err, "Luminescence: Jphot   .imag()=%e", Jphot   .imag());
					NEGF_FASSERT(fabs(Jphot_VB.imag()) < 100*constants::imag_err, "Luminescence: Jphot_VB.imag()=%e", Jphot_VB.imag());
					uint x1 = (xx>0)    ? xx-1 : xx;
					uint x2 = (xx<Nx-1) ? xx+1 : xx;
					double dx = 0.5 * (  xspace->get_vertex(xspace->get_global_vertex_index(x2))->get_coordinate(0)
						     	       - xspace->get_vertex(xspace->get_global_vertex_index(x1))->get_coordinate(0) );
					this->spectrum   (ww+1,xx+1) = 1.0/dx      * Jphot   .real() * hw;	// spectrally and spatially resolved recombination, [m-3s-1]
					this->spectrum_VB(ww+1,xx+1) = 1.0/dx      * Jphot_VB.real() * hw;
					this->output_spectrum   (2,ww+1) += this->area * Jphot   .real() * hw;	// spectrum, [s-1]
					this->output_spectrum_VB(2,ww+1) += this->area * Jphot_VB.real() * hw;
				}
			} else { 
				// -----------------
				// mode=lake
				// -----------------
				for (uint xx=0; xx<Nx-1; xx++) {
					for (uint ee=0; ee<NE; ee++) {
						double dE = energies->get_weight_from_global_idx(ee);
						
						this->jphot   (ww+1,xx+1) += spin * dE/(4.0*constants::pi*constants::pi*hbar) * hw * this->QQ_total   [xx][ee][ww];
						this->jphot_VB(ww+1,xx+1) += spin * dE/(4.0*constants::pi*constants::pi*hbar) * hw * this->QQ_total_VB[xx][ee][ww];
					}
				}
				
				// spatially resolved spectrum is the recombination spectrum multiplied by the energy
				// recombination is the divergence of the current
				// columns 1 and Nx will be empty
				for (uint xx=1; xx<Nx-1; xx++) 
				{
					double dx = 0.5*(  xspace->get_vertex(xspace->get_global_vertex_index(xx+1))->get_coordinate(0)
								     - xspace->get_vertex(xspace->get_global_vertex_index(xx-1))->get_coordinate(0));
					// photon current xx-1->xx is stored in jphot(xx)
					// photon current xx->xx+1 is stored in jphot(xx+1)
					cplx Jphot1 = this->jphot(ww+1,xx);
					cplx Jphot2 = this->jphot(ww+1,xx+1);
					NEGF_FASSERT(fabs(Jphot1.imag()) < 100*constants::imag_err, "Luminescence: Jphot1.imag()=%e", Jphot1.imag());
					NEGF_FASSERT(fabs(Jphot2.imag()) < 100*constants::imag_err, "Luminescence: Jphot2.imag()=%e", Jphot2.imag());
					this->spectrum(ww+1,xx+1) = hw * (Jphot2.real()-Jphot1.real()) / dx;
					cplx Jphot1_VB = this->jphot_VB(ww+1,xx);
					cplx Jphot2_VB = this->jphot_VB(ww+1,xx+1);
					NEGF_FASSERT(fabs(Jphot1_VB.imag()) < 100*constants::imag_err, "Luminescence: Jphot1_VB.imag()=%e", Jphot1_VB.imag());
					NEGF_FASSERT(fabs(Jphot2_VB.imag()) < 100*constants::imag_err, "Luminescence: Jphot2_VB.imag()=%e", Jphot2_VB.imag());
					this->spectrum_VB(ww+1,xx+1) = hw * (Jphot2_VB.real()-Jphot1_VB.real()) / dx;
					
					// integrate spatially resolved spectrum over space
					output_spectrum   (2,ww+1) += this->area * this->spectrum   (ww+1,xx+1) * dx;
					output_spectrum_VB(2,ww+1) += this->area * this->spectrum_VB(ww+1,xx+1) * dx;
				}
			}
			uint ww1 = (ww>0) ? ww-1 : ww;
			uint ww2 = (ww<Nhw-1) ? ww+1 : ww;
			double dhw = 0.5 * (this->hw_grid[ww2]-this->hw_grid[ww1]);
			this->output_power    += dhw * this->output_spectrum   (2,ww+1);
			this->output_power_VB += dhw * this->output_spectrum_VB(2,ww+1);
		}
		logmsg->emit(LOG_INFO,"Output power, calculated from radiative electron recombination (A=%.1ecm^2): %.5e[W]",
			this->area / constants::convert_from_SI(units::area, 1e-4),
			this->output_power    / constants::convert_from_SI(units::power, 1.0) );
		logmsg->emit(LOG_INFO,"Output power, calculated from radiative hole     recombination (A=%.1ecm^2): %.5e[W]",
			this->area / constants::convert_from_SI(units::area, 1e-4),
			this->output_power_VB / constants::convert_from_SI(units::power, 1.0) );
	}
	mpi->synchronize_processes();
);}


void Luminescence::write_recombination_to_file(const char * filename) throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(mpi->get_rank()==constants::mpi_master_rank, "spectrally resolved recombination is stored in master process only!");
	vector<double> xcoord;
	for (uint ii=0; ii<xspace->get_num_internal_vertices(); ii++) {
		xcoord.push_back(xspace->get_vertex(xspace->get_global_vertex_index(ii))->get_coordinate(0));
	}
	negf::write_xE_matrix(filename, this->spectrum, xcoord, this->hw_grid);
	string filename2(filename);
	filename2.append("_VB");
	negf::write_xE_matrix(filename2.c_str(), this->spectrum_VB, xcoord, this->hw_grid);
);}


void Luminescence::write_spectrum_to_file(const char * filename) throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(mpi->get_rank()==constants::mpi_master_rank, "luminescence spectrum is stored in master process only!");
	logmsg->emit(LOG_INFO, "Writing spectrum (in simulator units) to %s",filename);
	negf::write_matrix(filename, this->output_spectrum);
	
	string filename2(filename);
	filename2.append("_VB");
	logmsg->emit(LOG_INFO, "Writing spectrum (in simulator units) to %s",filename2.c_str());
	negf::write_matrix(filename2.c_str(), this->output_spectrum_VB);
);}


void Luminescence::snapshot(double voltage) throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(mpi->get_rank()==constants::mpi_master_rank, "snapshot should be taken in master process only!");
	
	if (this->power_snapshots.num_cols()!=3) {	
		// first voltage point
		this->power_snapshots = Matd(1,3);
		this->power_snapshots(1, 1) = voltage;
		this->power_snapshots(1, 2) = this->output_power    / constants::convert_from_SI(units::power, 1.0);
		this->power_snapshots(1, 3) = this->output_power_VB / constants::convert_from_SI(units::power, 1.0);
	} else {
		// enlarge snapshot matrix (very inefficient)
		uint nvoltages = power_snapshots.num_rows();
		Matd tmp(1,1); tmp = this->power_snapshots;
		this->power_snapshots = Matd(nvoltages+1,3);
		
		// fill previous points
		for (uint vv=1; vv<=nvoltages; vv++) {
			for (uint ii=1; ii<=3; ii++) {
				this->power_snapshots(vv,ii) = tmp(vv,ii);
			}
		}
		
		// fill new point
		this->power_snapshots(nvoltages+1, 1) = voltage;
		this->power_snapshots(nvoltages+1, 2) = this->output_power    / constants::convert_from_SI(units::power, 1.0) ; // storage in W/cm2!!!
		this->power_snapshots(nvoltages+1, 3) = this->output_power_VB / constants::convert_from_SI(units::power, 1.0) ; // storage in W/cm2!!!
	}
);}


void Luminescence::write_power_to_file(const char * filename) throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(mpi->get_rank()==constants::mpi_master_rank, "luminescence spectrum is stored in master process only!");
	string filename2(filename);
	filename2.append(".power");
	logmsg->emit(LOG_INFO, "Writing output power (in [W/cm2]) to %s",filename2.c_str());
	negf::write_matrix(filename2.c_str(), this->power_snapshots);
);}
