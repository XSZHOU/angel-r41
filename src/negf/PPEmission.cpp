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
#include "PPEmission.h"
using namespace negf;

PPEmission::PPEmission(const Overlap * ov_,
						const Geometry * xspace_, 
						const Kspace * kspace_, 
						const Energies * energies_, 
						const Options * options_,
						const Hamiltonian * ham_,
						const GreenFunctions * gf_,
						const SEPhotonSpontaneous * spont_):
	ov(ov_),
	xspace(xspace_),
	kspace(kspace_),
	energies(energies_),
	options(options_),
	ham(ham_),
	gf(gf_),
	spont(spont_),
	area(constants::convert_from_SI(units::area, 1e-4)) // 1cm2
{STACK_TRACE(
	NEGF_ASSERT(ov!=NULL && xspace!=NULL && kspace!=NULL && energies!=NULL && options!=NULL && ham!=NULL && spont!=NULL, "null pointer.");
	this->Nx = xspace->get_num_internal_vertices();
	this->Nn = options->get_num_degrees_of_freedom();
	this->Nk = kspace->get_number_of_points();
	this->NE = energies->get_number_of_points();
	this->myNE = energies->get_my_number_of_points();
	this->Nhw = 100;
	
	// --------------------------------
	// fix discretization of spectrum
	// --------------------------------
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
	this->QQ   .resize(Nx-1, tmp2); 
	this->QQ_VB.resize(Nx-1, tmp2); 
	this->jphot       = GEMatrix(Nhw, Nx-1);
	this->jphot_VB    = GEMatrix(Nhw, Nx-1);	
	this->spectrum    = DGEMatrix(Nhw, Nx); 	
	this->spectrum_VB = DGEMatrix(Nhw, Nx); 
);}



void PPEmission::determine_mpi_stuff()
{STACK_TRACE(
	logmsg->emit_small_header("Determining MPI communication stuff for spectrum");
		
	this->E0_idx = energies->get_start_global_idx(mpi->get_rank());
	double E0    = energies->get_energy_from_global_idx(E0_idx);
	this->E1_idx = energies->get_stop_global_idx(mpi->get_rank());
	double E1    = energies->get_energy_from_global_idx(E1_idx);
	
	// -----------------------------------------------------------------------------
	// compute the following quantities:
	// 1. E0_plus_Emin_idx - the energy just below E0+Emin, at E1
	// 2. E1_minus_Emin_idx - the energy just above E1+Emax, at most E_{n-1}
	// -----------------------------------------------------------------------------
	
	// determine lowest needed global energy index ABOVE interval
	this->E0_plus_Emin_idx = E1_idx;
	while (E0_plus_Emin_idx < NE-1 
		   && energies->get_energy_from_global_idx(E0_plus_Emin_idx) < E0+this->hw_min) {
		E0_plus_Emin_idx++;
	}
	if (E0_plus_Emin_idx>E1_idx) E0_plus_Emin_idx--;
	// determine highest needed global energy index ABOVE interval
	this->E1_plus_Emax_idx = E1_idx;
	while (E1_plus_Emax_idx < NE-1 
		   && energies->get_energy_from_global_idx(E1_plus_Emax_idx) < E1+this->hw_max) {
		E1_plus_Emax_idx++;
	}
	
	bool mpi_out = logmsg->get_master_thread_output_only();
	logmsg->set_master_thread_output_only(false);
	logmsg->emit(LOG_INFO,"p%d (energy indices %d...%d = %.3geV...%.3geV) needs indices %d(%.3geV)...%d(%.3geV)",
			mpi->get_rank(), E0_idx, E1_idx, E0, E1,
			E0_plus_Emin_idx,energies->get_energy_from_global_idx(E0_plus_Emin_idx),
			E1_plus_Emax_idx,energies->get_energy_from_global_idx(E1_plus_Emax_idx));
	logmsg->set_master_thread_output_only(mpi_out);
	
	// OBSOLETE: we want to get rid of the points E0_idx and Emax_idx in the auxiliary arrays because we don't need them	
	
	// allocate memory for pp: PP[xx][ee2][ee'] --> (Nx-1)*myNE*(E1_plus_Emax_idx-E0_plus_Emin_idx+1)
	vector<cplx> tmp1; tmp1.resize(E1_plus_Emax_idx-E0_plus_Emin_idx+1, 0.0);
	vector< vector<cplx> > tmp2; tmp2.resize(myNE, tmp1);
	this->PP   .clear(); this->PP   .resize(Nx-1, tmp2);
	this->PP_VB.clear(); this->PP_VB.resize(Nx-1, tmp2);
	
	mpi->synchronize_processes();
);}

void PPEmission::calculate()
{STACK_TRACE(
	this->determine_mpi_stuff();
	
	// ---------------------------------------------
	// initialize the matrix PP to zero
	// PP[xx][ee2][ee-this->E0_plus_Emin_idx]
	// ---------------------------------------------
	for (uint xx=0; xx<Nx-1; xx++) {
		NEGF_ASSERT(xx<PP.size(),"xx<PP.size()");
		for (uint ee2=0; ee2<myNE; ee2++) {
			NEGF_ASSERT(ee2<PP[xx].size(), "ee2<PP[xx].size()");
			for (uint ee=this->E0_plus_Emin_idx; ee<=this->E1_plus_Emax_idx; ee++) {
				NEGF_ASSERT(ee-this->E0_plus_Emin_idx < PP[xx][ee2].size(), "ee-E0_plus_Emin_idx < PP[xx][ee2].size()");
				PP[xx][ee2][ee-this->E0_plus_Emin_idx] = 0.0;
				PP_VB[xx][ee2][ee-this->E0_plus_Emin_idx] = 0.0;
			}
		}
	}
	
	this->communicate_As();
	
	this->compute_QQ_Jphot_spectrum();
);}


void PPEmission::determine_needed_processes(
					vector< vector<int> > & processes_needed,
					vector<int> & pp_E0_idx,
					vector<int> & pp_E1_idx,
					vector<int> & pp_E0_plus_Emin_idx,
					vector<int> & pp_E1_plus_Emax_idx) const
{STACK_TRACE(
	logmsg->emit(LOG_INFO,"Setting up processes_needed...");
	int root = constants::mpi_master_rank;
	int my_rank = mpi->get_rank();
	
	processes_needed.clear(); processes_needed.resize(mpi->get_num_procs());
	pp_E0_plus_Emin_idx.assign(mpi->get_num_procs(), -1);
	pp_E1_plus_Emax_idx.assign(mpi->get_num_procs(), -1);
	pp_E0_idx          .assign(mpi->get_num_procs(), -1);
	pp_E1_idx          .assign(mpi->get_num_procs(), -1);
	
	if (my_rank==root) 
	{
		for (int pp=0; pp < mpi->get_num_procs(); pp++) 
		{
			// process pp needs energies E0_minus_Emax_idx...E0_idx-1 and Emax_idx+1...E1_plus_Emax_idx
			
			// get these variables
			if (pp==my_rank) {
				pp_E0_plus_Emin_idx[pp]= this->E0_plus_Emin_idx;
				pp_E1_plus_Emax_idx[pp]= this->E1_plus_Emax_idx;
				pp_E0_idx[pp]           = this->E0_idx;
				pp_E1_idx[pp]           = this->E1_idx;
			} else {
				int source = pp;
				int tag = 1; mpi->recv(pp_E0_plus_Emin_idx[pp], source, tag);
				    tag = 2; mpi->recv(pp_E1_plus_Emax_idx[pp], source, tag);
				    tag = 3; mpi->recv(pp_E0_idx[pp]          , source, tag);
				    tag = 4; mpi->recv(pp_E1_idx[pp]          , source, tag);
			}
			
			// determine which processes compute the energies in question
			// and add them to the list
			for (int ee=pp_E0_plus_Emin_idx[pp]; ee <= pp_E1_plus_Emax_idx[pp]; ee++) 
			{
				int pp2 = energies->get_process_computing(uint(ee));
				if (pp2==pp) continue; // this is possible, ee could be E0_idx!
				bool pp2_already_in_list = false;
				for (uint ii=0; ii<processes_needed[pp].size(); ii++) {
					if (processes_needed[pp][ii]==pp2) {
						pp2_already_in_list = true;
						break;
					}
				}
				if (!pp2_already_in_list) {
					processes_needed[pp].push_back(pp2);
				}
			}
		}
		// screen output
		for (int pp=0; pp < mpi->get_num_procs(); pp++) {
			logmsg->emit_noendl(LOG_INFO_L3,"p%d needs ",pp);
			for (uint ii=0; ii < processes_needed[pp].size(); ii++) {
				logmsg->emit_noendl(LOG_INFO_L3,"p%d, ",processes_needed[pp][ii]);
			}
			logmsg->emit(LOG_INFO_L3,"");
		}
	} else {
		int E0_idx_copy           = this->E0_idx;
		int E1_idx_copy           = this->E1_idx;
		int E0_plus_Emin_idx_copy = this->E0_plus_Emin_idx;
		int E1_plus_Emax_idx_copy = this->E1_plus_Emax_idx;
		int tag = 1; mpi->send(E0_plus_Emin_idx_copy, root, tag);
			tag = 2; mpi->send(E1_plus_Emax_idx_copy, root, tag);
			tag = 3; mpi->send(E0_idx_copy          , root, tag);
			tag = 4; mpi->send(E1_idx_copy          , root, tag);
	}
	
	// -----------------------------------------------------
	// broadcast processed_needed and all the rest
	// -----------------------------------------------------
	logmsg->set_master_thread_output_only(true);
	logmsg->emit(LOG_INFO_L1,"Broadcasting processes_needed...");
	for (int pp=0; pp<mpi->get_num_procs(); pp++) {
		int num_processes_needed = 0;
		if (my_rank==root) {
			num_processes_needed = processes_needed[pp].size();
		}
		mpi->broadcast(num_processes_needed, root);
		if (my_rank!=root) {
			processes_needed[pp].resize(num_processes_needed);
		}
		mpi->broadcast(processes_needed[pp], root);
	}
	logmsg->emit(LOG_INFO_L1,"Broadcasting energy index arrays...");
	mpi->broadcast(pp_E0_plus_Emin_idx, root);
	mpi->broadcast(pp_E1_plus_Emax_idx, root);
	mpi->broadcast(pp_E0_idx          , root);
	mpi->broadcast(pp_E1_idx          , root);
	
	mpi->synchronize_processes();
);}

void PPEmission::communicate_As()
{STACK_TRACE(
	// ========================================================================
	// get AL from processes storing E0_plus_Emin_idx...E1_plus_Emax_idx
	// and compute the helper quantity "PP" (for an explanation-->report of S.Steiger)
	// PP = array of (Nx-1)*myNE*(E1_plus_Emax_idx-E0_plus_Emin_idx+1) complex numbers
	// ========================================================================
	int root = constants::mpi_master_rank;
	int my_rank = mpi->get_rank();
	bool mpi_out = logmsg->get_master_thread_output_only();
	
	// --------------------------------------------------------------------------------------------------------
	// set up an array which will store for every process all other processes it relies on (EXcluding itself)
	// computed only by master process
	// --------------------------------------------------------------------------------------------------------
	vector< vector<int> > processes_needed; 
	vector<int> pp_E0_plus_Emin_idx; 
	vector<int> pp_E1_plus_Emax_idx; 
	vector<int> pp_E0_idx;
	vector<int> pp_E1_idx;
	this->determine_needed_processes(processes_needed, pp_E0_idx, pp_E1_idx, pp_E0_plus_Emin_idx, pp_E1_plus_Emax_idx);
	
	// ----------------------------------------------------------------
	// zipped data arrays are already available in SEPhotonSpontaneous!
	// ----------------------------------------------------------------
	const vector<unsigned long>   & AL_real_comp_size  = spont->get_AL_real_comp_size();
	const vector<unsigned long>   & AL_imag_comp_size  = spont->get_AL_imag_comp_size();
	const vector<unsigned char *> & AL_real_compressed = spont->get_AL_real_compressed();
	const vector<unsigned char *> & AL_imag_compressed = spont->get_AL_imag_compressed();
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
		logmsg->set_master_thread_output_only(true);
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
		
		// no screen output! later on...	
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
						
			// sender sends his AL's to processes storing energies BELOW his energy interval
			for (uint ee=this->E0_idx; ee<=this->E1_idx; ee++) 
			{			
				// NO receivers BELOW current process
				// determine receiver ABOVE current process - if there is none, continue to next energy
				for (int pp=0; pp<mpi->get_num_procs(); pp++) 
				{
					bool pp_receives = receiver[pp] && pp_E0_plus_Emin_idx[pp]<=int(ee) && pp_E1_plus_Emax_idx[pp]>=int(ee);
					if (!pp_receives) continue;
					
					// so the current process needs to send AL[ee-this->E0_idx] plus the size information !
					buffer_size_needed += AL_real_comp_size[ee-this->E0_idx] + AL_imag_comp_size[ee-this->E0_idx] + 2*sizeof(unsigned long) + 4*MPI_BSEND_OVERHEAD;
				}	
			}
			} // if(i_am_sender)
		}
		logmsg->set_master_thread_output_only(mpi_out);
	
		// mark receivers as computed
		// needs to be performed in ALL threads (senders, receivers and those which are neither)
		for (int pp=0; pp < mpi->get_num_procs(); pp++) {
			if (receiver[pp]) {
				process_was_computed[pp] = true;
			}
		}
	} // while(true)
	logmsg->set_master_thread_output_only(false);
	logmsg->emit(LOG_INFO_L3,"p%d will need to send %d chars in total.", mpi->get_rank(), buffer_size_needed);
	logmsg->set_master_thread_output_only(mpi_out);
	mpi->synchronize_processes();
	
	// ------------------------------
	// allocate buffer!
	// ------------------------------
	logmsg->emit(LOG_INFO,"Allocating MPI buffer...");
	unsigned long buffersize_long = buffer_size_needed+1000;
	NEGF_ASSERT(buffersize_long < 2147483648, "Buffer size does not fit into an int!!@!");
	int buffersize = int(buffersize_long);
	char * buffer = new char[buffersize];
	int err = MPI_Buffer_attach(buffer, buffersize);
	NEGF_FASSERT(err==0, "MPI_Buffer_attach gave error %d",err);
	mpi->synchronize_processes();
	
	// ------------------------------
	// prepare Hamiltonian
	// ------------------------------
	logmsg->emit(LOG_INFO,"Preparing Hamiltonians...");
	const GEMatrix & M = ov->get_internal_overlap();
	NEGF_ASSERT(M.numCols()==(int)(Nx*Nn) && M.numRows()==M.numCols(), "wrong overlap matrix.");
	uint Nvert = xspace->get_num_vertices();
	GEMatrix H(Nvert*Nn,Nvert*Nn);
	GEMatrix Hsmall(Nx*Nn,Nx*Nn);
	GEMatrix HminusEM_tmp1(Nx*Nn,Nx*Nn);
	vector< GEMatrix >         HminusEM_tmp2; HminusEM_tmp2.resize(Nk, HminusEM_tmp1);
	vector< vector<GEMatrix> > HminusEM;      HminusEM.resize(myNE, HminusEM_tmp2);
	for (uint kk=0; kk<Nk; kk++) 
	{
		ham->get(kspace->get_point(kk), H);
		NEGF_ASSERT(H.numRows()==(int)(Nvert*Nn), "something went wrong.");
		
		// construct the Hamiltonian of the interior points only
		for (uint xx=1; xx<=Nx; xx++) {
			for (uint yy=1; yy<=Nx; yy++) {
				GEMatrix::View Hsmall_xy = Hsmall(_((xx-1)*Nn+1,xx*Nn),_((yy-1)*Nn+1,yy*Nn));
				uint gx = xspace->get_global_vertex_index(xx-1) + 1;
				uint gy = xspace->get_global_vertex_index(yy-1) + 1;
				const GEMatrix::ConstView H_xy = H(_((gx-1)*Nn+1,gx*Nn),_((gy-1)*Nn+1,gy*Nn));
				
				Hsmall_xy = H_xy;
			}
		}
	
		NEGF_ASSERT(xspace->get_dimension()==1, "only 1D is implemented!");
		for (uint ee2 = 0; ee2 < myNE; ee2++) 
		{
			uint ee = energies->get_global_index(ee2);
			const double E = energies->get_energy_from_global_idx(ee);
					
			GEMatrix & HmEM = HminusEM[ee2][kk];	
			HmEM = (-E) * M;
			HmEM = HmEM + Hsmall; // Hsmall+HminusEM would fail!
		}
	}	
	
	// ------------------------------
	// set up CB/VB distinction
	// ------------------------------
	vector<uint> cb_bands; options->get_conduction_degrees_of_freedom(cb_bands);
	vector<uint> vb_bands; options->get_valence_degrees_of_freedom(vb_bands);
	for (uint nn=0; nn<Nn; nn++) {
		bool cb = false;
		bool vb = false;
		for (uint ii=0; ii<cb_bands.size(); ii++) {	if (cb_bands[ii]==nn) {	cb = true; break; }	}
		for (uint ii=0; ii<vb_bands.size(); ii++) {	if (vb_bands[ii]==nn) {	vb = true; break; }	}
		NEGF_ASSERT((cb && !vb) || (!cb && vb), "band must be CB or VB!");
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
		logmsg->set_master_thread_output_only(true);
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
		logmsg->set_master_thread_output_only(false);
		if (my_rank==root) {
			logmsg->emit_noendl(LOG_INFO, "This time we have %d receivers: ",num_receivers);
			for (int pp=0; pp < mpi->get_num_procs(); pp++) {
				if (receiver[pp]) logmsg->emit_noendl(LOG_INFO, "%d ", pp);
			}
			logmsg->emit_noendl(LOG_INFO, " and %d senders: ",num_senders);
			for (int pp=0; pp < mpi->get_num_procs(); pp++) {
				if (sender[pp]) logmsg->emit_noendl(LOG_INFO, "%d ", pp);
			}
			logmsg->emit(LOG_INFO, "");
		}
		logmsg->set_master_thread_output_only(mpi_out);
		
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
			GEMatrix tmp(Nx*Nn,Nx*Nn);
			vector<GEMatrix> ALmat; ALmat.resize(Nk, tmp);
			
			double mpi_time = 0.0;
			double comp_time = 0.0;
			
			// receiver receives missing AL's ABOVE his energy interval
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
				logmsg->emit(LOG_INFO_L3,"p%d waits for ee=%d from p%d", my_rank, ee, sender_id);
				mpi->recv_antihermitians(ALmat, sender_id, tag);// get stuff from an energy below --> outscattering --> AL
				logmsg->emit(LOG_INFO_L3,"p%d just got ee=%d from p%d", my_rank, ee, sender_id);
				double t2 = MPI_Wtime();
				
				// ----------------------------------------------------------------------------
				// For all Ei (own energies) and k:
				// 1. Compute C(Ei,ee,k) = GR(Ei,k)*AL(ee,k)*GA(Ei,k)
				// 2. Compute t(x,x+1)*C(x+1,x) - t(x+1,x)*C(x,x+1) for all x
				// 3. Take trace, multiply with k-space weight, and store
				// ----------------------------------------------------------------------------
				GEMatrix CCtmp(Nx*Nn,Nx*Nn);
				for (uint ee2=0; ee2<myNE; ee2++) 
				{
					for (uint kk=0; kk<Nk; kk++) 
					{
						uint Ei_idx = energies->get_global_index(ee2);
						const GEMatrix & GR = gf->get_retarded(kk,Ei_idx);
						const GEMatrix & GA = gf->get_advanced(kk,Ei_idx);
						const GEMatrix & ALee = ALmat[kk];
						tmp = ALee * GA;
						CCtmp = GR * tmp;
						
						const GEMatrix & HmEM = HminusEM[ee2][kk];
							
						for (uint xx=0; xx<Nx-1; xx++) 
						{
							// get CC(xx,xx+1) and CC(xx+1,xx) and 
							const GEMatrix::ConstView CCtmp_part1 = CCtmp(_((xx+0)*Nn+1,(xx+1)*Nn),_((xx+1)*Nn+1,(xx+2)*Nn));
							const GEMatrix::ConstView CCtmp_part2 = CCtmp(_((xx+1)*Nn+1,(xx+2)*Nn),_((xx+0)*Nn+1,(xx+1)*Nn));
							
							// get the overlap-augmented hamiltonian matrix and extract the relevant parts
							const GEMatrix::ConstView HmEM_part1 = HmEM(_((xx+0)*Nn+1,(xx+1)*Nn),_((xx+1)*Nn+1,(xx+2)*Nn));	
							const GEMatrix::ConstView HmEM_part2 = HmEM(_((xx+1)*Nn+1,(xx+2)*Nn),_((xx+0)*Nn+1,(xx+1)*Nn));	
							
							// compute!
							GEMatrix PPP(Nn,Nn);
							PPP  = HmEM_part1 * CCtmp_part2;
							GEMatrix PPP2(Nn,Nn);
							PPP2 = HmEM_part2 * CCtmp_part1;
							PPP -= PPP2;				// PPP = t(x,x+1)*C(x+1,x) - t(x+1,x)*C(x,x+1)
							
							// take trace!
							/*cplx trace = 0.0;
							for (uint nn=1; nn<=Nn; nn++) {
								trace += PPP(nn,nn);
							}*/
							cplx trace_CB = 0.0;
							for (uint nn=0; nn<cb_bands.size(); nn++) {
								trace_CB += PPP(cb_bands[nn]+1,cb_bands[nn]+1);
							}
							cplx trace_VB = 0.0;
							for (uint nn=0; nn<vb_bands.size(); nn++) {
								trace_VB += PPP(vb_bands[nn]+1,vb_bands[nn]+1);
							}
							
							// multiply with k-weight and store
							double wk = kspace->get_point(kk).get_weight() / (constants::pi*constants::pi*4.0);
							NEGF_ASSERT(xx<PP.size(), "xx<PP.size() failed.");
							NEGF_ASSERT(ee2<PP[xx].size(), "ee2<PP[xx].size() failed.");
							NEGF_ASSERT(ee-this->E0_plus_Emin_idx < PP[xx][ee2].size(), "ee-this->E0_plus_Emin_idx<PP[xx][ee2].size() failed.");
							//this->PP[xx][ee2][ee-this->E0_plus_Emin_idx] += wk * trace; // PP was initialized to 0 at the beginning of this routine
							this->PP[xx][ee2][ee-this->E0_plus_Emin_idx] += wk * trace_CB; 
							this->PP_VB[xx][ee2][ee-this->E0_plus_Emin_idx] += wk * trace_VB;
						}
					}
				}
				
				double t3 = MPI_Wtime();
				mpi_time += t2-t1;
				comp_time += t3-t2;
				
			}
			logmsg->emit(LOG_ERROR,"receiver p%d needed %e for MPI and %e for computation",mpi->get_rank(),mpi_time,comp_time);
			
			
		} else 
		{
			if (!i_am_sender) { // this is possible!
			} else {
				
			// sender sends his AL's to processes storing energies BELOW his energy interval
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
				
					logmsg->emit(LOG_INFO_L3,"p%d AL sends ee=%d downwards to p%d", my_rank, ee, receiver_id);
					// BUFFERED SEND! mpi->send_antihermitians and mpi->recv_antihermitians send/receive the following 
					// quantities in order: real_char_size, imag_char_size, real_compressed, imag_compressed
					
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
					
					logmsg->emit(LOG_INFO_L3,"p%d AL just sent ee=%d to p%d", my_rank, ee, receiver_id);
				}	
			}
			//cout << "DONEs" << my_rank << endl;
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
	logmsg->set_master_thread_output_only(false);
	logmsg->emit_noendl(LOG_INFO,"p%d   ",mpi->get_rank());
	logmsg->set_master_thread_output_only(mpi_out);
	
	mpi->synchronize_processes();
);}

void PPEmission::compute_QQ_Jphot_spectrum()
{STACK_TRACE(
	// =====================================================================
	// COMPUTE QQ(xx,E_i,hw_j) from PP(xx,E_i,E_nu) by linear combination
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
			double E_plus_hw = E+hw;
			
			for (uint xx=0; xx<Nx-1; xx++) 
			{
				NEGF_ASSERT(xx<QQ.size(), "xx<QQ.size() failed.");
				NEGF_ASSERT(ee2<QQ[xx].size(), "ee2<QQ[xx].size() failed.");
				NEGF_ASSERT(ww<QQ[xx][ee2].size(), "ww<QQ[xx][ee2].size() failed.");
				this->QQ   [xx][ee2][ww] = 0.0;
				this->QQ_VB[xx][ee2][ww] = 0.0;
				
				if (E_plus_hw > energies->get_energy_from_global_idx(NE-1)) continue;
				
				// find energy grid points adjacent to E-hw				
				uint ee_lower = NE-1;
				while (ee_lower>0 && energies->get_energy_from_global_idx(ee_lower)>E_plus_hw) {
					ee_lower--;
				}
				if (ee_lower==NE-1) continue;
				uint ee_upper = ee_lower + 1;
				double Elow = energies->get_energy_from_global_idx(ee_lower);
				double Eupp = energies->get_energy_from_global_idx(ee_upper);
				double frac = (E_plus_hw - Elow) / (Eupp-Elow);
				NEGF_ASSERT(frac>=0.0 && frac<=1.0, "frac should be between 0 and 1.");
				
				// linear combination - frac=0 means E_minus_hw=Elow!
				NEGF_ASSERT(xx<PP.size(), "xx<PP.size() failed.");
				NEGF_ASSERT(ee2<PP[xx].size(), "ee2<PP[xx].size() failed.");
				NEGF_ASSERT(ee_lower>=this->E0_plus_Emin_idx, "ee_lower<this->E0_minus_Emax_idx");
				NEGF_ASSERT(ee_upper<=this->E1_plus_Emax_idx, "ee_upper>this->E1_minus_Emin_idx");
				NEGF_ASSERT(ee_lower-E0_plus_Emin_idx < PP[xx][ee2].size(), "ee_lower-E0_minus_Emax_idx<PP[xx][ee2].size() failed.");
				NEGF_ASSERT(ee_upper-E0_plus_Emin_idx < PP[xx][ee2].size(), "ee_upper-E0_minus_Emax_idx<PP[xx][ee2].size() failed.");
				this->QQ   [xx][ee2][ww] = (1.0-frac) * this->PP   [xx][ee2][ee_lower-E0_plus_Emin_idx] + frac * this->PP   [xx][ee2][ee_upper-E0_plus_Emin_idx];
				this->QQ_VB[xx][ee2][ww] = (1.0-frac) * this->PP_VB[xx][ee2][ee_lower-E0_plus_Emin_idx] + frac * this->PP_VB[xx][ee2][ee_upper-E0_plus_Emin_idx];				
			}
		}
	}
	
	mpi->synchronize_processes();
	
	// =====================================================================
	// COMPUTE Jphot, spectrum! (result is in master process only)
	// =====================================================================
	logmsg->emit(LOG_INFO,"Sending QQ to master thread");
	const double init_num = -8888.0;
	// send QQ-array to master thread
	if (mpi->get_rank()==constants::mpi_master_rank) 
	{
		vector<cplx> tmp1; tmp1.resize(Nhw, init_num);
		vector< vector<cplx> > tmp2; tmp2.resize(NE, tmp1);
		this->QQ_total.assign(Nx-1, tmp2); // array of size (Nx-1)*NE*Nhw
		this->QQ_total_VB.assign(Nx-1, tmp2); 
		
		// own contribution
		for (uint ee2=0; ee2<myNE; ee2++) {
			uint ee = energies->get_global_index(ee2);
			for (uint xx=0; xx<Nx-1; xx++) {
				for (uint ww=0; ww<Nhw; ww++) {
					NEGF_ASSERT(xx<QQ.size(), "xx<QQ.size() failed.");
					NEGF_ASSERT(ee2<QQ[xx].size(), "ee2<QQ[xx].size() failed.");
					NEGF_ASSERT(ww<QQ[xx][ee2].size(), "ww<QQ[xx][ee2].size() failed.");
					NEGF_ASSERT(xx<QQ_total.size(), "xx<QQ_tot.size() failed.");
					NEGF_ASSERT(ee<QQ_total[xx].size(), "ee<QQ_tot[xx].size() failed.");
					NEGF_ASSERT(ww<QQ_total[xx][ee].size(), "ww<QQ_tot[xx][ee].size() failed.");
					QQ_total   [xx][ee][ww] = this->QQ   [xx][ee2][ww];
					QQ_total_VB[xx][ee][ww] = this->QQ_VB[xx][ee2][ww];
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
			uint length = (Nx-1)*ppNE*Nhw;
			vector<double> tmp_real   ; tmp_real   .resize(length, init_num);
			vector<double> tmp_imag   ; tmp_imag   .resize(length, init_num);
			vector<double> tmp_real_VB; tmp_real_VB.resize(length, init_num);
			vector<double> tmp_imag_VB; tmp_imag_VB.resize(length, init_num);
			int tag = pp;
			mpi->recv(tmp_real, length, pp, tag);
			mpi->recv(tmp_imag, length, pp, tag);
			mpi->recv(tmp_real_VB, length, pp, tag);
			mpi->recv(tmp_imag_VB, length, pp, tag);
			
			// add to total matrix
			for (uint ee2=0; ee2<ppNE; ee2++) {
				uint ee = start_idx + ee2;
				for (uint xx=0; xx<Nx-1; xx++) {
					for (uint ww=0; ww<Nhw; ww++) {
						uint idx = xx*Nhw*ppNE + ee2*Nhw + ww;
						NEGF_ASSERT(idx<tmp_real.size() && idx<tmp_imag.size(), "tmp_real and/or tmp_imag has wrong size.");
						NEGF_ASSERT(xx<QQ_total.size(), "xx<QQ_tot.size() failed.");
						NEGF_ASSERT(ee<QQ_total[xx].size(), "ee<QQ_tot[xx].size() failed.");
						NEGF_ASSERT(ww<QQ_total[xx][ee].size(), "ww<QQ_tot[xx][ee].size() failed.");
						NEGF_ASSERT(abs(QQ_total[xx][ee][ww]-init_num) < 1e-14, "expected free place QQ_total[xx][ee][ww]");
						NEGF_ASSERT(fabs(tmp_real[idx]-init_num) > 1e-12, "expected filled tmp_real");
						NEGF_ASSERT(fabs(tmp_imag[idx]-init_num) > 1e-12, "expected filled tmp_imag");
						QQ_total   [xx][ee][ww] = tmp_real   [idx] + constants::imag_unit * tmp_imag   [idx];
						QQ_total_VB[xx][ee][ww] = tmp_real_VB[idx] + constants::imag_unit * tmp_imag_VB[idx];
					}
				}
			}
		}
	} else {		
		// prepare array to send
		vector<double> tmp_real;    tmp_real   .resize((Nx-1)*myNE*Nhw, init_num);
		vector<double> tmp_imag;    tmp_imag   .resize((Nx-1)*myNE*Nhw, init_num);
		vector<double> tmp_real_VB; tmp_real_VB.resize((Nx-1)*myNE*Nhw, init_num);
		vector<double> tmp_imag_VB; tmp_imag_VB.resize((Nx-1)*myNE*Nhw, init_num);
		for (uint ee2=0; ee2<myNE; ee2++) {
			for (uint xx=0; xx<Nx-1; xx++) {
				for (uint ww=0; ww<Nhw; ww++) {
					uint idx = xx*Nhw*myNE + ee2*Nhw + ww;
					NEGF_ASSERT(idx<tmp_real.size() && idx<tmp_imag.size(), "tmp_real and/or tmp_imag has wrong size.");
					NEGF_ASSERT(xx<QQ.size(), "xx<QQ>size() failed.");
					NEGF_ASSERT(ee2<QQ[xx].size(), "ee2<QQ[xx].size() failed.");
					NEGF_ASSERT(ww<QQ[xx][ee2].size(), "ww<QQ[xx][ee2].size() failed.");
					NEGF_FASSERT(fabs(tmp_real[idx]-init_num) < 1e-12, "expected free place in tmp_real: should be %e, is %e.",init_num,tmp_real[idx]);
					NEGF_ASSERT(fabs(tmp_imag[idx]-init_num) < 1e-12, "expected free place in tmp_imag.");
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
 
	// master thread performs energy integration
	logmsg->emit(LOG_INFO,"Computing Jphot, spectrum");
	if (mpi->get_rank()==constants::mpi_master_rank) 
	{
		const double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
		
		this->jphot       = GEMatrix(Nhw, Nx-1);	
		this->jphot_VB    = GEMatrix(Nhw, Nx-1);
		this->spectrum    = DGEMatrix(Nhw, Nx);
		this->spectrum_VB = DGEMatrix(Nhw, Nx);	
		this->output_spectrum   .assign(Nhw, 0.0);
		this->output_spectrum_VB.assign(Nhw, 0.0);
		
		for (uint ww=0; ww<Nhw; ww++) 
		{
			double hw = this->hw_grid[ww];
			
			for (uint xx=0; xx<Nx-1; xx++) {
				for (uint ee=0; ee<NE; ee++) {
					double dE = energies->get_weight_from_global_idx(ee);
					
					this->jphot   (ww+1,xx+1) += dE/(4.0*constants::pi*constants::pi*hbar) * hw * this->QQ_total   [xx][ee][ww];
					this->jphot_VB(ww+1,xx+1) += dE/(4.0*constants::pi*constants::pi*hbar) * hw * this->QQ_total_VB[xx][ee][ww];
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
				NEGF_FASSERT(fabs(Jphot1.imag()) < 100*constants::imag_err, "PPEmission: Jphot1.imag()=%e", Jphot1.imag());
				NEGF_FASSERT(fabs(Jphot2.imag()) < 100*constants::imag_err, "PPEmission: Jphot2.imag()=%e", Jphot2.imag());
				this->spectrum(ww+1,xx+1) = hw * (Jphot2.real()-Jphot1.real()) / dx;
				cplx Jphot1_VB = this->jphot_VB(ww+1,xx);
				cplx Jphot2_VB = this->jphot_VB(ww+1,xx+1);
				NEGF_FASSERT(fabs(Jphot1_VB.imag()) < 100*constants::imag_err, "PPEmission: Jphot1_VB.imag()=%e", Jphot1_VB.imag());
				NEGF_FASSERT(fabs(Jphot2_VB.imag()) < 100*constants::imag_err, "PPEmission: Jphot2_VB.imag()=%e", Jphot2_VB.imag());
				this->spectrum_VB(ww+1,xx+1) = hw * (Jphot2_VB.real()-Jphot1_VB.real()) / dx;
				
				// integrate spatially resolved spectrum over space
				output_spectrum   [ww] += this->area * this->spectrum   (ww+1,xx+1) * dx;
				output_spectrum_VB[ww] += this->area * this->spectrum_VB(ww+1,xx+1) * dx;
			}
		}
	}
	mpi->synchronize_processes();
);}


void PPEmission::write_recombination_to_file(const char * filename)
{STACK_TRACE(
	vector<double> xcoord;
	for (uint ii=0; ii<xspace->get_num_internal_vertices(); ii++) {
		xcoord.push_back(xspace->get_vertex(xspace->get_global_vertex_index(ii))->get_coordinate(0));
	}
	flens_io::write_xE_matrix(filename, this->spectrum, xcoord, this->hw_grid);
	string filename2(filename);
	filename2.append("VB");
	flens_io::write_xE_matrix(filename2.c_str(), this->spectrum_VB, xcoord, this->hw_grid);
)}


void PPEmission::write_spectrum_to_file(const char * filename)
{STACK_TRACE(
	ofstream fout(filename);
	NEGF_FASSERT(fout, "file %s could not be opened for write.", filename);	
	logmsg->emit(LOG_INFO_L2, "writing spectrum to %s",filename);
	fout.precision(12); 	
	for(uint ii = 0; ii < this->output_spectrum.size(); ii++) {
		fout << this->output_spectrum[ii] << "\n";
	}
	fout.close();
	
	string filename2(filename);
	filename2.append("VB");
	ofstream fout2(filename2.c_str());
	NEGF_FASSERT(fout2, "file %s could not be opened for write.", filename2.c_str());	
	logmsg->emit(LOG_INFO_L2, "writing spectrum to %s",filename2.c_str());
	fout2.precision(12); 	
	for(uint ii = 0; ii < this->output_spectrum_VB.size(); ii++) {
		fout2 << this->output_spectrum_VB[ii] << "\n";
	}
	fout2.close();
)}
