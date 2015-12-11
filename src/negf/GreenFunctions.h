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
#ifndef GREENFUNCTIONS_H_
#define GREENFUNCTIONS_H_

#include "all.h"

#include "NEGFObject.h"
#include "Overlap.h"
#include "Hamiltonian.h"
#include "SelfEnergy.h"

namespace negf {

#ifdef USE_BANDED
	typedef BandedNEGFObject GL_object_type;
#else
	typedef FullNEGFObject GL_object_type;
#endif
	
	/** Stores the electronic nonequilibrium Green functions G_retarded, G_advanced, G_lesser and G_greater 
	 *  Although GA and GG are easily computable elsewhere, they are stored anyway because they might
	 *  be needed at several other places and would hve to be computed several times otherwise */
	class GreenFunctions
	{
	public:
		GreenFunctions(const Options * options_, 
					   const Geometry * xspace_, 
					   const Kspace * kspace_, 
					   Energies * energies_,
					   const Hamiltonian * ham_,
					   const Overlap * ov_,
					   bool mangle_overlap_calc_ = true) throw (Exception *);
		~GreenFunctions();
		
		/** calculations involving Green functions only (no self-energies) and possibly overlap matrices 
		 *  these methods were moved here from class InnerLoop because interpolate_from_old_energies() relies on them */
		void 				calculate_retarded(const SelfEnergy * SEtot) throw (Exception *);
		void 				calculate_advanced() throw (Exception *);
		
		void 				calculate_lesser_greater(const SelfEnergy * SEtot) throw (Exception *);
		void 				calculate_lesser(const SelfEnergy * SEtot) throw (Exception *);
		void 				calculate_greater() throw (Exception *);
		void 				calculate_overlap_augmented_lesser() throw (Exception *);
		void 				calculate_overlap_augmented_greater() throw (Exception *);
		
		/** Interpolate GF from an old energy grid to the current one */
		void 				interpolate_from_old_energies() throw (Exception *);
		
		/** get interpolation coefficients c_i of 
		 *   
		 *    1/(Eupp-Elow) * int_Elow^Eupp f(E)dE    ~    sum_i c_i f(Ei) 
		 *
		 *  the resulting index array "indices" is 0-based (lowest energy has index 0)
		 */
		void 				get_interpol_coeffs(vector<uint> & indices, vector<double> & coeffs, const double & Elow, const double & Eupp) const;
		void 				get_interpol_coeffs(vector<uint> & indices, vector<double> & coeffs, const double & E) const; // simpler version
		
		/** Assign flat-band quasi-equilibrium values to GL and GG - only possible in a 2-band model */
		void 				assign_flat_band_lesser_greater(double EFn, double EFp, double Ec, double Ev) throw (Exception *);
		
		bool 				mangle_overlap_calculation() const { return mangle_overlap_calc; }
		
		/** functions to access entire NEGF objects (comprising all k, E) */ 
		FullNEGFObject * 	get_retarded() const { return this->GR; }
		FullNEGFObject * 	get_advanced() const { return this->GA; }
		GL_object_type * 	get_lesser()   const { return this->GL; }
		GL_object_type * 	get_greater()  const { return this->GG; }
		GL_object_type * 	get_overlap_augmented_lesser()  const { return this->MGLM; }
		GL_object_type * 	get_overlap_augmented_greater()  const { return this->MGGM; }
		
		/** functions to access locally stored matrices */
		Matc  & 			get_retarded(uint kidx, uint global_Eidx) const { return this->GR->get(kidx, global_Eidx); }
		Matc  & 			get_advanced(uint kidx, uint global_Eidx) const { return this->GA->get(kidx, global_Eidx); }
		GLMat & 			get_lesser  (uint kidx, uint global_Eidx) const { return this->GL->get(kidx, global_Eidx); }
		GLMat & 			get_greater (uint kidx, uint global_Eidx) const { return this->GG->get(kidx, global_Eidx); }
		GLMat & 			get_overlap_augmented_lesser  (uint kidx, uint global_Eidx) const { return this->MGLM->get(kidx, global_Eidx); }
		GLMat & 			get_overlap_augmented_greater (uint kidx, uint global_Eidx) const { return this->MGGM->get(kidx, global_Eidx); }
	
		/** other trivial access functions */
		const Geometry 	* 	get_xspace() 	const { return this->xspace; }
		const Kspace 	* 	get_kspace() 	const { return this->kspace; }
		const Options 	* 	get_options() 	const { return this->options; }
		Energies 		* 	get_energies()  const { return this->energies; }
		const Overlap   *   get_overlap()   const { return ov; }
		
	protected:
		
		// helper function; T will be BandedNEGFObject or FullNEGFObject, Mtype will be Matc or BMatc
		template<class T, class Mtype>
		void interpolate_from_old_energies(T * gf, const vector<int> & E1_index, const vector<double> & E1_coeff, const Mtype & empty_matrix); // called from interpolate_from_old_energies() only
		
		FullNEGFObject * GR;		// retarded GF 
		FullNEGFObject * GA;		// advanced GF (given by GA = GR+)
		GL_object_type * GL;		// lesser   GF
		GL_object_type * GG;		// greater  GF (given by GG = GR - GA + GL)
		GL_object_type * MGLM;		// = M*GL*M, needed for various self-energies
		GL_object_type * MGGM;		// = M*GG*M, needed for various self-energies
		
		const Options  * 	options;	// stores what scattering is used, etc.
		const Geometry * 	xspace;		// real-space grid
		const Kspace   * 	kspace;		// k-space grid (including integration weights for each point)
		Energies * 			energies;	// 1D energy space grid (including integration weights for each point)
		const Overlap * 	ov;
		const Hamiltonian * ham;
		uint 				Nx;
		uint 				NxNn;
		
		const bool 			security_checking;
		const bool 			mangle_overlap_calc;
	};
	
	
	// T will be Matc or BMatc
	// Gnew[Ei] = E1_coeff[Ei] * Gold[E1_index[Ei]] + (1-E1_coeff[Ei]) * Gold[E1_index[Ei]+1]
	template<class T, class Mtype>
	void GreenFunctions::interpolate_from_old_energies(T * gf, const vector<int> & E1_index, const vector<double> & E1_coeff, const Mtype & empty_matrix)
	{//STACK_TRACE(
		const vector<double> & new_grid = this->energies->get_energy_grid();
		const vector<double> & old_grid = this->energies->get_old_energy_grid();
		NEGF_ASSERT(old_grid.size()==new_grid.size(), "inconsistent sizes of old and new energy grid.");
		uint     nE = energies->get_number_of_points();
		uint   myNE = energies->get_my_number_of_points();
		NEGF_ASSERT(new_grid.size()==nE, "inconsistent number of energy points.");
		uint     nk = kspace->get_number_of_points();
		int    root = constants::mpi_master_rank;
		int my_rank = mpi->get_rank();
		
		// initialize temporary array
		//BMatc tmp = BMatc::create<offdiags>(NxNn, NxNn);
		vector<Mtype>           tmpvec;      tmpvec    .resize(nk, empty_matrix);
		vector<Mtype>           mpi_tmpvec;  mpi_tmpvec.resize(nk, empty_matrix);
		vector< vector<Mtype> > new_Gs;      new_Gs    .resize(myNE, tmpvec);
		
		// note: interpolation does not mix k's!
			
		// array which will store for every process all other processes it relies on (EXcluding itself)
		vector< vector<int> > processes_needed; processes_needed.resize(mpi->get_num_procs());
		
		// -------------------------------------------------------------
		// first do all possible local calculations in parallel
		// -------------------------------------------------------------
		if (my_rank==root) logmsg->emit_noendl(LOG_INFO,"local calculations...  ");
		vector<bool> computed;     computed.resize(nE, false);
		vector<bool> mpi_involved; mpi_involved.resize(nE, false);
		for (int pp=0; pp < mpi->get_num_procs(); pp++)
		{	
			for (uint ee = energies->get_start_global_idx(uint(pp)); ee <= energies->get_stop_global_idx(uint(pp)); ee++)
			{
				const uint ee2 = ee - energies->get_start_global_idx(uint(pp));
				int old_idx = E1_index[ee];
				
				// if new energy lies outside old energy interval, assign zero matrix
				if (old_idx==-1)
				{
					if (my_rank==pp) {
						NEGF_ASSERT(ee2 < new_Gs.size(),  "new_Gs is wrong.");
						NEGF_ASSERT(new_Gs[ee2].size()==nk, "new_Gs is wrong.");
						for (uint kk=0; kk < nk; kk++) {
							Mtype & G = new_Gs[ee2][kk];
							G = empty_matrix;
						}
					}
					computed[ee] = true;	// will be overwritten lots of times - doesn't matter. needs to be stored in all processes
					continue;
				}
				
				computed[ee] = true; // will be overwritten lots of times - doesn't matter. needs to be stored in all processes
				
				// check if MPI communication is necessary for new energy point ee 
				int proc1 = energies->get_process_computing(uint(old_idx));
				if (proc1!=pp && fabs(E1_coeff[ee])>1e-13) {
					computed[ee]          = false;
					mpi_involved[old_idx] = true;
					bool proc1_in_list    = false;
					for (uint ii=0; ii < processes_needed[pp].size(); ii++) {
						if (processes_needed[pp][ii]==proc1) proc1_in_list = true;
					}
					if (!proc1_in_list) processes_needed[pp].push_back(proc1);
				} else if (my_rank==pp && fabs(E1_coeff[ee])>1e-13) {
					for (uint kk=0; kk < nk; kk++) {
						Mtype & G = new_Gs[ee2][kk];
						add(gf->get(kk, old_idx), E1_coeff[ee], G); // G += E1_coeff[ee] * gf->get(kk, old_idx);
					}
				}
				
				int proc2 = energies->get_process_computing(uint(old_idx)+1);
				if (proc2!=pp && fabs(1.0-E1_coeff[ee])>1e-13) {
					computed[ee]            = false;
					mpi_involved[old_idx+1] = true;
					bool proc2_in_list      = false;
					for (uint ii=0; ii < processes_needed[pp].size(); ii++) {
						if (processes_needed[pp][ii]==proc2) proc2_in_list = true;
					}
					if (!proc2_in_list && proc2!=pp) processes_needed[pp].push_back(proc2);
				} else if (my_rank==pp && fabs(1.0-E1_coeff[ee])>1e-13) {
					for (uint kk=0; kk < nk; kk++) {
						Mtype & G = new_Gs[ee2][kk];
						add(gf->get(kk, old_idx+1),(1.0-E1_coeff[ee]), G); // G += (1.0-E1_coeff[ee]) * gf->get(kk, old_idx+1);
					}
				}
			}
		}
		// at this point computed[ee]=false means that new energy ee needs MPI communication 
		
		// compress the locally stored matrices
		if (my_rank==root) logmsg->emit_noendl(LOG_INFO,"compressing...  ");
		mpi->synchronize_processes();
		double *        G_real_data[myNE];			// will store copies of the complex G split into real and imag parts
		double *        G_imag_data[myNE];
		unsigned char * G_real_char[myNE];			// will store pointer to same memory as G_...._data, but different type
		unsigned char * G_imag_char[myNE];
		unsigned long   num_chars[myNE];			// will store how many chars the array of complex matrices corresponds to
		unsigned long   G_real_comp_size[myNE];		// will store the size of the zipped data
		unsigned long   G_imag_comp_size[myNE];
		unsigned char * G_real_compressed[myNE];	// will store pointers to the zipped data
		unsigned char * G_imag_compressed[myNE];
		for (uint ee2=0; ee2<myNE; ee2++) {		
			uint ee = energies->get_global_index(ee2);
			if (mpi_involved[ee]) {
				negf_math::do_compress<Mtype>(gf->get(ee), G_real_data[ee2], G_imag_data[ee2], G_real_char[ee2], G_imag_char[ee2], num_chars[ee2], 
							G_real_comp_size[ee2], G_imag_comp_size[ee2], G_real_compressed[ee2], G_imag_compressed[ee2]);	
			}
		}
		
		
		// --------------------------------------------------------------------------------------
		// determine how much data is sent in total
		// we do this by going through the same algorithm as when sending the data, just w/o sending
		// necessary for the amount of buffer that needs to be allocated for a safe operation
		// --------------------------------------------------------------------------------------
		logmsg->emit(LOG_INFO,"Computing total amount of data to be sent...    ");
		unsigned long buffer_size_needed = 0;
		vector<bool> process_was_computed; 
		process_was_computed.resize(mpi->get_num_procs(), false); 
		while (true)
		{
			vector<bool> receiver; receiver.resize(mpi->get_num_procs(), false);
			vector<bool> sender;   sender.resize(mpi->get_num_procs(), false);
			
			bool done = false;
			
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
			
			// no screen output! later on...	
			
			done = (num_receivers == 0);
			if (done) {
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
				// NO!!! if (!sender[my_rank]) continue; // the process was already computed and is not needed for any receiver at all
				
				// go through all receiver processes and look which of them need the current energy
				for (int pp=0; pp < mpi->get_num_procs(); pp++) 
				{	
					if (!receiver[pp]) continue;
					
					for (uint ee = energies->get_start_global_idx(uint(pp)); ee <= energies->get_stop_global_idx(uint(pp)); ee++)
					{
						if (computed[ee]) continue; // energy was already computed - MPI-receiver sets this to true only AFTER the communication was done
						
						int old_idx = E1_index[ee];
						NEGF_ASSERT(old_idx>=0, "something is wrong");
						
						// set up contribution from lower part
						int proc = energies->get_process_computing(uint(old_idx));
						if (my_rank==proc && fabs(E1_coeff[ee])>1e-13) 	//  old_idx is owned by the current process --> send it
						{
							const uint ee2 = old_idx - energies->get_start_global_idx();
							NEGF_FASSERT(ee2 < myNE, "(1) ee2=%d, myNE=%d",ee2,myNE);
							NEGF_ASSERT(mpi_involved[old_idx], "expected mpi_involved[old_idx]==true");
							
							buffer_size_needed += G_real_comp_size[ee2] + G_imag_comp_size[ee2] + 2*sizeof(unsigned long)+ 8*MPI_BSEND_OVERHEAD;
						}
						
						// set up contribution from upper part - if needed
						proc = energies->get_process_computing(uint(old_idx)+1);
						if (proc==my_rank && fabs(1.0-E1_coeff[ee])>1e-13) {	//  old_idx+1 is owned by the current process --> send it;
							const uint ee2 = old_idx - energies->get_start_global_idx() + 1;
							NEGF_FASSERT(ee2 < myNE, "(2) ee2=%d, myNE=%d",ee2,myNE);
							NEGF_ASSERT(mpi_involved[old_idx+1], "expected mpi_involved[old_idx+1]==true");
							
							buffer_size_needed += G_real_comp_size[ee2] + G_imag_comp_size[ee2] + 2*sizeof(unsigned long)+ 8*MPI_BSEND_OVERHEAD;
						}
					}
				}
					
				//} // if(i_am_sender)
			}
		
			// mark receivers as computed. needs to be performed in ALL threads (senders, receivers and those which are neither)
			for (int pp=0; pp < mpi->get_num_procs(); pp++) {
				if (receiver[pp]) {
					process_was_computed[pp] = true;
					for (uint ee = energies->get_start_global_idx(pp); ee <= energies->get_stop_global_idx(pp); ee++) {
						computed[ee] = true;
					}
				}
			}
		} // while(true)
		logmsg->emit_all(LOG_INFO_L3,"p%d will need to send %d chars in total.", my_rank, buffer_size_needed);
		mpi->synchronize_processes();
		
		// ------------------------------
		// allocate buffer!
		// ------------------------------
		logmsg->emit(LOG_INFO,"Allocating MPI buffer...   ");
		unsigned long buffersize_long = buffer_size_needed+1000;
		unsigned long max_array_size = 2147483648UL; // =2.1GB
		//NEGF_FASSERT(buffersize_long < max_array_size, "Buffer size (%ld) does not fit into an int!!@!",buffersize_long);
		const bool blocking_send = (buffersize_long >= max_array_size);
		int buffersize = 0;
		char * buffer  = 0;
		int err = 0;
		if (blocking_send) {
			logmsg->emit_all(LOG_WARN,"\n****WARNING****\np%d buffer size (%ld) does not fit into an int!!@! Falling back to blocking MPI_Send.\n\n",
					mpi->get_rank(), buffersize_long);
		} else {
			buffersize = int(buffersize_long);
			buffer     = new char[buffersize];
			err        = MPI_Buffer_attach(buffer, buffersize);
			NEGF_FASSERT(err==0, "MPI_Buffer_attach gave error %d",err);
		}
		mpi->synchronize_processes();
		
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
			
			bool done = false;
			
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
			
			done = (num_receivers == 0);
			if (done) {
				break;
			}
			
			// some screen output
			if (my_rank==root) {
				logmsg->emit_noendl_all(LOG_INFO, "This time we have %d receivers ",num_receivers);
				for (int pp=0; pp < mpi->get_num_procs(); pp++) {
					if (receiver[pp]) logmsg->emit_noendl_all(LOG_INFO_L2, "%d ", pp);
				}
				logmsg->emit_noendl_all(LOG_INFO_L3, " and %d senders ",num_senders);
				for (int pp=0; pp < mpi->get_num_procs(); pp++) {
					if (sender[pp]) logmsg->emit_noendl_all(LOG_INFO_L3, "%d ", pp);
				}
				logmsg->emit_all(LOG_INFO_L2, "");
			}			
					
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
				for (uint ee = energies->get_start_global_idx(); ee <= energies->get_stop_global_idx(); ee++)
				{
					if (computed[ee]) continue; // energy was already computed
					
					const uint ee2 = ee - energies->get_start_global_idx();
					int old_idx = E1_index[ee];
					NEGF_ASSERT(old_idx>=0, "something is wrong");
				
					// set up MPI-communicated contribution from lower part
					int proc = energies->get_process_computing(uint(old_idx));
					if (proc!=my_rank && fabs(E1_coeff[ee])>1e-13) { // receive from other process
						NEGF_ASSERT(mpi_involved[old_idx], "expected mpi_involved[old_idx]==true");
						
						int tag = my_rank;
						mpi->recv(mpi_tmpvec, proc, tag);	// MPI::recv(vector<Matc> & data, int & source, int & tag);
						
						for (uint kk=0; kk<nk; kk++) {
							Mtype & G = new_Gs[ee2][kk];
							add(mpi_tmpvec[kk], E1_coeff[ee], G); // G = E1_coeff[ee] * mpi_tmpvec[kk];
						}
					} 
					
					// set up MPI-communicated contribution from upper part
					proc = energies->get_process_computing(uint(old_idx)+1);
					if (proc!=my_rank && fabs(1.0-E1_coeff[ee])>1e-13) { // receive from other process
						NEGF_ASSERT(mpi_involved[old_idx+1], "expected mpi_involved[old_idx+1]==true");
						
						int tag = my_rank;
						mpi->recv(mpi_tmpvec, proc, tag);	// MPI::recv(vector<Mtype> & data, int & source, int & tag);
						
						for (uint kk=0; kk<nk; kk++) {
							Mtype & G = new_Gs[ee2][kk];
							add(mpi_tmpvec[kk], (1.0-E1_coeff[ee]), G); // G += (1.0-E1_coeff[ee]) * mpi_tmpvec[kk];
						}
					} 
					
					computed[ee] = true; 		
				}
			} else 
			{	
				// NO!!! if (!sender[my_rank]) continue; // the process was already computed and is not needed for any receiver at all
				
				// go through all receiver processes and look which of them need the current energy
				for (int pp=0; pp < mpi->get_num_procs(); pp++) 
				{	
					if (!receiver[pp]) continue;
					int receiver_id = pp;
					
					for (uint ee = energies->get_start_global_idx(uint(pp)); ee <= energies->get_stop_global_idx(uint(pp)); ee++)
					{
						if (computed[ee]) continue; // energy was already computed - MPI-receiver sets this to true only AFTER the communication was done
						
						int old_idx = E1_index[ee];
						NEGF_ASSERT(old_idx>=0, "something is wrong");
						
						// set up contribution from lower part
						int proc = energies->get_process_computing(uint(old_idx));
						int tag = pp;
						int tag2 = tag+1;
						if (my_rank==proc && fabs(E1_coeff[ee])>1e-13) 	//  old_idx is owned by the current process --> send it
						{				
							NEGF_ASSERT(mpi_involved[old_idx], "expected mpi_involved[old_idx]==true");
							
							const uint ee2 = old_idx - energies->get_start_global_idx();
							NEGF_FASSERT(ee2 < myNE, "(3) ee2=%d, myNE=%d",ee2,myNE);
							int real_char_size = int(G_real_comp_size[ee2]);
							int imag_char_size = int(G_imag_comp_size[ee2]);
							if (!blocking_send) {
								// ----------------------
								// buffered send
								// ----------------------
								
								// send array lengths
								err = MPI_Bsend(&real_char_size, 1, MPI_INT, receiver_id, tag, MPI_COMM_WORLD);
								NEGF_FASSERT(err==0, "MPI send of G_real_char_size failed: err=%d",err);
								err = MPI_Bsend(&imag_char_size, 1, MPI_INT, receiver_id, tag2, MPI_COMM_WORLD);
								NEGF_FASSERT(err==0, "MPI send of G_imag_char_size failed: err=%d",err);
					
								// send arrays
								err = MPI_Bsend(G_real_compressed[ee2], real_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag, MPI_COMM_WORLD);
								NEGF_FASSERT(err==0, "MPI send of G_real_compressed failed: err=%d",err);
								err = MPI_Bsend(G_imag_compressed[ee2], imag_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag2, MPI_COMM_WORLD);
								NEGF_FASSERT(err==0, "MPI send of G_imag_compressed failed: err=%d",err);
							} else {
								// ----------------------
								// blocking send
								// ----------------------
								
								// send array lengths
								err = MPI_Send(&real_char_size, 1, MPI_INT, receiver_id, tag, MPI_COMM_WORLD);
								NEGF_FASSERT(err==0, "MPI send of G_real_char_size failed: err=%d",err);
								err = MPI_Send(&imag_char_size, 1, MPI_INT, receiver_id, tag2, MPI_COMM_WORLD);
								NEGF_FASSERT(err==0, "MPI send of G_imag_char_size failed: err=%d",err);
					
								// send arrays
								err = MPI_Send(G_real_compressed[ee2], real_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag, MPI_COMM_WORLD);
								NEGF_FASSERT(err==0, "MPI send of G_real_compressed failed: err=%d",err);
								err = MPI_Send(G_imag_compressed[ee2], imag_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag2, MPI_COMM_WORLD);
								NEGF_FASSERT(err==0, "MPI send of G_imag_compressed failed: err=%d",err);
							}
						}
						
						// set up contribution from upper part 
						proc = energies->get_process_computing(uint(old_idx)+1);
						if (proc==my_rank && fabs(1.0-E1_coeff[ee])>1e-13) 	//  old_idx+1 is owned by the current process --> send it;
						{						
							NEGF_ASSERT(mpi_involved[old_idx+1], "expected mpi_involved[old_idx+1]==true");
							
							const uint ee2 = old_idx - energies->get_start_global_idx() + 1;
							NEGF_FASSERT(ee2 < myNE, "(4) ee2=%d, myNE=%d",ee2, myNE);
							int real_char_size = int(G_real_comp_size[ee2]); 
							int imag_char_size = int(G_imag_comp_size[ee2]);
							if (!blocking_send) {
								// ----------------------
								// buffered send
								// ----------------------
								
								// send array lengths
								err = MPI_Bsend(&real_char_size, 1, MPI_INT, receiver_id, tag, MPI_COMM_WORLD);
								NEGF_FASSERT(err==0, "MPI send of G_real_char_size failed: err=%d",err);
								err = MPI_Bsend(&imag_char_size, 1, MPI_INT, receiver_id, tag2, MPI_COMM_WORLD);
								NEGF_FASSERT(err==0, "MPI send of G_imag_char_size failed: err=%d",err);
					
								// send arrays
								err = MPI_Bsend(G_real_compressed[ee2], real_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag, MPI_COMM_WORLD);
								NEGF_FASSERT(err==0, "MPI send of G_real_compressed failed: err=%d",err);
								err = MPI_Bsend(G_imag_compressed[ee2], imag_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag2, MPI_COMM_WORLD);
								NEGF_FASSERT(err==0, "MPI send of G_imag_compressed failed: err=%d",err);
							} else {
								// ----------------------
								// blocking send
								// ----------------------
								
								// send array lengths
								err = MPI_Send(&real_char_size, 1, MPI_INT, receiver_id, tag, MPI_COMM_WORLD);
								NEGF_FASSERT(err==0, "MPI send of G_real_char_size failed: err=%d",err);
								err = MPI_Send(&imag_char_size, 1, MPI_INT, receiver_id, tag2, MPI_COMM_WORLD);
								NEGF_FASSERT(err==0, "MPI send of G_imag_char_size failed: err=%d",err);
					
								// send arrays
								err = MPI_Send(G_real_compressed[ee2], real_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag, MPI_COMM_WORLD);
								NEGF_FASSERT(err==0, "MPI send of G_real_compressed failed: err=%d",err);
								err = MPI_Send(G_imag_compressed[ee2], imag_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag2, MPI_COMM_WORLD);
								NEGF_FASSERT(err==0, "MPI send of G_imag_compressed failed: err=%d",err);
							}
						}
					}
				}
			}
		
			// mark receivers as finished and corresponding energies as computed
			// needs to be performed in ALL threads (senders, receivers and those which are neither)
			for (int pp=0; pp < mpi->get_num_procs(); pp++) {
				if (receiver[pp]) {
					process_was_computed[pp] = true;
					for (uint ee = energies->get_start_global_idx(pp); ee <= energies->get_stop_global_idx(pp); ee++) {
						computed[ee] = true;
					}
				}
			}
		}
		// security checks
		for (int pp=0; pp < mpi->get_num_procs(); pp++) {
			NEGF_ASSERT(process_was_computed[pp], "a process was not computed.");
		}
		for (uint ee=0; ee < nE; ee++) {
			NEGF_FASSERT(computed[ee], "an energy was not computed! ee=%d (p%d)",ee,energies->get_process_computing(ee));
		}
		
		// now the array new_Gs has been set up for every process and we may overwrite the old Green functions
		// this can be done in parallel
		for (uint ee = energies->get_start_global_idx(); ee <= energies->get_stop_global_idx(); ee++)
		{
			NEGF_ASSERT(computed[ee], "an energy was not computed!");
			uint ee2 = ee - energies->get_start_global_idx();
			for (uint kk=0; kk<nk; kk++) {
				Mtype & Gnew = new_Gs[ee2][kk];
				Mtype & Gold = gf->get(kk,ee);
				Gold = Gnew;
			}
		}
			
		if (my_rank==constants::mpi_master_rank) logmsg->emit(LOG_INFO,"");
		
		// ----------------------------------------------------------------
		// release buffer
		// ----------------------------------------------------------------
		if (!blocking_send) {
			logmsg->emit(LOG_INFO_L2,"Deallocating MPI buffer");
			err = MPI_Buffer_detach(&buffer, &buffersize);
			NEGF_FASSERT(err==0, "MPI_Buffer_attach gave error %d",err);
			delete [] buffer;
		}
		
		// ----------------------------------------------------------------
		// release zipped data
		// ----------------------------------------------------------------
		logmsg->emit(LOG_INFO_L2,"Deallocating zipped arrays");
		for (uint ee2=0; ee2<myNE; ee2++) {
			uint ee = energies->get_global_index(ee2);
			if (mpi_involved[ee]) {
				delete [] G_real_data[ee2];
				delete [] G_imag_data[ee2]; 
				delete [] G_real_compressed[ee2];
				delete [] G_imag_compressed[ee2];
			}
		}
		mpi->synchronize_processes();
	}//);}

	
} // end of namespace

#endif /*GREENFUNCTIONS_H_*/
