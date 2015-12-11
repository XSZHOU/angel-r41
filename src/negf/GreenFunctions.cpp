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
#include "GreenFunctions.h"
using namespace negf;

GreenFunctions::GreenFunctions(const Options * options_, 
							   const Geometry * xspace_, 
							   const Kspace * kspace_, 
							   Energies * energies_,
							   const Hamiltonian * ham_,
							   const Overlap * ov_,
							   bool mangle_overlap_calc_ /* true by default */) throw (Exception *):
	options(options_),
	xspace(xspace_),
	kspace(kspace_),
	energies(energies_),
	ov(ov_),
	ham(ham_),
	security_checking(true),
	mangle_overlap_calc(mangle_overlap_calc_)	// if true, GL, GG, MGGM and MGLM will all be computed in calculate_lesser_greater
{STACK_TRACE(
	logmsg->emit_header("setting up Green functions");
	NEGF_ASSERT(options!=NULL && xspace!=NULL && kspace!=NULL && energies!=NULL && ov!=NULL, "null pointer encountered.");
	
	this->Nx   = xspace->get_num_internal_vertices();
	this->NxNn = Nx*Nn;
	
	// create GF for the energies of the calling process
	GR = new FullNEGFObject(xspace, kspace, energies, options);
	GA = new FullNEGFObject(xspace, kspace, energies, options);
	GL = new BandedNEGFObject(xspace, kspace, energies, options, constants::odGL);
	GG = new BandedNEGFObject(xspace, kspace, energies, options, constants::odGL);
	
	MGLM = new BandedNEGFObject(xspace, kspace, energies, options, constants::odGL);
	MGGM = new BandedNEGFObject(xspace, kspace, energies, options, constants::odGL);
	
	mpi->synchronize_processes();
);}

							   
GreenFunctions::~GreenFunctions() 
{STACK_TRACE(
	delete GR;
	delete GA;
	delete GL;
	delete GG;
	delete MGLM;
	delete MGGM;
);}


/** calculate the retarded GF using the Dyson equation for each k- and energy point */
void GreenFunctions::calculate_retarded(const SelfEnergy/*<odSE>*/ * SEtot) throw (Exception *)
{STACK_TRACE(
	logmsg->emit_header("calculating GR");
	uint Nk    = kspace->get_number_of_points();
	uint myNE  = energies->get_my_number_of_points();
	
	const OVMat & M = ov->get_internal_overlap();
	NEGF_ASSERT(M.num_cols()==NxNn && M.num_rows()==M.num_cols(), "wrong overlap matrix.");
	
	// initialize helper matrices
	Matc  Hsmall(Nx*Nn,Nx*Nn);
	OVMat EM = OVMat_create(NxNn);
	
	// small numerical parameter to make E a little bit complex
	// if omitted, NaN's will appear sometimes
	const double eta = constants::convert_from_SI(units::energy, constants::dyson_eta * constants::SIec);
	
	for (uint ee2 = 0; ee2 < myNE; ee2++)
	{	
		uint ee = energies->get_global_index(ee2);
		mult(M, energies->get_energy_from_global_idx(ee)+eta, EM); // EM = (energies->get_energy_from_global_idx(ee)+eta) * M;
		
		if (ee % /*5*/10 == 0) logmsg->emit_noendl_all(LOG_INFO, "p%d: GR(E=%d,:)...   ",mpi->get_rank(),ee);
		for (uint kk = 0; kk < Nk; kk++)
		{			
			// get Hamiltonian of the interior points only, including electrostatic potential
			ham->get_internal(kspace->get_point(kk), Hsmall);
			NEGF_ASSERT(Hsmall.num_rows()==Nx*Nn, "something went wrong.");
					
			// get retarded GF and self-energy
			      Matc  & GRm    =  this->get_retarded(kk,ee);
			const SEMat & SigmaR = SEtot->get_retarded(kk,ee);
			
			// --------------------------
			// solve Dyson equation
			// --------------------------
			
			// assign GR = (E+eta)M - H - SigmaR
			GRm  = SigmaR;
			GRm *= -1.0;
			GRm -= Hsmall;
			GRm += EM;
			
			// invert
			invert(GRm);
		}
	}
	
	// timing
	double time = timer->get_time();
	uint slowest_process = constants::mpi_master_rank;
	string slowest_hostname = mpi->get_hostname();
	if (mpi->get_rank()==constants::mpi_master_rank) {
		double max_time = time;
		for (int proc=0; proc < mpi->get_num_procs(); proc++) {
			if (proc==constants::mpi_master_rank) continue;
			int source = proc;
			double proc_time = 0.0;
			string proc_hostname; uint num_chars = 100;
			
			int tag = 832;
			mpi->recv(proc_time, source, tag);
			tag = 833;
			mpi->recv(proc_hostname, num_chars, source, tag);
			
			if (proc_time > max_time) {
				max_time = proc_time;
				slowest_process = proc;
				slowest_hostname = proc_hostname;
			}
		}
		// strip white spaces from hostname
		while (slowest_hostname[slowest_hostname.length()-1]==' ') {
			slowest_hostname.erase(slowest_hostname.length()-1);
		}
	} else {
		int dest = constants::mpi_master_rank;
		int tag = 832;
		mpi->send(time, dest, tag);
		
		string host = mpi->get_hostname();
		host.resize(100, ' ');
		tag = 833;
		mpi->send(host, dest, tag);
	}
	mpi->synchronize_processes();
	logmsg->emit(LOG_INFO, "");
	if (mpi->get_rank()==constants::mpi_master_rank) logmsg->emit(LOG_INFO, "Process %d (%s) was slowest.",slowest_process, slowest_hostname.c_str());
);}


/** calculate the advanced GF for each k- and own energy point */
void GreenFunctions::calculate_advanced() throw (Exception *)
{STACK_TRACE(
	logmsg->emit_header("calculating GA");
	uint Nk    = kspace->get_number_of_points();
	uint myNE  = energies->get_my_number_of_points();
	
	for (uint ee2 = 0; ee2 < myNE; ee2++)
	{	
		uint ee = energies->get_global_index(ee2);
		
		for (uint kk = 0; kk < Nk; kk++)
		{
			const Matc & GRmat = this->get_retarded(kk,ee);
			Matc & GAmat = this->get_advanced(kk,ee);
			      
			conjtrans(GRmat, GAmat);
		}
	}
	mpi->synchronize_processes();
);}

/** calculate the lesser GF using the Keldysh equation for each k- and energy point */
void GreenFunctions::calculate_lesser_greater(const SelfEnergy/*<odSE>*/ * SEtot) throw (Exception *)
{STACK_TRACE(
	if (!mangle_overlap_calc) {
		logmsg->emit_header("calculating GL and GG");
	} else {
		logmsg->emit_header("calculating GL, M*GL*M, GG and M*GG*M");
	}
	uint Nk   = kspace->get_number_of_points();
	uint myNE = energies->get_my_number_of_points();
			
	const OVMat & M = ov->get_internal_overlap();
	NEGF_ASSERT(M.num_cols()==NxNn && M.num_rows()==M.num_cols(), "wrong overlap matrix.");
	
	// initialize helper matrices
	Matc SLGA(NxNn,NxNn);	// will store SL*GA
	Matc GLtmp(NxNn,NxNn);	// will store the full GL
	Matc GLM(NxNn,NxNn);	// will store GL*M
	Matc GGtmp(NxNn,NxNn);	// will store the full GG
	Matc GGM(NxNn,NxNn); 	// will store GG*M
	
	for (uint ee2 = 0; ee2 < myNE; ee2++)
	{
		uint ee = energies->get_global_index(ee2);
		if (Nk==1) {
			if (ee % 20 == 0) logmsg->emit_noendl_all(LOG_INFO, "p%d: GL(E=%d,:)...   ",mpi->get_rank(),ee);
		} else {
			if (ee % /*5*/10 == 0) logmsg->emit_noendl_all(LOG_INFO, "p%d: GL(E=%d,:)...   ",mpi->get_rank(),ee);
		}
		for (uint kk = 0; kk < Nk; kk++)
		{
			// get Green's functions
			const Matc & GRm = this->get_retarded(kk,ee);
			const Matc & GAm = this->get_advanced(kk,ee);
			     GLMat & GLm = this->get_lesser(kk,ee);
			     GLMat & GGm = this->get_greater(kk,ee);
						
			// get lesser self-energy
			const SEMat & SLm = SEtot->get_lesser(kk,ee);
			
			// solve Keldysh equation (double multiplication is not supported)
			mult(SLm, GAm, SLGA); // SLGA = SL * GA;	
			
			if (!mangle_overlap_calc) {
				mult(GRm, SLGA, GLm); // GL = GR * SLGA;
				
				GGm  = GRm;
				GGm -= GAm;
				GGm += GLm;
			} else {
				// store full matrix GL in GLtmp
				mult(GRm, SLGA, GLtmp);
				GLm = GLtmp;
				
				// calculate overlap-augmented M*GL*M
				GLMat & MGLMm = this->get_overlap_augmented_lesser(kk,ee);
				mult(GLtmp, M, GLM);
				mult(M, GLM, MGLMm);
				
				// calculate the full GG (using the full GL)
				GGtmp  = GRm;
				GGtmp -= GAm;
				GGtmp += GLtmp;
				GGm    = GGtmp;
				
				// calculate M*GG*M
				GLMat & MGGMm = this->get_overlap_augmented_greater(kk,ee);
				mult(GGtmp, M, GGM);
				mult(M, GGM, MGGMm);				
			}
			
			// anti-hermiticity check
			if (security_checking) {
				GLMat tmp = GLMat_create(NxNn); conjtrans(GLm, tmp); tmp += GLm;
				double tmp_norm = negf_math::matrix_norm(tmp);
				NEGF_FASSERT(tmp_norm < constants::antiherm_check, "InnerLoop: anti-hermiticity check of GL failed (delta=%e) : kk=%d, ee=%d=%.4g (ee2=%d)", tmp_norm, kk, ee, energies->get_energy_from_global_idx(ee), ee2);
			}
		}
	}
	mpi->synchronize_processes();
	logmsg->emit(LOG_INFO, "");
);}


/** calculate the lesser GF using the Keldysh equation for each k- and energy point */
void GreenFunctions::calculate_lesser(const SelfEnergy/*<odSE>*/ * SEtot) throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(!mangle_overlap_calc, "do not call this when GX and M*GX*M are calculated elsewhere");
	logmsg->emit_header("calculating GL");
	uint Nk   = kspace->get_number_of_points();
	uint myNE = energies->get_my_number_of_points();
	
	// initialize helper matrices
	Matc SLGA(NxNn,NxNn);
	
	for (uint ee2 = 0; ee2 < myNE; ee2++)
	{
		uint ee = energies->get_global_index(ee2);
		if (Nk==1) {
			if (ee % 20 == 0) logmsg->emit_noendl_all(LOG_INFO, "p%d: GL(E=%d,:)...   ",mpi->get_rank(),ee);
		} else {
			if (ee % 5 == 0) logmsg->emit_noendl_all(LOG_INFO, "p%d: GL(E=%d,:)...   ",mpi->get_rank(),ee);
		}
		for (uint kk = 0; kk < Nk; kk++)
		{			
			// get retarded and advanced GF
			const Matc & GRm = this->get_retarded(kk,ee);
			const Matc & GAm = this->get_advanced(kk,ee);
			
			// get lesser self-energy
			const SEMat & SLm = SEtot->get_lesser(kk,ee);
			
			// solve Keldysh equation (double multiplication is not supported)
			mult(SLm, GAm, SLGA);
			GLMat & GLm = this->get_lesser(kk,ee);
			mult(GRm, SLGA, GLm);
			
			// anti-hermiticity check
			if (security_checking) {
				GLMat tmp = GLMat_create(NxNn); conjtrans(GLm, tmp); tmp += GLm;
				double tmp_norm = negf_math::matrix_norm(tmp);
				NEGF_FASSERT(tmp_norm < constants::antiherm_check, "InnerLoop: anti-hermiticity check of GL failed (delta=%e) : kk=%d, ee=%d=%.4g (ee2=%d)", tmp_norm, kk, ee, energies->get_energy_from_global_idx(ee), ee2);
			}
		}
	}
	mpi->synchronize_processes();
	logmsg->emit(LOG_INFO, "");
);}


/** calculate the greater GF from GR, GR+ and G<: G> = GR - GA + G< */
void GreenFunctions::calculate_greater() throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(!mangle_overlap_calc, "do not call this when GX and M*GX*M are calculated elsewhere");
	logmsg->emit_header("calculating GG");
	
	uint Nk   = kspace->get_number_of_points();
	uint myNE = energies->get_my_number_of_points();
		
	for (uint ee2 = 0; ee2 < myNE; ee2++)
	{
		uint ee = energies->get_global_index(ee2);
		logmsg->emit_noendl_all(LOG_INFO_L2, "p%d: GG(E=%d,k=:)...   ",mpi->get_rank(),ee);
		for (uint kk = 0; kk < Nk; kk++)
		{			
			const Matc  & GRm = this->get_retarded(kk,ee);
			const Matc  & GAm = this->get_advanced(kk,ee);
			const GLMat & GLm = this->get_lesser(kk,ee); 
			      GLMat & GGm = this->get_greater(kk,ee);
			
			GGm  = GRm;
			GGm -= GAm;
			GGm += GLm;
		}
	}
	mpi->synchronize_processes();
);}


void GreenFunctions::calculate_overlap_augmented_lesser() throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(!mangle_overlap_calc, "do not call this when GX and M*GX*M are calculated elsewhere");
	logmsg->emit_header("calculating M*GL*M");
	uint Nk   = kspace->get_number_of_points();
	uint myNE = energies->get_my_number_of_points();
	
	const OVMat & M = ov->get_internal_overlap();
	NEGF_ASSERT(M.num_cols()==NxNn && M.num_rows()==M.num_cols(), "wrong overlap matrix.");

	Matc tmp(NxNn,NxNn);

	for (uint ee2 = 0; ee2 < myNE; ee2++)
	{
		uint ee = energies->get_global_index(ee2);
		for (uint kk = 0; kk < Nk; kk++)
		{			
			const GLMat & GLm   = this->get_lesser(kk,ee);
			      GLMat & MGLMm = this->get_overlap_augmented_lesser(kk,ee);
			
			// security check for anti-Hermiticity 
			if (security_checking) { 
				GLMat tmp2 = GLMat_create(NxNn);
				conjtrans(GLm,tmp2);
				tmp2 += GLm;
				double tmp_norm = negf_math::matrix_norm(tmp2);
				NEGF_FASSERT(tmp_norm < constants::antiherm_check, "GL was not anti-Hermitian: |GL - (-GL+)|=%e",tmp_norm);
			}
			
			mult(GLm, M, tmp); // tmp = GLmat * M;		mult(BMatc, BMatc,  Matc)
			mult(M,tmp,MGLMm); // MGLMmat = M * tmp;	mult(BMatc,  Matc, BMatc)
			
			// security check for anti-Hermiticity 
			if (security_checking) { 
				GLMat tmp2 = GLMat_create(NxNn);
				conjtrans(MGLMm, tmp2);
				tmp2 += MGLMm;
				double tmp_norm = negf_math::matrix_norm(tmp2);
				NEGF_FASSERT(tmp_norm < constants::antiherm_check, "MGLM was not anti-Hermitian: |MGLM - (-MGLM+)|=%e",tmp_norm);
			}
		}
	}
	mpi->synchronize_processes();
);}


void GreenFunctions::calculate_overlap_augmented_greater() throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(!mangle_overlap_calc, "do not call this when GX and M*GX*M are calculated elsewhere");
	logmsg->emit_header("calculating M*GG*M");
	uint Nk   = kspace->get_number_of_points();
	uint myNE = energies->get_my_number_of_points();
	
	const OVMat & M = ov->get_internal_overlap();
	NEGF_ASSERT(M.num_cols()==NxNn && M.num_rows()==M.num_cols(), "wrong overlap matrix.");
	
	Matc tmp(NxNn,NxNn);
	
	for (uint ee2 = 0; ee2 < myNE; ee2++)
	{
		uint ee = energies->get_global_index(ee2);
		for (uint kk = 0; kk < Nk; kk++)
		{			
			const GLMat & GGm   = this->get_greater(kk,ee);
			      GLMat & MGGMm = this->get_overlap_augmented_greater(kk,ee);
			
			mult(GGm, M, tmp);   // tmp = GGmat * M;	mult(BMatc, BMatc,  Matc)
			mult(M, tmp, MGGMm); // MGGMmat = M * tmp;	mult(BMatc,  Matc, BMatc)
		}
	}
	mpi->synchronize_processes();
);}



void GreenFunctions::interpolate_from_old_energies() throw (Exception *)
{STACK_TRACE(
	const vector<double> & new_grid = this->energies->get_energy_grid();
	const vector<double> & old_grid = this->energies->get_old_energy_grid();
	NEGF_ASSERT(old_grid.size()==new_grid.size(), "inconsistent sizes of old and new energy grid.");
	uint nE = energies->get_number_of_points();
	NEGF_ASSERT(new_grid.size()==nE, "inconsistent number of energy points.");
	
	// -------------------------------------------------------------------------------
	// compute interpolation coefficients old-->new:
	// E_new[ii] = E1_coeff[ii] * E1_index[ii] + (1-E1_coeff[ii]) * E1_index[ii+1]
	// -------------------------------------------------------------------------------
	vector<int>    E1_index; E1_index.resize(nE, -2);
	vector<double> E1_coeff; E1_coeff.resize(nE, 0.0);
	
	for (uint ii=0; ii < nE; ii++) {
		double E = new_grid[ii];
		if (E < old_grid[0] || E >= old_grid[nE-1]) {	// >= is necessary s.th. idx<nE later; topmost energy should not be populated anyway
			// new energy is outside the old grid --> will assign zero
			E1_index[ii] = -1;
			E1_coeff[ii] = 0.0;
			continue;
		}
		// determine index of old energy grid which is just ABOVE energy E
		uint idx = 0;
		while(old_grid[idx] <= E) {
			idx++;
			NEGF_ASSERT(idx < nE, "something went wrong (idx>=nE).");
		}
		NEGF_ASSERT(idx>0, "something went wrong (idx=0).");
		
		// assign E1_index, E1_coeff
		E1_index[ii] = idx-1;
		double E1 = old_grid[idx-1];
		double E2 = old_grid[idx];
		NEGF_ASSERT(E>=E1 && E<=E2, "something went wrong (E1<=E<=E2)");
		E1_coeff[ii] = (E2 - E) / (E2 - E1); // E = lambda*E1 + (1-lambda)*E2; lambda=0 --> E=E2, lambda=1 --> E=E1
	}
	
	// security check
	if (security_checking) {
		uint Nk = kspace->get_number_of_points();
		GLMat tmp = GLMat_create(NxNn);
		for (uint ee2=0; ee2<energies->get_my_number_of_points(); ee2++) {
			uint ee = energies->get_global_index(ee2);
			for (uint kk=0; kk < Nk; kk++) {
				const GLMat & GLm = this->get_lesser(kk,ee);
				conjtrans(GLm, tmp); // tmp = conjugateTranspose(GLm);
				tmp += GLm;
				double tmp_norm = negf_math::matrix_norm(tmp);
				NEGF_FASSERT(tmp_norm < constants::antiherm_check, "first anti-hermiticity check failed (delta=%e) : kk=%d, ee=%d=%.4g (ee2=%d)", tmp_norm, kk, ee, energies->get_energy_from_global_idx(ee), ee2);
			}
		}
	}
	
	// communicate GR and GL, then CALCULATE GA, GG, MGLM, MGGM!
	logmsg->emit_small_header("Interpolating GR");
	Matc empty_full_matrix(NxNn,NxNn);
	this->interpolate_from_old_energies<FullNEGFObject, Matc>(this->get_retarded(), E1_index, E1_coeff, empty_full_matrix);
	mpi->synchronize_processes();
	
	this->calculate_advanced();
	
	logmsg->emit_small_header("Interpolating GL");
	GLMat empty_GL_matrix = GLMat_create(NxNn);
	this->interpolate_from_old_energies<GL_object_type, GLMat>(this->get_lesser(), E1_index, E1_coeff, empty_GL_matrix);
	mpi->synchronize_processes();
	
	if (!mangle_overlap_calc) {
		this->calculate_greater();
		if (kspace->get_number_of_points()!=1) { // we don't need MGLM and MGGM in that case --> save time
			this->calculate_overlap_augmented_lesser();
			this->calculate_overlap_augmented_greater();
		}
	} else {
		logmsg->emit_small_header("Interpolating GG");
		this->interpolate_from_old_energies<GL_object_type, GLMat>(this->get_greater(), E1_index, E1_coeff, empty_GL_matrix);
		mpi->synchronize_processes();
		
		if (kspace->get_number_of_points()!=1) { // we don't need MGLM and MGGM in that case --> save time
			logmsg->emit_small_header("Interpolating M*GL*M");
			this->interpolate_from_old_energies<GL_object_type, GLMat>(this->get_overlap_augmented_lesser(), E1_index, E1_coeff, empty_GL_matrix);
			mpi->synchronize_processes();
			logmsg->emit_small_header("Interpolating M*GG*M");
			this->interpolate_from_old_energies<GL_object_type, GLMat>(this->get_overlap_augmented_greater(), E1_index, E1_coeff, empty_GL_matrix);
			mpi->synchronize_processes();
		}
	}
	
	// security check
	if (security_checking) {
		uint Nk = kspace->get_number_of_points();
		GLMat tmp = GLMat_create(NxNn);
		for (uint ee2=0; ee2<energies->get_my_number_of_points(); ee2++) {
			uint ee = energies->get_global_index(ee2);
			for (uint kk=0; kk < Nk; kk++) {
				const GLMat & GLm = this->get_lesser(kk,ee);
				conjtrans(GLm, tmp); // tmp = conjugateTranspose(GLm);
				tmp += GLm;
				double tmp_norm = negf_math::matrix_norm(tmp);
				NEGF_FASSERT(tmp_norm < constants::antiherm_check, "second anti-hermiticity check failed (delta=%e) : kk=%d, ee=%d=%.4g (ee2=%d)", tmp_norm, kk, ee, energies->get_energy_from_global_idx(ee), ee2);
			}
		}
	}
);}



/** determine the linear interpolation E = a*E_i + (1-a)*E_{i+1}, E_i<=E<E_i+1 */
void GreenFunctions::get_interpol_coeffs(vector<uint> & indices, vector<double> & coeffs, const double & E) const
{STACK_TRACE(
	uint nE = energies->get_number_of_points();
	indices.clear();
	coeffs.clear();
	
	// determine the highest index below Elow
	uint E0idx = nE-1;
	while (E0idx > 0 && energies->get_energy_from_global_idx(E0idx) > E) {
		E0idx--;
	}
	
	double E0 = energies->get_energy_from_global_idx(E0idx);
	if (E0idx==nE-1 || E0>E) {
		// E lies outside computed range
		return;
	}
	
	double E1 = energies->get_energy_from_global_idx(E0idx+1);
	
	
	indices.resize(2, 0);
	indices[0] = E0idx;
	indices[1] = E0idx + 1;
	coeffs.resize(2, 0.0);
	coeffs[0] = 1.0 - (E-E0)/(E1-E0);
	coeffs[1] =       (E-E0)/(E1-E0);	// E=E0 --> coeffs[0] = 1 , coeffs[1] = 0
	//coeffs[0] = 0.5;
	//coeffs[1] = 0.5;
	return;
);}


/** determine the coefficients c_i of 
 *
 *              1/(Eupp-Elow) * int_Elow^Eupp f(E)dE    ~    sum_i c_i f(Ei) 
 *
 *  to do so, we assume that f behaves linearly between the enbergy grid points
 * */
void GreenFunctions::get_interpol_coeffs(vector<uint> & indices, vector<double> & coeffs, const double & Elow, const double & Eupp) const
{STACK_TRACE(
	NEGF_FASSERT(fabs(Eupp-Elow)>1e-15, "Elow=Eupp=%.5e are identical! Will yield nan coefficients.", Elow);
	uint nE = energies->get_number_of_points();
	indices.clear();
	coeffs.clear();
	
	// determine the highest index below Elow
	uint E0idx = nE-1;
	while (E0idx > 0 && energies->get_energy_from_global_idx(E0idx) > Elow) {
		E0idx--;
	}
	
	// determine the lowest index above Eupp
	uint Enidx = 0;
	while (Enidx < nE-1 && energies->get_energy_from_global_idx(Enidx) < Eupp) {
		Enidx++;
	}
	
	// if Elow, Eupp are outside the energy range, we're done
	if (E0idx >= Enidx) {
		return;
	}
	
	// set up indices array (containing global energy indices)
	for (uint ii=E0idx; ii<=Enidx; ii++) {
		indices.push_back(ii);
	}
	
	// ----------------------
	// find coefficients
	// ----------------------
	const bool use_lumped_scheme = false;
	coeffs.resize(indices.size(), 0.0);
	
	uint NE = energies->get_number_of_points();
	for (uint ii=0; ii<indices.size(); ii++) {
		// boundaries of interval belonging to energy grid point indices[ii]
		//double Ilow = (indices[ii]>0)    ? energies->get_energy_from_global_idx(indices[ii]-1) : energies->get_energy_from_global_idx(indices[ii]);
		//double Iupp = (indices[ii]<NE-1) ? energies->get_energy_from_global_idx(indices[ii]+1) : energies->get_energy_from_global_idx(indices[ii]);
		// <ss 8.6.2010> was wrong!
		double E1 = energies->get_energy_from_global_idx(indices[ii]);
		double E0 = (indices[ii]>0   ) ? energies->get_energy_from_global_idx(indices[ii]-1) : E1;
		double E2 = (indices[ii]<NE-1) ? energies->get_energy_from_global_idx(indices[ii]+1) : E1;
		double Ilow = (E0+E1)/2;
		double Iupp = (E1+E2)/2;
		//if (indices[ii]<20) logmsg->emit(LOG_INFO,"indices[%d]=%d, E=%g", ii, indices[ii], E1);
		
		// determine the length of interval [Elow,Eupp] which is contained in [Ilow,Iupp]
		if (Ilow >= Eupp || Iupp <= Elow) {
			coeffs[ii] = 0.0;
		} else if (Ilow >= Elow && Iupp <= Eupp) {
			coeffs[ii] = Iupp - Ilow;
		} else if (Ilow <= Elow && Iupp >= Eupp) {
			coeffs[ii] = Eupp - Elow;
		} else if (Ilow >= Elow && Iupp >= Eupp) {
			coeffs[ii] = Eupp - Ilow;
		} else if (Ilow <= Elow && Iupp <= Eupp) {
			coeffs[ii] = Iupp - Elow;
		} else {
			NEGF_FEXCEPTION("We forgot a case: Elow=%e, Eupp=%e, Ilow=%e, Iupp=%e", Elow, Eupp, Ilow, Iupp);
		}
	}
	
	// take away points with coefficient ~0
	vector<uint> indices_old = indices;
	vector<double> coeffs_old = coeffs;
	indices.clear(); coeffs.clear();
	for (uint ii=0; ii<indices_old.size(); ii++) {
		if (fabs(coeffs_old[ii])>1e-14) {
			indices.push_back(indices_old[ii]);
			coeffs .push_back(coeffs_old[ii]);
		}
	}
	
	// DIVIDE ALL BY Eupp-Elow
	double simulated_interval = Eupp-Elow;
	for (uint ii=0; ii<coeffs.size(); ii++) {
		double E0 = energies->get_energy_from_global_idx(0);
		if (Elow<E0) { simulated_interval = Eupp-E0; }
		double EN = energies->get_energy_from_global_idx(energies->get_number_of_points()-1);
		if (Eupp>EN) { simulated_interval = EN-Elow; }
		
		coeffs[ii] = coeffs[ii] / simulated_interval;
		NEGF_FASSERT(!isnan(coeffs[ii]) && !isinf(coeffs[ii]), "Elow=%.3e, Eupp=%.3e: indices[%d]=%d, coeffs[%d]=%.3e", Elow,Eupp,ii,indices[ii],ii,coeffs[ii]);
	}
		
	// security checks
	// check that the function '1' interpolated gives 1
	if (security_checking) {
		double check = 0.0;
		for (uint ii=0; ii<coeffs.size(); ii++) {
			if (coeffs[ii]<-1e-14) {
				logmsg->emit_all(LOG_ERROR,"Elow=%e, Eupp=%e, E0idx=%d(E=%e), Enidx=%d(E=%e)",
						Elow,Eupp,E0idx,energies->get_energy_from_global_idx(E0idx),Enidx,energies->get_energy_from_global_idx(Enidx));
				for (uint jj=0; jj<coeffs.size(); jj++) {
					logmsg->emit_all(LOG_ERROR,"    coeffs[%d]=%.10e",jj,coeffs[jj]);
				}
				NEGF_FEXCEPTION("negative coefficient (%.10e) encountered.", coeffs[ii]);
			}
			check += coeffs[ii];
		}
		NEGF_FASSERT(fabs(check - 1.0) < 1e-12, "check failed: check=%e, Elow=%g, Eupp=%g, Eupp-Elow=%e, indices.size()=%d.",check,Elow,Eupp,Eupp-Elow, indices.size()); 
	}
	return;
);}


void GreenFunctions::assign_flat_band_lesser_greater(double EFn, double EFp, double Ec, double Ev) throw (Exception *)
{STACK_TRACE(
	uint Nk   = kspace->get_number_of_points();
	uint myNE = energies->get_my_number_of_points();
	NEGF_ASSERT(Nn==2, "need 2 DOFs!");
	NEGF_ASSERT(fabs(options->get("kp_method")-3.0)<1e-10, "need kpmethod=3!");
	
	const double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
	const double me   = constants::convert_from_SI(units::mass, 0.067 * constants::SIm0);
	const double mh   = constants::convert_from_SI(units::mass, 0.500 * constants::SIm0);
	const double kT   = constants::convert_from_SI(units::energy, constants::SIec * 0.02585);
	
	vector<double> V; V.resize(Nx, 0.0);
	double V_left = 0.0;
	double V_right = -0.3;
	for (uint xx=0; xx<Nx; xx++) {
		V[xx] = V_left + (V_right-V_left)/(Nx-1) * xx;
	}
	
	uint xx_left  =     Nx / 3; // integer division
	uint xx_right = 2 * Nx / 3; // integer division
	
	for (uint ee2 = 0; ee2 < myNE; ee2++)
	{
		uint ee = energies->get_global_index(ee2);
		double E = energies->get_energy_from_global_idx(ee);
		
		for (uint kk = 0; kk < Nk; kk++)
		{	
			double k = kspace->get_point(kk).get_coord_abs();
			double Eke = hbar*hbar*k*k / (2.0*me);
			double Ekh = hbar*hbar*k*k / (2.0*mh);
			
			GLMat & GLm = this->get_lesser(kk,ee);
			GLMat & GGm = this->get_greater(kk,ee);
			GLm = GLMat_create(NxNn);						// initialize to zero
			GGm = GLMat_create(NxNn);
			
			for (uint xx=xx_left; xx<=xx_right; xx++) 
			{
				double dx = 0.0;
				if (constants::old_orthogonal) {
					dx = 1.0;
				} else {
					const vector<Element *> elems_near_x = xspace->get_elems_near(xspace->get_vertex(xx-1));
					for (uint ii=0; ii<elems_near_x.size(); ii++) {
						dx += elems_near_x[ii]->get_edge(0)->get_length();
					}
				}
				double Ec_plus_V = Ec + V[xx];
				double Ev_plus_V = Ev + V[xx];
				
				double fn = 1.0 / (1.0 + negf_math::exp((E-Ec_plus_V)/kT));	// spatially varying Fermilevel
				double fp = 1.0 / (1.0 + negf_math::exp((E-Ev_plus_V)/kT));
				
				double rho_CB = 0.0;
				if (E > Ec_plus_V+Eke) { rho_CB = 1.0 / negf_math::sqrt(E - (Ec_plus_V+Eke));	}	// spatially varying DOS
				
				double rho_VB = 0.0;
				if (E < Ev_plus_V-Ekh) { rho_VB = 1.0 / negf_math::sqrt((Ev_plus_V-Ekh) - E); }
			
				Matc GLpart(2,2);
				GLpart(1,1) = constants::imag_unit * fn * rho_CB / dx;
				GLpart(2,2) = constants::imag_unit * fp * rho_VB / dx;
				GLm.fill_block(xx,xx, GLpart, Nx);
				Matc GGpart(2,2);
				GGpart(1,1) = -constants::imag_unit * (1.0-fn) * rho_CB / dx;
				GGpart(2,2) = -constants::imag_unit * (1.0-fp) * rho_VB / dx;
				GGm.fill_block(xx,xx, GGpart, Nx);
			}
			
		}
	}
	
	// security check
	if (security_checking) {
		GLMat tmp = GLMat_create(NxNn);
		for (uint ee2=0; ee2<energies->get_my_number_of_points(); ee2++) {
			uint ee = energies->get_global_index(ee2);
			for (uint kk=0; kk < Nk; kk++) {
				const GLMat & GLm = this->get_lesser(kk,ee);
				conjtrans(GLm, tmp);
				tmp += GLm;
				double tmp_norm = negf_math::matrix_norm(tmp);
				NEGF_FASSERT(tmp_norm < constants::antiherm_check, "anti-hermiticity check failed for GL (delta=%e) : kk=%d, ee=%d=%.4g (ee2=%d)", tmp_norm, kk, ee, energies->get_energy_from_global_idx(ee), ee2);
	
				const GLMat & GGm = this->get_greater(kk,ee);
				conjtrans(GGm, tmp);
				tmp += GGm;
				tmp_norm = negf_math::matrix_norm(tmp);
				NEGF_FASSERT(tmp_norm < constants::antiherm_check, "anti-hermiticity check failed for GG (delta=%e) : kk=%d, ee=%d=%.4g (ee2=%d)", tmp_norm, kk, ee, energies->get_energy_from_global_idx(ee), ee2);
			}
		}
	}
	mpi->synchronize_processes();
);}

