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
#include "InnerLoop.h"
using namespace negf;

InnerLoop::InnerLoop(Hamiltonian      * ham_, 
					 Overlap          * ov_, 
					 GreenFunctions   * gf_, 
					 SelfEnergies     * se_, 
					 PostProcessing   * pp_,
					 MaterialDatabase * material_) throw (Exception *):
	max_inner_iterations(constants::max_inner_iters),
	ham(ham_),
	ov(ov_),
	gf(gf_),
	se(se_),
	pp(pp_),
	material(material_),
	xspace(gf->get_xspace()),
	kspace(gf->get_kspace()),
	options(gf->get_options()),
	energies(gf->get_energies())
{STACK_TRACE(
	// default constructors for old_spectral_X --> 1x1-matrices
	logmsg->emit_header("setting up inner loop with Dyson and Keldysh equation");
	NEGF_ASSERT(ham!=NULL && gf!=NULL && se!=NULL && pp!=NULL && xspace!=NULL && kspace!=NULL && options!=NULL && energies !=NULL,
			"null pointer encountered.");
	
	// security checks
	NEGF_ASSERT(xspace==ham->get_xspace() && xspace==se->get_xspace() && xspace==pp->get_xspace(), "inconsistency in real space grid.");
	NEGF_ASSERT(kspace==ham->get_kspace() && kspace==se->get_kspace() && kspace==pp->get_kspace(), "inconsistency in k-space grid.");
	NEGF_ASSERT(options==ham->get_options() && options==se->get_options() && options==pp->get_options(), "inconsistency in used options.");
	NEGF_ASSERT(energies==se->get_energies() && energies==pp->get_energies(), "inconsistency in used energy grid.");
	
	uint   Nx = xspace->get_num_internal_vertices();
	uint myNE = energies->get_my_number_of_points();
	uint   Nk = kspace->get_number_of_points();
	uint   NE = energies->get_number_of_points();
	
	// set up Frey method for contact state broadening, renormalization
	if (options->exists("FreyModel") && options->get("FreyModel")==1) 
	{
		logmsg->emit_big_header("setting up Frey model for contact state broadening");
		
		for (uint cc=0; cc<2; cc++) {
			// ----------------------------------------------------------------------------------------------------------------------------
			// create trivial geometry - I think I need 2 vertices, otherwise there is no element, no region and consequently no material
			// ----------------------------------------------------------------------------------------------------------------------------
			logmsg->emit(LOG_INFO, "Contact %d: creating geometry",cc);
			Geometry * freyspace = new Geometry(2, 1, 0, 1);
			freyspace->set_dimension(1);
			
			double dx = this->se->get_contact_selfenergy()->get_dx(cc);
			freyspace->add_vertex( new Vertex(0, 0.0) ); freyspace->get_vertex(0)->set_index_external(0);
			freyspace->add_vertex( new Vertex(1,  dx) ); freyspace->get_vertex(1)->set_index_external(1);
			
			Edge * edg = new Edge(0,freyspace->get_vertex(0),freyspace->get_vertex(1)); edg->set_index_external(0);
			freyspace->add_edge(edg);
			
			Element * elem = new Element(0, element_type::interval); elem->set_index_external(0);
			elem->add_vertex(freyspace->get_vertex(0));	elem->add_vertex(freyspace->get_vertex(1));
			elem->add_edge(freyspace->get_edge(0));
			freyspace->add_element(elem);
			
			Region * reg = new Region("contact");
			reg->set_material_name(xspace->get_contact(cc)->get_adjacent_region()->get_material()->get_name().c_str());
			reg->set_material     (xspace->get_contact(cc)->get_adjacent_region()->get_material());
			freyspace->get_element(0)->set_region(reg);	
			reg->add_element(freyspace->get_element(0));	
			freyspace->add_region(reg);
			
			freyspace->set_num_dfise_elems(1);
			freyspace->set_num_dfise_regions(1);
			freyspace->prepare();
			freyspace->verify();
			
			mpi->synchronize_processes();
			
			logmsg->emit(LOG_INFO, "Contact %d: creating Hamiltonian", cc);
			// initialize Hamiltonian
			Hamiltonian * freyham = new Hamiltonian(freyspace, kspace, options, material, fnames->get_filename().c_str());
			//ham->set_electrostatic_potential(poiss->get_poisson_equation()->get_values());
			
			logmsg->emit(LOG_INFO, "Contact %d: creating Overlap", cc);
			// initialize overlap matrix
			Overlap * frey_overlap = new Overlap(freyham, freyspace);
			
			logmsg->emit(LOG_INFO, "Contact %d: creating Green's Functions", cc);
			// initialize Green's functions 
			GreenFunctions * frey_gf = new GreenFunctions(options, freyspace, kspace, energies, freyham, frey_overlap, false);
		
			logmsg->emit(LOG_INFO, "Contact %d: creating Self Energies", cc);
			// create self energies using that geometry
			// there is a problem: Nx is globally defined!!!
			SelfEnergies * frey_energies = new SelfEnergies(freyham, frey_overlap, freyspace, kspace, energies, options, frey_gf, material);
			
			if (cc==0) {
				this->left_freyspace      = freyspace;
				this->left_freyham        = freyham;
				this->left_frey_overlap   = frey_overlap;
				this->left_frey_gf        = frey_gf;
				this->left_frey_energies  = frey_energies;
			} else {
				this->right_freyspace     = freyspace;
				this->right_freyham       = freyham;
				this->right_frey_overlap  = frey_overlap;
				this->right_frey_gf       = frey_gf;
				this->right_frey_energies = frey_energies;			
			}
		}
		
		logmsg->emit_big_header("finished setting up Frey model");
	} else {
		this->left_freyspace      = 0;
		this->left_freyham        = 0;
		this->left_frey_overlap   = 0;
		this->left_frey_gf        = 0;
		this->left_frey_energies  = 0;
		this->right_freyspace     = 0;
		this->right_freyham       = 0;
		this->right_frey_overlap  = 0;
		this->right_frey_gf       = 0;
		this->right_frey_energies = 0;	
	}
	
	if (options->exists("FreyModel") && options->get("FreyModel")>=1) {
		if (mpi->get_rank()==0) {
			vector< cplx >         tmp1; tmp1.resize(  Nn, 0.0);
			vector< vector<cplx> > tmp2; tmp2.resize(  Nk, tmp1);
			this->left_total_frey_broadening.resize(NE, tmp2);
			this->right_total_frey_broadening.resize(NE, tmp2);
		}		
	}
	
	
	// set up contact 0 fermilevel
	this->contact_0_fermi = new ContactFermilevel(energies, ham, ov, se->get_contact_selfenergy());
	
	// allocate memory
	this->old_spectral_edensity = Matd(myNE,Nx);
	this->old_spectral_hdensity = Matd(myNE,Nx);
	this->old_spectral_ecurrent = Matd(myNE,Nx-1);
	this->old_spectral_hcurrent = Matd(myNE,Nx-1);
	this->old_edensity.resize(Nx, 0.0);
	this->old_hdensity.resize(Nx, 0.0);
	this->old_ecurrent.resize(Nx-1,0.0);
	this->old_hcurrent.resize(Nx-1,0.0);
	mpi->synchronize_processes();
);}


void InnerLoop::set_max_inner_iterations(uint new_number) 
{
	NEGF_ASSERT(new_number>0 && new_number<1000, "invalid maximum number of inner iterations.");
	logmsg->emit(LOG_INFO,"Maximum number of inner iterations is now %d.", new_number);
	this->max_inner_iterations = new_number;
}


bool InnerLoop::perform(uint outer_iteration, const double err_crit) throw (Exception *)
{STACK_TRACE(
	
	// determine contact broadening by Frey method (if wanted)
	if (options->exists("FreyModel") && options->get("FreyModel")==1) {
		this->perform_frey_iteration();
	}

	char buf[1000];
	const uint max_num_inner_iters_in_first_outer_iter = constants::max_inner_iters_start;
	uint iter;
	int mpi_master = constants::mpi_master_rank; 
	for (iter = 1; iter <= this->max_inner_iterations; iter++)
	{
		logmsg->emit_big_header("Inner iteration %d...", iter);
		
		// perform iteration (calculation of GR,GL,GG and then the new SigmaR, SigmaL, SigmaG) on every process
		this->iterate();
		
		// wait until every process has finished.
		mpi->synchronize_processes();	
		
		// test convergence of own energy points
		int conv = (this->converged(err_crit)) ? 1 : 0; // also test convergence in first loop because density etc is computed
		if (se->is_ballistic()) {
			logmsg->emit(LOG_INFO, "Only one inner iteration is performed because simulation is ballistic.");
			conv = 1;
		} else if (iter==1) {
			logmsg->emit(LOG_INFO, "No inner convergence because first inner loop.");
			conv = 0;
		}
	
		// DEBUG
		bool lumi_debug = false;
		if (se->has_self_energy(SEtype_spont_photon) && lumi_debug) {
			se->get_spontaneous_photon_selfenergy()->output_debug_info();
		}
		bool POPdebug = false;
		if (this->se->has_self_energy(SEtype_optical_phonon) && POPdebug) {
			this->pp->compute_scattering_current(SEtype_optical_phonon);
			se->get_optical_phonon_selfenergy()->output_debug_info();
		}
		
		//pp->check_conservation(SEtype_all);
		if (!this->se->is_ballistic())                            pp->check_conservation(SEtype_contact);
		if (this->se->has_self_energy(SEtype_optical_phonon))     pp->check_conservation(SEtype_optical_phonon);
		if (this->se->has_self_energy(SEtype_acoustic_phonon))    pp->check_conservation(SEtype_acoustic_phonon);
		if (this->se->has_self_energy(SEtype_spont_photon))       pp->check_conservation(SEtype_spont_photon);
		if (this->se->has_self_energy(SEtype_golizadeh_momentum)) pp->check_conservation(SEtype_golizadeh_momentum);
		if (this->se->has_self_energy(SEtype_golizadeh_phase))    pp->check_conservation(SEtype_golizadeh_phase);
		if (this->se->has_self_energy(SEtype_ion_imp))   		  pp->check_conservation(SEtype_ion_imp);

        if (options->exists("Transmission") && options->get("Transmission")==1) {
            // call it by all MPI processes
            pp->compute_transmission();
            mpi->synchronize_processes();
        }

		if (options->exists("QuasiFermilevels") && options->get("QuasiFermilevels")==1) {
		    // call it by all MPI processes
		    pp->compute_quasi_fermilevel();
		    mpi->synchronize_processes();
		}


		if(mpi->get_rank()==mpi_master && !se->is_ballistic()) // in a ballistic simulation there is just one inner loop 
		{		
			const vector<double> & edens = this->get_post_processing()->get_edensity();
			logmsg->emit_noendl(LOG_INFO_L2,"edensity=");
			for (uint ii=0; ii < edens.size(); ii++) {
				logmsg->emit_noendl(LOG_INFO_L2,"%.2e   ",edens[ii]);
			}
			logmsg->emit(LOG_INFO_L2,"");
			
			if (options->get("kp_method")>1.000001) {
				const vector<double> & hdens = this->get_post_processing()->get_hdensity();
				logmsg->emit_noendl(LOG_INFO_L2,"hdensity=");
				for (uint ii=0; ii < hdens.size(); ii++) {
					logmsg->emit_noendl(LOG_INFO_L2,"%.2e   ",hdens[ii]);
				}
				logmsg->emit(LOG_INFO_L2,"");
			}
			
			const vector<double> & ecurrent = this->get_post_processing()->get_ecurrent();
			logmsg->emit_noendl(LOG_INFO_L2,"ecurrent=");
			for (uint ii=0; ii < ecurrent.size(); ii++) {
				logmsg->emit_noendl(LOG_INFO_L2,"%.2e   ",ecurrent[ii]);
			}
			logmsg->emit(LOG_INFO_L2,"");
						
			string output_filename = fnames->get_outfiledirectory() + "inner";
			InputParser parser;
			sprintf(buf,"%s_step%d_nE",output_filename.c_str(),iter);
			parser.write_xE_matrix(buf, this->get_post_processing()->get_entire_spectral_edensity(), xspace, energies->get_energy_grid());
			sprintf(buf,"%s_step%d_pE",output_filename.c_str(),iter);
			parser.write_xE_matrix(buf, this->get_post_processing()->get_entire_spectral_hdensity(), xspace, energies->get_energy_grid());
			sprintf(buf,"%s_step%d_JEe",output_filename.c_str(),iter);
			parser.write_current_matrix(buf, this->get_post_processing()->get_entire_spectral_ecurrent(), xspace, energies->get_energy_grid());
			//sprintf(buf,"%s_step%d_JEe2",output_filename.c_str(),iter);
			//parser.write_xE_matrix(buf, this->get_post_processing()->get_entire_spectral_ecurrent2(), xspace, energies->get_energy_grid());
			sprintf(buf,"%s_step%d_JEh",output_filename.c_str(),iter);
			parser.write_current_matrix(buf, this->get_post_processing()->get_entire_spectral_hcurrent(), xspace, energies->get_energy_grid());
			//sprintf(buf,"%s_step%d_JEh2",output_filename.c_str(),iter);
			//parser.write_xE_matrix(buf, this->get_post_processing()->get_entire_spectral_hcurrent2(), xspace, energies->get_energy_grid());
			sprintf(buf,"%s_step%d_LDOS",output_filename.c_str(),iter);
			parser.write_xE_matrix(buf, this->get_post_processing()->get_entire_local_dos(), xspace, energies->get_energy_grid());
			//sprintf(buf,"%s_step%d_LDOS_k0",output_filename.c_str(),iter);
			//parser.write_xE_matrix(buf, this->get_post_processing()->get_entire_local_dos_k0(), xspace, energies->get_energy_grid());
			//if (fabs(Nn-2.0) < 1e-10) {
			//	sprintf(buf,"%s_step%d_LDOS_VB",output_filename.c_str(),iter);
			//	parser.write_xE_matrix(buf, this->get_post_processing()->get_entire_local_dos_VB(), xspace, energies->get_energy_grid());
			//}
			if (this->se->has_self_energy(SEtype_optical_phonon) && POPdebug) {
				sprintf(buf,"%s_step%d_Jephon",output_filename.c_str(),iter);
				parser.write_xE_matrix(buf, pp->get_entire_scattering_ecurrent(SEtype_optical_phonon), xspace, energies->get_energy_grid());
				sprintf(buf,"%s_step%d_Jhphon",output_filename.c_str(),iter);
				parser.write_xE_matrix(buf, pp->get_entire_scattering_hcurrent(SEtype_optical_phonon), xspace, energies->get_energy_grid());
			}
		}
		
		// "conv" of master process shall be minimum of "conv" of all processes (0 or 1)
		// note that when the convergence of total current and densities is tested (which only the master thread does),
		// the other threads return true
		int tag = 777;
		if(mpi->get_rank()==mpi_master) {
			for (int cpu = 0; cpu < mpi->get_num_procs(); cpu++) {
				if (cpu==mpi_master) continue;
				int conv_part;
				mpi->recv(conv_part, cpu, tag);
				NEGF_ASSERT(conv_part==0 || conv_part==1, "unexpected value received.");
				conv = min(conv,conv_part);
			}
		} else {
			mpi->send(conv,mpi_master,tag);
		}
		
		mpi->broadcast(conv, mpi_master);	// performed by all processes, conv is send for mpi_master, recv for rest
		if (conv==1) {
			// all processes are converged
			logmsg->emit(LOG_INFO, "Inner loop converged in %d iterations!",iter);
			break;
		}
		if (outer_iteration==1 && iter==min(this->max_inner_iterations, max_num_inner_iters_in_first_outer_iter)) {
			logmsg->emit(LOG_WARN,"%d inner iterations reached without convergence.", min(this->max_inner_iterations, max_num_inner_iters_in_first_outer_iter));
			logmsg->emit(LOG_WARN,"We continue because it is the first outer iteration.");
			break;
		}
	}
	if (outer_iteration==1) {
		return true;
	}
	if (iter==max_inner_iterations+1) {
		//NEGF_FEXCEPTION("Maximum number of inner iterations (%d) reached without convergence.",max_inner_iterations);
		logmsg->emit(LOG_WARN,"Maximum number of inner iterations (%d) reached without convergence.",max_inner_iterations);
		return false;
	} else {
		logmsg->emit(LOG_INFO, "Inner Convergence reached in %d iterations", iter);
		return true;
	}
);}


void InnerLoop::initial_guess() throw(Exception)
{STACK_TRACE(
	logmsg->emit_header("Initial guess for retarded and advanced Green functions");
	gf->calculate_retarded(se);
	gf->calculate_advanced();
);}


void InnerLoop::iterate()
{STACK_TRACE(
	if (options->exists("FreyModel") && options->get("FreyModel")==2) {
		this->broaden_contact_states();
	}
	
	// ----------------------------
	// calculate self energies
	// ----------------------------
	// note: in the first iteration, all but the contact self-energies will be zero
	// because they depend on the GF (which are zero in the beginning).
	se->calculate();
	
	// ----------------------------
	// some checks, info
	// ----------------------------
	bool debug = false;
	if (debug) {
		if (this->se->has_self_energy(SEtype_optical_phonon))     pp->check_conservation(SEtype_optical_phonon);
		if (this->se->has_self_energy(SEtype_acoustic_phonon))    pp->check_conservation(SEtype_acoustic_phonon);
		if (this->se->has_self_energy(SEtype_spont_photon))       pp->check_conservation(SEtype_spont_photon);
		if (this->se->has_self_energy(SEtype_golizadeh_momentum)) pp->check_conservation(SEtype_golizadeh_momentum);
		if (this->se->has_self_energy(SEtype_golizadeh_phase))    pp->check_conservation(SEtype_golizadeh_phase);
		if (this->se->has_self_energy(SEtype_ion_imp))   		  pp->check_conservation(SEtype_ion_imp);
	}
	//pp->compute_spectral_current();	// computes both GH-GH and SLGG-SGGL; also recomputes contact current
	
	// ----------------------------
	// calculate Green functions
	// ----------------------------
	gf->calculate_retarded(se);
	gf->calculate_advanced();
	
	if (gf->mangle_overlap_calculation()) {
		gf->calculate_lesser_greater(se);
	} else {
		gf->calculate_lesser(se);
		gf->calculate_greater();
		logmsg->emit(LOG_INFO,"Nk=%d", kspace->get_number_of_points());
		if (kspace->get_number_of_points()!=1) { // we don't need MGLM and MGGM in that case --> save time
			gf->calculate_overlap_augmented_lesser();
			gf->calculate_overlap_augmented_greater();
		}
	}
);}


/** our convergence criterion: spectral density and current do not change anymore */
bool InnerLoop::converged(const double err_crit)
{STACK_TRACE(
	logmsg->emit_header("testing convergence of inner loop");
	
	uint   Nx = xspace->get_num_internal_vertices();
	uint NxNn = Nx*Nn;
	uint myNE = energies->get_my_number_of_points();
	
	pp->compute_local_dos();
	pp->compute_spectral_edensity();
	pp->compute_spectral_hdensity();
	pp->compute_spectral_current();	// computes both GH-GH and SLGG-SGGL; also recomputes contact current
	mpi->synchronize_processes();
	pp->compute_edensity();
	pp->compute_hdensity();
	pp->compute_current();
	mpi->synchronize_processes();
	
	// check
	if (mpi->get_rank()==constants::mpi_master_rank) {
		Matd n_plus_p(NxNn,NxNn);
		n_plus_p  = pp->get_entire_spectral_edensity();
		n_plus_p += pp->get_entire_spectral_hdensity();
		
		Matd LDOS_2pi(NxNn,NxNn);
		mult(pp->get_entire_local_dos(), 2.0*constants::pi, LDOS_2pi); // LDOS_2pi = 2.0*constants::pi * pp->get_entire_local_dos();
		n_plus_p -= LDOS_2pi;
		double diff_norm = negf_math::matrix_norm(n_plus_p);
		
		double LDOS_norm = negf_math::matrix_norm(LDOS_2pi);
		
		if (diff_norm / LDOS_norm >= 1e-3) {
			//NEGF_FEXCEPTION("n+p did not correspond to 2pi*LDOS! |n+p-LDOS|=%e, |LDOS|=%e",diff_norm,LDOS_norm);
			logmsg->emit(LOG_INFO,"n+p did not correspond to 2pi*LDOS! |n+p-LDOS|=%e, |LDOS|=%e",diff_norm,LDOS_norm);
		}
	}
		
	// switch to test spectrally resolved quantities or just the energy-integrated densities/current
	bool test_matrix_norms = false;
	if (test_matrix_norms) 
	{
		Matd spectral_ecurrent_difference(myNE,Nx-1);
		Matd spectral_hcurrent_difference(myNE,Nx-1);
		Matd spectral_edensity_difference(myNE,Nx);
		Matd spectral_hdensity_difference(myNE,Nx);
		sub(this->old_spectral_ecurrent, pp->get_spectral_ecurrent(), spectral_ecurrent_difference); // sub(A,B,C) means C=A-B
		sub(this->old_spectral_hcurrent, pp->get_spectral_hcurrent(), spectral_hcurrent_difference); 
		sub(this->old_spectral_edensity, pp->get_spectral_edensity(), spectral_edensity_difference); 
		sub(this->old_spectral_hdensity, pp->get_spectral_hdensity(), spectral_hdensity_difference); 
		this->old_spectral_ecurrent = pp->get_spectral_ecurrent();
		this->old_spectral_hcurrent = pp->get_spectral_hcurrent();
		this->old_spectral_edensity = pp->get_spectral_edensity();
		this->old_spectral_hdensity = pp->get_spectral_hdensity();
		
		// get the Frobenius norm |A|=sqrt(sum_ij |A_ij|^2) of the spectral current and densities
		// this norm scales with the matrix dimensions: |A|~sqrt(m*n) if A is an n*m-matrix
		double  Jenorm = negf_math::matrix_norm(pp->get_spectral_ecurrent());
		double dJenorm = negf_math::matrix_norm(spectral_ecurrent_difference);
		double  Jhnorm = negf_math::matrix_norm(pp->get_spectral_hcurrent());
		double dJhnorm = negf_math::matrix_norm(spectral_hcurrent_difference);
		double   Nnorm = negf_math::matrix_norm(pp->get_spectral_edensity());
		double  dNnorm = negf_math::matrix_norm(spectral_edensity_difference);
		double   Pnorm = negf_math::matrix_norm(pp->get_spectral_hdensity());
		double  dPnorm = negf_math::matrix_norm(spectral_hdensity_difference);
		
		NEGF_FASSERT(!isnan(Jenorm) && !isnan(dJenorm), "p%d: NaN Jenorm or dJenorm!",mpi->get_rank());
		NEGF_FASSERT(!isnan(Jhnorm) && !isnan(dJhnorm), "p%d: NaN Jhnorm or dJhnorm!",mpi->get_rank());
		NEGF_FASSERT(!isnan(Nnorm)  && !isnan(dNnorm),  "p%d: NaN  Nnorm or  dNnorm!",mpi->get_rank());
		NEGF_FASSERT(!isnan(Pnorm)  && !isnan(dPnorm),  "p%d: NaN  Pnorm or  dPnorm!",mpi->get_rank());
		
		// an entry in J(x,E) has the units charge/(energy*time*area)
		// --> typical quantity is 1e2 A cm^-2 / 10 meV
		double Jref = constants::convert_from_SI(units::electrical_current_density_3d, constants::refcurr)
					/ constants::convert_from_SI(units::energy, 0.01*constants::SIec);
		//          ~ 1e6*1e-11 / 0.01 = 1e-3 for ec-nm-ps
		double Jnormref = negf_math::sqrt(Nx*myNE) * Jref;
		// an entry in n(x,E) has the unit density/energy --> typical quantity is 1e16 cm^-3 / 10meV
		double Nref = constants::convert_from_SI(units::density_3d, constants::refdens)
					/ constants::convert_from_SI(units::energy, 0.01*constants::SIec);
		//          ~ 1e22*1e-27 / 0.01 = 1e-3 for ec-nm-ps
		double Nnormref = negf_math::sqrt(Nx*myNE) * Nref;
		
		logmsg->emit_all(LOG_INFO_L2,"p%d: |Je(x,E)|=%.3e, rel_change=%.3e; |n(x,E)|=%.3e, rel_change=%.3e; |p(x,E)|=%.3e, rel_change=%.3e",
				mpi->get_rank(), Jenorm, dJenorm, Nnorm, dNnorm, Pnorm, dPnorm);
				
		if (   dJenorm / (Jenorm + Jnormref) < err_crit 
			&& dJhnorm / (Jhnorm + Jnormref) < err_crit 
			&&  dNnorm / ( Nnorm + Nnormref) < err_crit 
			&&  dPnorm / ( Pnorm + Nnormref) < err_crit) {
			logmsg->emit_all(LOG_INFO_L2, "Inner convergence achieved in process %d! eps=%.2e",mpi->get_rank(), err_crit);
			logmsg->emit_all(LOG_INFO_L2, "    |n(x,E)|=%.2e, rel_change=%.2e; |p(x,E)|=%.2e, rel_change=%.2e;",
					Nnorm, dNnorm/(Nnorm+Nnormref), Pnorm, dPnorm/(Pnorm+Nnormref));
			logmsg->emit_all(LOG_INFO_L2, "    |Je(x,E)|=%.2e, rel_change=%.2e; |Jh(x,E)|=%.2e, rel_change=%.2e;",
					Jenorm, dJenorm / (Jenorm+Jnormref), Jhnorm, dJhnorm / (Jhnorm+Jnormref));
			return true;
		} else {
			logmsg->emit_all(LOG_INFO, "Inner convergence failed in process %d: eps=%.2e",mpi->get_rank(), err_crit);
			logmsg->emit_all(LOG_INFO_L2, "    |n(x,E)|=%.2e, rel_change=%.2e; |p(x,E)|=%.2e, rel_change=%.2e;",
					Nnorm, dNnorm/(Nnorm+Nnormref), Pnorm, dPnorm/(Pnorm+Nnormref));
			logmsg->emit_all(LOG_INFO_L2, "     |Je(x,E)|=%.2e, rel_change=%.2e; |Jh(x,E)|=%.2e, rel_change=%.2e;",
					Jenorm, dJenorm / (Jenorm+Jnormref), Jhnorm, dJhnorm / (Jhnorm+Jnormref));
			return false;
		}
	} else 
	{
		if (mpi->get_rank()==constants::mpi_master_rank)
		{
			vector<double> edensity_change = old_edensity;
			vector<double> hdensity_change = old_hdensity;
			vector<double> ecurrent_change = old_ecurrent;
			vector<double> hcurrent_change = old_hcurrent;
			for (uint ii=0; ii<Nx; ii++) {
				edensity_change[ii] -= (pp->get_edensity())[ii];
				hdensity_change[ii] -= (pp->get_hdensity())[ii];
				if (ii<Nx-1) {
					ecurrent_change[ii] -= (pp->get_ecurrent())[ii];
					hcurrent_change[ii] -= (pp->get_hcurrent())[ii];
				}
			}
			this->old_edensity = pp->get_edensity();
			this->old_hdensity = pp->get_hdensity();
			this->old_ecurrent = pp->get_ecurrent();
			this->old_hcurrent = pp->get_hcurrent();
			
			double  Jenorm = negf_math::vector_norm(pp->get_ecurrent());
			double dJenorm = negf_math::vector_norm(ecurrent_change);
			double  Jhnorm = negf_math::vector_norm(pp->get_hcurrent());
			double dJhnorm = negf_math::vector_norm(hcurrent_change);
			double   Nnorm = negf_math::vector_norm(pp->get_edensity());
			double  dNnorm = negf_math::vector_norm(edensity_change);
			double   Pnorm = negf_math::vector_norm(pp->get_hdensity());
			double  dPnorm = negf_math::vector_norm(hdensity_change);
			
			NEGF_FASSERT(!isnan(Jenorm) && !isnan(dJenorm), "p%d: NaN Jenorm or dJenorm!",mpi->get_rank());
			NEGF_FASSERT(!isnan(Jhnorm) && !isnan(dJhnorm), "p%d: NaN Jhnorm or dJhnorm!",mpi->get_rank());
			NEGF_FASSERT(!isnan( Nnorm) && !isnan( dNnorm), "p%d: NaN  Nnorm or  dNnorm!",mpi->get_rank());
			NEGF_FASSERT(!isnan( Pnorm) && !isnan( dPnorm), "p%d: NaN  Pnorm or  dPnorm!",mpi->get_rank());
			
			double Jref = constants::convert_from_SI(units::electrical_current_density_3d, constants::refcurr);
			//          ~ 1e6*1e-11 = 1e-5 for ec-nm-ps
			double Nref = constants::convert_from_SI(units::density_3d, constants::refdens);
			//          ~ 1e22*1e-27 = 1e-5 for ec-nm-ps
			
			logmsg->emit_all(LOG_INFO_L2,"|Je(x)|=%.2e, change=%.2e; |Jh(x)=%.2e, change=%.2e; |n(x)|=%.2e, change=%.2e; |p(x)|=%.2e, change=%.2e",
					Jenorm, dJenorm, Jhnorm, dJhnorm, Nnorm, dNnorm, Pnorm, dPnorm);
					
			if (   dJenorm / (Jenorm + (Nx-1)*Jref) < err_crit 
				&& dJhnorm / (Jhnorm + (Nx-1)*Jref) < err_crit 
				&&  dNnorm / ( Nnorm + Nx    *Nref) < err_crit 
				&&  dPnorm / ( Pnorm + Nx    *Nref) < err_crit) 
			{
				logmsg->emit(LOG_INFO, "Inner convergence achieved! eps=%.2e", err_crit);
				logmsg->emit(LOG_INFO, "      |n(x)|=%.2e, rel_change=%.2e;  |p(x)|=%.2e, rel_change=%.2e;",
					Nnorm, dNnorm/(Nnorm+Nx*Nref), Pnorm, dPnorm/(Pnorm+Nx*Nref));
				logmsg->emit(LOG_INFO, "     |Je(x)|=%.2e, rel_change=%.2e; |Jh(x)|=%.2e, rel_change=%.2e;",
					Jenorm, dJenorm / (Jenorm+(Nx-1)*Jref), Jhnorm, dJhnorm / (Jhnorm+(Nx-1)*Jref));
				return true;
			} else {
				logmsg->emit(LOG_INFO, "Inner convergence failed: eps=%.1e", err_crit);
				logmsg->emit(LOG_INFO, "      |n(x)|=%.2e, rel_change=%.2e;  |p(x)|=%.2e, rel_change=%.2e.", 
					Nnorm, dNnorm/(Nnorm+Nx*Nref), Pnorm, dPnorm/(Pnorm+Nx*Nref));
				logmsg->emit(LOG_INFO, "     |Je(x)|=%.2e, rel_change=%.2e; |Jh(x)|=%.2e, rel_change=%.2e.",
					Jenorm, dJenorm / (Jenorm+(Nx-1)*Jref), Jhnorm, dJhnorm / (Jhnorm+(Nx-1)*Jref));
				return false;
			}
		} else {
			return true;
		}
	}
);}


// execute this by all threads
void InnerLoop::perform_frey_iteration()
{STACK_TRACE(
	logmsg->emit_header("Performing Frey iteration");
	NEGF_ASSERT(options->exists("FreyModel") && options->get("FreyModel")==1 && left_freyham!=0 && right_freyham!=0, "Cannot do Frey iteration.");
	NEGF_ASSERT(options->get("kp_method") <= 3.00001, "Frey model is at the moment only applicable to effective mass band structure!"); 
	
	// assumptions electrostatic potential has already been assigned to left and right contact Hamiltonians
	const uint max_freyiter = 100;	// max. # iterations between GF approximation and SE calculation
	const uint freycrit     = 1e-5; // convercence criterion in eV for each (k,E)-point in the diagnoal of the contact scattering self-energy
	
	const uint     NE = energies->get_number_of_points();
	const uint   myNE = energies->get_my_number_of_points();
	const uint     Nk = kspace->get_number_of_points();
	const double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
	const double   ec = constants::convert_from_SI(units::charge, constants::SIec);
	
	 // set up array for diagonal part of retarded contact scattering self-energy belonging to this process
	logmsg->emit(LOG_INFO,"setting up scattering self-energy array");
	vector<cplx>                                          tmp1;            tmp1.resize(   2, 0.0);
	vector< vector<cplx> >                                tmp2;            tmp2.resize(  Nn, tmp1);
	vector< vector< vector<cplx> > >                      tmp3;            tmp3.resize(  Nk, tmp2);
	vector< vector< vector< vector<cplx> > > > frey_broadening; frey_broadening.resize(myNE, tmp3);
	
	for (uint cc=0; cc<2; cc++) 	// loop over both contacts
	{
		logmsg->emit(LOG_INFO,"CONTACT %d:",cc);
		Geometry       * freyspace     = (cc==0) ? this->left_freyspace     : this->right_freyspace;
		Hamiltonian    * freyham       = (cc==0) ? this->left_freyham       : this->right_freyham;
		//Overlap        * frey_overlap  = (cc==0) ? this->left_frey_overlap  : this->right_frey_overlap;
		GreenFunctions * frey_gf       = (cc==0) ? this->left_frey_gf       : this->right_frey_gf;
		SelfEnergies   * frey_energies = (cc==0) ? this->left_frey_energies : this->right_frey_energies;
		vector< vector< vector<cplx> > > & total_frey_broadening = (cc==0) ? this->left_total_frey_broadening : this->right_total_frey_broadening; // only used in master process
		vector<double> masses; masses.resize(Nn, 0.0);
		vector<double> bandedges; bandedges.resize(Nn, 0.0);
		const PropertyContainer<double> * mat =  freyspace->get_region(0)->get_material();
		masses[0] = constants::convert_from_SI(units::mass, constants::SIm0 * mat->get("electron_effective_mass"));
		bandedges[0] = TdkpInfoDesk::get_cbedge(mat, options->get("temperature"), material);
		if (Nn==2) {
		masses[1] = constants::convert_from_SI(units::mass, constants::SIm0 * mat->get("hole_effective_mass"));
		bandedges[1] = mat->get("valence_band_edge");
		}
		
		// ---------------------------------------------------------
		// assign correct electrostatic potential to Hamiltonian
		// ---------------------------------------------------------
		const vector<double> & pot = this->ham->get_electrostatic_potential();
		NEGF_ASSERT(pot.size()==this->xspace->get_num_vertices(), "inconsistent num. vertices");
		double phi = pot[this->xspace->get_contact(cc)->get_contact_vertex(0)->get_index_global()];
		vector<double> contact_potential(2, phi);
		freyham->set_electrostatic_potential(contact_potential);
		
		// ---------------------------------------------------------
		// do Gummel iteration!
		// ---------------------------------------------------------
		uint iter = 0;
		for (iter=0; iter<max_freyiter; iter++) 
		{
			logmsg->emit(LOG_INFO,"*** iteration %d ***",iter);
			logmsg->emit(LOG_INFO,"---------> setting up GF");
			// ---------------------------------------------------------------------
			// set up Green's functions
			// ---------------------------------------------------------------------
			// RETARDED, ADVANCED AND LESSER
			for (uint ee2=0; ee2<myNE; ee2++) {
				uint ee = energies->get_global_index(ee2);
				double E = energies->get_energy_from_global_idx(ee);
				for (uint kk=0; kk<Nk; kk++) {
					// --------------------------------------------------
					// find own longitudinal k-vectors, band dependent
					// --------------------------------------------------
					vector<cplx> kvectors; kvectors.resize(Nn, 0.0);
					NEGF_ASSERT(options->get("kp_method")<=3+1.0e-8, "Right now the kx calculation only works w/ effective mass.");
					double kval = kspace->get_point(kk).get_coord_abs();
					for (uint nn=0; nn<Nn; nn++) {
						cplx SigmaR = frey_broadening[ee2][kk][nn][cc];
						if (nn==0) {
							double Ek = hbar*hbar * kval*kval / (2.0 * masses[0]);
							double Ec = bandedges[0];
							kvectors[0] = 1.0/hbar * sqrt(2.0*masses[0]*(E - Ek - (Ec-ec*phi) - SigmaR));
						} else {
							double Ek = hbar*hbar * kval*kval / (2.0 * masses[1]);
							double Ev = bandedges[1];
							kvectors[1] = 1.0/hbar * sqrt(2.0*masses[1]*(E - Ek - (Ev-ec*phi) - SigmaR));
						}
					}
					
					// --------------------------------------------------
					// find own Green's function approximation
					// --------------------------------------------------
					// RETARDED
					Matc & GRmat = frey_gf->get_retarded(kk,ee);
					NEGF_ASSERT(GRmat.num_rows()==2*Nn && GRmat.num_cols()==2*Nn, "strange matrix.");
					GRmat = Matc(2*Nn,2*Nn);
					for (uint nn=1; nn<=Nn; nn++) {
						cplx   kn = kvectors[nn-1];
						double mn = masses[nn-1];
						for (uint xx=1; xx<=2; xx++) {
							for (uint yy=1; yy<=2; yy++) {
								uint idx1 = get_mat_idx(xx,nn,2); // get_mat_idx(xx,nn,Nx), xx and nn 1-based
								uint idx2 = get_mat_idx(yy,nn,2); 
								double x = freyspace->get_vertex(xx-1)->get_coordinate(0);
								double y = freyspace->get_vertex(yy-1)->get_coordinate(0);
								
								// travelling GF
								GRmat(idx1,idx2) = -constants::imag_unit * mn / (hbar*hbar * kn) * exp(-constants::imag_unit * kn * abs(x-y));
								if (kk==10) {
									cout << "GRmat[ee=" << ee << "][kk=10][cc=" << cc << "](" << idx1 << ", " << idx2 << ")=" << GRmat(idx1,idx2) << endl;
								}
								NEGF_FASSERT(!isnan(GRmat(idx1,idx2).real()) && !isnan(GRmat(idx1,idx2).imag()), 
										"encountered a NaN in GR[ee=%d][kk=%d](%d,%d).", ee, kk, idx1, idx2);
							}
						}
					}
			
				}
			}
			// ADVANCED
			frey_gf->calculate_advanced();
			// LESSER by fluctuation-dissipation GL = -f*(GR-GA)
			Matc tmp(2*Nn,2*Nn);
			for (uint ee2=0; ee2<myNE; ee2++) {
				uint ee = energies->get_global_index(ee2);
				for (uint kk=0; kk<Nk; kk++) {
					Matc & GRmat = frey_gf->get_retarded(kk,ee);
					Matc & GAmat = frey_gf->get_advanced(kk,ee);
					tmp  = GRmat;
					tmp -= GAmat;
					BMatc & GLmat = frey_gf->get_lesser(kk,ee);
					double f = se->get_contact_selfenergy()->get_fermi(cc, ee); // 1/(1+exp(...))
					mult(tmp, -f, GLmat);
				}
			}			
			// GREATER
			frey_gf->calculate_greater();
			// OVERLAP-AUGMENTED
			frey_gf->calculate_overlap_augmented_lesser();
			frey_gf->calculate_overlap_augmented_greater();
			mpi->synchronize_processes();
		
			// ---------------------------------------------------------------------
			// calculate self energies
			// ---------------------------------------------------------------------
			logmsg->emit(LOG_INFO,"---------> calculating self-energies");
			frey_energies->calculate();
			mpi->synchronize_processes();
			
			// ---------------------------------------------------------------------
			// extract diagonal part of the retarded total scattering self-energy
			// ---------------------------------------------------------------------
			logmsg->emit(LOG_INFO,"---------> extracting diagonal part of retarded self-energy");
			for (uint ee2=0; ee2<myNE; ee2++) {
				uint ee = energies->get_global_index(ee2);
				for (uint kk=0; kk<Nk; kk++) {
					const BMatc & SRmat = frey_energies->get_retarded(kk,ee);
					NEGF_ASSERT(SRmat.num_rows()==2*Nn && SRmat.num_cols()==2*Nn, "strange matrix.");
					for (uint nn=0; nn<Nn; nn++) {
						cplx val1 = SRmat(get_mat_idx(1,nn+1,2),get_mat_idx(1,nn+1,2));	// get_mat_idx(xx,nn,Nx), xx and nn 1-based
						cplx val2 = SRmat(get_mat_idx(2,nn+1,2),get_mat_idx(2,nn+1,2));
						NEGF_ASSERT(!isnan(val1.real()) && !isnan(val2.real()) && !isnan(val1.imag()) && !isnan(val2.imag()), "encountered a NaN.");
						if (abs(val1-val2)/(abs(val1)+abs(val2)+1e-10) > 1e-5) {
							NEGF_FEXCEPTION("Dunno what to take: val1=(%e,%e), val2=(%e,%e)", val1.real(), val1.imag(), val2.real(), val2.imag());
						}
						frey_broadening[ee2][kk][nn][cc] = 0.5*(val1+val2);
						// WHAT ABOUT SPACE DISCRETIZATION!!!!! SE HAS UNITS ENERGY*SPACE!
					}
				}
			}
			
			// ---------------------------------------------------------------------
			// communicate diagonal part to master process
			// master process tests convergence
			// ---------------------------------------------------------------------
			logmsg->emit(LOG_INFO,"---------> MPI communication");
			bool frey_converged = true;
			if (mpi->get_rank()==0) {
				// copy old values to separate array
				vector< vector< vector<cplx> > > old_total_frey_broadening = total_frey_broadening;
				
				// own part
				for (uint ee2=0; ee2<myNE; ee2++) {
					uint ee = energies->get_global_index(ee2);
					for (uint kk=0; kk<Nk; kk++) {
						for (uint nn=0; nn<Nn; nn++) {
							total_frey_broadening[ee][kk][nn] = frey_broadening[ee2][kk][nn][cc];
						}
					}
				}
				
				// part from other processes
				uint array_size = Nn*Nk;
				vector<double> tmp_real; tmp_real.resize(array_size, 0.0); // storage: nn*Nn+kk
				vector<double> tmp_imag; tmp_imag.resize(array_size, 0.0); 
				for (int proc=0; proc<mpi->get_num_procs(); proc++) {
					if (proc==constants::mpi_master_rank) continue;
					
					uint start_idx = energies->get_start_global_idx(proc);
					uint stop_idx  = energies->get_stop_global_idx(proc);
				
					int tag = proc;
					for (uint ee=start_idx; ee<=stop_idx; ee++) {
						// receive from other process
						mpi->recv(tmp_real, array_size, proc, tag);
						mpi->recv(tmp_imag, array_size, proc, tag);
						// add to total matrix
						for (uint kk=0; kk<Nk; kk++) {
							for (uint nn=0; nn<Nn; nn++) {
								total_frey_broadening[ee][kk][nn] = tmp_real[nn*Nn+kk] + constants::imag_unit * tmp_imag[nn*Nn+kk];
							}
						}
					}
				}
				
				// TEST CONVERGENCE - relative change of each entry must not exceed freycrit
				logmsg->emit(LOG_INFO," master process tests convergence");
				for (uint ee=0; ee<NE; ee++) {
					for (uint kk=0; kk<Nk; kk++) {
						for (uint nn=0; nn<Nn; nn++) {
							if (abs(old_total_frey_broadening[ee][kk][nn] - total_frey_broadening[ee][kk][nn]) > freycrit) {
								frey_converged = false;
								break;
							}
						}
					}
					if (!frey_converged) break;
				}
				
				// some screen output
				//if (frey_converged) {
					logmsg->emit(LOG_INFO,"Displaying n=0 self-energy:");
					NEGF_ASSERT(total_frey_broadening.size()==NE, "something is wrong (NE).");
					for (uint ee=0; ee<NE; ee++) {
						NEGF_ASSERT(total_frey_broadening[ee].size()==Nk, "something is wrong (Nk).");
						logmsg->emit(LOG_INFO,"E=%7.3g:  ",energies->get_energy_from_global_idx(ee));
						for (uint kk=0; kk<Nk; kk++) {
							NEGF_ASSERT(total_frey_broadening[ee][kk].size()>=1, "something is wrong (Nn).");
							logmsg->emit_noendl(LOG_INFO,"(%7.3e,%7.3e)   ", total_frey_broadening[ee][kk][0].real(), total_frey_broadening[ee][kk][0].imag());
						}
						logmsg->emit(LOG_INFO,"");
					}
				//}
					
				// communicate convergence
				logmsg->emit(LOG_INFO," master process communicates convergence");
				int root = 0;
				mpi->broadcast(frey_converged, root);
				
			} else {
				// send own part to master process
				int dest = constants::mpi_master_rank;
				int tag = mpi->get_rank();
				vector<double> tmp_real; tmp_real.resize(Nn*Nk, 0.0); // storage: nn*Nn+kk
				vector<double> tmp_imag; tmp_imag.resize(Nn*Nk, 0.0); 
				for (uint ee2=0; ee2<myNE; ee2++) {
					for (uint kk=0; kk<Nk; kk++) {
						for (uint nn=0; nn<Nn; nn++) {
							tmp_real[nn*Nn+kk] = frey_broadening[ee2][kk][nn][cc].real();
							tmp_imag[nn*Nn+kk] = frey_broadening[ee2][kk][nn][cc].imag();
						}
					}
					mpi->send(tmp_real, dest, tag);
					mpi->send(tmp_imag, dest, tag);
				}
				
				// get boolean if converged or not
				int root = 0;
				mpi->broadcast(frey_converged, root);
			}
			
			if (frey_converged) {
				logmsg->emit(LOG_INFO,"Convergence of self-energy iteration in contact %d achieved after %d iterations!", cc, iter+1);
				break;
			}
		}
		if (iter>=max_freyiter-1) {
			NEGF_FEXCEPTION("No convergence of self-energy iteration in contact %d achieved after %d iterations - aborting.", cc, max_freyiter);
		}
	}	
			
	// assign numbers to contact self-energy
	logmsg->emit(LOG_INFO,"assigning k- and E-dependent broadening to contact self-energies");
	this->se->get_contact_selfenergy()->assign_broadening(frey_broadening);
);}


void InnerLoop::broaden_contact_states()
{STACK_TRACE(
	const uint     NE = energies->get_number_of_points();
	const uint   myNE = energies->get_my_number_of_points();
	const uint     Nk = kspace->get_number_of_points();
	const uint     Nx = xspace->get_num_internal_vertices();
	
	 // set up array for diagonal part of retarded contact scattering self-energy belonging to this process
	logmsg->emit(LOG_INFO,"setting up scattering self-energy array");
	vector<cplx>                                          tmp1;            tmp1.resize(   2, 0.0);
	vector< vector<cplx> >                                tmp2;            tmp2.resize(  Nn, tmp1);
	vector< vector< vector<cplx> > >                      tmp3;            tmp3.resize(  Nk, tmp2);
	vector< vector< vector< vector<cplx> > > > frey_broadening; frey_broadening.resize(myNE, tmp3);
	
	for (uint cc=0; cc<2; cc++) 	// loop over both contacts
	{
		// ---------------------------------------------------------------------
		// get internal vertex number of vertex at which SR is taken
		// NOT first/last vertex because there the contact SE is added
		// ---------------------------------------------------------------------
		//uint vidx = (cc==0) ? 5 : Nx-4;
		uint vidx = (cc==0) ? 4 : Nx-3;
		logmsg->emit(LOG_INFO,"Taking vertex at x=%.3g to be near contact", xspace->get_vertex(xspace->get_global_vertex_index(vidx-1))->get_coordinate(0));
		
		// ---------------------------------------------------------------------
		// extract diagonal part of the retarded total scattering self-energy
		// and add to frey_broadening array
		// ---------------------------------------------------------------------
		logmsg->emit(LOG_INFO,"---------> extracting diagonal part of retarded self-energy");
		for (uint ee2=0; ee2<myNE; ee2++) {
			uint ee = energies->get_global_index(ee2);
			for (uint kk=0; kk<Nk; kk++) {
				const BMatc & SRmat  = se->get_retarded(kk,ee);
				const BMatc & SRcont = se->get_contact_selfenergy()->get_retarded(kk,ee);
				// ATTENTION: SRmat has units energy*length!
				// we assume here that the grid spacing near the second vertex is the same as at the contact-device interface
				for (uint nn=0; nn<Nn; nn++) {
					uint idx = get_mat_idx(vidx,nn+1,Nx);
					frey_broadening[ee2][kk][nn][cc] = (-5.0) * (SRmat(idx,idx) - SRcont(idx,idx)); // MINUS SIGN
					// get_mat_idx(xx,nn,Nx), xx and nn 1-based
				}
			}
		}
	
		// ---------------------------------------------------------------------
		// communicate diagonal part to master process
		// master process tests convergence
		// ---------------------------------------------------------------------
		logmsg->emit(LOG_INFO,"---------> MPI communication");
		vector< vector< vector<cplx> > > & total_frey_broadening = (cc==0) ? this->left_total_frey_broadening : this->right_total_frey_broadening; // only used in master process
		
		if (mpi->get_rank()==0) {
			// copy old values to separate array
			vector< vector< vector<cplx> > > old_total_frey_broadening = total_frey_broadening;
			
			// own part
			for (uint ee2=0; ee2<myNE; ee2++) {
				uint ee = energies->get_global_index(ee2);
				for (uint kk=0; kk<Nk; kk++) {
					for (uint nn=0; nn<Nn; nn++) {
						total_frey_broadening[ee][kk][nn] = frey_broadening[ee2][kk][nn][cc];
					}
				}
			}
			
			// part from other processes
			uint array_size = Nn*Nk;
			vector<double> tmp_real; tmp_real.resize(array_size, 0.0); // storage: nn*Nn+kk
			vector<double> tmp_imag; tmp_imag.resize(array_size, 0.0); 
			for (int proc=0; proc<mpi->get_num_procs(); proc++) {
				if (proc==constants::mpi_master_rank) continue;
				
				uint start_idx = energies->get_start_global_idx(proc);
				uint stop_idx  = energies->get_stop_global_idx(proc);
			
				int tag = proc;
				for (uint ee=start_idx; ee<=stop_idx; ee++) {
					// receive from other process
					mpi->recv(tmp_real, array_size, proc, tag);
					mpi->recv(tmp_imag, array_size, proc, tag);
					// add to total matrix
					for (uint kk=0; kk<Nk; kk++) {
						for (uint nn=0; nn<Nn; nn++) {
							total_frey_broadening[ee][kk][nn] = tmp_real[nn*Nn+kk] + constants::imag_unit * tmp_imag[nn*Nn+kk];
						}
					}
				}
			}
					
			// some screen output
			//if (frey_converged) {
				logmsg->emit(LOG_INFO_L3,"Displaying n=0 self-energy for every 5th k-point (units energy*length):");
				NEGF_ASSERT(total_frey_broadening.size()==NE, "something is wrong (NE).");
				for (uint ee=0; ee<NE; ee++) {
					NEGF_ASSERT(total_frey_broadening[ee].size()==Nk, "something is wrong (Nk).");
					logmsg->emit(LOG_INFO_L3,"E=%7.3g:  ",energies->get_energy_from_global_idx(ee));
					for (uint kk=0; kk<Nk; kk+=5) {
						NEGF_ASSERT(total_frey_broadening[ee][kk].size()>=1, "something is wrong (Nn).");
						logmsg->emit_noendl(LOG_INFO_L3,"(%7.1e,%7.1e)   ", total_frey_broadening[ee][kk][0].real(), total_frey_broadening[ee][kk][0].imag());
					}
					logmsg->emit(LOG_INFO_L3,"");
				}
			//}
						
		} else {
			// send own part to master process
			int dest = constants::mpi_master_rank;
			int tag = mpi->get_rank();
			vector<double> tmp_real; tmp_real.resize(Nn*Nk, 0.0); // storage: nn*Nn+kk
			vector<double> tmp_imag; tmp_imag.resize(Nn*Nk, 0.0); 
			for (uint ee2=0; ee2<myNE; ee2++) {
				for (uint kk=0; kk<Nk; kk++) {
					for (uint nn=0; nn<Nn; nn++) {
						tmp_real[nn*Nn+kk] = frey_broadening[ee2][kk][nn][cc].real();
						tmp_imag[nn*Nn+kk] = frey_broadening[ee2][kk][nn][cc].imag();
					}
				}
				mpi->send(tmp_real, dest, tag);
				mpi->send(tmp_imag, dest, tag);
			}
		}
	
	} // contact loop
	
	this->se->get_contact_selfenergy()->assign_broadening(frey_broadening); // for both contacts
);}

