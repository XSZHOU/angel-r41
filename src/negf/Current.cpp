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
#include "Current.h"
using namespace negf;

Current::Current(const Hamiltonian * ham_, const Overlap * ov_, const GreenFunctions * gf_, 
		const SelfEnergies * se_, quantities::PhysicalQuantity e_or_h_) throw (Exception *):
	e_or_h(e_or_h_),
	ham(ham_),
	ov(ov_),
	gf(gf_),
	se(se_),
	xspace(gf->get_xspace()),
	kspace(gf->get_kspace()),
	options(gf->get_options()),
	energies(gf->get_energies()),
	Nx((xspace==0) ? 0 : xspace->get_num_internal_vertices()),
	NxNn(Nx*Nn)
{STACK_TRACE(
	NEGF_ASSERT(ham!=0 && ov!=0 && gf!=0 && xspace!=0 && kspace!=0 && options!=0 && energies!=0, "null pointer encountered.");
	NEGF_ASSERT(e_or_h==quantities::electron_density || e_or_h==quantities::hole_density, "expected electrons or holes!");
	
	this->i_am_master = (mpi->get_rank()==constants::mpi_master_rank) ? true : false;
	
	uint myNE = energies->get_my_number_of_points();
	this->spectral_current = Matd(myNE,Nx-1);  // (ee,xx) stores current between vertices xx,xx+1, energy ee
	this->entire_spectral_current = Matd(1,1); // only used in the master process
	this->current.clear();				 		   // only used in master process
	
	this->spectral_current2 = Matd(myNE,Nx);    // DIVERGENCE of SLGG-SGGL
	this->entire_spectral_current2 = Matd(1,1); // only used in the master process
	
	if (i_am_master) {
		//this->contact_current = new ContactCurrent(xspace,energies,&this->entire_spectral_current);
		// need to initialize the matrices because of possible output
		uint NE  = energies->get_number_of_points();
		uint num_x_points = Nx; 
		this->Jleft      = 0.0;
		this->Jright     = 0.0;
		this->Jcoh       = Matd(NE,num_x_points);
		this->Jbuettiker = Matd(NE,num_x_points);
		this->Jgolim     = Matd(NE,num_x_points);
		this->Jgolip     = Matd(NE,num_x_points);
		this->Jpop       = Matd(NE,num_x_points);
		this->Jac        = Matd(NE,num_x_points);
		this->Jspont     = Matd(NE,num_x_points);
		this->Jionimp    = Matd(NE,num_x_points);
	} else{
		//this->contact_current = NULL;
	}
);}


// compute J(x,E) = 1/hbar * Tr(H_{i,i+1}*GL_{i+1,i} - GL_{i,i+1}*H_{i+1,i}) 
void Current::compute_spectral_current() throw (Exception *)
{STACK_TRACE(
	//logmsg->emit_small_header("computing spectral %s current (HG-GH)", ((e_or_h==quantities::electron_density) ? "electron" : "hole"));
	logmsg->emit(LOG_INFO,"Computing spectral %s current (HG-GH)...", ((e_or_h==quantities::electron_density) ? "electron" : "    hole"));
	uint Nk    = kspace->get_number_of_points();
	uint myNE  = energies->get_my_number_of_points();
	double kspace_factor = 1.0 / (negf_math::pow(2.0*constants::pi, kspace->get_dimension()));
	const double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
	// note: 1/(L^d) (see report) compensates with L^d from change sum_k --> int dk
	// the computed quantity has units length^{-2} energy^{-1} time^{-1}
	
	vector<uint> bands;
	double spin = 0.0;
	if (this->e_or_h==quantities::electron_density) {
		options->get_conduction_degrees_of_freedom(bands);
		spin = get_spin_degeneracy(options->get("kp_method"), quantities::electron_density); 
	} else {
		options->get_valence_degrees_of_freedom(bands);
		spin = get_spin_degeneracy(options->get("kp_method"), quantities::hole_density);
	}
	
	const OVMat & M = ov->get_internal_overlap();
	NEGF_ASSERT(M.num_cols()==NxNn && M.num_rows()==M.num_cols(), "wrong overlap matrix.");
	
	// re-initialize to zero
	this->spectral_current = Matd(myNE,Nx-1); // zeros; spectral_current(ee,xx) gives current between node xx and node xx+1


	// initialize helper matrices
	Matc Hsmall(NxNn,NxNn);
	Matc HminusEM(NxNn,NxNn);
	Matc tmp1(Nn,Nn);
	Matc tmp2(Nn,Nn);
	Matc product(Nn,Nn);
	
	/*if (mpi->get_rank()==0) {
        ofstream fout("hamiltonian");
        ham->get_internal(kspace->get_point(0), Hsmall);
        for (uint ii=1; ii<=NxNn; ii++) {
            for (uint jj=1; jj<=NxNn; jj++) {
                fout << Hsmall(ii,jj) << "   ";
            }
            fout << "\n";
        }
        fout.close();
	}*/

	this->J_FisherLee = 0.0;
	
	// integrate k-space
	for (uint kk=0; kk<Nk; kk++) 
	{
		if (kk % 10 == 0) logmsg->emit_noendl_all(LOG_INFO_L2, "p%d: J(kk=%d,:)...   ",mpi->get_rank(),kk);
		double wk = kspace->get_point(kk).get_weight() * kspace_factor;

        // if there is only 1 k-point, assume ballistic calculation and parabolic bands
        // --> use Fermi-integral of order 0 with different EF's for the 2 contacts instead of k-weight
        // that was added to SigmaL
        if (Nk==1) {
            NEGF_ASSERT(Nn==1, "Nk=1 is possible only for single-band effective mass model.");
            wk = 0.5; // spin is added later on
        }

		// construct the Hamiltonian of the interior points only
		ham->get_internal(kspace->get_point(kk), Hsmall);
		NEGF_ASSERT(Hsmall.num_rows()==NxNn, "something went wrong.");
	
		NEGF_ASSERT(xspace->get_dimension()==1, "only 1D current is implemented!");
		for (uint ee2 = 0; ee2 < myNE; ee2++) 
		{
			uint ee = energies->get_global_index(ee2);
			const double E = energies->get_energy_from_global_idx(ee);
					
			const GLMat & GLG = (e_or_h==quantities::electron_density) ? gf->get_lesser(kk, ee) : gf->get_greater(kk, ee);
			
			mult(M, -E, HminusEM); // HminusEM = (-E) * M;
			HminusEM += Hsmall; // HminusEM = HminusEM + Hsmall; // Hsmall+HminusEM would fail!
			for (uint xx=1; xx<=Nx-1; xx++) 
			{
				vector<cplx> tmp3; tmp3.resize(Nn,0.0);
				vector<cplx> tmp4; tmp4.resize(Nn,0.0);
				for (uint mm=1; mm<=Nn; mm++) {
					for (uint nn=1; nn<=Nn; nn++) {
					    // "normal" current - Nx-1 values
						tmp3[mm-1]  += HminusEM(get_mat_idx(xx  ,mm,Nx),get_mat_idx(xx+1,nn,Nx)) * GLG(get_mat_idx(xx+1,nn,Nx),get_mat_idx(xx  ,mm,Nx));
						tmp4[mm-1]  += HminusEM(get_mat_idx(xx+1,mm,Nx),get_mat_idx(xx  ,nn,Nx)) * GLG(get_mat_idx(xx  ,nn,Nx),get_mat_idx(xx+1,mm,Nx));

/*					    // "symmetrized" current - Nx-2 values (kubis)
					    bool xx1 = false;
					    if (xx==1) { xx1 = true; xx=2; }
					    tmp3[mm-1]  += 0.5 * HminusEM(get_mat_idx(xx  ,mm,Nx),get_mat_idx(xx+1,nn,Nx)) * GLG(get_mat_idx(xx+1,nn,Nx),get_mat_idx(xx  ,mm,Nx));
					    tmp4[mm-1]  += 0.5 * HminusEM(get_mat_idx(xx+1,mm,Nx),get_mat_idx(xx  ,nn,Nx)) * GLG(get_mat_idx(xx  ,nn,Nx),get_mat_idx(xx+1,mm,Nx));
                        tmp3[mm-1]  += 0.5 * HminusEM(get_mat_idx(xx-1,mm,Nx),get_mat_idx(xx  ,nn,Nx)) * GLG(get_mat_idx(xx  ,nn,Nx),get_mat_idx(xx-1,mm,Nx));
                        tmp4[mm-1]  += 0.5 * HminusEM(get_mat_idx(xx  ,mm,Nx),get_mat_idx(xx-1,nn,Nx)) * GLG(get_mat_idx(xx-1,nn,Nx),get_mat_idx(xx  ,mm,Nx));
                        if (xx1==true) xx=1;
*/
                        // "effective mass" current, WRONG UNITS!
//					    double tmp1 = abs(HminusEM(get_mat_idx(1  ,1,Nx),get_mat_idx(2,1,Nx)) * GLG(get_mat_idx(2,1,Nx),get_mat_idx(1  ,1,Nx)));
//					    double tmp2 = abs(GLG(get_mat_idx(2,1,Nx),get_mat_idx(1,1,Nx)));
//					    tmp3[mm-1]  += tmp1/(tmp2+1e-20) * GLG(get_mat_idx(xx+1,nn,Nx),get_mat_idx(xx  ,mm,Nx));
//					    tmp4[mm-1]  += tmp1/(tmp2+1e-20) * GLG(get_mat_idx(xx  ,nn,Nx),get_mat_idx(xx+1,mm,Nx));
					}
				}
				
				cplx trace = 0.0;
				for (uint nn=0; nn<bands.size(); nn++) { 
					trace += spin * (tmp3[bands[nn]] - tmp4[bands[nn]]); // bands[nn] starts at 0
				}
				NEGF_FASSERT(fabs(trace.imag()) < constants::imag_err, "encountered complex trace: (%.4e,%.4e)", trace.real(), trace.imag());
				double trace_dbl = trace.real();
				
				spectral_current(ee2+1,xx) += wk * 1.0/hbar * trace_dbl;
			}
			
			// ---------------
			// FISHER-LEE
			// ---------------
			if (false) {
				//const Matc & GR = gf->get_retarded(kk, ee);	
				//const Matc & SR = se->get_retarded(kk, ee);	
				/*	
				const Matc & GA = gf->get_advanced(kk, ee);	
				cplx Gamma_11 = -2.0*SR(1,1).imag(); 
				cplx Gamma_NN = -2.0*SR(Nx,Nx).imag(); 
				double kT = constants::convert_from_SI(units::energy, constants::SIkb * options->get("temperature"));
				double EF1 = 1.437478e+00;
				double EFN = EF1 - 0.050;
				double f1 = 1.0/(1.0+negf_math::exp((E - EF1)/kT));
				double fN = 1.0/(1.0+negf_math::exp((E - EFN)/kT));
				this->J_FisherLee += Gamma_11 * GR(1,Nx) * Gamma_NN * GA(Nx,1) * (f1-fN) 
						* wk * energies->get_point_from_local_idx(ee2).get_weight()/(2.0*constants::pi*hbar)
						* 2.0; // spin
				*//*
				const Matc & SL = se->get_lesser(kk, ee);
				this->J_FisherLee += (SL(1,1) * 2.0*constants::imag_unit*GR(1,1).imag() - GL(1,1) * 2.0*constants::imag_unit*SR(1,1).imag()) 
						* wk * energies->get_point_from_local_idx(ee2).get_weight()/(2.0*constants::pi*hbar);*/
				/*
				cplx HmEM_12 = Hsmall(1,2)-E*M(1,2);
				cplx HmEM_21 = Hsmall(2,1)-E*M(2,1);
				this->J_FisherLee += (HmEM_12*GL(2,1) - HmEM_21*GL(1,2)) 
						* wk * energies->get_point_from_local_idx(ee2).get_weight()/(2.0*constants::pi*hbar);
				*/
			}
		}
	}

	// set up entire_spectral_current of master process
	if (i_am_master) {
		uint NE  = energies->get_number_of_points();
		this->entire_spectral_current = Matd(NE, Nx-1);
		
		// own part
		for (uint ee2 = 0; ee2 < myNE; ee2++) {
			uint ee = energies->get_global_index(ee2);
			for (uint jj=1; jj<=Nx-1; jj++) {
				this->entire_spectral_current(ee+1,jj) = this->spectral_current(ee2+1, jj);
			}
		}
		
		// collect the pieces
		for (int pp=0; pp<mpi->get_num_procs(); pp++) 
		{
			if (pp==constants::mpi_master_rank) continue;
			uint num_energy_points = energies->get_number_of_points(pp);
			
			// receive from other process
			Matd tmp_JE(num_energy_points, Nx-1);
			int tag = pp;
			mpi->recv(tmp_JE, pp, tag);
			
			// add to total matrix
			uint start_idx = energies->get_start_global_idx(pp) + 1;
			uint stop_idx = energies->get_stop_global_idx(pp) + 1;
			for (uint ii=start_idx; ii<=stop_idx; ii++) {
				for (uint jj=1; jj<=Nx-1; jj++) {
					this->entire_spectral_current(ii,jj) = tmp_JE(ii-start_idx+1, jj);
				}
			}
		}
		
		// computation of contact current is now done in PostProcessing::compute_spectral_current()
	} else {
		// send to master process
		int dest = constants::mpi_master_rank;
		int tag = mpi->get_rank();
		mpi->send(this->spectral_current, dest, tag);
	}
	
	mpi->synchronize_processes();
);}


// compute dJ/dx(x,E) = 1/hbar * Tr_perp(SigL*GG - GL*SigG)_ii 
void Current::compute_spectral_current2() throw (Exception *)
{STACK_TRACE(
	//logmsg->emit_small_header("computing spectral %s current (divergence form SL*GG-SG*GL)", ((e_or_h==quantities::electron_density) ? "electron" : "hole"));
	logmsg->emit(LOG_INFO,"Computing spectral %s current (divergence SL*GG-SG*GL)...", ((e_or_h==quantities::electron_density) ? "electron" : "    hole"));
	uint Nk    = kspace->get_number_of_points();
	uint myNE  = energies->get_my_number_of_points();
	double kspace_factor = 1.0 / (negf_math::pow(2.0*constants::pi, kspace->get_dimension()));
	const double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
	// note: 1/(L^d) (see report) compensates with L^d from change sum_k --> int dk
	// the computed quantity has units length^{-2} energy^{-1} time^{-1}
	
	vector<uint> bands;
	double spin = 0.0;
	if (this->e_or_h==quantities::electron_density) {
		options->get_conduction_degrees_of_freedom(bands);
		spin = get_spin_degeneracy(options->get("kp_method"), quantities::electron_density); 
	} else {
		options->get_valence_degrees_of_freedom(bands);
		spin = get_spin_degeneracy(options->get("kp_method"), quantities::hole_density);
	}

	logmsg->emit(LOG_INFO,"1...");
	// re-initialize to zero
	this->spectral_current2 = Matd(myNE,Nx); // zeros
	
	// integrate k-space
	for (uint kk=0; kk<Nk; kk++) 
	{
		if (kk % 10 == 0) logmsg->emit_noendl_all(LOG_INFO_L2, "p%d: J(kk=%d,:)...   ",mpi->get_rank(),kk);
		double wk = kspace->get_point(kk).get_weight() * kspace_factor;

        // if there is only 1 k-point, assume ballistic calculation and parabolic bands
        // --> use Fermi-integral of order 0 with different EF's for the 2 contacts instead of k-weight
        // that was added to SigmaL
        if (Nk==1) {
            NEGF_ASSERT(Nn==1, "Nk=1 is possible only for single-band effective mass model.");
            wk = 0.5; // spin is added later on
        }

		NEGF_ASSERT(xspace->get_dimension()==1, "only 1D current is implemented!");
		for (uint ee2 = 0; ee2 < myNE; ee2++) 
		{
			uint ee = energies->get_global_index(ee2);
					
			const GLMat & GL = gf->get_lesser(kk, ee);
			const GLMat & GG = gf->get_greater(kk, ee);	
			const SEMat & SL = se->get_lesser(kk, ee);
			const SEMat & SG = se->get_greater(kk, ee);
			
			vector<cplx> SLGG; get_diag(SL, GG, SLGG);
			vector<cplx> SGGL; get_diag(SG, GL, SGGL);
			vector<cplx> GGSL; get_diag(GG, SL, GGSL);
			vector<cplx> GLSG; get_diag(GL, SG, GLSG);
			
			for (uint xx=1; xx<=Nx; xx++) 
			{
				cplx trace = 0.0;
				for (uint nn=0; nn<bands.size(); nn++) {					
					uint idx = get_mat_idx(xx,bands[nn]+1,Nx) - 1; // vector<cplx> indices of course start at 0
					// idx is indeed 0-based
					
					//NEGF_ASSERT(SLGG.size()>idx && SGGL.size()>idx && GGSL.size()>idx && GLSG.size()>idx, "invalid idx.");
					trace += spin * 0.5 * (SLGG[idx] - SGGL[idx] + GGSL[idx] - GLSG[idx]);
					//trace += spin * (SLGG[idx] - SGGL[idx]);
				}
				NEGF_FASSERT(fabs(trace.imag()) < 1000*constants::imag_err, "encountered complex trace (2): (%.4e,%.4e)", trace.real(), trace.imag());
				double trace_dbl = trace.real();
				
				double x_cell_length = 0.0;
				Vertex * v = xspace->get_vertex(xspace->get_global_vertex_index(xx-1));
				for (uint ii=0; ii < xspace->get_edges_near(v).size(); ii++) {
					x_cell_length += 0.5 * xspace->get_edges_near(v)[ii]->get_length();
				}
				
				spectral_current2(ee2+1,xx) += wk * 1.0/hbar * 1.0/x_cell_length * trace_dbl;
			}
		}
	}

    logmsg->emit(LOG_INFO,"2...");
	// set up entire_spectral_current of master process
	if (i_am_master) {
		uint NE  = energies->get_number_of_points();
		this->entire_spectral_current2 = Matd(NE, Nx);

        // own part
        for (uint ee2 = 0; ee2 < myNE; ee2++) {
            uint ee = energies->get_global_index(ee2);
            for (uint jj=1; jj<=Nx; jj++) {
                this->entire_spectral_current2(ee+1,jj) = this->spectral_current2(ee2+1, jj);
            }
        }

		// collect the pieces
		for (int pp=0; pp<mpi->get_num_procs(); pp++) 
		{
			if (pp==constants::mpi_master_rank) continue;
			uint num_energy_points = energies->get_number_of_points(pp);
			
			// receive from other process
			Matd tmp_JE(num_energy_points, Nx);
			int tag = pp;
			mpi->recv(tmp_JE, pp, tag);
			
			// add to total matrix
			uint start_idx = energies->get_start_global_idx(pp) + 1;
			uint stop_idx = energies->get_stop_global_idx(pp) + 1;
			for (uint ii=start_idx; ii<=stop_idx; ii++) {
				for (uint jj=1; jj<=Nx; jj++) {
					this->entire_spectral_current2(ii,jj) = tmp_JE(ii-start_idx+1, jj);
				}
			}
		}
	
		// integrate for screen output
		bool screen = false;
		if (screen) {
			double left_curr = 0.0;
			double right_curr = 0.0;
			for (uint ee=0; ee<NE; ee++) {
				double dE = energies->get_weight_from_global_idx(ee) / (2.0*constants::pi);
				left_curr  += this->entire_spectral_current2(ee+1, 1) * dE;
				right_curr += this->entire_spectral_current2(ee+1,Nx) * dE;
			}
			double dxl = 2 * 0.5 * (  xspace->get_vertex(xspace->get_global_vertex_index(   1))->get_coordinate(0)
								    - xspace->get_vertex(xspace->get_global_vertex_index(   0))->get_coordinate(0) );
			double dxr = 2 * 0.5 * (  xspace->get_vertex(xspace->get_global_vertex_index(Nx-1))->get_coordinate(0)
								    - xspace->get_vertex(xspace->get_global_vertex_index(Nx-2))->get_coordinate(0) );
			logmsg->emit(LOG_INFO,"%s current from dx*dJ/dx: left %.3e [A/cm2], right %.3e [A/cm2]",
				((e_or_h==quantities::electron_density) ? "Electron" : "Hole"),
				dxl *  left_curr / constants::convert_from_SI(units::electrical_current_density_3d,1.0)*1e-4,
				dxr * right_curr / constants::convert_from_SI(units::electrical_current_density_3d,1.0)*1e-4);
		}
	} else {		
		// send to master process
		int dest = constants::mpi_master_rank;
		int tag = mpi->get_rank();
		mpi->send(this->spectral_current2, dest, tag);
	}

    logmsg->emit(LOG_INFO,"3...");
	mpi->synchronize_processes();
);}


/* compute dJ_{scatt}/dx(x,E) = 1/hbar * Tr_perp(SigLs*GG - GL*SigGs)_ii, where 
 * SigLGs is the self-energy of the scattering mechanism only and GLG is the normal lesser/greater GF. */
void Current::compute_scattering_current(SelfEnergyType scattering_type) throw (Exception *)
{STACK_TRACE(	
	char buf[1000]; sprintf(buf,"%s",((e_or_h==quantities::electron_density) ? "electron" : "hole"));
	switch (scattering_type) {
	case SEtype_contact:			logmsg->emit_small_header("computing coherent %s current",buf); break;
	case SEtype_buettiker:			logmsg->emit_small_header("computing Buettiker scattering %s current",buf); break;
	case SEtype_golizadeh_momentum:	logmsg->emit_small_header("computing Golizadeh momentum scattering %s current",buf); break;
	case SEtype_golizadeh_phase:	logmsg->emit_small_header("computing Golizadeh phase scattering %s current",buf); break;
	case SEtype_optical_phonon:		logmsg->emit_small_header("computing optical phonon scattering %s current",buf); break;
	case SEtype_acoustic_phonon:	logmsg->emit_small_header("computing acoustic phonon %s current",buf); break;
	case SEtype_spont_photon:		logmsg->emit_small_header("computing spontaneous photon scattering %s current",buf); break;
	case SEtype_ion_imp:			logmsg->emit_small_header("computing ionized impurity scattering %s current",buf); break;
	default: 						logmsg->emit_small_header("computing unknown scattering %s current",buf); break;
	}
	NEGF_ASSERT(se->has_self_energy(scattering_type), "Need to have this scattering type turned on!");
	
	uint Nk    = kspace->get_number_of_points();
	uint myNE  = energies->get_my_number_of_points();
	double kspace_factor = 1.0 / (negf_math::pow(2.0*constants::pi, kspace->get_dimension()));
	const double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
		
	vector<uint> bands;
	double spin = 0.0;
	if (this->e_or_h==quantities::electron_density) {
		options->get_conduction_degrees_of_freedom(bands);
		spin = get_spin_degeneracy(options->get("kp_method"), quantities::electron_density); 
	} else {
		options->get_valence_degrees_of_freedom(bands);
		spin = get_spin_degeneracy(options->get("kp_method"), quantities::hole_density);
	}
	
	const OVMat & M = ov->get_internal_overlap();
	NEGF_ASSERT(M.num_cols()==NxNn && M.num_rows()==M.num_cols(), "wrong overlap matrix.");
	
	// re-initialize to zero
	Matd Jscatt(myNE,Nx);
		
	// integrate k-space
	for (uint kk=0; kk<Nk; kk++) 
	{
		if (kk % 10 == 0) logmsg->emit_noendl_all(LOG_INFO_L2, "p%d: J(kk=%d,:)...   ",mpi->get_rank(),kk);
		double wk = kspace->get_point(kk).get_weight() * kspace_factor;

        // if there is only 1 k-point, assume ballistic calculation and parabolic bands
        // --> use Fermi-integral of order 0 with different EF's for the 2 contacts instead of k-weight
        // that was added to SigmaL
        if (Nk==1) {
            NEGF_ASSERT(Nn==1, "Nk=1 is possible only for single-band effective mass model.");
            wk = 0.5; // spin is added later on
        }

		for (uint ee2 = 0; ee2 < myNE; ee2++) 
		{
			uint ee = energies->get_global_index(ee2);
					
			const GLMat & GL = gf->get_lesser(kk,ee);
			const GLMat & GG = gf->get_greater(kk,ee);
			const SEMat & SL = se->get_selfenergy(scattering_type)->get_lesser(kk,ee);
			const SEMat & SG = se->get_selfenergy(scattering_type)->get_greater(kk,ee);
			
			vector<cplx> SLGG; get_diag(SL, GG, SLGG);
			vector<cplx> SGGL; get_diag(SG, GL, SGGL);
			vector<cplx> GGSL; get_diag(GG, SL, GGSL);
			vector<cplx> GLSG; get_diag(GL, SG, GLSG);
			
			for (uint xx=1; xx<=Nx; xx++) 
			{
				cplx trace = 0.0;
				for (uint nn=0; nn<bands.size(); nn++) {
					uint idx = get_mat_idx(xx,bands[nn]+1,Nx) - 1;
					trace += spin * 0.5 * (SLGG[idx] - SGGL[idx] + GGSL[idx] - GLSG[idx]);
				}
				NEGF_FASSERT(fabs(trace.imag()) < 10000*constants::imag_err, "encountered complex trace: (%.4e,%.4e)", trace.real(), trace.imag());
				double trace_dbl = trace.real();
				
				// we have calculated the divergence of the current, or rather div(J)*dx!
				
				uint x1 = (xx>1)  ? xx-1 : xx;
				uint x2 = (xx<Nx) ? xx+1 : xx;
				double dx = 0.5 * (  xspace->get_vertex(xspace->get_global_vertex_index(x2-1))->get_coordinate(0)
					     	       - xspace->get_vertex(xspace->get_global_vertex_index(x1-1))->get_coordinate(0) );
				
				// for xx=1 and nn=Nx, dx must be doubled because the self-energy assumes that the points are connected to xx=0 and xx=Nx+1 
				// if you don't believe this, compare the ballistic current calculated like this with the result of compute_spectral_current!
				if (xx==1 || xx==Nx) {
					dx *= 2;
				}
				
				Jscatt(ee2+1,xx) += wk * 1.0/hbar * trace_dbl / dx; // units m-3s-1
			}
		}
	}
	
	// set up entire_Jscatt of master process
	if (i_am_master) {
		uint NE  = energies->get_number_of_points();
		Matd * Jscatt_tot = 0;
		switch (scattering_type) {
			case SEtype_contact:			Jscatt_tot = &this->Jcoh; break;
			case SEtype_buettiker:			Jscatt_tot = &this->Jbuettiker; break;
			case SEtype_golizadeh_momentum:	Jscatt_tot = &this->Jgolim; break;
			case SEtype_golizadeh_phase:	Jscatt_tot = &this->Jgolip; break;
			case SEtype_optical_phonon:		Jscatt_tot = &this->Jpop; break;
			case SEtype_acoustic_phonon:	Jscatt_tot = &this->Jac; break;
			case SEtype_spont_photon:		Jscatt_tot = &this->Jspont; break;
			case SEtype_ion_imp:			Jscatt_tot = &this->Jionimp; break;
			default: 						NEGF_EXCEPTION("Unknwon scattering type.");; break;
		}
		*(Jscatt_tot) = Matd(NE, Nx);

        // own part
        for (uint ee2 = 0; ee2 < myNE; ee2++) {
            uint ee = energies->get_global_index(ee2);
            for (uint jj=1; jj<=Nx; jj++) {
                (*Jscatt_tot)(ee+1,jj) = Jscatt(ee2+1, jj);
            }
        }
						
		// collect the pieces
		for (int pp=0; pp<mpi->get_num_procs(); pp++) 
		{
			if (pp==constants::mpi_master_rank) continue;
			uint num_energy_points = energies->get_number_of_points(pp);
			
			// receive from other process
			Matd tmp_Jscatt(num_energy_points, Nx);
			int tag = pp;
			mpi->recv(tmp_Jscatt, pp, tag);
			
			// add to total matrix
			uint start_idx = energies->get_start_global_idx(pp) + 1;
			uint stop_idx = energies->get_stop_global_idx(pp) + 1;
			for (uint ii=start_idx; ii<=stop_idx; ii++) {
				for (uint jj=1; jj<=Nx; jj++) {
					(*Jscatt_tot)(ii,jj) = tmp_Jscatt(ii-start_idx+1,jj);
				}
			}
		}
		
		// in case of contact self-energy, integrate over energy and multiply by dx to get contact current!
		// (we have computed dJ/dx --> multiply by dx to get dJ = J1-J0 = J1)
		if (scattering_type==SEtype_contact) {
			const double ec = constants::convert_from_SI(units::charge, constants::SIec);
			
			double dxleft  = 2 * 0.5 * (  xspace->get_vertex(xspace->get_global_vertex_index(1))->get_coordinate(0)
					     	            - xspace->get_vertex(xspace->get_global_vertex_index(0))->get_coordinate(0) );
			double dxright = 2 * 0.5 * (  xspace->get_vertex(xspace->get_global_vertex_index(Nx-1))->get_coordinate(0)
					     	            - xspace->get_vertex(xspace->get_global_vertex_index(Nx-2))->get_coordinate(0) );
				
			this->Jleft = 0.0;
			this->Jright = 0.0;
			for (uint ee=0; ee<NE; ee++) {
				double dE = energies->get_weight_from_global_idx(ee);
				
				Jleft  += dE/(2.0*constants::pi) * ec * dxleft  * (*Jscatt_tot)(ee+1, 1);
				Jright += dE/(2.0*constants::pi) * ec * dxright * (*Jscatt_tot)(ee+1,Nx);
			}
			
			logmsg->emit(LOG_INFO,"Left  contact current, calculated from div. of coherent c.: %.3e[A/cm2]",
					Jleft  / constants::convert_from_SI(units::electrical_current_density_3d,1.0)*1e-4);
			logmsg->emit(LOG_INFO,"Right contact current, calculated from div. of coherent c.: %.3e[A/cm2]",
					Jright / constants::convert_from_SI(units::electrical_current_density_3d,1.0)*1e-4);
		}
	} else {		
		// send to master process
		int dest = constants::mpi_master_rank;
		int tag = mpi->get_rank();
		mpi->send(Jscatt, dest, tag);
	}
	
	mpi->synchronize_processes();
	logmsg->emit(LOG_INFO, "");
);}


void Current::compute_current() throw (Exception *)
{STACK_TRACE(
	//logmsg->emit_small_header("Integrating over E for %s current (HG-GH)",((e_or_h==quantities::electron_density) ? "electron" : "hole"));
	logmsg->emit_noendl(LOG_INFO,"Integrating over E for %s current (HG-GH)...",((e_or_h==quantities::electron_density) ? "electron" : "    hole"));
	this->current.clear();
	uint myNE = energies->get_my_number_of_points();
	
	// integrate own part
	current.resize(Nx-1, 0.0);
	for (uint ee2 = 0; ee2 < myNE; ee2++) {
		uint ee = energies->get_global_index(ee2);
		double dE = energies->get_weight_from_global_idx(ee)/(2.0*constants::pi);
		for (uint xx = 0; xx < Nx-1; xx++) {
			current[xx] += dE * spectral_current(ee2+1,xx+1);
		}
	}
	
	if (i_am_master) {
		
		// collect the pieces
		for (int pp=0; pp<mpi->get_num_procs(); pp++) {
			if (pp==constants::mpi_master_rank) continue;
			vector<double> tmp_curr;
			int tag = pp;
			uint siz = current.size();
			mpi->recv(tmp_curr, siz, pp, tag);
			NEGF_ASSERT(tmp_curr.size()==Nx-1, "something went wrong.");
			// add to master current
			for (uint xx = 0; xx < Nx-1; xx++) {
				current[xx] += tmp_curr[xx];
			}
		}
		//logmsg->emit(LOG_INFO,"%s current (from HG-GH): %7.4g[A/cm2] @ left contact, %7.4g[A/cm2] @ right contact",
		//		((e_or_h==quantities::electron_density) ? "electron" : "hole"),
		/*logmsg->emit(LOG_INFO," left: %9.3e[A/cm2], right: %9.3e[A/cm2]",
				current   [0]/constants::convert_from_SI(units::electrical_current_density_3d,1.0)*1e-4,
				current[Nx-2]/constants::convert_from_SI(units::electrical_current_density_3d,1.0)*1e-4);*/
		logmsg->emit(LOG_INFO,"");
	} else {
		/*double J_norm = 0.0;
		for (uint xx=0; xx<Nx-1; xx++) {
			J_norm += current[xx]*current[xx];
		}
		logmsg->emit(LOG_INFO,"p%d: |J|=%e",mpi->get_rank(),sqrt(J_norm));*/
		
		// send to master thread
		int dest = constants::mpi_master_rank;
		int tag = mpi->get_rank();
		mpi->send(current, dest, tag);
		
		current.clear();
	}
	mpi->synchronize_processes();
	
	/*
	// communicate fisher-lee current
	if (i_am_master) {
		for (int pp=0; pp<mpi->get_num_procs(); pp++) {
			if (pp==constants::mpi_master_rank) continue;			
			double J_FL_tmp = 0.0;
			int tag = pp;
			mpi->recv(J_FL_tmp, pp, tag);
			this->J_FisherLee += J_FL_tmp;
		}
		
		logmsg->emit(LOG_INFO,"************ FISHER-LEE, 50mV: %e[A/cm^2]",
				J_FisherLee.real()/constants::convert_from_SI(units::electrical_current_density_3d, 1.0)*1e-4);
	} else {
		// send to master thread
		int dest = constants::mpi_master_rank;
		int tag = mpi->get_rank();		
		double JFL_tmp = this->J_FisherLee.real();
		mpi->send(JFL_tmp, dest, tag);
	}
	mpi->synchronize_processes();
	*/
);}


/** get the spectrally resolved current for all energies */
const Matd & Current::get_entire_spectral_current() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(i_am_master, "This routine may be called in the master process only!");
	return this->entire_spectral_current;
);}


/** get the spectrally resolved current for all energies */
const Matd & Current::get_entire_spectral_current2() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(i_am_master, "This routine may be called in the master process only!");
	return this->entire_spectral_current2;
);}


const Matd & Current::get_entire_scattering_current(SelfEnergyType scattering_type) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(i_am_master, "This routine may be called in the master process only!");
	switch(scattering_type) {
	case SEtype_contact:			return this->Jcoh; break;
	case SEtype_buettiker:			return this->Jbuettiker; break;
	case SEtype_golizadeh_momentum:	return this->Jgolim; break;
	case SEtype_golizadeh_phase:	return this->Jgolip; break;
	case SEtype_optical_phonon:		return this->Jpop; break;
	case SEtype_acoustic_phonon:	return this->Jac; break;
	case SEtype_spont_photon:		return this->Jspont; break;
	case SEtype_ion_imp:			return this->Jionimp; break;
	default: 						NEGF_EXCEPTION("Unknwon scattering type.");; break;
	}
	return this->Jcoh;
);}


/** get the total current */
const vector<double> & Current::get_current() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(i_am_master, "This routine may be called in the master process only!");
	return this->current;
);}


ContactCurrent::ContactCurrent(const Geometry * xspace_, const Energies * energies_):
	xspace(xspace_),
	energies(energies_),
	entire_spectral_ecurrent(NULL),
	entire_spectral_hcurrent(NULL)
{STACK_TRACE(	
	NEGF_ASSERT(xspace!=0 && energies!=0, "null pointer enecountered.");
		
	this->voltages.clear();
	this->contact_currents.clear();
	
	// standard eqn stuff
	//this->dependencies.clear(); // no dependencies - update will be called from Current class
	//this->its_type = quantities::electrical_current_density;
	this->number_of_variables = xspace->get_num_contacts();
	this->timestamp = 0;
	this->current_variable_values.clear();
);}


void ContactCurrent::set_entire_spectral_ecurrent(const Matd * entire_spectral_ecurrent_)
{STACK_TRACE(
	NEGF_ASSERT(entire_spectral_ecurrent_!=NULL, "null pointer encountered.");
	this->entire_spectral_ecurrent = entire_spectral_ecurrent_;
);}


void ContactCurrent::set_entire_spectral_hcurrent(const Matd * entire_spectral_hcurrent_)
{STACK_TRACE(
	NEGF_ASSERT(entire_spectral_hcurrent_!=NULL, "null pointer encountered.");
	this->entire_spectral_hcurrent = entire_spectral_hcurrent_;
);}

void ContactCurrent::compute_values(uint new_timestamp)
{STACK_TRACE(
    if (this->current_variable_values.size()!=number_of_variables) {
        this->current_variable_values.resize(number_of_variables, 0.0);
    }
    for (uint ii=0; ii<number_of_variables; ii++) {
        this->current_variable_values[ii] = this->compute_value(ii);
    }
    this->timestamp = new_timestamp;
);}


double ContactCurrent::compute_value(uint line) const
{STACK_TRACE(
	NEGF_ASSERT(xspace->get_dimension()==1, "only 1D real-space is implemented at the moment.");
	NEGF_ASSERT(entire_spectral_ecurrent!=NULL && entire_spectral_hcurrent!=NULL, "null pointer encountered.");
	NEGF_ASSERT(entire_spectral_ecurrent->num_rows()==energies->get_number_of_points()
			 && entire_spectral_ecurrent->num_cols()==xspace->get_num_internal_vertices()-1,
			 "unexpected spectral ecurrent matrix size.");
	NEGF_ASSERT(entire_spectral_hcurrent->num_rows()==energies->get_number_of_points()
			 && entire_spectral_hcurrent->num_cols()==xspace->get_num_internal_vertices()-1,
			 "unexpected spectral hcurrent matrix size.");
	
	Contact * contact = xspace->get_contact(line);
	// try to find a contact vertex which is adjacent to another region
	Vertex * thevertex = 0;
	for (uint ii=0; ii<contact->get_num_contact_vertices(); ii++) {
		Vertex * v = contact->get_contact_vertex(ii);
		const vector<Region *> & regs_near_v = xspace->get_regions_near(v);
		if (regs_near_v.size()>1) {
			NEGF_ASSERT(thevertex==0, "found more than one vertex adjacent to >1 regions.");
			thevertex = v;
		}
	}
	NEGF_ASSERT(thevertex!=0, "could not find any contact vertex directly adjacent to the actual device.");
	
	// get adjacent vertex which is inside the device
	const vector<Edge *> & edges_near_vertex = xspace->get_edges_near(thevertex);
	Vertex * v1 = 0;
	for (uint ii=0; ii < edges_near_vertex.size() /*1 or 2*/; ii++) {
		Vertex * v = edges_near_vertex[ii]->get_lower_vertex()==thevertex ? edges_near_vertex[ii]->get_upper_vertex() : edges_near_vertex[ii]->get_lower_vertex();
		if (!v->is_at_contact()) {
			NEGF_ASSERT(v1==0, "encountered more than 1 device-interior vertex adjacent to the contact.");
			v1 = v;
		}
	}
	NEGF_ASSERT(v1!=0, "could not find device-interior vertex.");
	int internal_idx = v1->get_index_internal();
	NEGF_ASSERT(internal_idx>=0, "expected nonnegative internal vertex index.");
	
	// get adjacent vertex of vertex which is even more inside
	const vector<Edge *> & edges_near_vertex2 = xspace->get_edges_near(v1);
	Vertex * v2 = 0;
	for (uint ii=0; ii < edges_near_vertex.size() /*1 or 2*/; ii++) {
		Vertex * v = edges_near_vertex2[ii]->get_lower_vertex()==v1 ? edges_near_vertex2[ii]->get_upper_vertex() : edges_near_vertex2[ii]->get_lower_vertex();
		if (v!=thevertex) {
			v2 = v;
		}
	}
	NEGF_ASSERT(v2!=0, "could not find device-interior vertex.");
	int internal_idx2 = v2->get_index_internal();
	NEGF_ASSERT(internal_idx2>=0, "expected nonnegative internal vertex index 2.");
	NEGF_ASSERT(fabs(internal_idx2-internal_idx)==1.0, "expected consecutive internal vertex indices.");
	
	uint bigger_index = max(internal_idx,internal_idx2);
	uint current_index = bigger_index-1;
	
	// now integrate over energy!
	double eresult = 0.0;
	double hresult = 0.0;
	uint nE = energies->get_number_of_points();
	const double ec = constants::convert_from_SI(units::charge, constants::SIec);
	for (uint ee=0; ee<nE; ee++) {
		double wE = energies->get_weight_from_global_idx(ee) / (2.0*constants::pi);
		
		cplx jEe = (*entire_spectral_ecurrent)(ee+1,current_index+1);
		NEGF_FASSERT(fabs(jEe.imag())<constants::imag_err, "encountered imaginary je(E)=(%e,%e)",jEe.real(),jEe.imag());
		eresult += jEe.real() * wE * ec;
		
		cplx jEh = (*entire_spectral_hcurrent)(ee+1,current_index+1);
		NEGF_FASSERT(fabs(jEh.imag())<constants::imag_err, "encountered imaginary jh(E)=(%e,%e)",jEh.real(),jEh.imag());
		hresult += jEh.real() * wE * ec;
	}
	
	logmsg->emit(LOG_INFO,"Contact %d: ecurrent %10.3e[A/cm2], hcurrent %10.3e[A/cm2], total %10.3e[A/cm2]", line, 
			eresult/constants::convert_from_SI(units::electrical_current_density_3d,1.0)*1e-4,
			hresult/constants::convert_from_SI(units::electrical_current_density_3d,1.0)*1e-4,
			(eresult+hresult)/constants::convert_from_SI(units::electrical_current_density_3d,1.0)*1e-4);
	
	return eresult+hresult;
);}	


void ContactCurrent::snapshot(double voltage) 
{STACK_TRACE(
	if (mpi->get_rank()!=constants::mpi_master_rank) return;
	
	this->voltages.push_back(voltage);
	this->Jlefts  .push_back(fabs(*(this->eJleft )) + fabs(*(this->hJleft )));
	this->Jrights .push_back(fabs(*(this->eJright)) + fabs(*(this->hJright)));
	vector<double> contact_values;
	NEGF_ASSERT(this->get_num_variables()==xspace->get_num_contacts(), "expected num_contacts values.");
	NEGF_ASSERT(this->current_variable_values.size()==xspace->get_num_contacts(), "values were not yet computed.");
	this->contact_currents.push_back(current_variable_values);
);}


void ContactCurrent::save_to_file(const char * filename)
{STACK_TRACE(
	if (mpi->get_rank()!=constants::mpi_master_rank) return;
	
	NEGF_ASSERT(this->contact_currents.size()>0, "no voltages to save!");
	uint num_voltages = contact_currents.size();
	Matd tmp(num_voltages, 1+2+xspace->get_num_contacts()); // +2 for Jleft, Jright
	for (uint ii=1; ii<=num_voltages; ii++) {
		tmp(ii,1) = this->voltages[ii-1];
		tmp(ii,2) = this->Jlefts  [ii-1] / constants::convert_from_SI(units::electrical_current_density_3d, 1.0) * 1e-4;
		tmp(ii,3) = this->Jrights [ii-1] / constants::convert_from_SI(units::electrical_current_density_3d, 1.0) * 1e-4;
		NEGF_ASSERT(contact_currents[ii-1].size()==xspace->get_num_contacts(), "something went wrong.");
		for (uint jj=1; jj<=xspace->get_num_contacts(); jj++) {
			tmp(ii,3+jj) = contact_currents[ii-1][jj-1] / constants::convert_from_SI(units::electrical_current_density_3d, 1.0) * 1e-4;
		}
	}
	char buf[1000];
	sprintf(buf,"%s.ccurrent",filename);
	logmsg->emit(LOG_INFO,"Saving contact currents (in [A/cm2]) to %s ...",buf);
	string description = "Voltage [V], Jleft (SLGG-SGGL) [A/cm2], Jright (SLGG-SGGL) [A/cm2], Jleft (HG-GH) [A/cm2], Jright (HG-GH) [A/cm2]";
	write_matrix(buf, tmp, description);
);}

