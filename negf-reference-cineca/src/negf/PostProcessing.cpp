/*
Copyright (c) 2010 Sebastian Steiger, Integrated Systems Laboratory, ETH Zurich.
Comments, suggestions, criticism or bug reports are welcome: steiger@purdue.edu. 

This file is part of ANGEL, a simulator for LEDs based on the NEGF formalism.
The software is distributed under the Lesser GNU General Public License (LGPL).
ANGEL is free software: you can redistribute it and/or modify it under the terms 
of the Lesser GNU General Public License v3 or later. ANGEL is distributed
without any warranty; without even the implied warranty of merchantability or 
fitness for a particular purpose. See also <http://www.gnu.org/licenses/>.
*/
#include "PostProcessing.h"
using namespace negf;

PostProcessing::PostProcessing( const Hamiltonian * ham_,
								const Overlap * ov_,
								const GreenFunctions * gf_, 
								const SelfEnergies * se_, 
								const Options * opts_,
								const Geometry * xspace_,
								const Kspace * kspace_, 
								const Energies * energies_) throw (Exception *):
	xspace(xspace_),
	kspace(kspace_),
	energies(energies_),
	opts(opts_),
	ham(ham_),
	ov(ov_),
	gf(gf_),
	se(se_),
	Nx((xspace==0) ? 0 : xspace->get_num_internal_vertices()),
	NxNn(Nx*Nn),
	ecurrent(NULL),
	hcurrent(NULL),
	contact_current(NULL),
	transmission(NULL),
	qfl(NULL),
	luminescence(NULL),
	luminescence2(NULL)
{STACK_TRACE(
	logmsg->emit_header("setting up post-processing quantities");
	NEGF_ASSERT(xspace!=NULL && kspace!=NULL && energies!=NULL && opts!=NULL && ham!=NULL && ov!=NULL && gf!=NULL && se!=NULL, "null pointer encountered.");
	
	this->i_am_master = (mpi->get_rank()==constants::mpi_master_rank) ? true : false;
		
	// allocate memory for LDOS and spectral density for the energy points of the calling process
	uint myNE               = energies->get_my_number_of_points();
	this->LDOS              = Matd(myNE,Nx);
	this->LDOS_k0           = Matd(myNE,Nx);
	this->LDOS_VB           = Matd(myNE,Nx);
	this->spectral_edensity = Matd(myNE,Nx);
	this->spectral_hdensity = Matd(myNE,Nx);
	this->ecurrent          = new Current(ham,ov,gf,se,quantities::electron_density);
	this->hcurrent          = new Current(ham,ov,gf,se,quantities::hole_density);
	if (mpi->get_rank()==constants::mpi_master_rank) {
		this->contact_current = new ContactCurrent(xspace,energies);
		this->contact_current->set_entire_spectral_ecurrent(&this->ecurrent->get_entire_spectral_current());
		this->contact_current->set_entire_spectral_hcurrent(&this->hcurrent->get_entire_spectral_current());
		this->contact_current->set_eJleft_eJright(&this->ecurrent->get_Jleft(), &this->ecurrent->get_Jright());
		this->contact_current->set_hJleft_hJright(&this->hcurrent->get_Jleft(), &this->hcurrent->get_Jright());
	}
	if (opts->exists("Transmission") && opts->get("Transmission")==1) {
		this->transmission = new Transmission(xspace, kspace, energies, gf, se->get_contact_selfenergy());
	}
    if (opts->exists("QuasiFermilevels") && opts->get("QuasiFermilevels")==1) {
        this->qfl = new QuasiFermilevel(xspace, kspace, energies, gf);
    }
			
	// the following will only be used in the master process
	this->entire_spectral_edensity = Matd(1,1);
	this->entire_spectral_hdensity = Matd(1,1);
	this->entire_LDOS              = Matd(1,1);
	this->entire_LDOS_k0           = Matd(1,1);
	this->entire_LDOS_VB           = Matd(1,1);
	this->edensity.clear();
	this->hdensity.clear();
);}
								
								
PostProcessing::~PostProcessing()
{STACK_TRACE(
	delete this->ecurrent;
	delete this->hcurrent;
	if (this->contact_current != NULL) {
		delete this->contact_current;
	}
	if (this->luminescence != NULL) {
		delete this->luminescence;
	}
	if (this->luminescence2 != NULL) {
		delete this->luminescence2;
	}
);}


void PostProcessing::compute_local_dos() throw (Exception *)
{STACK_TRACE(
	//logmsg->emit_small_header("computing local DOS");
	logmsg->emit(LOG_INFO,"Computing local DOS...");
	uint Nk = kspace->get_number_of_points();
	uint myNE = energies->get_my_number_of_points();
	double kspace_factor = 1.0 / (negf_math::pow(2.0*constants::pi, kspace->get_dimension()));
	// note: 1/(L^d) (see report) compensates with L^d from change sum_k --> int dk
	// the computed quantity has units length^{-3} energy^{-3}
	vector<uint> vb_bands;
	opts->get_valence_degrees_of_freedom(vb_bands);
	
	this->LDOS    = Matd(myNE,Nx); // zeros
	this->LDOS_k0 = Matd(myNE,Nx); // zeros
	this->LDOS_VB = Matd(myNE,Nx); // zeros
	uint negative_count = 0;
	for (uint ee2 = 0; ee2 < myNE; ee2++) 
	{
		uint ee = energies->get_global_index(ee2);
		logmsg->emit_noendl_all(LOG_INFO_L2, "p%d: LDOS(E=%d,:)...   ",mpi->get_rank(),ee);
		
		for (uint kk = 0; kk < Nk; kk++) 
		{
			double wk = kspace->get_point(kk).get_weight() * kspace_factor;
			const Matc & GR = gf->get_retarded(kk, ee);

	        // if there is only 1 k-point, assume ballistic calculation and parabolic bands
			// --> use Fermi-integral of order 0 with different EF's for the 2 contacts instead of k-weight
			// that was added to SigmaL
			if (Nk==1) {
			    NEGF_ASSERT(Nn==1, "Nk=1 is possible only for single-band effective mass model.");
			    wk = 0.5; // spin is added later on
			}
		
			for (uint xx = 1; xx <= Nx; xx++) 
			{
				Matc GRxx(Nn,Nn); GR.get_block(xx,xx, GRxx, Nx);
				for (uint nn = 1; nn <= Nn; nn++) {
					double spin = get_spin_degeneracy(opts->get("kp_method"), quantities::electron_density); 
					
					// LDOS = i/2pi * (GR-GA)
					double contrib = (-1.0/constants::pi) * spin * GRxx(nn,nn).imag();
					LDOS(ee2+1,xx) += wk*contrib;
					if (kk==0) {
						LDOS_k0(ee2+1,xx) += contrib;
					}
					if (Nn==2 && nn==2) {
						LDOS_VB(ee2+1,xx) += wk*contrib;
					}
				}
				if (LDOS(ee2+1,xx) < -1e-9) {
					negative_count++;
					if (negative_count<100) {
						string buf;
						if (Nn==2) {
							if (fabs(GRxx(1,1).imag()) < 1e-20) { buf = " (from VB)"; }
							if (fabs(GRxx(2,2).imag()) < 1e-20) { buf = " (from CB)"; }
						}
						logmsg->emit_noendl_all(LOG_INFO,"Negative LDOS at xx=%d, E=%.2e%s: %.2e       ",xx,energies->get_energy_from_global_idx(ee),buf.c_str(),LDOS(ee2+1,xx));
					} else if (negative_count==100) {
						logmsg->emit_noendl_all(LOG_INFO,"Skipping further messages about negative LDOS (>100)...      ");
					}
				}
				if (LDOS_VB(ee2+1,xx) < -1e-9) {
					logmsg->emit(LOG_INFO_L2,"Negative VB-LDOS at xx=%d, E=%.2e: %.2e",xx,energies->get_energy_from_global_idx(ee),LDOS_VB(ee2+1,xx));
				}
			}
		}
		
		// correct with geometrical factor in the case of dumb orthogonal-basis-hamiltonian
		// in that case, GF have units energy^{-1} and not  energy^{-1}*length^{-1}
		if (constants::old_orthogonal 
			&& (fabs(opts->get("kp_method") - 0.0) < 1e-14 || fabs(opts->get("kp_method") - 3.0) < 1e-14)) {
			for (uint xx = 1; xx <= Nx; xx++) {	
				double x_cell_length = 0.0;
				Vertex * v = xspace->get_vertex(xspace->get_global_vertex_index(xx-1));
				for (uint ii=0; ii < xspace->get_edges_near(v).size(); ii++) {
					x_cell_length += 0.5 * xspace->get_edges_near(v)[ii]->get_length();
				}
				LDOS   (ee2+1,xx) = LDOS   (ee2+1,xx) / x_cell_length;
				LDOS_k0(ee2+1,xx) = LDOS_k0(ee2+1,xx) / x_cell_length;
				LDOS_VB(ee2+1,xx) = LDOS_VB(ee2+1,xx) / x_cell_length;
			}
		}
	}
	
	// set up entire_LDOS of master process
	mpi->synchronize_processes();
	if (i_am_master) {
		uint NE  = energies->get_number_of_points();
		this->entire_LDOS    = Matd(NE, Nx);
		this->entire_LDOS_k0 = Matd(NE, Nx);
		this->entire_LDOS_VB = Matd(NE, Nx);
		
		// own part
		for (uint ee2 = 0; ee2 < myNE; ee2++) {
			uint ee = energies->get_global_index(ee2);
			for (uint xx=1; xx<=Nx; xx++) {
				this->entire_LDOS   (ee+1,xx) = LDOS   (ee2+1,xx);
				this->entire_LDOS_k0(ee+1,xx) = LDOS_k0(ee2+1,xx);
				this->entire_LDOS_VB(ee+1,xx) = LDOS_VB(ee2+1,xx);
			}
		}
		
		// collect the pieces
		for (int pp=0; pp<mpi->get_num_procs(); pp++) 
		{
			if (pp==constants::mpi_master_rank) continue;
			uint num_energy_points = energies->get_number_of_points(pp);
			
			// receive from other process
			int tag = pp;
			Matd tmp_dos(num_energy_points, Nx);
			Matd tmp_dos_k0(num_energy_points, Nx);
			Matd tmp_dos_VB(num_energy_points, Nx);
			mpi->recv(tmp_dos, pp, tag);
			mpi->recv(tmp_dos_k0, pp, tag);
			mpi->recv(tmp_dos_VB, pp, tag);
			
			// add to total matrix
			uint start_idx = energies->get_start_global_idx(pp) + 1;
			uint stop_idx  = energies->get_stop_global_idx(pp)  + 1;
			for (uint ee=start_idx; ee<=stop_idx; ee++) {
				for (uint xx=1; xx<=Nx; xx++) {
					this->entire_LDOS   (ee,xx) = tmp_dos   (ee-start_idx+1, xx);
					this->entire_LDOS_k0(ee,xx) = tmp_dos_k0(ee-start_idx+1, xx);
					this->entire_LDOS_VB(ee,xx) = tmp_dos_VB(ee-start_idx+1, xx);
				}
			}
		}
	} else {
		// send to master process
		int dest = constants::mpi_master_rank;
		int tag = mpi->get_rank();
		mpi->send(LDOS, dest, tag);
		
		// send k=0 to master process
		mpi->send(LDOS_k0, dest, tag);
		
		// send VB to master process
		mpi->send(LDOS_VB, dest, tag);
	}
	
	mpi->synchronize_processes();
);}

void PostProcessing::compute_spectral_edensity() throw (Exception *)
{STACK_TRACE(
	vector<uint> cb_bands; opts->get_conduction_degrees_of_freedom(cb_bands);
	this->compute_spectral_xdensity(cb_bands, false);
);}

void PostProcessing::compute_spectral_hdensity() throw (Exception *)
{STACK_TRACE(
	vector<uint> vb_bands; opts->get_valence_degrees_of_freedom(vb_bands);
	this->compute_spectral_xdensity(vb_bands, true);
);}

void PostProcessing::compute_spectral_xdensity(const vector<uint> & bands, bool e_or_h /* true --> holes, false --> electrons */)
{STACK_TRACE(
	//logmsg->emit_small_header("computing spectral %s density", ((e_or_h==false) ? "electron" : "hole"));
	logmsg->emit(LOG_INFO,"Computing spectral %s density...", ((e_or_h==false) ? "electron" : "    hole"));
	uint Nk = kspace->get_number_of_points();
	uint myNE = energies->get_my_number_of_points();
	uint NE = energies->get_number_of_points();
	double kspace_factor = 1.0 / (negf_math::pow(2.0*constants::pi, kspace->get_dimension()));
	// note: 1/(L^d) (see report) compensates with L^d from change sum_k --> int dk
	// the computed quantity has units length^{-3} energy^{-3}
	
	double spin = get_spin_degeneracy(opts->get("kp_method"), (e_or_h) ? quantities::hole_density : quantities::electron_density);
	
	Matd & spectral_density = (e_or_h) ? this->spectral_hdensity : this->spectral_edensity;
	// re-initialize to zero!!
	spectral_density = Matd(myNE,Nx);
	
	// compute!
	for (uint ee2 = 0; ee2 < myNE; ee2++) 
	{
		uint ee = energies->get_global_index(ee2);
		if (ee % 5 == 0) logmsg->emit_noendl_all(LOG_INFO_L2, "p%d: %s(E=%d,:)...   ",	mpi->get_rank(),(e_or_h ? "p" : "n"), ee);
		// FLENS indices start with 1!
		for (uint xx=1; xx<=Nx; xx++) 
		{
			cplx result = 0.0;
			// integrate!
			for (uint kk=0; kk<Nk; kk++) 
			{
				double wk = kspace->get_point(kk).get_weight() * kspace_factor;
				const GLMat & GLG = (e_or_h==false) ? gf->get_lesser(kk, ee) : gf->get_greater(kk, ee);
				NEGF_FASSERT(GLG.num_rows()==NxNn && GLG.num_cols()==NxNn, 
						"GL or GG is %dx%d instead of %dx%d", GLG.num_rows(),GLG.num_cols(), NxNn,NxNn);

	            // if there is only 1 k-point, assume ballistic calculation and parabolic bands
	            // --> use Fermi-integral of order 0 with different EF's for the 2 contacts instead of k-weight
	            // that was added to SigmaL
	            if (Nk==1) {
	                NEGF_ASSERT(Nn==1, "Nk=1 is possible only for single-band effective mass model.");
			    	wk = 0.5; // spin is added later on
	            }

				
				Matc GLGxx(Nn,Nn); GLG.get_block(xx,xx, GLGxx, Nx);
				for (uint nidx = 0; nidx < bands.size()/*Nn*/; nidx++) 
				{
					uint nn = bands[nidx]/*nidx*/ + 1;
					cplx tmp = GLGxx(nn,nn);
					
					// check for NaN's
					if (isnan(tmp.real()) || isnan(tmp.imag())) 
					{
						logmsg->emit_noendl_all(LOG_ERROR,"GLG(xx=%d, nn=%d; kk=%d, E=%.3e) = (%.3e, %.3e)",
								xx,nn,kk,energies->get_energy_from_global_idx(ee),
								tmp.real(), tmp.imag());
						// recovery attempt: if neighbouring k-values at same energy are ~0, it's still OK
						bool recover = true;
						if (kk>0 && kk<Nk-1) {
							const GLMat & GLG2 = (e_or_h==false) ? gf->get_lesser(kk-1, ee) : gf->get_greater(kk-1, ee);
							Matc GLG2xx(Nn,Nn); GLG2.get_block(xx,xx, GLG2xx, Nx);
							cplx tmp2 = GLG2xx(nn,nn);
							if (isnan(tmp2.real()) || isnan(tmp2.imag()) || abs(tmp2) > 1e-20) {
								recover = false;
							}
							
							const GLMat & GLG3 = (e_or_h==false) ? gf->get_lesser(kk+1, ee) : gf->get_greater(kk+1, ee);
							Matc GLG3xx(Nn,Nn); GLG3.get_block(xx,xx, GLG3xx, Nx);
							cplx tmp3 = GLG3xx(nn,nn);
							if (isnan(tmp3.real()) || isnan(tmp3.imag()) || abs(tmp3) > 1e-20) {
								recover = false;
							}
						} else {
							recover = false;
						}
						if (recover) {
							logmsg->emit_all(LOG_ERROR," recovering...");
							tmp = 0.0;
						} else {
							NEGF_EXCEPTION("NaN in GL or GG encountered!");
						}
					}
					
					result += wk * tmp;
				}
			}
			result = 
				((e_or_h==false) ? -1.0 : 1.0)/*-1.0*/ * constants::imag_unit // +- i * ...
				 * spin	// YES! NO!!!! spin was already treated in contact fermilevel setting and lesser contact self-energy
				* result;
			if (isnan(result.real()) || isnan(result.imag())) {
				NEGF_EXCEPTION("NaN in spectral density encountered!");
			}
			NEGF_FASSERT(fabs(result.imag()) < constants::imag_err, "imaginary spectral density encountered! (%e, %e)",
						result.real(), result.imag());
			if (result.real() <= constants::convert_from_SI(units::density_3d, 1e13)) {
				if (result.real() < -constants::convert_from_SI(units::density_3d, 1e13)) {
					//NEGF_FEXCEPTION("negative spectral density encountered! (%e, %e) E=%.3g, xx=%d",
					logmsg->emit_noendl_all(LOG_ERROR,"negative %s(xx=%d,E=%.2e): %.2g         ",
							((e_or_h==false) ? "n" : "p"), xx, energies->get_energy_from_global_idx(ee), result.real()/*, result.imag()*/);
				}
				result = constants::convert_from_SI(units::density_3d, 1e13);
			}
			
			// correct with geometrical factor in the case of dumb orthogonal-basis-hamiltonian
			// in that case, GF have units energy^{-1} and not  energy^{-1}*length^{-1}
			if (constants::old_orthogonal 
				&& (fabs(opts->get("kp_method") - 0.0) < 1e-14 || fabs(opts->get("kp_method") - 3.0) < 1e-14)) {
				double x_cell_length = 0.0;
				Vertex * v = xspace->get_vertex(xspace->get_global_vertex_index(xx-1));
				for (uint ii=0; ii < xspace->get_edges_near(v).size(); ii++) {
					x_cell_length += 0.5 * xspace->get_edges_near(v)[ii]->get_length();
				}
				result = result / x_cell_length;
			}
			
			spectral_density(ee2+1,xx) = result.real();
		}
	}
	
	// set up entire_spectral_xdensity of master process
	mpi->synchronize_processes();
	if (i_am_master) {
		Matd & entire_spectral_xdensity = (e_or_h) ? this->entire_spectral_hdensity : this->entire_spectral_edensity;
		entire_spectral_xdensity = Matd(NE, Nx);
		
		// own contribution
		for (uint ee2=0; ee2<myNE; ee2++) {
			uint ee = energies->get_global_index(ee2);
			for (uint xx=1; xx<=Nx; xx++) {
				entire_spectral_xdensity(ee+1,xx) = spectral_density(ee2+1,xx);
			}
		}		

		// collect the pieces
		for (int pp=0; pp<mpi->get_num_procs(); pp++) 
		{
			if (pp==constants::mpi_master_rank) continue;
			uint num_energy_points = energies->get_number_of_points(pp);
			
			// receive from other process
			Matd tmp_np(num_energy_points, Nx);
			int tag = pp;
			mpi->recv(tmp_np, pp, tag);
			
			// add to total matrix
			uint start_idx = energies->get_start_global_idx(pp) + 1;
			uint stop_idx  = energies->get_stop_global_idx(pp)  + 1;
			for (uint ee=start_idx; ee<=stop_idx; ee++) {
				for (uint xx=1; xx<=Nx; xx++) {
					entire_spectral_xdensity(ee,xx) = tmp_np(ee-start_idx+1, xx);
				}
			}
		}
	} else {
		// send to master process
		int dest = constants::mpi_master_rank;
		int tag = mpi->get_rank();
		mpi->send(spectral_density, dest, tag);
	}
	mpi->synchronize_processes();
);}


void PostProcessing::compute_edensity() throw (Exception *)
{STACK_TRACE(
	//logmsg->emit_small_header("computing electron density");
	logmsg->emit(LOG_INFO,"Computing electron density...");
	this->compute_xdensity(this->edensity, this->spectral_edensity);	
);}

void PostProcessing::compute_hdensity() throw (Exception *)
{STACK_TRACE(
	//logmsg->emit_small_header("computing hole density");
	logmsg->emit(LOG_INFO,"Computing     hole density...");
	this->compute_xdensity(this->hdensity, this->spectral_hdensity);	
);}

void PostProcessing::compute_xdensity(vector<double> & density, const Matd & spectral_density)
{STACK_TRACE(
	density.clear();
	
	if (i_am_master) {
		density.resize(Nx, 0.0);
		
		// integrate own part
		uint myNE = energies->get_my_number_of_points();
		// FLENS indices start with 1, but density array starts with 0!
		for (uint xx = 0; xx < Nx; xx++) {
			for (uint ee2 = 0; ee2 < myNE; ee2++) {
				uint ee = energies->get_global_index(ee2);
				NEGF_ASSERT(energies->get_energy_from_local_idx(ee2)==energies->get_energy_from_global_idx(ee),	"inconsistency.");
				double tmp = energies->get_weight_from_global_idx(ee)/(2.0*constants::pi)
							* spectral_density(ee2+1,xx+1);
				density[xx] += tmp;
			}
		}
		
		// collect the pieces
		for (int pp=0; pp<mpi->get_num_procs(); pp++) {
			if (pp==constants::mpi_master_rank) continue;
			vector<double> tmp_dens;
			int tag = pp;
			uint siz = density.size();
			mpi->recv(tmp_dens, siz, pp, tag);
			NEGF_ASSERT(tmp_dens.size()==Nx, "something went wrong.");
			// add to master density
			for (uint xx = 0; xx < Nx; xx++) {
				density[xx] += tmp_dens[xx];
			}
		}
	} else {
		density.resize(Nx, 0.0);
		
		// integrate own part
		// FLENS indices start with 1, but density array starts with 0!
		uint myNE = energies->get_my_number_of_points();
		for (uint xx = 0; xx < Nx; xx++) {
			for (uint ee2 = 0; ee2 < myNE; ee2++) {
				uint ee = energies->get_global_index(ee2);
				double tmp = energies->get_weight_from_global_idx(ee)/(2.0*constants::pi) * spectral_density(ee2+1,xx+1);
				density[xx] += tmp;
			}
		}
		double n_norm = 0.0;
		for (uint xx=0; xx<Nx; xx++) {
			n_norm += density[xx]*density[xx];
		}
		logmsg->emit_noendl_all(LOG_INFO_L3,"p%d: |n|=%e       ",mpi->get_rank(),sqrt(n_norm));
		
		// send to master thread
		int dest = constants::mpi_master_rank;
		int tag = mpi->get_rank();
		mpi->send(density, dest, tag);
		
		density.clear();
	}
	mpi->synchronize_processes();
);}


const vector<double> & PostProcessing::get_edensity() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(i_am_master, "This routine may be called in the master process only!");
	return this->edensity;
);}

const vector<double> & PostProcessing::get_hdensity() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(i_am_master, "This routine may be called in the master process only!");
	return this->hdensity;
);}


/** get the spectrally resolved edensity for all energies */
const Matd & PostProcessing::get_entire_spectral_edensity() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(i_am_master, "This routine may be called in the master process only!");
	return this->entire_spectral_edensity;
);}


/** get the spectrally resolved hdensity for all energies */
const Matd & PostProcessing::get_entire_spectral_hdensity() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(i_am_master, "This routine may be called in the master process only!");
	return this->entire_spectral_hdensity;
);}


/** get the spectrally resolved local density of states for all energies */
const Matd & PostProcessing::get_entire_local_dos() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(i_am_master, "This routine may be called in the master process only!");
	return this->entire_LDOS;
);}


/** get the spectrally resolved local density of states for all energies, k=0 only b(--> different units) */
const Matd & PostProcessing::get_entire_local_dos_k0() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(i_am_master, "This routine may be called in the master process only!");
	return this->entire_LDOS_k0;
);}

/** get the spectrally resolved local density of states for all energies, k=0 only b(--> different units) */
const Matd & PostProcessing::get_entire_local_dos_VB() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(i_am_master, "This routine may be called in the master process only!");
	return this->entire_LDOS_VB;
);}

void PostProcessing::set_up_luminescence(SEPhotonSpontaneous * spont) throw (Exception *)
{STACK_TRACE(
	logmsg->emit_header("Setting up luminescence calculator");
	this->luminescence  = new Luminescence(ov,xspace,kspace,energies,opts,ham,gf,spont,lake);
	this->luminescence2 = new Luminescence(ov,xspace,kspace,energies,opts,ham,gf,spont,galerpin);
);}

void PostProcessing::compute_luminescence() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->luminescence!=NULL, "Luminescence calculator is not set up!");
	//this->luminescence->calculate();
	logmsg->emit(LOG_INFO,"****** --> Computation of Lake-style luminescence is skipped! <-- ******");
	NEGF_ASSERT(this->luminescence2!=NULL, "Luminescence calculator is not set up!");
	this->luminescence2->calculate();
);}

Luminescence * PostProcessing::get_luminescence() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->luminescence!=NULL, "Luminescence calculator is not set up!");
	return this->luminescence;
);}

Luminescence * PostProcessing::get_luminescence2() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->luminescence2!=NULL, "Luminescence2 calculator is not set up!");
	return this->luminescence2;
);}


void PostProcessing::compute_spectral_current() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->ecurrent!=NULL && this->hcurrent!=NULL, "null pointer encountered.");
	this->ecurrent->compute_spectral_current();
	this->hcurrent->compute_spectral_current();
	this->ecurrent->compute_spectral_current2();
	this->hcurrent->compute_spectral_current2();
	
	// integrate divergence current for screen output
	if (mpi->get_rank()==constants::mpi_master_rank) 
	{
		double left_ecurr  = 0.0;
		double right_ecurr = 0.0;
		double left_hcurr  = 0.0;
		double right_hcurr = 0.0;
		for (uint ee=0; ee<energies->get_number_of_points(); ee++) {
			double dE = energies->get_weight_from_global_idx(ee) / (2.0*constants::pi);
			 left_ecurr += this->ecurrent->get_entire_spectral_current2()(ee+1, 1) * dE;
			right_ecurr += this->ecurrent->get_entire_spectral_current2()(ee+1,Nx) * dE;
			 left_hcurr += this->hcurrent->get_entire_spectral_current2()(ee+1, 1) * dE;
			right_hcurr += this->hcurrent->get_entire_spectral_current2()(ee+1,Nx) * dE;
		}
		double dxl = 2 * 0.5 * (  xspace->get_vertex(xspace->get_global_vertex_index(   1))->get_coordinate(0)
							    - xspace->get_vertex(xspace->get_global_vertex_index(   0))->get_coordinate(0) );
		double dxr = 2 * 0.5 * (  xspace->get_vertex(xspace->get_global_vertex_index(Nx-1))->get_coordinate(0)
							    - xspace->get_vertex(xspace->get_global_vertex_index(Nx-2))->get_coordinate(0) );
		
		const double   ec = constants::convert_from_SI(units::charge, constants::SIec);
		const double conv = constants::convert_from_SI(units::electrical_current_density_3d,1.0)*1e4; // 1e+4, NOT 1e-4
		logmsg->emit(LOG_INFO,"Left  contact (SLGG-SGGL): ecurrent %10.3e[A/cm2], hcurrent %10.3e[A/cm2], total %10.3e[A/cm2]",
				dxl * ec *  left_ecurr / conv, dxl * ec *  left_hcurr / conv, dxl * ec * (left_ecurr + left_hcurr) / conv);
		logmsg->emit(LOG_INFO,"Right contact (SLGG-SGGL): ecurrent %10.3e[A/cm2], hcurrent %10.3e[A/cm2], total %10.3e[A/cm2]",
				dxr * ec * right_ecurr / conv, dxr * ec * right_hcurr / conv, dxr * ec * (right_ecurr + right_hcurr) / conv);
	}
	
	logmsg->emit(LOG_INFO,"Computing contact current (from HG-GH):");
	if (mpi->get_rank()==constants::mpi_master_rank) {
		this->contact_current->set_timestamp(9999);
		this->contact_current->compute_values(10000);	// contact_current took HG-GH - calculated current
	}
);}

void PostProcessing::compute_current() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->ecurrent!=NULL && this->hcurrent!=NULL, "null pointer encountered.");
	this->ecurrent->compute_current();
	this->hcurrent->compute_current();
);}

void PostProcessing::compute_scattering_current(SelfEnergyType type) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->ecurrent!=NULL && this->hcurrent!=NULL, "null pointer encountered.");
	this->ecurrent->compute_scattering_current(type);
	this->hcurrent->compute_scattering_current(type);
);}

ContactCurrent * PostProcessing::get_contact_current() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(mpi->get_rank()==constants::mpi_master_rank && this->contact_current!=NULL, "not master thread or ccurrent not set up.");
	return this->contact_current;
);}

Transmission * PostProcessing::get_transmission() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(mpi->get_rank()==constants::mpi_master_rank && this->transmission!=NULL, "not master thread or transmission not set up.");
	return this->transmission;
);}

void PostProcessing::compute_transmission() const throw (Exception *) // needs to be called by all threads
{STACK_TRACE(
	NEGF_ASSERT(this->transmission!=NULL, "Transmission not set up.");
	this->transmission->calculate();
);}

QuasiFermilevel * PostProcessing::get_quasi_fermilevel() const throw (Exception *)
{STACK_TRACE(
    NEGF_ASSERT(mpi->get_rank()==constants::mpi_master_rank && this->qfl!=NULL, "not master thread or QFL not set up.");
    return this->qfl;
);}

void PostProcessing::compute_quasi_fermilevel() const throw (Exception *) // needs to be called by all threads
{STACK_TRACE(
    NEGF_ASSERT(this->qfl!=NULL, "QFL not set up.");
    if (mpi->get_rank()==constants::mpi_master_rank) {
        this->qfl->set_densities(this->get_edensity(), this->get_hdensity());
    }
    this->qfl->calculate();
);}


void PostProcessing::check_conservation(SelfEnergyType type) const throw (Exception *)
{STACK_TRACE(
	// screen output
	/*switch (type) {
	case SEtype_contact:			logmsg->emit_small_header("checking coherent particle conservation"); break;
	case SEtype_buettiker:			logmsg->emit_small_header("checking Buettiker scattering particle conservation"); break;
	case SEtype_golizadeh_momentum:	logmsg->emit_small_header("checking Golizadeh momentum scattering particle conservation"); break;
	case SEtype_golizadeh_phase:	logmsg->emit_small_header("checking Golizadeh phase scattering particle conservation"); break;
	case SEtype_optical_phonon:		logmsg->emit_small_header("checking optical phonon scattering particle conservation"); break;
	case SEtype_acoustic_phonon:	logmsg->emit_small_header("checking acoustic phonon particle conservation"); break;
	case SEtype_spont_photon:		logmsg->emit_small_header("checking spontaneous photon scattering particle conservation"); break;
	case SEtype_ion_imp:			logmsg->emit_small_header("checking ionized impurity scattering particle conservation"); break;
	case SEtype_all:				logmsg->emit_small_header("checking overall particle conservation"); break;
	default: 						logmsg->emit_small_header("checking unknown scattering particle conservation"); break;
	}*/
	switch (type) {
	case SEtype_contact:			logmsg->emit_noendl(LOG_INFO,"Coherent particle conservation:           "); break;
	case SEtype_buettiker:			logmsg->emit_noendl(LOG_INFO,"Buettiker particle conservation:          "); break;
	case SEtype_golizadeh_momentum:	logmsg->emit_noendl(LOG_INFO,"Golizadeh momentum particle conservation: "); break;
	case SEtype_golizadeh_phase:	logmsg->emit_noendl(LOG_INFO,"Golizadeh phase particle conservation:    "); break;
	case SEtype_optical_phonon:		logmsg->emit_noendl(LOG_INFO,"Optical phonon particle conservation:     "); break;
	case SEtype_acoustic_phonon:	logmsg->emit_noendl(LOG_INFO,"Acoustic phonon particle conservation:    "); break;
	case SEtype_spont_photon:		logmsg->emit_noendl(LOG_INFO,"Spontaneous photon particle conservation: "); break;
	case SEtype_ion_imp:			logmsg->emit_noendl(LOG_INFO,"Ionized impurity particle conservation:   "); break;
	case SEtype_all:				logmsg->emit_noendl(LOG_INFO,"Overall particle conservation:            "); break;
	default: 						logmsg->emit_noendl(LOG_INFO,"Unknown scattering particle conservation: "); break;
	}
		
	uint Nk    = kspace->get_number_of_points();
	uint myNE  = energies->get_my_number_of_points();
	double kspace_factor = 1.0 / (negf_math::pow(2.0*constants::pi, kspace->get_dimension()));
	const double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);

	double espin = get_spin_degeneracy(opts->get("kp_method"), quantities::electron_density);
	double hspin = get_spin_degeneracy(opts->get("kp_method"), quantities::hole_density );
	vector<uint> ebands;
	opts->get_conduction_degrees_of_freedom(ebands);
	vector<uint> hbands;
	opts->get_valence_degrees_of_freedom(hbands);
	
	// own energies
	vector<cplx> trace; trace.resize(Nn, 0.0);
	for (uint ee2=0; ee2< myNE; ee2++) 
	{
		uint ee = energies->get_global_index(ee2);
		double dE = energies->get_weight_from_global_idx(ee)/(2.0*constants::pi);
		
		for (uint kk=0; kk<Nk; kk++) 
		{			
			double wk = kspace->get_point(kk).get_weight() * kspace_factor;

            // if there is only 1 k-point, assume ballistic calculation and parabolic bands
            // --> use Fermi-integral of order 0 with different EF's for the 2 contacts instead of k-weight
            // that was added to SigmaL
            if (Nk==1) {
                NEGF_ASSERT(Nn==1, "Nk=1 is possible only for single-band effective mass model.");
			    wk = 0.5; // spin is added later on
            }

			const GLMat & GL = gf->get_lesser(kk,ee);
			const GLMat & GG = gf->get_greater(kk,ee);
			const SEMat & SL = (type==SEtype_all) ? se->get_lesser(kk,ee) :  se->get_selfenergy(type)->get_lesser(kk,ee);
			const SEMat & SG = (type==SEtype_all) ? se->get_greater(kk,ee) : se->get_selfenergy(type)->get_greater(kk,ee);
		
			vector<cplx> SLGG; get_diag(SL, GG, SLGG);
			vector<cplx> SGGL; get_diag(SG, GL, SGGL);
			
			for (uint xx=1; xx<=Nx; xx++) {
				for (uint nn=1; nn<=Nn; nn++) {
					uint idx = get_mat_idx(xx,nn,Nx) - 1;
					trace[nn-1] += dE * wk * 1/hbar * (SLGG[idx] - SGGL[idx]);
				}
			}
		}
	}
	vector<double> trace_dbl; trace_dbl.resize(Nn,0.0);
	for (uint nn=0; nn<Nn; nn++) {
		NEGF_FASSERT(fabs(trace[nn].imag()) < constants::imag_err, "encountered complex trace: (%.4e,%.4e)", trace[nn].real(), trace[nn].imag());
		trace_dbl[nn] = trace[nn].real();
	}
	
	// communicate to master process
	if (i_am_master) {
		
		// collect the pieces
		for (int pp=0; pp<mpi->get_num_procs(); pp++) {
			if (pp==constants::mpi_master_rank) continue;
			vector<double> tmp_trace;
			int tag = pp;
			uint siz = Nn;
			mpi->recv(tmp_trace, siz, pp, tag);
			NEGF_ASSERT(tmp_trace.size()==Nn, "something went wrong.");
			// add to master density
			for (uint nn = 0; nn < Nn; nn++) {
				trace_dbl[nn] += tmp_trace[nn];
			}
		}
		
		// screen output
		//const double time_conv = constants::convert_from_SI(units::time, 1.0);   // 1s
		//const double      area = constants::convert_from_SI(units::area, 1e-18); // 1nm^2
		const double      coul = constants::convert_from_SI(units::charge, 1.0);   // 1C
		const double time_conv = constants::convert_from_SI(units::time, 1.0);     // 1s
		const double      area = constants::convert_from_SI(units::area, 1e-4);    // 1cm^2
		for (uint nn=0; nn<Nn; nn++) {
			double spin = 0.0;
			for (uint nn2=0; nn2<ebands.size(); nn2++) {
				if (ebands[nn2]==nn) {
					spin = espin;
					break;
				}
			}
			for (uint nn2=0; nn2<hbands.size(); nn2++) {
				if (hbands[nn2]==nn) {
					spin = hspin;
					break;
				}
			}
			NEGF_ASSERT(spin!=0.0, "could not determine spin.");
			//logmsg->emit(LOG_INFO,"n=%d: sum_k int dE 1/hbar sum_x (SLGG - SGGL)_xx = %e [s-1nm-2]", nn, trace_dbl[nn] * time_conv * area);
			//logmsg->emit_noendl(LOG_INFO,"%10.3e [s-1nm-2,n=%d]  ", trace_dbl[nn] * time_conv * area, nn);
			logmsg->emit_noendl(LOG_INFO,"%10.3e [A/cm2,n=%d]  ", spin * trace_dbl[nn] / coul * time_conv * area, nn);
		}
		logmsg->emit(LOG_INFO,"");
	} else {
		// send to master thread
		int dest = constants::mpi_master_rank;
		int tag = mpi->get_rank();
		mpi->send(trace_dbl, dest, tag);
	}
	mpi->synchronize_processes();			
	
);}


void PostProcessing::output_selfenergies(double Etarget, const char * filename)
{STACK_TRACE(
	// find energy index closest to E
	double dist = 1e100;
	uint Eidx = 0;
	for (uint ee=0; ee<energies->get_number_of_points(); ee++) {
		double E = energies->get_energy_from_global_idx(ee);
		if (fabs(E-Etarget) < dist) {
			dist = fabs(E-Etarget);
			Eidx = ee;
		}
	}
	logmsg->emit(LOG_INFO,"Output of self-energies associated to ee=%d (E=%.3f) which is closest to target %.3f", 
			Eidx, energies->get_energy_from_global_idx(Eidx), Etarget);
	
	int pp=energies->get_process_computing(Eidx);
	if (pp!=mpi->get_rank()) {
		return; // only the process computing energy index Eidx has to do something
	}
	
	uint Nk = kspace->get_number_of_points();
	double kspace_factor = 1.0/(4.0*constants::pi*constants::pi);
	char buf[1000];
	
	// ---------------------------------------
	// output contact self-energy
	// ---------------------------------------
	Matc SLcont(NxNn,NxNn);
	for (uint kk=0; kk<Nk; kk++) {
		add(se->get_selfenergy(SEtype_contact)->get_lesser(kk,Eidx), kspace->get_point(kk).get_weight()*kspace_factor, SLcont);
	}
	
	sprintf(buf,"%s.SLcont",filename);
	logmsg->emit(LOG_INFO,"   contact SL: %s",buf);
	string description1 = "Contact self-energy integrated over k (size NxNn*NxNn)";
	write_matrix(buf, SLcont, description1);
	
	// --------------------------------------
	// output acoustic phonon self-energy
	// --------------------------------------
	if (se->has_self_energy(SEtype_acoustic_phonon)) {
		Matc SLac(NxNn,NxNn);
		for (uint kk=0; kk<Nk; kk++) {
			add(se->get_selfenergy(SEtype_acoustic_phonon)->get_lesser(kk,Eidx), kspace->get_point(kk).get_weight()*kspace_factor, SLac);
		}
		
		sprintf(buf,"%s.SLac",filename);
		logmsg->emit_all(LOG_INFO,"   AC SL: %s",buf);
		string description = "AC self-energy integrated over k (size NxNn*NxNn)";
		write_matrix(buf, SLac, description);
	}
	
	// --------------------------------------
	// output optical phonon self-energy
	// --------------------------------------
	if (se->has_self_energy(SEtype_optical_phonon)) {
		Matc SLopt(NxNn,NxNn);
		for (uint kk=0; kk<Nk; kk++) {
			add(se->get_selfenergy(SEtype_optical_phonon)->get_lesser(kk,Eidx), kspace->get_point(kk).get_weight()*kspace_factor, SLopt);
		}
		
		sprintf(buf,"%s.SLopt",filename);
		logmsg->emit_all(LOG_INFO,"   POP SL: %s",buf);
		string description = "POP self-energy integrated over k (size NxNn*NxNn)";
		write_matrix(buf, SLopt, description);
	}
	
	// --------------------------------------
	// output Golizadeh momentum self-energy
	// --------------------------------------
	if (se->has_self_energy(SEtype_golizadeh_momentum)) {
		Matc SLgolim(NxNn,NxNn);
		for (uint kk=0; kk<Nk; kk++) {
			add(se->get_selfenergy(SEtype_golizadeh_momentum)->get_lesser(kk,Eidx), kspace->get_point(kk).get_weight()*kspace_factor, SLgolim);
		}
		
		sprintf(buf,"%s.SLgolim",filename);
		logmsg->emit_all(LOG_INFO,"   Golizadeh momentum SL: %s",buf);
		string description = "Golizadeh momentum self-energy integrated over k (size NxNn*NxNn)";
		write_matrix(buf, SLgolim, description);
	}
	
	// --------------------------------------
	// output Golizadeh momentum self-energy
	// --------------------------------------
	if (se->has_self_energy(SEtype_spont_photon)) {
		Matc SLphot(NxNn,NxNn);
		for (uint kk=0; kk<Nk; kk++) {
			add(se->get_selfenergy(SEtype_spont_photon)->get_lesser(kk,Eidx), kspace->get_point(kk).get_weight()*kspace_factor, SLphot);
		}
		
		sprintf(buf,"%s.SLphot",filename);
		logmsg->emit_all(LOG_INFO,"   Spontaneous photon SL: %s",buf);
		string description = "Spontaneous photon self-energy integrated over k (size NxNn*NxNn)";
		write_matrix(buf, SLphot, description);
	}
	
	// --------------------------------------
	// output ionized impurity self-energy
	// --------------------------------------
	if (se->has_self_energy(SEtype_ion_imp)) {
		Matc SLion(NxNn,NxNn);
		for (uint kk=0; kk<Nk; kk++) {
			add(se->get_selfenergy(SEtype_ion_imp)->get_lesser(kk,Eidx), kspace->get_point(kk).get_weight()*kspace_factor, SLion);
		}
		
		sprintf(buf,"%s.SLion",filename);
		logmsg->emit_all(LOG_INFO,"   Ionized impurity SL: %s",buf);
		string description = "Ionized impurity self-energy integrated over k (size NxNn*NxNn)";
		write_matrix(buf, SLion, description);
	}
);}


