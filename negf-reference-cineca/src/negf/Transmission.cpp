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
#include "Transmission.h"
using namespace negf;


Transmission::Transmission(const Geometry * xspace_, const Kspace * kspace_, 
			 const Energies * energies_, const GreenFunctions * gf_,
			 const SelfEnergy * se_contact_) throw (Exception *):
	xspace(xspace_),
	kspace(kspace_),
	energies(energies_),
	gf(gf_),
	se_contact(se_contact_)
{STACK_TRACE(
	NEGF_ASSERT(xspace!=0 && kspace!=0 && energies!=0 && gf!=0 && se_contact!=0, "null pointer encountered.");
	this->Nx   = xspace->get_num_internal_vertices();
	this->NxNn = Nx*Nn;
	this->Nk   = kspace->get_number_of_points();
	this->NE   = energies->get_number_of_points();
	this->myNE = energies->get_my_number_of_points();
	
	if (mpi->get_rank()==constants::mpi_master_rank) {
		this->transmission_spectrum  .resize(NE, 0.0);
		this->transmission_spectrum_k0.resize(NE, 0.0);
	} else {
		this->transmission_spectrum  .clear();
		this->transmission_spectrum_k0.clear();
	}
);}


void Transmission::calculate() throw (Exception *)
{STACK_TRACE(
	logmsg->emit_small_header("Calculating transmission T(E)");
	
	// calculate own points
	SEMat        SAmat = SCMat_create(NxNn);
	SEMat        Gamma = SCMat_create(NxNn);
	SEMat       Gamma1 = SCMat_create(NxNn);
	SEMat       Gamma2 = SCMat_create(NxNn);
	Matc        Gam1GA(NxNn,NxNn);
	Matc      GRGam1GA(NxNn,NxNn);
	Matc  Gam2GRGam1GA(NxNn,NxNn);
	Matc        result(NxNn,NxNn);
	Matc     result_k0(NxNn,NxNn);
	vector<double>    my_T(myNE, 0.0);
	vector<double> my_T_k0(myNE, 0.0);
	for (uint ee2=0; ee2<myNE; ee2++) 
	{
		double ee = energies->get_global_index(ee2);
		
		result = Matc(NxNn,NxNn);	// init to 0
		for (uint kk=0; kk<Nk; kk++) 
		{
			double wk = kspace->get_point(kk).get_weight(); // includes 2pi from angular integration but NOT 1/(2pi)^2
			wk = wk / (4.0*constants::pi*constants::pi);
			
			const Matc  & GRmat = gf->get_retarded(kk,ee);
			const Matc  & GAmat = gf->get_advanced(kk,ee);
			const SEMat & SRmat = se_contact->get_retarded(kk,ee);
			conjtrans(SRmat, SAmat);
			
			Gamma  = SRmat;
			Gamma -= SAmat;
			Gamma *= constants::imag_unit;
			
			// it should be possible to split Gamma into contributions from the left and right contact
			for (uint mm=1; mm<=Nn; mm++) {
				for (uint nn=1; nn<=Nn; nn++) {
					for (uint xx=1; xx<=Nx; xx++) {
						for (uint yy=1; yy<=Nx; yy++) {
							uint ii = get_mat_idx(xx,mm,Nx);
							uint jj = get_mat_idx(yy,nn,Nx);
#ifdef USE_BANDED
							if (fabs(ii-jj) > Gamma.num_offdiags+1e-8) continue;
#endif
							if (xx<Nx/2 && yy<Nx/2) {  // belongs to left contact
								Gamma1(ii,jj) = Gamma(ii,jj);
							} else if (xx>Nx/2 && yy>Nx/2) {  // belongs to right contact
								Gamma2(ii,jj) = Gamma(ii,jj);
							} else {
								// not clear to which contact it belongs, but should be 0 anyway
								NEGF_ASSERT(abs(Gamma(ii,jj))<1e-15, "entry should be zero");
							}
						}
					}
				}
			}
			
			// calculate transmission 
			mult(Gamma1,    GAmat,       Gam1GA);
			mult( GRmat,   Gam1GA,     GRGam1GA);
			mult(Gamma2, GRGam1GA, Gam2GRGam1GA);
			
			// k=0 only: without wk
			if (kk==0) {
				result_k0 = Gam2GRGam1GA;
			}
			
			// add to integral
			add(Gam2GRGam1GA, wk, result);
		}
		
		// take trace over space and band and add to array of transmission coefficients of local energies
		for (uint nn=1; nn<=Nn; nn++) {
			for (uint xx=1; xx<=Nx; xx++) {
				uint ii = get_mat_idx(xx,nn,Nx);
				
				cplx tmp = result(ii,ii);
				NEGF_FASSERT(fabs(tmp.imag()) < 1e-10, "imaginary transmission encountered: (%e,%e)", tmp.real(), tmp.imag());
				my_T[ee2] += tmp.real();
				
				cplx tmp_k0 = result_k0(ii,ii);
				NEGF_FASSERT(fabs(tmp_k0.imag()) < 1e-10, "imaginary transmission encountered (k=0): (%e,%e)", tmp_k0.real(), tmp_k0.imag());
				my_T_k0[ee2] += tmp_k0.real();
			}
		}
	}
	
	// communicate to master process
	mpi->synchronize_processes();
	if (mpi->get_rank()==constants::mpi_master_rank) 
	{
		this->transmission_spectrum.assign(NE, 0.0);
		
		// own part
		for (uint ee2 = 0; ee2 < myNE; ee2++) {
			uint ee = energies->get_global_index(ee2);
			this->transmission_spectrum   [ee] = my_T   [ee2];
			this->transmission_spectrum_k0[ee] = my_T_k0[ee2];
		}
		
		// collect the pieces
		for (int pp=0; pp<mpi->get_num_procs(); pp++) 
		{
			if (pp==constants::mpi_master_rank) continue;
			
			// receive from other process
			uint num_energy_points_pp = energies->get_number_of_points(pp);
			int tag = pp;
			vector<double> tmp_array(num_energy_points_pp, 0.0); 
			mpi->recv(tmp_array, num_energy_points_pp, pp, tag);
			vector<double> tmp_array_k0(num_energy_points_pp, 0.0); 
			mpi->recv(tmp_array_k0, num_energy_points_pp, pp, tag);
			
			// add to total matrix
			uint start_idx = energies->get_start_global_idx(pp) + 1;
			uint stop_idx  = energies->get_stop_global_idx(pp)  + 1;
			for (uint ee=start_idx; ee<=stop_idx; ee++) {
				this->transmission_spectrum   [ee] = tmp_array   [ee-start_idx];
				this->transmission_spectrum_k0[ee] = tmp_array_k0[ee-start_idx];
			}
		}
	} else {
		// send to master process
		int dest = constants::mpi_master_rank;
		int tag = mpi->get_rank();
		mpi->send(my_T   , dest, tag);
		mpi->send(my_T_k0, dest, tag);
	}
	mpi->synchronize_processes();
);}


void Transmission::write_to_file(const char * filename) throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(mpi->get_rank()==constants::mpi_master_rank, "transmission spectrum is stored in master process only!");
	const vector<double> & egrid = energies->get_energy_grid();
	Matd output(3, NE);
	for (uint ee=0; ee<NE; ee++) {
		output(1, ee+1) = egrid[ee];
		output(2, ee+1) = this->transmission_spectrum   [ee];
		output(3, ee+1) = this->transmission_spectrum_k0[ee];
	}
	negf::write_matrix(filename, output);
);}



