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
#include "ContactFermilevel.h"
using namespace negf;

ContactFermilevel::ContactFermilevel(const Energies * energies_, Hamiltonian * ham_, const Overlap * ov_,
									SEContacts * se_contacts_) throw (Exception *):
	ham(ham_),
	xspace(ham->get_xspace()),
	kspace(ham->get_kspace()),
	options(ham->get_options()),
	energies(energies_),
	ov(ov_),
	se_contacts(se_contacts_),
	Nx((xspace==0) ? 0 : xspace->get_num_internal_vertices()),
	NxNn(Nx*Nn)
{STACK_TRACE(
	NEGF_ASSERT(xspace!=NULL && kspace!=NULL && options!=NULL && energies!=NULL && ov!=NULL, "null pointer encountered.");
	NEGF_ASSERT(se_contacts!=NULL, "null pointer encountered.");
		
	// calculate
	//this->calculate_dos_contact0();
	this->calculate_contact0_dos_from_device_GR();
);}


// phi=0 in contact0 is assumed because the contact DOS is constructed from H(k) without phi
// most of this is copy-paste from SEContacts.cpp :-(
void ContactFermilevel::calculate_dos_contact0()
{STACK_TRACE(
	logmsg->emit_header("Calculating contact 0 density of states from Hamiltonian");
	uint Nvert = xspace->get_num_vertices();
	uint nE    = energies->get_number_of_points();
	uint myNE  = energies->get_my_number_of_points();
	uint Nk    = kspace->get_number_of_points();
	
	const OVMat & M = ov->get_overlap(); // defined on all indices (not just internal), (Nvert*Nn)^2
	NEGF_ASSERT(M.num_rows()==Nvert*Nn && M.num_cols()==Nvert*Nn, "inconsistent overlap matrix size.");
	
	// get the band indices in the Hamiltonian corresponding to conduction and valence bands
	vector<uint> cb_dofs;
	options->get_conduction_degrees_of_freedom(cb_dofs);
	vector<uint> vb_dofs;
	options->get_valence_degrees_of_freedom(vb_dofs);
	
	const vector<double> elstat_pot = ham->get_electrostatic_potential();
	vector<double> zero_potential;
	zero_potential.resize(Nvert, 0.0);
	ham->set_electrostatic_potential(zero_potential);
	
	this->contact_0_dos_n.clear();
	this->contact_0_dos_p.clear();
	
	// --------------------------------------
	// find DOS for own energy points
	// --------------------------------------
	logmsg->emit(LOG_INFO,"Computing DOS for individual energies...");
	vector<double> my_dos_n; my_dos_n.resize(myNE, 0.0);
	vector<double> my_dos_p; my_dos_p.resize(myNE, 0.0);
	
	for (uint ee2 = 0; ee2 < myNE; ee2++)
	{	
		uint ee = energies->get_global_index(ee2);
		logmsg->emit_noendl_all(LOG_INFO_L2, "p%d: SR_cont(E=%d,k=:)...  ",mpi->get_rank(),ee);
		for (uint kk = 0; kk < Nk; kk++)
		{	
			// get Hamiltonian including energy from transversal k-vector BUT excluding any potential!
			Matc H(Nvert*Nn,Nvert*Nn);
			ham->get(kspace->get_point(kk), H);	
			NEGF_ASSERT(H.num_rows()==Nvert*Nn && H.num_cols()==Nvert*Nn, "inconsistent Hamiltonian size.");
			double wk = kspace->get_point(kk).get_weight();

            // if there is only 1 k-point, assume ballistic calculation and parabolic bands
            // --> use Fermi-integral of order 0 with different EF's for the 2 contacts instead of k-weight
            // that was added to SigmaL
            if (Nk==1) {
                NEGF_ASSERT(Nn==1, "Nk=1 is possible only for single-band effective mass model.");
                wk = 0.5; // spin is added later on
            }

			// contact 0 only!
			
			// ---------------------------------------------------------------------------------------------------------
			// construct the matrix Hc of size 2NcNn*2NcNn, Nc = number of vertices at the interface contact-device
			// Hc(      1:NcNn,      1:NcNn) will be denoted H00
			// Hc(NcNn+1:2NcNn,      1:NcNn) will be denoted H10
			// Hc(      1:NcNn,NcNn+1:2NcNn) will be denoted H01
			// Hc(NcNn+1:2NcNn,NcNn+1:2NcNn) will be denoted H11
			// ---------------------------------------------------------------------------------------------------------
			uint Nc = se_contacts->get_interface_vertices(0).size();
			Matc Hc(2*Nc*Nn, 2*Nc*Nn);
			Matc Mc(2*Nc*Nn, 2*Nc*Nn);
			for (uint ii=1; ii<=Nc; ii++) 
			{
				uint gii1 = se_contacts->get_interface_vertices(0)[ii-1]->get_index_global() + 1;
				uint gii2 = se_contacts->get_second_row_vertices(0)[ii-1]->get_index_global() + 1;
				
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
				/*
				Hc.fill_block(   ii,    ii, H, gii1, gii1, 2*Nc);	// 2*Nc --> target matrix has 2*Nc entries, not Nx
				Hc.fill_block(   ii, Nc+ii, H, gii1, gii2, 2*Nc);
				Hc.fill_block(Nc+ii,    ii, H, gii2, gii1, 2*Nc);
				Hc.fill_block(Nc+ii, Nc+ii, H, gii2, gii2, 2*Nc);
				
				Mc.fill_block(   ii,    ii, M, gii1, gii1, 2*Nc);
				Mc.fill_block(   ii, Nc+ii, M, gii1, gii2, 2*Nc);
				Mc.fill_block(Nc+ii,    ii, M, gii2, gii1, 2*Nc);
				Mc.fill_block(Nc+ii, Nc+ii, M, gii2, gii2, 2*Nc);*/
			}			
			Matc Hc_00(Nc*Nn,Nc*Nn); Hc.get_submatrix(      1,   Nc*Nn,       1,   Nc*Nn, Hc_00);	
			Matc Hc_01(Nc*Nn,Nc*Nn); Hc.get_submatrix(      1,   Nc*Nn, Nc*Nn+1, 2*Nc*Nn, Hc_01);
			Matc Hc_10(Nc*Nn,Nc*Nn); Hc.get_submatrix(Nc*Nn+1, 2*Nc*Nn,       1,   Nc*Nn, Hc_10);
			Matc Hc_11(Nc*Nn,Nc*Nn); Hc.get_submatrix(Nc*Nn+1, 2*Nc*Nn, Nc*Nn+1, 2*Nc*Nn, Hc_11);
			
			Matc Mc_00(Nc*Nn,Nc*Nn); Mc.get_submatrix(      1,   Nc*Nn,       1,   Nc*Nn, Mc_00);	
			Matc Mc_01(Nc*Nn,Nc*Nn); Mc.get_submatrix(      1,   Nc*Nn, Nc*Nn+1, 2*Nc*Nn, Mc_01);
			Matc Mc_10(Nc*Nn,Nc*Nn); Mc.get_submatrix(Nc*Nn+1, 2*Nc*Nn,       1,   Nc*Nn, Mc_10);
			
			double dx = se_contacts->get_dx(0);
			
			// ----------------------------------------------
			// construct D  = E*Mc_00-Hc_00
			//           Tu = E*Mc_01-Hc_01
			//           Tl = E*Mc_10-Hc_10
			// ----------------------------------------------
			double E = energies->get_energy_from_global_idx(ee);
			Matc D(Nc*Nn, Nc*Nn);
			mult(Mc_00, E, D); // D  = E*Mc_00;
			D -= Hc_00;
			Matc Tu(Nc*Nn, Nc*Nn);
			mult(Mc_01, E, Tu); // Tu  = E*Mc_01;
			Tu -= Hc_01;
			Matc Tl(Nc*Nn, Nc*Nn);
			mult(Mc_10, E, Tl); // Tl  = E*Mc_10;
			Tl -= Hc_10;
			
			// ----------------------------------------------
			// construct A = -Tu^-1*D
			// ----------------------------------------------
			// find tmp=Tu^-1
			Matc tmp(Nc*Nn, Nc*Nn);
			tmp = Tu;
			invert(Tu);
			
			// multiply, negative
			Matc A(Nc*Nn, Nc*Nn);
			mult(tmp, D, A); // A = tmp * tmp2_c;
			A *= -1;
			
			// ----------------------------------------------
			// solve eigenvalue problem A*psi=lambda*psi
			// ----------------------------------------------
			uint n = A.num_rows();
			cplx AA[n*n];				// SOOOOOOO INEFFICIENT
			for (uint ii=1; ii <= n; ii++) {
				for (uint jj=1; jj <= n; jj++) {
					AA[(jj-1)*n+(ii-1)] = A(ii,jj);
				}
			}
			cplx lambda[n];			// eigenvalues
			cplx vr[n*n];
			negf_math::zgeev(n, AA, lambda, vr); 
			vector< vector<cplx> > psi;
			psi.resize(2*n);
			for (uint ii=0; ii<n; ii++) {
				psi[2*ii]  .resize(n);				// psi[2*ii]   stores eigenvector ii      --> +k
				psi[2*ii+1].resize(n);				// psi[2*ii+1] also stores eigenvector ii --> -k
				for (uint jj=0; jj<n; jj++) {
					psi[2*ii]  [jj] = vr[ii*n+jj];	// vr(jj,ii)
					psi[2*ii+1][jj] = vr[ii*n+jj];	// vr(jj,ii)
				}
			}
			
			// ------------------------------------------------------------------------
			// get the k-vectors (in transport direction) from the eigenvalues lambda
			// using analytic continuation of cosine
			// there are always 2 solutions: k and -k!
			// ------------------------------------------------------------------------
			vector<cplx> k;
			k.resize(2*n, 0.0); 
			for (uint ii=0; ii < n; ii++) 
			{
				cplx z = lambda[ii]/2.0;
				k[ii*2]   = +1.0/dx * negf_math::acos(z);	// dx>0
				k[ii*2+1] = -1.0/dx * negf_math::acos(z);
			}
			
			// ------------------------------------------------------------------------
			// get dE/dk (in transport direction) 
			// (H_00-E(k)M_00 - 2cos(k*dx) T_u) psi = 0, T_u = EM_01 - H_01
			//  --> (-dE/dk M_00 + 2*dx*sin(k*dx) T_u - 2cos(k*dx) dE/dk M_01) psi = 0
			//  -->   dE/dk (M_00 + 2cos(k*dx) M_01) = 2*dx*sin(k*dx) * T_u
			//  -->   dE/dk = (M_00 + 2cos(k*dx) M_01)^-1 * (2*dx*sin(k*dx) * T_u)
			//              = dEdk_helper2                * dEdk_helper1
			// ------------------------------------------------------------------------
			
			vector<cplx> dEdk; dEdk.resize(k.size(), 0.0);
			Matc dEdk_helper1(Nc*Nn, Nc*Nn);
			Matc dEdk_helper2(Nc*Nn, Nc*Nn);
			Matc dEdk_helper3(Nc*Nn, Nc*Nn);
			for (uint ii=0; ii < k.size(); ii++) 
			{
				dEdk_helper1 = Tu;
				dEdk_helper1 *= (+2.0*dx*sin(k[ii]*dx));
				
				dEdk_helper2 = Mc_01;
				dEdk_helper2 *= (+2.0*cos(k[ii]*dx));
				dEdk_helper2 += Mc_00;
				invert(dEdk_helper2);
				
				mult(dEdk_helper2, dEdk_helper1, dEdk_helper3); // h3 = h2*h1
								
				// now we have dEdk_helper3 * psi = dE/dk * psi
				// we find the scalar dE/dk by calculating the LHS and comparing it to the RHS whenever psi_i!=0
				vector<cplx> dEdk_psi;
				dEdk_psi.resize(Nc*Nn, 0.0);
				for (uint jj=0; jj<Nc*Nn; jj++) {
					for (uint ll=0; ll<Nc*Nn; ll++) {
						dEdk_psi[jj] += dEdk_helper3(jj+1,ll+1) * psi[ii][ll];
					}
				}
				uint num_candidates = 0;
				//if (kk==0 && cc==0 && ee%20==0) cout << "   candidates for dE/dk: ";
				for (uint jj=0; jj<Nc*Nn; jj++) {
					if (abs(psi[ii][jj]) > 1e-13) {
						dEdk[ii] = dEdk_psi[jj] / psi[ii][jj];
						//if (kk==0 && cc==0 && ee%20==0) cout << dEdk[ii] << "   ";
						//if (kk==0 && ee%20==0) cout << "E=" << E.real() << ", kt=0, cc=" << cc << ", n=" << ii << ", k=" << k[ii] << ", dEdk=" << dEdk[ii] << endl;
						num_candidates++;
					}
				}
				//if (kk==0 && cc==0 && ee%20==0) cout << endl;
				NEGF_ASSERT(num_candidates>0, "is the eigenvector entirely 0???");
				NEGF_ASSERT(num_candidates==1, "there seemed to be more than one possibility for dE/dk!");
			}
			
			// ----------------------------------------------
			// filter out unwanted k-vectors
			// ----------------------------------------------
			vector<cplx> k_filtered;
			vector< vector<cplx> > psi_filtered;
			const double eps = 1e-13;
			for (uint ii=0; ii < k.size(); ii++) {
				const cplx & num = /*k[ii]*/ dEdk[ii];
				bool filter_criterion = false;
				if (true) {
					filter_criterion =  /* (fabs(num.real())<eps &&      num.imag() <-eps)  // state decays to the right
								 	  ||*/ (     num.real() >eps && fabs(num.imag())<eps); // state propagates to the right
				}
				if (false) {
					filter_criterion =   (fabs(num.real())<eps &&      num.imag() >eps)  // state decays to the left
								 	  || (     num.real() <eps && fabs(num.imag())<eps); // state propagates to the left
				}
				//if (fabs(num.imag())<eps) filter_criterion = true; // state must be propagating (for DOS)
				if (filter_criterion) {
					k_filtered.push_back(k[ii]);
					psi_filtered.push_back(psi[ii]);
					//if (kk==0 && ee%20==0) cout << "E=" << E.real() << ", kt=0, cc=" << cc << ", n=" << ii << " filtered: k=" << k[ii] << ", dEdk=" << dEdk[ii] << endl;
				} else {
					// throw away
				}
			}
			uint Nkx = k_filtered.size(); 
			if (kk==0) {
				logmsg->emit_all(LOG_INFO/*_L2*/,"E=%g: %d out of %d k-vectors remain after filtering.", E, Nkx, k.size());
			}
			if (Nkx==0) {
				continue;
			}
			
			// -----------------------------------------------
			// construct the matrix P(-1)
			// -----------------------------------------------
			Matc Pm1(Nkx,Nkx);
			for (uint ii=0; ii<Nkx; ii++) {
				Pm1(ii+1,ii+1) = exp(-constants::imag_unit*k_filtered[ii]*(-dx)); 
			}
			
			// -----------------------------------------------
			// construct the matrix P(1)
			// -----------------------------------------------
			Matc Pp1(Nkx,Nkx);
			for (uint ii=0; ii<Nkx; ii++) {
				Pp1(ii+1,ii+1) = exp(+constants::imag_unit*k_filtered[ii]*(-dx));
			}
			
			// -----------------------------------------------------------
			// construct the matrix Phi: columns of Phi are eigenvectors
			// -----------------------------------------------------------
			Matc Phi(Nc*Nn,Nkx);
			for (uint ii=0; ii<Nc*Nn; ii++) {
				for (uint ll=0; ll<Nkx; ll++) {
					Phi(ii+1,ll+1) = psi_filtered[ll][ii];
				}
			}
		
			Matc PhiT(Nkx,Nc*Nn);
			trans(Phi, PhiT);
			
			// ---------------------------------------------------------------------------
			// construct the matrix g_tilde^-1 = PhiT*(D*Phi + Tu*Phi*P(-1) + Tl*Phi*P(+1))
			//                                 = PhiT*(D*Phi + Tu*tmp3      + Tl*tmp20)
			//                                 = PhiT*(D*Phi + tmp4         + tmp21)
			//                                 = PhiT*(D*Phi + tmp4         + tmp21)
			//                                 = PhiT*tmp5
			// ---------------------------------------------------------------------------
			
			Matc tmp3(Nkx,Nkx);
			mult(Phi,Pm1,tmp3); // tmp3 = Phi*Pm1;
			
			Matc tmp4(Nkx,Nkx);
			mult(Tu,tmp3,tmp4); // tmp4 = Tu*tmp3;
			
			Matc tmp20(Nkx,Nkx);
			mult(Phi,Pp1,tmp20); // tmp20 = Phi * Pp1;
			
			Matc tmp21(Nkx,Nkx);
			mult(Tl,tmp20,tmp21); // tmp21 = Tl * tmp20;
						
			Matc tmp5(Nkx,Nkx);
			mult(D,Phi,tmp5); // tmp5 = D*Phi;
			tmp5 += tmp4; 
			//tmp5 += tmp21; 			// NEW
			
			Matc g_tilde(Nkx,Nkx);
			mult(PhiT,tmp5,g_tilde); // g_tilde = PhiT * tmp5;
			
			// -----------------
			// invert
			// -----------------
			invert(g_tilde);
			
			// -----------------------------------------------
			// construct the contact GF
			// -----------------------------------------------
			// g00 = Phi * g_tilde * PhiT 
			//     = Phi * g_tilde * tmp6
			//     = Phi * tmp7	
			
			Matc tmp6(Nc*Nn,Nc*Nn);
			tmp6 = PhiT;
			
			Matc tmp7(Nc*Nn,Nc*Nn);	
			mult(g_tilde, tmp6, tmp7); // tmp7 = g_tilde * tmp6;
			
			Matc g00(Nc*Nn,Nc*Nn);
			mult(Phi, tmp7, g00); // g00 = Phi_c * tmp7;
			
			// -----------------------------------------------
			// add diagonal to the density of states (GR-GA)
			// copy-paste from PostProcessing.cpp
			// -----------------------------------------------
			for (uint nn=1; nn<=Nn; nn++) 
			{
				bool electron_band = false;
				for (uint jj=0; jj<cb_dofs.size(); jj++) {
					if(cb_dofs[jj]+1==nn) {
						electron_band = true;
						break;
					}
				}
				double spin = 0.0;
				if (electron_band) {
					spin = get_spin_degeneracy(options->get("kp_method"), quantities::electron_density);
				} else {
					spin = get_spin_degeneracy(options->get("kp_method"), quantities::hole_density);
				}
				for (uint jj = 1; jj <= Nc; jj++) {
					// ordering in g00 is always of type (xx-1)*Nn+nn!
					uint idx = (nn-1)*Nc+jj;
					int sign = +1;
					// LDOS = i/2pi * (GR-GA)
					if (electron_band) {
						my_dos_n[ee2] += sign * wk * (-1.0/constants::pi) * spin * g00(idx,idx).imag();
						logmsg->emit_all(LOG_INFO_L3,"E=%g, kk=%d: g00.imag=%e, wk=%e, spin=%g, my_dos[%d]=%e",E,kk,g00(idx,idx).imag(),wk,spin,ee2,my_dos_n[ee2]);
					} else {
						my_dos_p[ee2] += sign * wk * (-1.0/constants::pi) * spin * g00(idx,idx).imag();
					}
				}
			}
		}
	}
	mpi->synchronize_processes();
	
	// -----------------------------------------------------------------------
	// communicate own DOS to master thread, which puts everything together
	// -----------------------------------------------------------------------
	logmsg->emit(LOG_INFO_L1,"Aggregating total DOS...");
	if (mpi->get_rank()==constants::mpi_master_rank)
	{
		const double factor = 2.0;	// +k and -k
		this->contact_0_dos_n.resize(nE, -1.0);
		this->contact_0_dos_p.resize(nE, -1.0);
		for (int pp=0; pp<mpi->get_num_procs(); pp++) {
			if (pp==mpi->get_rank()) {
				for (uint ee2=0; ee2 < myNE; ee2++) {
					uint ee = energies->get_global_index(ee2);
					this->contact_0_dos_n[ee] = factor * my_dos_n[ee2];
					this->contact_0_dos_p[ee] = factor * my_dos_p[ee2];
				}
			} else {
				uint num_energy_points_pp = energies->get_number_of_points(pp);
				vector<double> ndos_pp;	ndos_pp.resize(num_energy_points_pp, 0.0);
				vector<double> pdos_pp;	pdos_pp.resize(num_energy_points_pp, 0.0);
				int tag = 987;
				int source = pp;
				mpi->recv(ndos_pp, num_energy_points_pp, source, tag);
				mpi->recv(pdos_pp, num_energy_points_pp, source, tag);
				NEGF_ASSERT(num_energy_points_pp==energies->get_stop_global_idx(pp)-energies->get_start_global_idx(pp)+1, "inconsistent number of points.");
				for (uint ee=energies->get_start_global_idx(pp); ee<=energies->get_stop_global_idx(pp); ee++) {
					this->contact_0_dos_n[ee] = factor * ndos_pp[ee-energies->get_start_global_idx(pp)];
					this->contact_0_dos_p[ee] = factor * pdos_pp[ee-energies->get_start_global_idx(pp)];
				}
			}
		}
		// check if every point was received
		
		logmsg->emit(LOG_INFO_L1,"TOTAL DOS [ContactFermilevel::calculate_dos_contact0()]:");;
		for (uint ee=0; ee < contact_0_dos_n.size(); ee++) {
			NEGF_ASSERT(contact_0_dos_n[ee]!=-1.0, "did not find an energy point.");
			NEGF_ASSERT(contact_0_dos_p[ee]!=-1.0, "did not find an energy point.");
			logmsg->emit_noendl(LOG_INFO_L1,"%g   ",contact_0_dos_n[ee]);
		}
		logmsg->emit_noendl(LOG_INFO_L1,"");
		
		/*
		cout << "\nANALYTICAL: " << endl;
		const PropertyContainer<double> * mat = xspace->get_contact(0)->get_adjacent_region()->get_material();
		const double me = constants::convert_from_SI(units::mass, mat->get("electron_effective_mass") * constants::SIm0);
		const double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
		const double T    = this->options->get("temperature");
		const double Ec   = TdkpInfoDesk::get_cbedge(mat, T, ham->get_material_db());
		//const double kmax = kspace->get_kmax();
		//const double Ekmax = hbar*hbar*kmax*kmax / (2.0*me);
		for (uint ee=0; ee < contact_0_dos_n.size(); ee++) {
			double E = energies->get_energy_from_global_idx(ee);
			double Elong = E - Ec;
			double DOS_factor = 2.0 * negf_math::pow(me, 1.5) / (sqrt(2.0)*constants::pi*constants::pi*hbar*hbar);
			if (E<Ec) {
				cout << 0 << "   ";
			} else {
				cout << DOS_factor * sqrt(Elong) << "   ";
				//cout << DOS_factor * 2.0*constants::pi*2.0*me/(3.0*hbar*hbar) * (negf_math::pow(Ekmax + Elong, 1.5) - negf_math::pow(Elong, 1.5)) << "   ";
			}
		}
		cout << endl;*/
	} else {
		int tag = 987;
		int dest = constants::mpi_master_rank;
		mpi->send(my_dos_n, dest, tag);
		mpi->send(my_dos_p, dest, tag);
		this->contact_0_dos_n.resize(nE, -1.0);
		this->contact_0_dos_p.resize(nE, -1.0);
	}
	mpi->synchronize_processes();
	
	// -----------------------------------------------------------------------
	// communicate total DOS to all processes
	// -----------------------------------------------------------------------
	logmsg->emit(LOG_INFO_L1,"Communicating total DOS...");
	int root = constants::mpi_master_rank;
	mpi->broadcast(contact_0_dos_n, root);
	mpi->broadcast(contact_0_dos_p, root);
	
	mpi->synchronize_processes();
	
	if (elstat_pot.size()>0) {
		ham->set_electrostatic_potential(elstat_pot);	// restore original potential
	}
);}


void ContactFermilevel::calculate_contact0_dos_from_device_GR()
{STACK_TRACE(
	logmsg->emit_small_header("Calculating contact 0 density of states from GR");
	uint Nvert = xspace->get_num_vertices();
	uint nE    = energies->get_number_of_points();
	uint myNE  = energies->get_my_number_of_points();
	uint Nk    = kspace->get_number_of_points();
	double kspace_factor = 1.0 / (negf_math::pow(2.0*constants::pi, kspace->get_dimension()));
	
	const OVMat & M = ov->get_internal_overlap();
	NEGF_ASSERT(M.num_cols()==NxNn && M.num_rows()==M.num_cols(), "wrong overlap matrix.");
	
	// get the band indices in the Hamiltonian corresponding to conduction and valence bands
	vector<uint> cb_dofs;
	options->get_conduction_degrees_of_freedom(cb_dofs);
	vector<uint> vb_dofs;
	options->get_valence_degrees_of_freedom(vb_dofs);
	
	const vector<double> elstat_pot = ham->get_electrostatic_potential();
	vector<double> zero_potential;
	zero_potential.resize(Nvert, 0.0);
	ham->set_electrostatic_potential(zero_potential);
	
	//---------------------------------------------
	// calculate contact SE w/o potential
	// --------------------------------------------
	logmsg->emit(LOG_INFO,"Calculating retarded contact self-energy w/o potential");
	se_contacts->calculate_retarded();	
        std::cout<<"before synchronization...S.Z."<<std::endl;
	mpi->synchronize_processes();
	
	//---------------------------------------------
	// calculate GR w/o potential
	// --------------------------------------------
	logmsg->emit(LOG_INFO,"Calculating GR w/o potential");
	
	// initialize helper matrices
	Matc  Hsmall(NxNn,NxNn);
	OVMat EM = OVMat_create(NxNn);
	
	// small numerical parameter to make E a little bit complex
	// if omitted, NaN's will appear sometimes
	const double eta = constants::convert_from_SI(units::energy, constants::dyson_eta * constants::SIec);
	
	vector<double> my_dos_n; my_dos_n.resize(myNE, 0.0);
	vector<double> my_dos_p; my_dos_p.resize(myNE, 0.0);
	
	for (uint ee2 = 0; ee2 < myNE; ee2++)
	{	
		uint ee = energies->get_global_index(ee2);
		double E = energies->get_energy_from_global_idx(ee);
		mult(M, E+eta, EM); // EM = (E+eta) * M;
		
		if (ee % /*5*/10 == 0) logmsg->emit_noendl_all(LOG_INFO, "p%d: GR(E=%d,:)...   ",mpi->get_rank(),ee);
		for (uint kk = 0; kk < Nk; kk++)
		{				
			// get Hamiltonian of the interior points only
			ham->get_internal(kspace->get_point(kk), Hsmall);
			NEGF_ASSERT(Hsmall.num_rows()==Nx*Nn, "something went wrong.");
			
			// get retarded GF and self-energy
			Matc GR(NxNn,NxNn);
			const SEMat & SigmaR = se_contacts->get_retarded(kk,ee);
			
			// --------------------------
			// solve Dyson equation
			// --------------------------
			
			// assign GR = (E+eta)M - H - SigmaR
			mult(SigmaR, -1.0, GR); // GR = -SigmaR;
			GR -= Hsmall;
			GR += EM;
			
			// invert
			invert(GR);			
			
			// -----------------------------------------------
			// add diagonal to the density of states (GR-GA)
			// copy-paste from PostProcessing.cpp
			// -----------------------------------------------
			double wk = kspace->get_point(kk).get_weight() * kspace_factor;
			for (uint nn=1; nn<=Nn; nn++) 
			{
				bool electron_band = false;
				for (uint jj=0; jj<cb_dofs.size(); jj++) {
					if(cb_dofs[jj]+1==nn) {
						electron_band = true;
						break;
					}
				}
				double spin = 0.0;
				if (electron_band) {
					spin = get_spin_degeneracy(options->get("kp_method"), quantities::electron_density);
				} else {
					spin = get_spin_degeneracy(options->get("kp_method"), quantities::hole_density);
				}
				uint xx = 1;
				uint idx = get_mat_idx(xx,nn,Nx); // (xx-1)*Nn+nn;
				// LDOS = i/2pi * (GR-GA)
				if (electron_band) {
					my_dos_n[ee2] += wk * (-1.0/constants::pi) * spin * GR(idx,idx).imag();
					if (my_dos_n[ee2]>0.05 && fabs(GR(idx,idx).imag()) > 1e-3) {
						logmsg->emit_all(LOG_INFO_L3,"E=%g, kk=%d: GR.imag=%e, wk=%e, spin=%g, my_dos[%d]=%e",E,kk,GR(idx,idx).imag(),wk,spin,ee2,my_dos_n[ee2]);
					}
				} else {
					my_dos_p[ee2] += wk * (-1.0/constants::pi) * spin * GR(idx,idx).imag();
				}
			}
		}
		
		// correct with geometrical factor in the case of dumb orthogonal-basis-hamiltonian
		if (constants::old_orthogonal
			&& (fabs(options->get("kp_method") - 0.0) < 1e-14 || fabs(options->get("kp_method") - 3.0) < 1e-14)) {
			uint xx = 1;
			double x_cell_length = 0.0;
			Vertex * v = xspace->get_vertex(xspace->get_global_vertex_index(xx-1));
			for (uint ii=0; ii < xspace->get_edges_near(v).size(); ii++) {
				x_cell_length += 0.5 * xspace->get_edges_near(v)[ii]->get_length();
			}
			my_dos_n[ee2] = my_dos_n[ee2] / x_cell_length;
			my_dos_p[ee2] = my_dos_p[ee2] / x_cell_length;
		}
	}
	mpi->synchronize_processes();
	
	this->contact_0_dos_n.clear();
	this->contact_0_dos_p.clear();
	
	// -----------------------------------------------------------------------
	// communicate own DOS to master thread, which puts everything together
	// -----------------------------------------------------------------------
	logmsg->emit(LOG_INFO,"Aggregating total DOS...");
	if (mpi->get_rank()==constants::mpi_master_rank)
	{
		this->contact_0_dos_n.resize(nE, -1.0);
		this->contact_0_dos_p.resize(nE, -1.0);
		for (int pp=0; pp<mpi->get_num_procs(); pp++) {
			if (pp==mpi->get_rank()) {
				for (uint ee2=0; ee2 < myNE; ee2++) {
					uint ee = energies->get_global_index(ee2);
					this->contact_0_dos_n[ee] = my_dos_n[ee2];
					this->contact_0_dos_p[ee] = my_dos_p[ee2];
				}
			} else {
				uint num_energy_points_pp = energies->get_number_of_points(pp);
				vector<double> ndos_pp;	ndos_pp.resize(num_energy_points_pp, 0.0);
				vector<double> pdos_pp;	pdos_pp.resize(num_energy_points_pp, 0.0);
				int tag = 987;
				int source = pp;
				mpi->recv(ndos_pp, num_energy_points_pp, source, tag);
				mpi->recv(pdos_pp, num_energy_points_pp, source, tag);
				NEGF_ASSERT(num_energy_points_pp==energies->get_stop_global_idx(pp)-energies->get_start_global_idx(pp)+1, "inconsistent number of points.");
				for (uint ee=energies->get_start_global_idx(pp); ee<=energies->get_stop_global_idx(pp); ee++) {
					this->contact_0_dos_n[ee] = ndos_pp[ee-energies->get_start_global_idx(pp)];
					this->contact_0_dos_p[ee] = pdos_pp[ee-energies->get_start_global_idx(pp)];
				}
			}
		}
		// check if every point was received
		logmsg->emit(LOG_INFO_L2,"TOTAL n-DOS [ContactFermilevel::calculate_contact0_dos_from_device_GR()]: ");
		for (uint ee=0; ee < contact_0_dos_n.size(); ee++) {
			NEGF_ASSERT(contact_0_dos_n[ee]!=-1.0, "did not find an energy point.");
			logmsg->emit_noendl(LOG_INFO_L2,"%g   ",contact_0_dos_n[ee]);
		}
		logmsg->emit(LOG_INFO_L2,"");
		if (Nn>=2) {
			logmsg->emit(LOG_INFO_L2,"TOTAL p-DOS [ContactFermilevel::calculate_contact0_dos_from_device_GR()]: ");
			for (uint ee=0; ee < contact_0_dos_p.size(); ee++) {
				NEGF_ASSERT(contact_0_dos_p[ee]!=-1.0, "did not find an energy point.");
				logmsg->emit_noendl(LOG_INFO_L2,"%g   ",contact_0_dos_p[ee]);
			}
			logmsg->emit(LOG_INFO_L2,"");
		}
	} else {
		int tag = 987;
		int dest = constants::mpi_master_rank;
		mpi->send(my_dos_n, dest, tag);
		mpi->send(my_dos_p, dest, tag);
		this->contact_0_dos_n.resize(nE, -1.0);
		this->contact_0_dos_p.resize(nE, -1.0);
	}
	mpi->synchronize_processes();
	
	// -----------------------------------------------------------------------
	// communicate total DOS to all processes
	// -----------------------------------------------------------------------
	logmsg->emit(LOG_INFO_L1,"Communicating total DOS...");
	int root = constants::mpi_master_rank;
	mpi->broadcast(contact_0_dos_n, root);
	mpi->broadcast(contact_0_dos_p, root);
	
	mpi->synchronize_processes();
	
	if (elstat_pot.size()>0) {
		ham->set_electrostatic_potential(elstat_pot);	// restore original potential
	}	
);}


void ContactFermilevel::compute_contact_0_fermilevel(double doping) throw (Exception *)
{STACK_TRACE(
	char buf[1000]; sprintf(buf,"Calculating contact 0 fermilevel (doping=%.3e)",doping);
	logmsg->emit_small_header(buf);
	
	// --------------------------------------------------------------------------
	// find parabolic single band result
	// --------------------------------------------------------------------------
	Region * reg = xspace->get_contact(0)->get_adjacent_region();
	/*for (uint ii=0; ii< xspace->get_contact(0)->get_contact_vertices().size(); ii++) {
		if (xspace->get_regions_near(xspace->get_contact(0)->get_contact_vertex(ii)).size()==1) {
			reg = xspace->get_regions_near(xspace->get_contact(0)->get_contact_vertex(ii))[0];
		}
	}*/
	const PropertyContainer<double> * mat = reg->get_material();
	double T    = this->options->get("temperature");
	double kT   = constants::convert_from_SI(units::energy, constants::SIkb * T);
	double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
	double Ec = TdkpInfoDesk::get_cbedge(mat, T, ham->get_material_db());
	double Ev = constants::convert_from_SI(units::energy, mat->get("valence_band_edge") * constants::SIec);
	double me = constants::convert_from_SI(units::mass, mat->get("electron_effective_mass") * constants::SIm0);
	double Nc = 2.0 * negf_math::pow(me*kT/(2.0*constants::pi*hbar*hbar), 1.5);
	double mh = constants::convert_from_SI(units::mass, mat->get("hole_effective_mass") * constants::SIm0);
	double Nv = 2.0 * negf_math::pow(mh*kT/(2.0*constants::pi*hbar*hbar), 1.5);
	double mu = 0.0;
	if (doping > 0.0) {
		// ND ~ Nc * F_{0.5}((mu-Ec)/kT)  or   mu ~ Ec + kT * F_{0.5}^{-1}(ND/Nc)      (phi=0)
	   	mu = Ec + kT * negf_math::fermihalf_inverse(doping / Nc);
	} else {
		// NA ~ Nc * F_{0.5}((Ev-mu)/kT)  or   mu ~ Ev - kT * F_{0.5}^{-1}(NA/Nv)      (phi=0)
		mu = Ev - kT * negf_math::fermihalf_inverse(-doping / Nv);
	}
		
	// --------------------------------------------------------------------------
	// if kpmethod=0/1 (single-band effective mass) use analytical formula
	// --------------------------------------------------------------------------
	double kpmethod = options->get("kp_method");
	if ( ((fabs(kpmethod - 0.0) < 1e-14	|| fabs(kpmethod - 1.0) < 1e-14) && false)
	     || kspace->get_number_of_points()==1)
	{
		NEGF_ASSERT(doping>0.0, "must have positive doping in single-band effective mass case.");
		NEGF_ASSERT(get_spin_degeneracy(kpmethod,quantities::hole_density)==2.0, "hole band degeneracy of 2 expected.");
		NEGF_ASSERT(get_spin_degeneracy(kpmethod,quantities::electron_density)==2.0, "electron band degeneracy of 2 expected.");
		
		logmsg->emit(LOG_INFO,"Fermilevel at contact 0 lies at %.6e (Ec=%.6e, kT=%.3e, dop=%.2e, Nc=%.3e, F_0.5^-1()=%.3e",
				mu, Ec, kT, doping, Nc, negf_math::fermihalf_inverse(doping / Nc));
		
		this->contact_0_fermilevel = mu;
		return;
	} else 
	{		
		// -------------------------
		// small newton solver
		// -------------------------
		vector<double> n_p;
		n_p.resize(2, 0.0);
		
		this->calculate_contact_density_from_fermilevel(mu, false, n_p); // mu contains initial guess
		logmsg->emit(LOG_INFO,"   initial guess: mu=%.6e, n=%.3e, p=%.3e (Ec=%.6e, Ev=%.6e, dop=%.2e)", 
				mu, n_p[0], n_p[1], Ec, Ev, doping);
		
		// iterate!
		// find root of F(mu) = n(mu) - p(mu) - doping
		// --> mu(i+1) = mu(i) - dF_dmu(mu(i)) * F(mu(i))
		uint max_iterations = 100;
		uint iter;
		for (iter = 1; iter <= max_iterations; iter++)
		{
			// find F(mu)
			this->calculate_contact_density_from_fermilevel(mu, false, n_p);
			
			double F = n_p[0] - n_p[1] - doping;
			logmsg->emit(LOG_INFO," Iteration %d: F = %e(n) - %e(p) - %e(dop) = %e", iter, n_p[0], n_p[1], doping, F);
			// find dF_dmu(mu)
			this->calculate_contact_density_from_fermilevel(mu, true, n_p);
			double dF_dmu = n_p[0] - n_p[1];
			logmsg->emit(LOG_INFO_L2," Iteration %d: dF_dmu = dn_dmu - dp_dmu = %e - %e = %e", iter, n_p[0], n_p[1], dF_dmu);
			
			// find update
			double update = -1.0/dF_dmu * F;
			
			// limit update!
			double max_update = constants::convert_from_SI(units::energy, constants::max_kpquasi_change_V * constants::SIec);
			if (fabs(update) > max_update) {
				update = negf_math::sign(update) * max_update;
				logmsg->emit(LOG_INFO," Iteration %d: Limiting to %e", iter, update);
			}
			
			// peform update
			mu += update;
			
			// convergence check
			if (fabs(update) < constants::convert_from_SI(units::energy, 1e-9*constants::SIec)) {
				logmsg->emit(LOG_INFO,"Fermilevel convergence reached after %d iterations: mu=%.6g",iter, mu);
				break;
			} else {
				this->calculate_contact_density_from_fermilevel(mu, false, n_p);
				logmsg->emit(LOG_INFO," Iteration %d: mu_new=%e, n_new=%.3e, p_new=%.3e, dop=%.2e", iter, mu, n_p[0], n_p[1], doping);
			}
		}
		if (iter==max_iterations+1) {
			NEGF_FEXCEPTION("Convergence not reached after %d iterations: mu=%.3e, n=%.3e, p=%.3e", max_iterations, mu, n_p[0], n_p[1]);
		}
		
		this->contact_0_fermilevel = mu;
		return;
	}
);}


// attention: contact_0_dos contain the DOS for phi=0 in contact 0!
void ContactFermilevel::calculate_contact_density_from_fermilevel(double EF, bool calculate_derivative, vector<double> & n_p) const
{STACK_TRACE(
	n_p.resize(2);
	double & n = n_p[0];
	double & p = n_p[1];
	n = 0.0;
	p = 0.0;
	
	// get the band indices in the Hamiltonian corresponding to conduction and valence bands
	vector<uint> cb_dofs;
	options->get_conduction_degrees_of_freedom(cb_dofs);
	vector<uint> vb_dofs;
	options->get_valence_degrees_of_freedom(vb_dofs);
	
	double T = options->get("temperature");
	double kT   = constants::convert_from_SI(units::energy, constants::SIkb * T);
	double kpmethod = options->get("kp_method");
	
	uint nE = energies->get_number_of_points(); 
	for (uint ee=0; ee<nE; ee++) 
	{
		double E = energies->get_energy_from_global_idx(ee);
		double dE = energies->get_weight_from_global_idx(ee);
		double nu = (E - EF) / kT;

		// contribution to n
		//double spin = get_spin_degeneracy(kpmethod, quantities::electron_density); // electron band degeneracy
		double spin = 1.0;
		if (calculate_derivative) {
			double f = 1.0 / (2.0 + negf_math::exp(nu) + negf_math::exp(-nu)) * 1.0/kT;
			n += spin*f * contact_0_dos_n[ee] * dE;
		} else {
			double f = 1.0 / (1.0 + negf_math::exp(nu));
			n += spin*f * contact_0_dos_n[ee] * dE;
		}

		// contribution to p
		spin = get_spin_degeneracy(kpmethod, quantities::hole_density); // hole band degeneracy
		if (calculate_derivative) {
			double f = 1.0 / (2.0 + negf_math::exp(nu) + negf_math::exp(-nu)) * 1.0/kT;
			p += -spin*f * contact_0_dos_p[ee] * dE;
		} else {
			double f = 1.0 / (1.0 + negf_math::exp(nu));
			p += spin*(1.0-f) * contact_0_dos_p[ee] * dE;
		}
	}
	return;
);}

