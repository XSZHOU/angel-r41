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
#include "SEIonizedImpurities.h"
using namespace negf;

SEIonizedImpurities::SEIonizedImpurities(const Overlap * ov_,
					const Geometry * xspace_, 
					const Kspace * kspace_, 
					const Energies * energies_, 
					const Options * options_,
					const GreenFunctions * gf_) throw (Exception *):
	SelfEnergy(xspace_,kspace_,energies_,options_, constants::odSII),
	ov(ov_),
	gf(gf_),
	prepared(false),
	scaling(1.0),
	F_neglect(1e-10)
{STACK_TRACE(
	NEGF_ASSERT(ov!=NULL && xspace!=NULL && kspace!=NULL && energies!=NULL && options!=NULL && gf!=NULL, "null pointer encountered.");
	
	double Nl_dbl = Nx*(Nx+1.0)/2.0;
	NEGF_ASSERT(fabs(ceil(Nl_dbl)-Nl_dbl)<1e-14 && fabs(floor(Nl_dbl)-Nl_dbl)<1e-14, "Nl must be integer");
	this->Nl = uint(Nl_dbl);
	
	mpi->synchronize_processes();
);}


void SEIonizedImpurities::set_scaling(double new_scaling)
{STACK_TRACE(
	NEGF_ASSERT(new_scaling>=0.0 && new_scaling<=1.0, "bad scaling factor"); 
	logmsg->emit(LOG_INFO,"Setting scaling of ionized impurities scattering to %.4g", new_scaling);
	this->scaling = new_scaling;
);}


void SEIonizedImpurities::set_up(const vector<double> & doping,
                                 const vector<double> & dielectric,
                                 const vector<double> & Nc,
                                 const vector<double> & Nv) throw (Exception *)
{STACK_TRACE(
	logmsg->emit(LOG_INFO,"Setting up ionized impurity self-energy");
    NEGF_ASSERT(    doping.size()==xspace->get_num_vertices() &&
                dielectric.size()==xspace->get_num_elements() &&
                        Nc.size()==xspace->get_num_vertices() &&
                        Nv.size()==xspace->get_num_vertices(), "invalid number of variables.");
	
	vector<double>     dopingvec;     dopingvec.resize(xspace->get_num_internal_vertices(), 0.0);
	vector<double> dielectricvec; dielectricvec.resize(xspace->get_num_internal_vertices(), 0.0);
	vector<double>            nu;            nu.resize(xspace->get_num_internal_vertices(), 0.0);
	logmsg->emit_noendl(LOG_INFO_L3,"doping: ");
	for (uint xx=1; xx<=Nx; xx++) 
	{
		uint gx = xspace->get_global_vertex_index(xx-1);
		
		//dopingvec[xx-1] = doping->get_value(gx);
        dopingvec[xx-1] = doping[gx];
		logmsg->emit_noendl(LOG_INFO_L3,"%.1e    ",dopingvec[xx-1] * 1e21);
		    
		const vector<Element *> & elems_near_gx = xspace->get_elems_near(xspace->get_vertex(gx));
		double total_meas = 0.0;
		for (uint ii=0; ii<elems_near_gx.size(); ii++) {
			dielectricvec[xx-1] += elems_near_gx[ii]->get_edge(0)->get_length()
			     //   * dielectric->get_value(elems_near_gx[ii]->get_index_global());
                      * dielectric[elems_near_gx[ii]->get_index_global()];
			total_meas += elems_near_gx[ii]->get_edge(0)->get_length();
		}
		NEGF_ASSERT(total_meas > 0.0, "zero measure encountered.");
		dielectricvec[xx-1] /= total_meas;	
		//logmsg->emit_noendl(LOG_INFO,"%.1e    ",dielectricvec[xx-1]/constants::convert_from_SI(units::dielectric,constants::SIeps0));
		
		//if (doping->get_value(gx) > 0.0) {  //n-doping
        if (doping[gx] > 0.0) {  //n-doping
			//nu[xx-1] = fabs(doping->get_value(gx)) / Nc->get_value(gx);
		    nu[xx-1] = fabs(doping[gx]) / Nc[gx];
		} else { 							// p-doping
		    //nu[xx-1] = fabs(doping->get_value(gx)) / Nv->get_value(gx);
		    nu[xx-1] = fabs(doping[gx]) / Nv[gx];
		}
		logmsg->emit_noendl(LOG_INFO_L3,"(nu=%g)    ",nu[xx-1]);
	}
	logmsg->emit(LOG_INFO_L3,"");
	
	// Nl was set in constructor already
	this->calculate_F(dopingvec, dielectricvec, nu);
	
	this->prepared = true;
);}


/* Need to store ii entries in line ii. total required entries: \sum_{i=1}^{Nx} i = Nx*(Nx+1)/2
 * storage scheme: 
 * (1,1) (2,1) (2,2) (3,1) ... (m,n)
 * 0     1     2     3         m(m-1)/2 + n - 1
 * note that (m,n) and (n,m) entries are the same */
uint SEIonizedImpurities::find_F_index(uint xx, uint yy) const		// result is 0-based
{STACK_TRACE(	
	uint smaller = negf_math::min(xx,yy);
	uint larger  = negf_math::max(xx,yy);

	double l_dbl = larger*(larger-1.0)/2.0 + smaller - 1.0;	// -1.0 --> 0-based ll
	NEGF_FASSERT(fabs(ceil(l_dbl)-l_dbl)<1e-14 && fabs(floor(l_dbl)-l_dbl)<1e-14, "l must be integer, instead l=%.14e",l_dbl);
	uint ll = uint(l_dbl);
	NEGF_FASSERT(ll<Nl, "wrong l (ll=%d, Nl=%d; xx=%d, yy=%d, Nx=%d).",ll,Nl,xx,yy,Nx);
	return ll;
);}


void SEIonizedImpurities::calculate_F(const vector<double> & doping, 
									  const vector<double> & dielectric, 
									  const vector<double> & nu)
{STACK_TRACE(
	mpi->synchronize_processes();
	logmsg->emit_small_header("Calculating F");
	NEGF_ASSERT(    doping.size()==xspace->get_num_internal_vertices() && 
				dielectric.size()==xspace->get_num_internal_vertices(), "wrong handed over arrays doping/dielectric");
	const double ec = constants::convert_from_SI(units::charge, constants::SIec);
	
	
	// small test
	vector<bool> test; test.resize(Nl, true);
	for (uint xx=1; xx<=Nx; xx++) {
		for (uint yy=1; yy<=xx; yy++) {
			uint ll = find_F_index(xx,yy);
			NEGF_ASSERT(test[ll], "ll was already used");
			test[ll] = false;
		}
	}
	for (uint ll=0; ll<Nl; ll++) {
		NEGF_FASSERT(!test[ll], "did not use ll=%d.",ll);
	}
	
	// -----------------------------------------
	// precompute some stuff
	// -----------------------------------------
	logmsg->emit(LOG_INFO,"Precomputations...");
	
	const uint Nphi = 100;
	
	// set up Debye inverse screening length - is function of density, which is approximated as doping
	const double q0_min = constants::convert_from_SI(units::density_1d, constants::imp_q0min); // minimum inverse screening length!
	vector<double> q0; q0.resize(Nx,0.0);
	const double kT = constants::convert_from_SI(units::energy, constants::SIkb * options->get("temperature"));
	logmsg->emit_noendl(LOG_INFO,"q0: ");
	for (uint xx=1; xx<=Nx; xx++) {
		double dop =     doping[xx-1];
		double eps = dielectric[xx-1];
		
		// Fermi correction nu / F_{-0.5}(F_{0.5}^{-1}(nu))
		//double corr = 1.0;	// nondegenerate result overestimates the screening length at high densities - scattering is weaker including the fermi-correction
		double corr = nu[xx-1] / negf_math::fermi_int(-0.5, negf_math::fermihalf_inverse(nu[xx-1]));
		
		double debye_screening_length = negf_math::sqrt(eps*kT/(ec*ec*fabs(dop)*corr));	
		
		//q0[xx-1] = constants::convert_from_SI(units::density_1d, 5e8); // ad hoc!
		q0[xx-1] = max(q0_min, 1.0 / debye_screening_length);
		logmsg->emit_noendl(LOG_INFO,"%g    ",q0[xx-1]);
	}
	logmsg->emit(LOG_INFO,"");
	
	// set up length for which each vertex stands
	vector<double> dx; dx.resize(Nx, 0.0);
	logmsg->emit_noendl(LOG_INFO_L3,"dx: ");
	for (uint xx=1; xx<=Nx; xx++) {
		uint  gx = xspace->get_global_vertex_index(xx-1);
		uint gx0 = (gx>0) ? gx-1 : gx;
		uint gx2 = (gx<xspace->get_num_vertices()-1) ? gx+1 : gx;
		dx[xx-1] = 0.5*xspace->get_distance(xspace->get_vertex(gx0),xspace->get_vertex(gx2));
		logmsg->emit_noendl(LOG_INFO_L3,"%g    ",dx[xx-1]);
	}
	logmsg->emit(LOG_INFO_L3,"");
	
	// set up prefactor (in which F is linear)
	double cheat = 1.0;
	if (options->exists("IonizedImpurityCheatFactor")) {
		cheat = options->get("IonizedImpurityCheatFactor");
		NEGF_ASSERT(cheat>=0.0 && cheat <= 100000.0, "invalid ionized impurity cheat factor");
	}
	logmsg->emit_noendl(LOG_INFO_L3,"prefactor: ");
	vector<double> fac; fac.resize(Nx, 0.0);
	for (uint xx=1; xx<=Nx; xx++) {
		double dop =     doping[xx-1];
		double eps = dielectric[xx-1];
		
		fac[xx-1] = ec*ec*ec*ec/(32.0*constants::pi*constants::pi*eps*eps) * fabs(dop) * dx[xx-1] * cheat;
		logmsg->emit_noendl(LOG_INFO_L3,"%g    ",fac[xx-1]);
	}
	logmsg->emit(LOG_INFO_L3,"");
	
	// set up array storing distance term involving xi, xj, xc
	vector<double> tmp0; tmp0.resize(Nx, 0.0);
	vector< vector<double> > tmp1; tmp1.resize(Nx,tmp0);
	vector< vector< vector<double> > > gx_gy_gc_dist;  gx_gy_gc_dist.resize(Nx,tmp1);
	for (uint xx=1; xx<=Nx; xx++) {
		uint gx = xspace->get_global_vertex_index(xx-1);
		for (uint yy=1; yy<=Nx; yy++) {
			uint gy = xspace->get_global_vertex_index(yy-1);
			for (uint cc=1; cc<=Nx; cc++) {
				uint gc = xspace->get_global_vertex_index(cc-1);
				
				double gx_gc_dist = xspace->get_distance(xspace->get_vertex(gx),xspace->get_vertex(gc));
				double gy_gc_dist = xspace->get_distance(xspace->get_vertex(gy),xspace->get_vertex(gc));
				
				gx_gy_gc_dist[xx-1][yy-1][cc-1] = gx_gc_dist + gy_gc_dist;
			}
		}
	}
	
	// ----------------------------------------------
	// determine which processes compute which kk,qq
	// ----------------------------------------------
	uint num_computations = 0;
	for (uint kk=0; kk<Nk; kk++) {
		for (uint qq=0; qq<=kk; qq++) {
			num_computations++;
		}
	}
	logmsg->emit(LOG_INFO,"Approximate number of [kk][qq]-pairs per process: %g",double(num_computations)/double(mpi->get_num_procs()));
	// divide count evenly amongst processes
	vector<int> computing_process; computing_process.resize(num_computations, 0);
	for (uint ii=0; ii<num_computations; ii++) {
		computing_process[ii] = floor(double(ii)/double(num_computations) * double(mpi->get_num_procs()));
		NEGF_ASSERT(computing_process[ii] < mpi->get_num_procs(), "invalid determination of computing process");
	}
	
	// -----------------------
	// compute F
	// -----------------------
	logmsg->emit(LOG_INFO,"Distributed calculation...");
	const double check = 1e100;
	const double exp_neglect = 15.0;  // exponents bigger than -exp_neglect are neglected
	vector<double> tmpvec; tmpvec.resize(Nl,check);
	vector< vector<double> > tmpvec2; tmpvec2.resize(Nk,tmpvec);
	this->F.resize(Nk, tmpvec2);
	uint kqcount = 0; uint num_kq_calc = 0;
	for (uint kk=0; kk<Nk; kk++) 
	{
		for (uint qq=0; qq<=kk; qq++) 							// F[kk][qq][.] = F[qq][kk][.] !
		{
			kqcount++;
			if (computing_process[kqcount-1]!=mpi->get_rank()) {
				continue;
			}
			logmsg->emit_noendl(LOG_INFO,"kk=%d qq=%d...   ",kk,qq);
			//logmsg->emit_noendl_all(LOG_INFO_L3,"p%d:kk=%d qq=%d...   ",mpi->get_rank(),kk,qq);
			num_kq_calc++;
			for (uint ll=0; ll<Nl; ll++) {
				F[kk][qq][ll] = 0.0;
			}
			
			double k = kspace->get_point(kk).get_coord_abs();
			double q = kspace->get_point(qq).get_coord_abs();
			vector<double> q02q2k2; q02q2k2.resize(Nx, 0.0);
			for (uint xx=1; xx<=Nx; xx++) {
				q02q2k2[xx-1] = q0[xx-1]*q0[xx-1] + q*q + k*k;
			}
			
			for (uint pp=0; pp<Nphi; pp++) 
			{
				double          dphi =      2.0*constants::pi / (Nphi-1);
				double           phi = pp * dphi;
				double       cos_phi = negf_math::cos(phi);
				
				vector<double>      phi_term;      phi_term.resize(Nx, 0.0);
				vector<double> sqrt_phi_term; sqrt_phi_term.resize(Nx, 0.0);
				for (uint xx=1; xx<=Nx; xx++) {
					     phi_term[xx-1] = q02q2k2[xx-1] - 2.0*k*q*cos_phi;
					sqrt_phi_term[xx-1] = negf_math::sqrt(phi_term[xx-1]);
				}
				
				for (uint xx=1; xx<=Nx; xx++)
				{		
					for (uint yy=1; yy<=xx; yy++) 
					{						
						// find index in array corresponding to (i,j)
						uint ll = this->find_F_index(xx,yy);
						
						double tmp = 0.0;
						for (uint cc=1; cc<=Nx; cc++) {	
							if (sqrt_phi_term[cc-1] * gx_gy_gc_dist[xx-1][yy-1][cc-1] /*positive*/ < exp_neglect) { 
								double exp_term = negf_math::exp( -sqrt_phi_term[cc-1] * gx_gy_gc_dist[xx-1][yy-1][cc-1] );
								tmp += fac[cc-1] * exp_term / phi_term[cc-1];
							}
						}
						
						this->F[kk][qq][ll] += tmp * dphi;
					}
				}
			}
		}
	}
	logmsg->emit_all(LOG_INFO_L3,"p%d finished after calculating %d [kk][qq]-pairs",mpi->get_rank(),num_kq_calc);
	logmsg->emit_noendl_all(LOG_INFO,"p%d  ",mpi->get_rank());
	mpi->synchronize_processes();
	
	// MPI communication
	logmsg->emit(LOG_INFO,"MPI communication of calculated data");
	kqcount = 0;
	for (uint kk=0; kk<Nk; kk++) {
		logmsg->emit_noendl(LOG_INFO,"kk=%d... ",kk);
		for (uint qq=0; qq<=kk; qq++) {
			int daprocess = computing_process[kqcount];
			vector<double> tmparray;
			tmparray = F[kk][qq]; // copy!
			mpi->broadcast(tmparray, daprocess);
			F[kk][qq] = tmparray; // copy!
			
			kqcount++;
		}
	}
	mpi->synchronize_processes();
	
	// compute F[kk][qq>kk][.]
	logmsg->emit(LOG_INFO,"Filling in F[kk][qq>kk][.]");
	for (uint kk=0; kk<Nk; kk++) {
		for (uint qq=0; qq<=kk; qq++) {
			for (uint ll=0; ll<Nl; ll++) {
				NEGF_ASSERT(fabs(F[kk][qq][ll] - check) > check/1000.0, "F[kk][qq][ll] did not seem to be computed.");
			}
		}
		for (uint qq=kk+1; qq<Nk; qq++) {
			for (uint ll=0; ll<Nl; ll++) {
				F[kk][qq][ll] = F[qq][kk][ll];
			}
		}
	}
	
	// compute Fmax
	logmsg->emit_noendl(LOG_INFO,"Calculating maximum F... ");
	this->Fmax = 0.0;
	for (uint kk=0; kk<Nk; kk++) {
		for (uint qq=0; qq<Nk; qq++) {
			for (uint ll=0; ll<Nl; ll++) {
				Fmax = max(Fmax, F[kk][qq][ll]);
			}
		}
	}
	logmsg->emit(LOG_INFO,"Fmax=%.3e",Fmax);
	
	// output
	if (mpi->get_rank()==0) {
		// output F[kk=5][.][.] for xx=20
		const uint kk = 5;
		const uint xx = 20;
		Matd out(Nk,xx);
		for (uint qq=1; qq<=Nk; qq++) {
			for (uint yy=1; yy<=xx; yy++) {
				uint ll = this->find_F_index(xx,yy);
				
				out(qq,xx-yy+1) = this->F[kk][qq-1][ll];
			}
		}
		 write_matrix("impurity_F_matrix", out, 
				 "Impurity F-matrix, kk=5, xx=20; first column is yy=20, second column is yy=19, etc.; rows are qq");
	}
	
	mpi->synchronize_processes();
);}


void SEIonizedImpurities::calculate()
{STACK_TRACE(
	logmsg->emit_small_header("calculating ionized impurities self-energy");
	NEGF_ASSERT(this->prepared, "call set_up(...) first!");
	const OVMat & M = ov->get_internal_overlap();
	NEGF_ASSERT(M.num_cols()==NxNn && M.num_rows()==M.num_cols(), "wrong overlap matrix.");
	
	const bool diagonals_only = true; // set to true if GF are diagonal in the bands
	
	Matc tmp(NxNn,NxNn);
	vector<Matc> MGRMs; MGRMs.resize(Nk,tmp);
	for (uint ee2 = 0; ee2 < myNE; ee2++)
	{
		uint ee = energies->get_global_index(ee2);
		
		// precompute MGRM for this energy and all k-vectors
		for (uint kk = 0; kk < Nk; kk++) {
			const Matc & GRmat = gf->get_retarded(kk,ee);
			mult(GRmat,   M,       tmp);
			mult(    M, tmp, MGRMs[kk]);
		}
		
		// prepare matrix norms --> will skip one if too small
		vector<double> MGLM_norms; MGLM_norms.resize(Nk, 0.0);
		vector<double> MGGM_norms; MGGM_norms.resize(Nk, 0.0);
		vector<double> MGRM_norms; MGRM_norms.resize(Nk, 0.0);
		for (uint kk = 0; kk < Nk; kk++) {
			MGLM_norms[kk] = negf_math::matrix_norm(gf->get_overlap_augmented_lesser(kk,ee));
			MGGM_norms[kk] = negf_math::matrix_norm(gf->get_overlap_augmented_greater(kk,ee));
			MGRM_norms[kk] = negf_math::matrix_norm(MGRMs[kk]);
		}
		
		for (uint kk = 0; kk < Nk; kk++)
		{
			SEMat & SLmat = this->get_lesser(kk,ee);
			SEMat & SGmat = this->get_greater(kk,ee);
			SEMat & SRmat = this->get_retarded(kk,ee);
			
			// re-initialize to 0
			SLmat = SIIMat_create(NxNn);
			SGmat = SIIMat_create(NxNn);
			SRmat = SIIMat_create(NxNn);
			
			if (this->scaling==0.0) {				// shortcut
				continue;
			}
			
			for (uint qq=0; qq<Nk; qq++) 
			{
				double wq = kspace->get_point(qq).get_weight() / (constants::pi*2.0); // get_weight() includes 2*pi!!!
			
				const GLMat & MGLM = gf->get_overlap_augmented_lesser (qq,ee);
				const GLMat & MGGM = gf->get_overlap_augmented_greater(qq,ee);
				const  Matc & MGRM = MGRMs[qq];
				if (   MGLM_norms[qq] < constants::GXnorm_neglect_ion 
					&& MGGM_norms[qq] < constants::GXnorm_neglect_ion 
					&& MGRM_norms[qq] < constants::GXnorm_neglect_ion) {
					continue;
				}
				
				for (uint xx=1; xx<=Nx; xx++) {
					for (uint yy=1; yy<=Nx; yy++) {
						uint ll = this->find_F_index(xx,yy);
						double Fkql = this->F[kk][qq][ll]; 
						if (Fkql < this->F_neglect*this->Fmax) continue;	// don't bother...
						
						double tmpfactor = Fkql * wq * this->scaling;
						for (uint mm=1; mm<=Nn; mm++) {
							if (diagonals_only) {
								// find matrix entry corresponding to ((x,m),(y,m))
								int ii = get_mat_idx(xx,mm,Nx);
								int jj = get_mat_idx(yy,mm,Nx);
#ifdef USE_BANDED
								if (fabs(ii-jj) > SLmat.num_offdiags+1e-8) continue;
#endif
								SLmat(ii,jj) += tmpfactor * MGLM(ii,jj);
								SGmat(ii,jj) += tmpfactor * MGGM(ii,jj);
								SRmat(ii,jj) += tmpfactor * MGRM(ii,jj);
							} else {
								for (uint nn=1; nn<=Nn; nn++) {
									// find matrix entry corresponding to ((x,m),(y,n))
									int ii = get_mat_idx(xx,mm,Nx);
									int jj = get_mat_idx(yy,nn,Nx);
#ifdef USE_BANDED
									if (fabs(ii-jj) > SLmat.num_offdiags+1e-8) continue;
#endif
									SLmat(ii,jj) += tmpfactor * MGLM(ii,jj);
									SGmat(ii,jj) += tmpfactor * MGGM(ii,jj);
									SRmat(ii,jj) += tmpfactor * MGRM(ii,jj);
								}
							}
						}
					}
				}
			}
		}
	}
	mpi->synchronize_processes();
);}

