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
#include "SEOpticalPhonon.h"
using namespace negf;

// the ANTIHERM flag defines whether special MPI routines for antihermitian routines are used or not.
//#define ANTIHERM

SEOpticalPhonon::SEOpticalPhonon(const Overlap * ov_,
					const Geometry * xspace_, 
					const Kspace * kspace_, 
					const Energies * energies_, 
					const Options * options_,
					const GreenFunctions * gf_,
					const MaterialDatabase * db):
	SelfEnergy(xspace_,kspace_,energies_,options_, constants::odSPOP),
	ov(ov_),
	gf(gf_),
	scaling(1.0),
	complicated_retarded((options->exists("LuisierSRpop") && options->get("LuisierSRpop")==1) ? true : false),
	nemo_retarded       ((options->exists("LuisierSRpop") && options->get("LuisierSRpop")==2) ? true : false),
	security_checking(true),
	diagonals_only(true),	// true-> SE is diagonal IN BANDS (NOT in space!)
	q0(constants::convert_from_SI(units::density_1d, constants::POP_q0)),
	F_neglect(1e-20),
	shift(0)
{STACK_TRACE(
	NEGF_ASSERT(ov!=NULL && xspace!=NULL && kspace!=NULL && energies!=NULL && options!=NULL && gf!=NULL,
					"null pointer encountered.");
	
	// find the bulk material, which is characterized by having the biggest volume
	vector<double> material_volumes;
	material_volumes.resize(db->get_num_materials(), 0.0);
	for (uint ii=0; ii<xspace->get_num_regions(); ii++) {
		const vector<Element *> & region_elems = xspace->get_region(ii)->get_elements();
		double region_volume = 0.0;
		for (uint jj=0; jj < region_elems.size(); jj++) {
			region_volume += region_elems[jj]->get_volume();
		}
		
		uint material_id = xspace->get_region(ii)->get_material()->get_id();
		NEGF_ASSERT(material_id < material_volumes.size(), "invalid material id.");
		material_volumes[material_id] += region_volume;
	}
	uint biggest_material_id = 0;
	for (uint ii=1; ii < material_volumes.size(); ii++) {
		if (material_volumes[ii] > material_volumes[biggest_material_id]) {
			biggest_material_id = ii;
		}
	}
	const PropertyContainer<double> * biggest_mat = db->get_material(biggest_material_id);
	logmsg->emit(LOG_INFO,"Material \"%s\" is taken for the bulk optical phonon properties.",biggest_mat->get_name().c_str());
		
	// find some material parameters
	// assume storage of lattice_constant in Angstrom, storage of optical_phonon_energy in eV
	bool is_wurtzite = (biggest_mat->get_name().find("GaN")!=string::npos);
	if (is_wurtzite) {
		logmsg->emit(LOG_INFO,"Taking wurtzite material parameters.");
		this->lattice_constant = constants::convert_from_SI(units::length, 1e-10*biggest_mat->get("lattice_constant_c"));
	} else {
		this->lattice_constant = constants::convert_from_SI(units::length, 1e-10*biggest_mat->get("lattice_constant"));
	}
	this->hw          = constants::convert_from_SI(units::energy, constants::SIec * biggest_mat->get("LO_phonon_energy"));
	const double eps0 = constants::convert_from_SI(units::dielectric, constants::SIeps0);
	const double static_dielectric = biggest_mat->get("static_dielectric_constant") * eps0;
	const double optic_dielectric  = biggest_mat->get("optic_dielectric_constant")  * eps0;
	
	logmsg->emit(LOG_INFO, "Static dielectric constant: %.5g", biggest_mat->get("static_dielectric_constant"));
	logmsg->emit(LOG_INFO, "Optic  dielectric constant: %.5g",  biggest_mat->get("optic_dielectric_constant"));
	logmsg->emit(LOG_INFO, "Lattice constant: %.4g[A]",        lattice_constant / constants::convert_from_SI(units::length, 1.0) * 1e10);
	logmsg->emit(LOG_INFO, "LO phonon energy: %.10g[meV]", 1000*biggest_mat->get("LO_phonon_energy"));
	logmsg->emit(LOG_INFO, "Screening length: %.10g[nm]", 1.0/this->q0 / constants::convert_from_SI(units::length, 1.0) * 1e9);
	
	NEGF_ASSERT(fabs(static_dielectric-optic_dielectric) > 1e-14, 
		"The static and dielectric constants are the same and therefore POP scattering will be zero. \nPlease turn off the option to save computation time.");
	
	// construct Froehlich prefactor
	const double ec = constants::convert_from_SI(units::charge, constants::SIec);
	this->prefactor = ec*ec * hw / (4.0*constants::pi*constants::pi)
					* (1.0/optic_dielectric - 1.0/static_dielectric);
	
	// multiply with factor for debugging
	if (options->exists("PhononCheatFactor")) {
		logmsg->emit(LOG_WARN,"WARNING: Artificially multiplying optical phonon strength by %g (%g-->%g)",
				options->get("PhononCheatFactor"),this->prefactor,this->prefactor * options->get("PhononCheatFactor"));
		this->prefactor *= options->get("PhononCheatFactor");
	}
	
	// calculate Bose-Einstein factor
	double kT = constants::convert_from_SI(units::energy, constants::SIkb * options->get("temperature"));
	this->Nphonon = 1.0 / (negf_math::exp(hw/kT) - 1.0);
	logmsg->emit(LOG_INFO,"Bose-Einstein distribution factor at LO phonon energy: %.3e",this->Nphonon);
		
	// set up matrices containing interpolation coefficients for coupling to E+-hw
	// G(E_i-hw) = sum_j downcoupling(i,j)*G(E_j): see line ~1000
	// G(E_i+hw) = sum_j   upcoupling(i,j)*G(E_j): see line ~1100
	this->determine_interpolations();
	
	this->determine_mpi_stuff();
		
	// determine number of differences xi-xj
	double Nl_dbl = 1 + Nx*(Nx-1)/2.0;
	NEGF_ASSERT(fabs(ceil(Nl_dbl)-Nl_dbl)<1e-14 && fabs(floor(Nl_dbl)-Nl_dbl)<1e-14, "Nl must be integer");
	this->Nl = uint(Nl_dbl);
	
	this->calculate_F();
	mpi->synchronize_processes();
);}


void SEOpticalPhonon::determine_interpolations()
{STACK_TRACE(
	// G(E_i-hw) = sum_j downcoupling(i+1,j+1)*G(E_j): see line ~1000. downcoupling is FLENS matrix --> indices start with 1
    // the contents of downcoupling and upcoupling are 0-based
	this->downcoupling = Matd(NE,NE);
    double max_norm_downcoupling = -1e10;
    double min_norm_downcoupling =  1e10;
	for(uint ee=0; ee<NE; ee++)
	{
		double  E = energies->get_energy_from_global_idx(ee);
		double E0 = (ee>0)            ? energies->get_energy_from_global_idx(ee-1) : E;
		double E1 = (ee<(this->NE-1)) ? energies->get_energy_from_global_idx(ee+1) : E;
				
		double Elow = (E+E0)/2.0 - hw;
		double Eupp = (E+E1)/2.0 - hw;
		if (Elow<energies->get_energy_from_global_idx(0)) {	
			// this is very well possible; however, Eupp<Emin was excluded by E_minus_hw_is_computed
			Elow=energies->get_energy_from_global_idx(0);
			if (fabs(Eupp-Elow)<1e-13) {
				// Elow=Eupp=lowest computed energy. in that case indices and coeffs can just remain empty.
				continue;
			}
		}
		NEGF_FASSERT(fabs(Eupp-Elow)>1e-15, "Elow=Eupp=%.5e are identical! Will yield nan coefficients. ee=%d, E=%.3e, E0=%.3e, E1=%.3e", ee, Elow, E, E0, E1);
		vector<uint> indices;
		vector<double> coeffs;
		gf->get_interpol_coeffs(indices, coeffs, Elow, Eupp); // indices-array is 0-based
		NEGF_ASSERT(indices.size()==coeffs.size(), "inconsistent indices/coeffs arrays.");
		
		for (uint ii=0; ii<indices.size(); ii++) {
			//if (indices[ii]==0) { logmsg->emit(LOG_INFO, "ee=%d couples to 0 w/ coefficient %.3e", ee, coeffs[ii]); } 
			if (ee<40 && ee>=30) {
				logmsg->emit(LOG_INFO_L1, "ee=%d (E=%g, Elow=%g, Eupp=%g) couples to %d(%g) w/ coefficient %.3e", ee, E, Elow, Eupp,
					indices[ii], energies->get_energy_from_global_idx(indices[ii]), coeffs[ii]);
			}
			NEGF_FASSERT(!isnan(coeffs[ii]) && !isinf(coeffs[ii]), "ee=%d: indices[%d]=%d, coeffs[%d]=%e", ee, ii, indices[ii], ii, coeffs[ii]);
			
			this->downcoupling(ee+1,indices[ii]+1) = coeffs[ii];// ee, indices-array are 0-based
		}

        // check
		if (indices.size()>0) {
            double check = 0.0;
            for (uint ff=0; ff<NE; ff++) {
                check += downcoupling(ee+1,ff+1);
            }
            //NEGF_FASSERT(fabs(check - 1.0) < 1e-8, "error in downcoupling: sum of coefficients was %g instead of 1. ee=%d.", check, ee);
            max_norm_downcoupling = negf_math::max(max_norm_downcoupling, check);
            min_norm_downcoupling = negf_math::min(min_norm_downcoupling, check);
		}
	}
	
	mpi->synchronize_processes();
	// G(E_i+hw) = sum_j upcoupling(i,j)*G(E_j): see line ~1100
	this->upcoupling = Matd(NE,NE);
    double max_norm_upcoupling = -1e10;
    double min_norm_upcoupling =  1e10;
	for(uint ee=0; ee<NE; ee++) {
		double dEi = energies->get_weight_from_global_idx(ee);
		NEGF_FASSERT(dEi>0.0, "ee=%d: dE=%e", ee, dEi);
		for(uint jj=0; jj<NE; jj++) {
			double dEj = energies->get_weight_from_global_idx(jj);
			this->upcoupling(ee+1,jj+1) = dEj/dEi * downcoupling(jj+1,ee+1);	// condition to fulfil particle conservation
		}

        // check
        double check = 0.0;
        for (uint ff=0; ff<NE; ff++) {
            check += upcoupling(ee+1,ff+1);
        }
        if (check>0.0) {
            // NEGF_FASSERT(fabs(check - 1.0) < 1e-8, "error in upcoupling: sum of coefficients was %g instead of 1. ee=%d.", check, ee);
            max_norm_upcoupling = negf_math::max(max_norm_upcoupling, check);
            min_norm_upcoupling = negf_math::min(min_norm_upcoupling, check);
        }
	}

	logmsg->emit(LOG_INFO,"Interpolation coefficients (should ideally sum up to 1):");
	logmsg->emit(LOG_INFO,"    downcoupling: min=%g, max=%g", min_norm_downcoupling, max_norm_downcoupling);
    logmsg->emit(LOG_INFO,"    upcoupling:   min=%g, max=%g", min_norm_upcoupling,   max_norm_upcoupling);
);}


void SEOpticalPhonon::determine_mpi_stuff()
{STACK_TRACE(
	logmsg->emit(LOG_INFO_L1,"Determining MPI communication stuff for optical phonons");
	
	// MPI: set up the array of matrices "A" which is needed to calculate the own self-energies
	// SG(k,E) = (Nq+1) AG(k,E-hw) + Nq AG(k,E+hw)
	// SL(k,E) = (Nq+1) AL(k,E+hw) + Nq AL(k,E-hw)
	
	this->myNE     = energies->get_my_number_of_points();
	this->E0_idx   = energies->get_start_global_idx(mpi->get_rank());
	this->Emax_idx = energies->get_stop_global_idx(mpi->get_rank());
				
	// -------------------------------------------------
	// needed energies in new interpolation
	// -------------------------------------------------
	vector<uint> needed_energies_below; // 0-based
	vector<uint> needed_energies_above; // 0-based
	for (uint ee=E0_idx; ee<=Emax_idx; ee++) {
		for (uint ff=0; ff<E0_idx; ff++) {
			if (downcoupling(ee+1,ff+1)!=0.0) {
				bool found = false;
				for (uint jj=0; jj<needed_energies_below.size(); jj++) {
					if(needed_energies_below[jj]==ff) {
						found = true;
						break;
					}
				}
				if (!found) needed_energies_below.push_back(ff);
			}
		}
		for (uint ff=Emax_idx+1; ff<NE; ff++) {
			if (upcoupling(ee+1,ff+1)!=0.0) {
				bool found = false;
				for (uint jj=0; jj<needed_energies_above.size(); jj++) {
					if(needed_energies_above[jj]==ff) {
						found = true;
						break;
					}
				}
				if (!found) needed_energies_above.push_back(ff);
			}
		}
	}
	
	uint below_min = NE-1;
	uint below_max = 0;
	for (uint ii=0; ii<needed_energies_below.size(); ii++) {
		below_min = negf_math::min(below_min, needed_energies_below[ii]);
		below_max = negf_math::max(below_max, needed_energies_below[ii]);
	}
	if (needed_energies_below.size()==0) below_min = 0;
	uint above_min = NE-1;
	uint above_max = 0;
	for (uint ii=0; ii<needed_energies_above.size(); ii++) {
		above_min = negf_math::min(above_min, needed_energies_above[ii]);
		above_max = negf_math::max(above_max, needed_energies_above[ii]);
	}
	if (needed_energies_above.size()==0) above_max = NE-1;
	logmsg->emit_all(LOG_INFO_L3,"p%d: below_min=%d, below_max=%d, above_min=%d, above_max=%d.",mpi->get_rank(),below_min,below_max,above_min,above_max);
	
	this->E0_minus_hw_idx   = below_min;
	this->Emax_minus_hw_idx = below_max;
	this->E0_plus_hw_idx    = above_min;
	this->Emax_plus_hw_idx  = above_max;
	
	// number of energy points below, above, and own energy points
	this->nE_below = Emax_minus_hw_idx - E0_minus_hw_idx + 1;
	this->nE_above = Emax_plus_hw_idx  - E0_plus_hw_idx  + 1;
	this->nE_self  = Emax_idx - E0_idx + 1;
	NEGF_ASSERT(nE_self == myNE, "something went wrong.");
	
	// now the possibility exists that Emax_minus_hw_idx==E0_idx or E0_plus_hw_idx==Emax_idx
	// this happens when the shift hw is small, so the own process is also involved in the interpolation of E+-hw.
	// in that case we want to get rid of the points E0_idx and Emax_idx in the auxiliary arrays because we don't need them
	if (Emax_minus_hw_idx==E0_idx) {
		nE_below--;
	}
	if (E0_plus_hw_idx==Emax_idx) {
		nE_above--;
	}
	// these cases need to be checked manually later on
	
	// allocate space; AL[...][kj] will store the lesser help matrices at point (...,kj)
	SEMat tmp = SPOPMat_create(NxNn);
	vector<SEMat> Atmp; Atmp.resize(Nk, tmp);
	this->AL.clear(); this->AL.resize(nE_below + nE_above + nE_self, Atmp);
	this->AG.clear(); this->AG.resize(nE_below + nE_above + nE_self, Atmp);
	if (complicated_retarded || nemo_retarded) {
		this->AR.clear(); this->AR.resize(nE_below + nE_above + nE_self, Atmp);
	}
	// set up control array which keeps track of the stored energy points
	this->control.assign(nE_below + nE_above + nE_self, 0);
	for (uint ii=0; ii<nE_below; ii++) { this->control[           ii] = E0_minus_hw_idx+ii; }
	for (uint ii=0; ii<nE_self ; ii++) { this->control[nE_below + ii] = E0_idx+ii; }

	// if E0_plus_hw_idx==Emax_idx, the third part of the control array stores energy indices which
	// start at E0_plus_hw_idx+1, NOT E0_plus_hw_idx
	this->shift = (E0_plus_hw_idx==Emax_idx) ? 1 : 0;
	for (uint ii=0; ii<nE_above; ii++) { this->control[nE_below + nE_self + ii] = E0_plus_hw_idx+shift + ii; }
	// if (shift==1) logmsg->emit_all(LOG_INFO,"p%d: shift=1!", mpi->get_rank()); // typically the last MPI process has shift=1
	
	mpi->synchronize_processes();
);}


/* 
how many Delta_l? Delta_l = |x_i-x_j|   i=1...Nx, j=1...Nx
since |x_i-x_j|=|x_j-x_i|, we need to save entries (i,1)...(i,i-1) AND i=j (Delta_ij=0), the latter only once
hence l has 
1 + \sum_{i=1}^{Nx} \sum_{j=1}^{i} 1 = 1 + \sum_{i=1}^{Nx} i = 1 + \sum_{i=0}^{Nx-1} i = 1 + Nx*(Nx-1)/2 entries

Storage scheme:
(i,j)     i=j   (2,1)  (3,1)  (3,2)  (4,1)  (4,2)  ...
l         0     1      2      3      4      5      ...

QUESTION: Given (i,j) where i>j, what is l?
ANSWER: 1. How many entries up to and including i-1?  
			--> 1 + \sum_{n=1}^{i-1} \sum_{m=1}^{n-1} 1 
			  = 1 + \sum_{n=1}^{i-1} n-1 
			  = 1 + \sum_{n=0}^{i-2} n 
			  = 1 + (i-1)(i-2)/2
		2. so entry (i,j), i>j, is at
			   1 + (i-1)(i-2)/2 + (j-1)       (j-1 because j is 1-based)
			= (i-1)(i-2)/2 + j.
		   entry i=j is at l=0.
		3. test, i=4, j=2: 1+3*2/2+1=5: correct!
*/
uint SEOpticalPhonon::find_F_index(uint xx, uint yy) const		// result is 0-based
{STACK_TRACE(
	if (xx==yy) {
		return 0;
	}
	
	uint smaller = negf_math::min(xx,yy);
	uint larger  = negf_math::max(xx,yy);

	double l_dbl = (larger-1.0)*(larger-2.0)/2.0 + smaller;
	NEGF_FASSERT(fabs(ceil(l_dbl)-l_dbl)<1e-14 && fabs(floor(l_dbl)-l_dbl)<1e-14, "l must be integer, instead l=%.14e",l_dbl);
	uint ll = uint(l_dbl);
	NEGF_FASSERT(ll<Nl, "wrong l (ll=%d, Nl=%d).",ll,Nl);
	return ll;
);}


/** set up F[k_i,qpar_j,Delta_l]
 *  = int_0^q_{xmax} dq_x \cos(q_x*Delta_l) * [\sqrt{(qx^2+qpar_j^2+k_i^2+q_0^2)^2
 * 			 - 4k_i^2qpar_j^2}^{-1} - q_0^2*(q_x^2+qpar_j^2+k_i^2+q_0^2)*\sqrt{(qx^2+qpar_j^2+k_i^2+q_0^2)^2 - 4k_i^2qpar_j^2}^{-3/2}]
 */
void SEOpticalPhonon::calculate_F()
{STACK_TRACE(
	logmsg->emit_header("Setting up optical-phonon F");
	
	// --------------------------------------------------------
	// initialize memory for F 
	// memory consumption will be approx. Nk*Nk*Nx*Nx/2 doubles
	// Nk=Nx=100, sizeof(double)=8 --> 8e8 bytes ~ 380MB 
	// --------------------------------------------------------
	vector<double> tmp; tmp.resize(this->Nl,0.0);
	
	// k ranges from 1 to Nk
	vector< vector<double> > tmp2; tmp2.resize(Nk,tmp);
	
	// qpar also ranges from 1 to Nk!
	//this->F.reserve(Nk, tmp2);
	//NEGF_ASSERT(this->F.capacity() >= Nk, "not enough memory for F");
	this->F.resize(Nk,tmp2);
	
	// ============================
	// calculate F!
	// ============================
		
	// set up some parameters for the integral
	const uint Nqx = 200; 
	const double qxmax = constants::pi / this->lattice_constant;
	
	// -------------------------------
	// precompute some quantities
	// -------------------------------
	vector< vector<double> > cos_qx_delta;
	cos_qx_delta.resize(Nqx);
	for (uint qqxx=0; qqxx < Nqx; qqxx++) 
	{
		cos_qx_delta[qqxx].resize(Nl, 0.0);
		
		double dqx = qxmax/(Nqx-1);		// interval length (equally spaced)
		double  qx = qqxx * dqx;			// midpoint of interval
		
		vector<bool> ll_used; ll_used.resize(Nl, false);
		for (int xx=1; xx<=int(Nx); xx++)
		{
			cos_qx_delta[qqxx][0] = 1.0;  //xx==yy

			for (int yy=1; yy<xx; yy++) 
			{
				// find |xi-xj|
				uint gx = xspace->get_global_vertex_index(xx-1);
				uint gy = xspace->get_global_vertex_index(yy-1);
				double delta = xspace->get_distance(xspace->get_vertex(gx),xspace->get_vertex(gy));
				
				// find index in array corresponding to (i,j)
				uint ll = this->find_F_index(xx,yy);
				
				// check 
				uint ll2 = this->find_F_index(yy,xx);
				NEGF_ASSERT(ll==ll2, "something went wrong (find_F_index)");
				NEGF_ASSERT(!ll_used[ll], "ll_used[ll]");
				ll_used[ll] = true;
				
				cos_qx_delta[qqxx][ll] = negf_math::cos(qx*delta);
			}
		}
	}
	
	int next = 0; int count = 0;
	logmsg->init_progress_bar(LOG_INFO, "   progress of master thread, Nk=", Nk);
	// compute "lower half" and diagonal of F[k][q]
	for (uint kk=0; kk<Nk; kk++)
	{
		double k = kspace->get_point(kk).get_coord_abs();
		for (uint qq=0; qq/*<Nk*/<=kk; qq++)
		{
			double q = kspace->get_point(qq).get_coord_abs();
			
			double sum1 = k*k + q*q + q0*q0;
			double fac1 = 4.0 * k*k * q*q;
			
			for (uint qqxx=0; qqxx<Nqx; qqxx++)
			{
				double dqx = qxmax/(Nqx-1);		// interval length (equally spaced)
				double qx = qqxx * dqx;			// midpoint of interval
				if (qqxx==0 || qqxx==Nqx-1) {	// adjust interval length for first and last point
					dqx = dqx / 2.0;
				}
				
				double sum2 = sum1 + qx*qx;
				
				double denominator = negf_math::sqrt(sum2*sum2-fac1);	// sqrt((qx^2+qpar^2+k^2+q0^2)^2-4k^2qpar^2)
				
				double term1 = 1.0 / denominator;
				double term2 = q0*q0*sum2 / (denominator*denominator*denominator);
				
				double dqx_term12 = dqx * (term1-term2);
					
				// l=0 (xi=xj)
				F[kk][qq][0] += 1.0 * dqx_term12;
				
				// l>0
				for (uint xx=1; xx<=Nx; xx++) {
					for (uint yy=1; yy<xx; yy++) {
						uint ll = this->find_F_index(xx,yy);
						
						F[kk][qq][ll] += cos_qx_delta[qqxx][ll] * dqx_term12;
						// this is done Nk*Nk*Nqx*Nx*(Nx-1)/2 times. Nk=100, Nqx=200, Nx=100 --> 10^10 times!!!!!
					}
				}
			}
		}
		if(next == count++) next = logmsg->set_progress_bar(count, Nk);
	}
	logmsg->end_progress_bar();
	
	// fill "upper half" (in k,q)
	for (uint kk=0; kk<Nk; kk++) {
		for (uint qq=kk+1; qq<Nk; qq++) {
			for (uint ll=0; ll<Nl; ll++) {
				F[kk][qq][ll] = F[qq][kk][ll];
			}
		}
	}
	
	// check for NaN's of Inf's
	for (uint kk=0; kk<Nk; kk++) {
		for (uint qq=0; qq<Nk; qq++) {
			for (uint ll=0; ll<Nl; ll++) {
				NEGF_FASSERT(!isnan(F[kk][qq][ll]) && !isinf(F[kk][qq][ll]), "NaN or Inf F-function encountered at kk=%d, qq=%d, ll=%d",kk,qq,ll);
			}
		}
	}
		
	// output F[kk=5][qq][ll] for all qq and xx=20,yy=1...20
	if (mpi->get_rank()==0) {
		const uint kk = 5;
		const uint xx = 20;
		Matd out(Nk,xx);
		for (uint qq=1; qq<=Nk; qq++) {
			for (uint yy=1; yy<=xx; yy++) {
				uint ll = this->find_F_index(xx,yy);
				
				out(qq,xx-yy+1) = this->F[kk][qq-1][ll];
			}
		}
		 write_matrix("poloptphon_F_matrix", out, 
				 "Polar optical phonon F-matrix, kk=5, xx=20; first column is yy=20, second column is yy=19, etc.; rows are qq");
	}

	/** HACK for DEBUG: F=0 except for k=q, there F[k][q][l] = F[k][q][0]! */
	/*for (uint kk=0; kk<Nk; kk++) {
		for (uint qq=0; qq<Nk; qq++) {
			for (uint ll=0; ll<Nl; ll++) {
				if (qq==kk) {
					//F[kk][qq][ll] = F[kk][qq][0];
					//F[kk][qq][ll] *= 100.0;
				} else {
					//F[kk][qq][ll] = 0.0;
				}
			    //F[kk][qq][ll] = F[kk][qq][0];
				//if (ll!=0) F[kk][qq][ll] = 0;
			}
		}
	}*/
		
	// now multiply F[k,qpar,dx] already with qpar*dqpar since it will only used in conjunction with that!
	logmsg->emit(LOG_INFO,"Multiplying F(q,k,dx) with q*dq...");
	for (uint qq=0; qq<Nk; qq++) {
		//double qpar_dqpar = kspace->get_point(qq).get_weight() / (2.0*constants::pi); 
		// p.get_weight() = q*dq*2pi - don't want angle integration 2pi!
		
		// <ss 8.6.10> maybe SE was 2pi too high!
		double kspace_factor = 1.0 / (negf_math::pow(2.0*constants::pi, kspace->get_dimension()));
		double qpar_dqpar = kspace->get_point(qq).get_weight() * kspace_factor; // get_weight() includes 2*pi
		
		for (uint kk=0; kk<Nk; kk++) {
			for (uint ll=0; ll<Nl; ll++) {
				F[kk][qq][ll] = F[kk][qq][ll] * qpar_dqpar;
			}
		}
	}
	
	
	logmsg->emit_noendl(LOG_INFO,"Calculating maximum F*q*dq... ");
	this->Fmax = 0.0;
	for (uint kk=0; kk<Nk; kk++) {
		for (uint qq=0; qq<Nk; qq++) {
			for (uint ll=0; ll<Nl; ll++) {
				Fmax = max(Fmax, F[kk][qq][ll]);
			}
		}
	}
	logmsg->emit(LOG_INFO,"Fmax=%.3e",Fmax);
	
	
	// small test
	/*double Fsum1 = 0.0;
	for (uint kk=0; kk<Nk; kk++) {
		double wk = kspace->get_point(kk).get_weight();
		for (uint qq=0; qq<Nk; qq++) {
			Fsum1 += wk * F[kk][qq][1];
		}
	}
	double Fsum2 = 0.0;
	for (uint kk=0; kk<Nk; kk++) {
		for (uint qq=0; qq<Nk; qq++) {
			double wq = kspace->get_point(qq).get_weight();
			Fsum2 += wq * F[qq][kk][1];
		}
	}
	logmsg->emit(LOG_INFO,"Fsum1=%e, Fsum2=%e",Fsum1,Fsum2);*/
	
	// another small test
	/*for (uint kk=0; kk<Nk; kk++) {
		double wk = kspace->get_point(kk).get_weight();
		for (uint qq=0; qq<Nk; qq++) {
			double wq = kspace->get_point(kk).get_weight();
			for (uint ll=0; ll<Nl; ll++) {
				NEGF_FASSERT(fabs(F[kk][qq][ll] - F[qq][kk][ll]) < 1e-10, 
					"something went wrong: kk=%d, qq=%d, ll=%d, wk*Fkql=%e, wq*Fqkl=%e",kk,qq,ll,wk*F[kk][qq][ll],wq*F[qq][kk][ll]);
			}
		}
	}*/
);}


/** Test the following:
	Given that (MGLM)_ij(k,E) = imag_unit*(i+j+E)     for all k 
		   and (MGGM)_ij(k,E) = imag_unit*(i+j+E+8.8) for all k 
		   and F[k][q][l] = 1,
	test that the interpolation SG/SL also equals i+j+E+-hw (times prefactor an Nk). */
void SEOpticalPhonon::test_operation() throw (Exception *)
{STACK_TRACE(
	logmsg->emit_header("Testing correct operation of POP self-energy");
	mpi->synchronize_processes();

	// test spatial symmetry of F
    logmsg->emit(LOG_INFO,"Testing symmetry of find_F_index");
    GLMat & arb_mat = gf->get_lesser(0,0);
	for (uint xx=1; xx<=Nx; xx++) {
	    for (uint yy=1; yy<=xx; yy++) {
#ifdef USE_BANDED
             if (fabs(xx-yy) > arb_mat.num_offdiags+1e-8) continue;
#endif
             uint l1 = this->find_F_index(xx,yy);
             uint l2 = this->find_F_index(yy,xx);
             NEGF_FASSERT(fabs((double)l1-(double)l2) < 1e-8, "xx=%d, yy=%d: l1=%d, l2=%d", xx, yy, l1, l2);
	    }
	}

	// -------------------------------------------------------------
	// Assign (MG.M)_ij(E) = imag_unit*(i+j+E) for all k (needs to be anti-Hermitian)
	// MG.M IS RENDERED UNUSABLE FOR LATER BY THIS!
	// -------------------------------------------------------------
	logmsg->emit(LOG_INFO,"Assigning MG.M");
	for (uint ee2=0; ee2<this->myNE; ee2++) {
		uint   ee = energies->get_global_index(ee2);
		double E  = energies->get_energy_from_global_idx(ee);
		
		for (uint kk=0; kk < Nk; kk++) {
			GLMat & MGLM = gf->get_overlap_augmented_lesser(kk,ee);
			GLMat & MGGM = gf->get_overlap_augmented_greater(kk,ee);

			for (int ii=1; ii<=int(NxNn); ii++) {
				for (int jj=1; jj<=int(NxNn); jj++) {
#ifdef USE_BANDED
					if (fabs(ii-jj) > MGLM.num_offdiags+1e-8) continue;
#endif
					double real_part = 0;
					if (ii<jj) real_part = E;
					if (ii>jj) real_part = -E;
					MGLM(ii,jj) = real_part + constants::imag_unit * (ii + jj + E);
					MGGM(ii,jj) = real_part + constants::imag_unit * (ii + jj + E + 8.8);
				}
			}
		}
	}

	// -------------------------------------------------------------
	// Assign F[k][q][l] = 1 and prefactor = 1
	// F AND prefactor ARE RENDERED UNUSABLE FOR LATER BY THIS!
	// -------------------------------------------------------------
	logmsg->emit(LOG_INFO,"Assigning F=1 and prefactor=1");
	mpi->synchronize_processes();
	for (uint kk=0; kk<Nk; kk++) {
		for (uint qq=0; qq<Nk; qq++) {
			for (uint ll=0; ll<Nl; ll++) {
				F[kk][qq][ll] = 1.0;
			}
		}
	}
	//this->prefactor = 1.0;
	//this->Nphonon = 0.0;
	//this->diagonals_only = false;

	// --------------------------------------------------------------
	// Calculate lesser and greater self-energy
	// includes determination of interpolation and MPI stuff
	// --------------------------------------------------------------
	logmsg->emit(LOG_INFO,"Calculating lesser and greater (and retarded)");
	logmsg->set_level(LOG_INFO_L3);
	this->calculate();
    mpi->synchronize_processes();

	// --------------------------------------------------------------
	// Test 
	// AL_ij[E][k] = \int dqpar*qpar*F[k][qpar][dij] * (M*GL[E][qpar]*M)_ij,
	// AG: same
	// but in this implementation  dqpar*qpar is already included in F, which has been set to 1
	// hence AL_ij[E][k] = Nk*(M*GL[E][qpar]*M)_ij.
	//
	// SG(E) = (N+1)*p*AG(E-hw) + N*p*AG(E+hw)
	// SL(E) = (N+1)*p*AL(E+hw) + N*p*AL(E-hw)
	// --------------------------------------------------------------
	logmsg->emit(LOG_INFO,"Testing");
	for (uint ee2=0; ee2<this->myNE; ee2++) {
		uint   ee = energies->get_global_index(ee2);
		double E  = energies->get_energy_from_global_idx(ee);
		double Em = (ee>0)    ? 0.5*(E+energies->get_energy_from_global_idx(ee-1)) : E;
        double Ep = (ee<NE-1) ? 0.5*(E+energies->get_energy_from_global_idx(ee+1)) : E;

		double Emin = energies->get_energy_from_global_idx(0);
		double Emax = energies->get_energy_from_global_idx(this->NE-1);
		bool E_minus_hw_is_computed = (Ep-hw >= Emin) ? true : false;
		bool E_plus_hw_is_computed  = (Em+hw <= Emax) ? true : false;

		for (uint kk=0; kk<Nk; kk++) {
			SEMat & SL = this->get_lesser(kk,ee);
			SEMat & SG = this->get_greater(kk,ee);
			
			for (int ii=1; ii<=int(NxNn); ii++) {
				for (int jj=1; jj<=int(NxNn); jj++) {
#ifdef USE_BANDED
					if (fabs(ii-jj) > SL.num_offdiags+1e-8) continue;
#endif
					// should-be values are independent of kk
					double real_part_minus_hw = 0; if (ii<jj) real_part_minus_hw = E-hw; if (ii>jj) real_part_minus_hw = -(E-hw);
					double real_part_plus_hw = 0;  if (ii<jj) real_part_plus_hw  = E+hw; if (ii>jj) real_part_plus_hw  = -(E+hw);
					cplx should_be_SG = 0.0;
					if (E_minus_hw_is_computed) should_be_SG += (Nphonon+1)*1*Nk*prefactor * (real_part_minus_hw + (ii+jj+E-hw+8.8)*constants::imag_unit);
					if (E_plus_hw_is_computed)  should_be_SG +=     Nphonon*1*Nk*prefactor * (real_part_plus_hw  + (ii+jj+E+hw+8.8)*constants::imag_unit);	
					cplx should_be_SL = 0.0;
					if (E_plus_hw_is_computed)  should_be_SL += (Nphonon+1)*1*Nk*prefactor * (real_part_plus_hw  + (ii+jj+E+hw)*constants::imag_unit);
					if (E_minus_hw_is_computed) should_be_SL +=     Nphonon*1*Nk*prefactor * (real_part_minus_hw + (ii+jj+E-hw)*constants::imag_unit);

					if (!E_minus_hw_is_computed || !E_plus_hw_is_computed) continue; // test fails for first few points...
					
					NEGF_FASSERT(abs(SL(ii,jj)-should_be_SL)/(abs(should_be_SL)+1e-15) < 1e-10, "ee=%d, kk=%d, ii=%d, jj=%d: SL=(%.5e,%.5e), should be (%.5e,%.5e), AL_ij[E][k]=(%.e, %.e).", ee, kk, ii, jj, SL(ii,jj).real(), SL(ii,jj).imag(), should_be_SL.real(), should_be_SL.imag(), AL[ee-E0_idx+nE_below][kk](ii,jj).real(), AL[ee-E0_idx+nE_below][kk](ii,jj).imag());
					NEGF_FASSERT(abs(SG(ii,jj)-should_be_SG)/(abs(should_be_SG)+1e-15) < 1e-10, "ee=%d, kk=%d, ii=%d, jj=%d: SG=(%.5e,%.5e), should be (%.5e,%.5e), AG_ij[E][k]=(%.e, %.e).", ee, kk, ii, jj, SG(ii,jj).real(), SG(ii,jj).imag(), should_be_SG.real(), should_be_SG.imag(), AG[ee-E0_idx+nE_below][kk](ii,jj).real(), AG[ee-E0_idx+nE_below][kk](ii,jj).imag());
				}
			}
		}
		logmsg->emit_noendl_all(LOG_INFO,"ee=%d is ok!    ",ee);
	}

	mpi->synchronize_processes();
	logmsg->emit(LOG_INFO, "\nPOP test succeeded.");
	NEGF_EXCEPTION("Must abort after test (set F to an artificial value).");
);}

	
void SEOpticalPhonon::calculate()
{STACK_TRACE(
	this->determine_interpolations();
	this->determine_mpi_stuff();
	this->calculate_lesser_greater();
	if (!complicated_retarded && !nemo_retarded) {	// otherwise the retarded self-energy was computed in calculate_lesser_greater
		this->calculate_retarded();
	}
	
	// screen output of renormalization energy. note: has units M*Sig*M
	for (uint ee2=0; ee2<this->myNE; ee2++) {
		uint   ee = energies->get_global_index(ee2);
		double E  = energies->get_energy_from_global_idx(ee);
		if (E>1.1 || E<0.4) { continue; }
		
		uint xx = 22;
		char buf[2000];
		sprintf(buf," ");
		bool enter = false;
		for (uint nn=1; nn<=Nn; nn++) {
			double re = (this->get_retarded(0,ee))(get_mat_idx(xx,nn,Nx),get_mat_idx(xx,nn,Nx)).real();
			double im = (this->get_retarded(0,ee))(get_mat_idx(xx,nn,Nx),get_mat_idx(xx,nn,Nx)).imag();
			if (fabs(re)>1e-10 || fabs(im)>1e-10) {
				sprintf(buf,"%sSR(E=%.3g,k=0,xx=%d,xx,n=%d,n): Re=%.2e, Im=%.2e     ",buf, E, xx, nn, re, im);
				enter = true;
			}
		}
		if (enter) logmsg->emit_all(LOG_INFO_L2,buf);
	}
);}


/** this method involves a lot of complicated MPI communication 
 *  please refer to the NEGF report of S. Steiger for more details */
void SEOpticalPhonon::calculate_lesser_greater()
{STACK_TRACE(
	if (complicated_retarded) {
		logmsg->emit_small_header("calculating SL,SG and Luisier-SR of optical phonon self-energy");
	} else if (nemo_retarded) {
		logmsg->emit_small_header("calculating SL,SG and Nemo1D-SR of optical phonon self-energy");
	} else {
		logmsg->emit_small_header("calculating lesser and greater optical phonon self-energy");
	}
	
	if (this->scaling==0.0) {	// shortcut...
		logmsg->emit(LOG_INFO,"Skipping detailed computation because scaling=0");
		for (uint ee2=0; ee2<this->myNE; ee2++) {
			uint ee=energies->get_global_index(ee2);
			for (uint kk=0; kk<Nk; kk++) {
				SEMat & SL = this->get_lesser(kk,ee);
				SEMat & SG = this->get_greater(kk,ee);
				SEMat & SR = this->get_retarded(kk,ee);
				SL = SPOPMat_create(NxNn);
				SG = SPOPMat_create(NxNn);
				SR = SPOPMat_create(NxNn);
			}
		}
		return;
	}
		
	for (uint aa=0; aa < this->AL.size(); aa++) {
		for (uint kk=0; kk<Nk; kk++) {
			(this->AL[aa][kk])(1,1) = 888.888;
			(this->AG[aa][kk])(1,1) = 888.888;
			if (complicated_retarded || nemo_retarded) {
				(this->AR[aa][kk])(1,1) = 888.888;
			}
		}
	}
	
	// -------------------------------------------------------------------
	// calculate AL and AG of own energies
	// this effectively means performing the integral over q-space
	// result: AL and AG arrays with indices nE_below...nE_below+myNE-1
	// -------------------------------------------------------------------
	logmsg->emit(LOG_INFO,"Own energies...");
	const OVMat & M = ov->get_internal_overlap();
	NEGF_ASSERT(M.num_cols()==NxNn && M.num_rows()==M.num_cols(), "wrong overlap matrix.");
	Matc  GRM(NxNn,NxNn); 						// helper
	vector<Matc> MGRMs; MGRMs.resize(Nk, GRM);  // helper
	for (uint ee2=0; ee2<this->myNE; ee2++)
	{
		uint ee=energies->get_global_index(ee2);
		
		if (complicated_retarded || nemo_retarded) {
			// precompute M*GR*M for later use
			for (uint qq=0; qq < Nk; qq++) {
				
				const Matc & GR = gf->get_retarded(qq,ee);
				mult(GR, M, GRM); 
				mult( M, GRM, MGRMs[qq]); // MGRM[qq] =  M * GR[qq][ee] * M;
				
#ifdef USE_BANDED
				if (!gf->mangle_overlap_calculation()) {
					// in this case M*GL*M was calculated from the banded GL
					// hence we need to crop GR to banded form, like GL and GG, before multiplying with Ms!
					GLMat GR2 = GLMat_create(NxNn); GR2 = GR;
					mult(GR2, M, GRM); 
					mult( M, GRM, MGRMs[qq]);
				}
#endif
			}
		}
		
		for (uint kk=0; kk<Nk; kk++)
		{
			NEGF_ASSERT(AL.size()>ee-E0_idx+nE_below, "invalid AL size");
			NEGF_FASSERT(control[ee-this->E0_idx + nE_below]==ee, "control failed (1): control[ee-this->E0_idx + nE_below]=%d, ee=%d.",control[ee-this->E0_idx + nE_below],ee);
			SEMat & ALmat = this->AL[ee-this->E0_idx + nE_below][kk];
			SEMat & AGmat = this->AG[ee-this->E0_idx + nE_below][kk];
			
			// re-init to zero
			ALmat = SPOPMat_create(NxNn);
			AGmat = SPOPMat_create(NxNn);
			
			// ---------------------------------------------------------------------
			// lower trigonal part of AL, AG
			// EDIT: currently, total matrix
			// AL_ij[E][k] = int_dqpar qpar*F[k][qpar][dij] * (M*GL[E][qpar]*M)_ij
			// ---------------------------------------------------------------------
							
			// integrate over q_parallel
			for (uint qq=0; qq < Nk; qq++) 
			{			
				const GLMat & MGLM = gf->get_overlap_augmented_lesser(qq,ee);
				const GLMat & MGGM = gf->get_overlap_augmented_greater(qq,ee);

				for (int xx=1; xx<=int(Nx); xx++) 
				{
					for (int yy=1; yy<=/*xx*/int(Nx); yy++) 
					{
						// find index in F-array corresponding to (i,j)
						uint ll = this->find_F_index(xx,yy);

						const double Fkql = this->F[kk][qq][ll]; 			// already includes qpar*dqpar
						if (Fkql < this->F_neglect*this->Fmax) continue;	// don't bother...
						
						for (uint mm=1; mm<=Nn; mm++) 
						{
							if (diagonals_only) {
								// find matrix entry corresponding to ((x,m),(y,m))
								int ii = get_mat_idx(xx,mm,Nx);
								int jj = get_mat_idx(yy,mm,Nx);
#ifdef USE_BANDED
								if (fabs(ii-jj) > ALmat.num_offdiags+1e-8) continue;
#endif
								ALmat(ii,jj) += Fkql * MGLM(ii,jj);
								AGmat(ii,jj) += Fkql * MGGM(ii,jj);
							} else {
								for (uint nn=1; nn<=Nn/*mm*/; nn++) 
								{
									// find matrix entry corresponding to ((x,m),(y,n))
									int ii = get_mat_idx(xx,mm,Nx);
									int jj = get_mat_idx(yy,nn,Nx);
#ifdef USE_BANDED
									if (fabs(ii-jj) > ALmat.num_offdiags+1e-8) continue;
#endif
									ALmat(ii,jj) += Fkql * MGLM(ii,jj);
									AGmat(ii,jj) += Fkql * MGGM(ii,jj);
								}
							}
						}
					}
				}
			}
			// -------------------------------------------
			// fill upper part of AL,AG accordingly
			// -------------------------------------------
			/*for (int ii=1; ii<=int(NxNn); ii++) {
				for (int jj=ii+1; jj<=int(NxNn); jj++) {
#ifdef USE_BANDED
					if (fabs(ii-jj) > ALmat.num_offdiags+1e-8) continue;
#endif
					ALmat(ii,jj) = -ALmat(jj,ii).real() + constants::imag_unit * ALmat(jj,ii).imag();
					AGmat(ii,jj) = -AGmat(jj,ii).real() + constants::imag_unit * AGmat(jj,ii).imag();
				}
			}*/
						
			if (complicated_retarded || nemo_retarded) 
			{
				NEGF_ASSERT(control[ee-this->E0_idx + nE_below]==ee, "control failed (1b).");
				SEMat & ARmat = this->AR[ee-this->E0_idx + nE_below][kk];
				ARmat = SPOPMat_create(NxNn); // re-init to 0
				// integrate over q_parallel
				for (uint qq=0; qq < Nk; qq++) 
				{
					const  Matc & MGRM = MGRMs[qq];
					
					for (int xx=1; xx<=int(Nx); xx++) 
					{
						for (int yy=1; yy<=int(Nx); yy++) 
						{
							// find index in F-array corresponding to (i,j)
							uint ll = this->find_F_index(xx,yy);

							const double Fkql = this->F[kk][qq][ll]; // already includes qpar*dqpar
							if (Fkql < this->F_neglect*this->Fmax) continue;	     // don't bother...
							
							for (uint mm=1; mm<=Nn; mm++) 
							{
								for (uint nn=1; nn<=Nn; nn++) 
								{
									// find matrix entry corresponding to ((x,m),(y,n))
									int ii = get_mat_idx(xx,mm,Nx);
									int jj = get_mat_idx(yy,nn,Nx);
	#ifdef USE_BANDED
									if (fabs(ii-jj) > ARmat.num_offdiags+1e-8) continue;
	#endif
									ARmat(ii,jj) += Fkql * MGRM(ii,jj);
								}
							}
						}
					}
				}
			} // complicated_retarded || nemo_retarded
		}
	}
	
	// security check: check Sigma = -Sigma+
	if (security_checking) {
		SEMat ALtmp = SPOPMat_create(NxNn);
		SEMat AGtmp = SPOPMat_create(NxNn);
		for (uint ee2=0; ee2<myNE; ee2++)
		{
			uint ee=energies->get_global_index(ee2);
			for (uint kk=0; kk<Nk; kk++) {
				const SEMat & ALmat = this->AL[ee-this->E0_idx + nE_below][kk];
				const SEMat & AGmat = this->AG[ee-this->E0_idx + nE_below][kk];
				
				conjtrans(ALmat, ALtmp); // ALtmp = conjugateTranspose(ALmat); 
				ALtmp += ALmat;
				conjtrans(AGmat, AGtmp); // AGtmp = conjugateTranspose(AGmat); 
				AGtmp += AGmat;
				double ALnorm = negf_math::matrix_norm(ALtmp);
				double AGnorm = negf_math::matrix_norm(AGtmp);
				NEGF_FASSERT(ALnorm < constants::antiherm_check, "AL was not anti-hermitian: |AL-(-AL+)|=%e, ee=%d=%.4g, kk=%d", ALnorm, ee, energies->get_energy_from_global_idx(ee), kk);
				NEGF_FASSERT(AGnorm < constants::antiherm_check, "AG was not anti-hermitian: |AG-(-AG+)|=%e, ee=%d=%.4g, kk=%d", AGnorm, ee, energies->get_energy_from_global_idx(ee), kk);
			}
		}
	}
	
	logmsg->emit_noendl_all(LOG_INFO_L2,"p%d  ",mpi->get_rank());
	mpi->synchronize_processes();
	
	// -------------------------------------------------------------------------
	// MPI communication of AL, AG such that in the end the
	// whole arrays of AL and AG are filled for every process
	// -------------------------------------------------------------------------
	double t1 = MPI_Wtime();
	this->communicate_As();
	logmsg->emit_noendl_all(LOG_INFO_L2,"p%d  ",mpi->get_rank());
	
	// security check
	if (security_checking) {
		bool ok = true;
		for (uint aa=0; aa < this->AL.size(); aa++) {
			for (uint kk=0; kk<Nk; kk++) {
				uint ee = aa+E0_minus_hw_idx;
				if (aa>=nE_below)         {	ee = aa-nE_self          + this->E0_idx;	}
				if (aa>=nE_self+nE_below) {	ee = aa-nE_self-nE_below + this->E0_plus_hw_idx; }
				if ((this->AL[aa][kk])(1,1) == 888.888) {
					logmsg->emit_noendl_all(LOG_WARN,"p%d,ee=%d,kk=%d: AL(1,1)=888.888    ", mpi->get_rank(), ee, kk);
					ok = false;
				}
				if ((this->AG[aa][kk])(1,1) == 888.888) {
					logmsg->emit_noendl_all(LOG_WARN,"p%d,ee=%d,kk=%d: AG(1,1)=888.888    ", mpi->get_rank(), ee, kk);
					ok = false;
				}
				if ((complicated_retarded || nemo_retarded) && ((this->AR[aa][kk])(1,1) == 888.888)) {
					logmsg->emit_noendl_all(LOG_WARN,"p%d,ee=%d,kk=%d: AR(1,1)=888.888    ", mpi->get_rank(), ee, kk);
					ok = false;
				}
			}
		}
		//if (!ok) NEGF_EXCEPTION("Aborting.");
	}
	mpi->synchronize_processes();
	double t2 = MPI_Wtime();
	logmsg->emit(LOG_INFO_L2, "Time needed for MPI communication: %.5g[s].   ",t2-t1);
	
	
	logmsg->emit(LOG_INFO,"Assigning self-energies   ");
	// ---------------------------------------------------------
	// compute lesser and greater self-energies for given k!
	// ---------------------------------------------------------
	for (uint ee2=0; ee2<this->myNE; ee2++)
	{
		uint   ee         = energies->get_global_index(ee2);
		double E          = energies->get_energy_from_global_idx(ee);

		for (uint kk=0; kk<Nk; kk++) 
		{
			// ------------------------------------
			// DO IT!
			// SG = (N+1)*p*AG(E-hw) + N*p*AG(E+hw)
			// SL = (N+1)*p*AL(E+hw) + N*p*AL(E-hw)
			// ------------------------------------
			SEMat & SL = this->get_lesser(kk,ee);
			SEMat & SG = this->get_greater(kk,ee);
			SEMat & SR = this->get_retarded(kk,ee);
			SL = SPOPMat_create(NxNn);
			SG = SPOPMat_create(NxNn);
			if (complicated_retarded || nemo_retarded) {
				SR = SPOPMat_create(NxNn);
			}
			
            // --------------------------------------------------------------------------
            // find indices of energy points between which E-hw lies
            // linear interpolation between two computed Green functions on energy axis
            // --------------------------------------------------------------------------
            vector<uint> indices;
            vector<double> coeffs;
            for (uint ii=0; ii<NE; ii++) {
                if (this->downcoupling(ee+1,ii+1)!=0.0) {
                    indices.push_back(ii);
                    coeffs.push_back(this->downcoupling(ee+1,ii+1));
                }
            }
            if (kk==0 && E>1.45 && E<1.5) {
                logmsg->emit_all(LOG_INFO_L3,"E=%.5g (ee=%d) (-hw):",E,ee);
                for (uint ii=0; ii<indices.size(); ii++) {
                    logmsg->emit_all(LOG_INFO_L3,"    relies on E=%.5g (ee=%d) with fraction %.5g",
                            energies->get_energy_from_global_idx(indices[ii]),indices[ii], coeffs[ii]);
                }
            }
            // indices[ii] stores global energy indices
            // such an index can be either within E0_idx...Emax_idx --> stored at AX[nE_below + idx-E0_idx]
            //                          or below  E0_idx            --> stored at AX[idx-E0_minus_hw_idx]
            for (uint ii=0; ii<indices.size(); ii++)
            {
                int local_idx = indices[ii] - this->E0_minus_hw_idx;
                if (indices[ii]>=this->E0_idx) {
                    // energy index index[ii] belongs to the own MPI process
                    NEGF_ASSERT(indices[ii]<=this->Emax_idx, "E-hw should not couple to indices above Emax!");
                    local_idx = nE_below + indices[ii]-this->E0_idx;
                }
                NEGF_FASSERT(local_idx < int(AL.size()), "invalid local_idx (-hw): local_idx=%d, Al.size()=%d", local_idx, AL.size());
                NEGF_FASSERT(control[local_idx]==indices[ii], "p%d (ee=%d...%d) control failed (2): control[%d]=%d, indices[%d]=%d, nE_below=%d.",
                        mpi->get_rank(),E0_idx,Emax_idx,local_idx,control[local_idx],ii,indices[ii],nE_below);

                add(this->AL[local_idx][kk], (coeffs[ii] * (this->Nphonon)   * this->prefactor * this->scaling), SL);
                // SL +=   ...(N) * AL[local_idx][kk];
                add(this->AG[local_idx][kk], (coeffs[ii] * (this->Nphonon+1) * this->prefactor * this->scaling), SG);
                // SG += ...(N+1) * AG[local_idx][kk];

                if (complicated_retarded) {
                    add(this->AR[local_idx][kk], (coeffs[ii] * (this->Nphonon+1) * this->prefactor * this->scaling), SR);
                    // SR += ...(N+1) * AR[local_idx][kk];
                    add(this->AL[local_idx][kk], (coeffs[ii] * 0.5               * this->prefactor * this->scaling), SR);
                    // SR += ...(0.5) * AL[local_idx][kk];
                }
                if (nemo_retarded) {
                    add(this->AR[local_idx][kk], (coeffs[ii] * (this->Nphonon+1) * this->prefactor * this->scaling), SR);
                    // SR += ...(N+1) * AR[local_idx][kk];
                }
            }
			
            // --------------------------------------------------------------------------
            // find indices of energy points between which E+hw lies
            // linear interpolation between two computed Green functions on energy axis
            // --------------------------------------------------------------------------
            indices.clear();
            coeffs.clear();
            for (uint ii=0; ii<NE; ii++) {
                if (this->upcoupling(ee+1,ii+1)!=0.0) {
                    indices.push_back(ii);
                    coeffs.push_back(this->upcoupling(ee+1,ii+1));
                }
            }
            if (kk==0 && E>1.4 && E<1.45) {
                logmsg->emit_all(LOG_INFO_L3,"E=%.5g (ee=%d) (+hw):",E,ee);
                for (uint ii=0; ii<indices.size(); ii++) {
                    logmsg->emit_all(LOG_INFO_L3,"    relies on E=%.5g (ee=%d) with fraction %.5g",
                            energies->get_energy_from_global_idx(indices[ii]),indices[ii], coeffs[ii]);
                }
            }
            // now "indices" and "coeffs" store the information of 1/(Eupp-Elow) * int_Elow^Eupp f(E)dE ~ sum_i c_i f(Ei)

            // indices[ii] stores global energy indices
            // such an index can be either within E0_idx...Emax_idx --> stored at AX[nE_below + idx-E0_idx]
            //                          or above  Emax_idx          --> stored at AX[nE_below + nE_self + idx-E0_plus_hw_idx (+shift)]
            for (uint ii=0; ii<indices.size(); ii++)
            {
                NEGF_FASSERT(indices[ii]>=this->E0_idx, "E+hw should not couple to indices below E0! indices[%d]=%d, coeffs[%d]=%.3e, E0_idx=%d",
                        ii, indices[ii], coeffs[ii], E0_idx);
                NEGF_FASSERT(indices[ii]<=this->Emax_plus_hw_idx, "p%d: index %d was not computed! (Emax_plus_hw_idx=%d)",mpi->get_rank(),indices[ii],Emax_plus_hw_idx);

                int local_idx = indices[ii] - (this->E0_plus_hw_idx+shift) + nE_self + nE_below;
                if (indices[ii]<=this->Emax_idx) {
                    // energy index index[ii] belongs to the own MPI process
                    // local_idx = nE_below + indices[ii]-this->E0_idx;
                    NEGF_ASSERT(local_idx == int(nE_below + indices[ii]-this->E0_idx), "index mess.");
                }
                NEGF_FASSERT(local_idx < int(AL.size()), "invalid local_idx (+hw): local_idx=%d, AL.size()=%d", local_idx, AL.size());
                NEGF_FASSERT(control[local_idx]==indices[ii], "p%d control failed (3): control[local_idx]=%d, indices[%d]=%d.",mpi->get_rank(),control[local_idx],ii,indices[ii]);

                add(this->AL[local_idx][kk], (coeffs[ii] * (this->Nphonon+1) * this->prefactor * this->scaling), SL);
                // SL +=  ...(N+1) * AL[local_idx][kk];
                add(this->AG[local_idx][kk], (coeffs[ii] * (this->Nphonon)   * this->prefactor * this->scaling), SG);
                // SG +=    ...(N) * AG[local_idx][kk];

                if (complicated_retarded) {
                    add(this->AR[local_idx][kk], (coeffs[ii] * this->Nphonon * this->prefactor * this->scaling), SR);
                    // SR +=    ...(N) * AR[local_idx][kk];
                    add(this->AL[local_idx][kk], -(coeffs[ii] * 0.5          * this->prefactor * this->scaling), SR); // minus sign!
                    // SR -=  ...(0.5) * AL[local_idx][kk];
                }
                if (nemo_retarded) {
                    add(this->AR[local_idx][kk], (coeffs[ii] * this->Nphonon * this->prefactor * this->scaling), SR);
                    // SR +=    ...(N) * AR[local_idx][kk];
                }
			}
		}
	}
	mpi->synchronize_processes();
	
	// security check: check Sigma = -Sigma+
	if (security_checking) {
		logmsg->emit(LOG_INFO_L2,"Security check...");
		SEMat SLtmp = SPOPMat_create(NxNn);
		SEMat SGtmp = SPOPMat_create(NxNn);
		for (uint ee2=0; ee2<this->myNE; ee2++)
		{
			uint ee  = energies->get_global_index(ee2);
			for (uint kk=0; kk<Nk; kk++)
			{
				const SEMat & SL = this->get_lesser(kk,ee);
				const SEMat & SG = this->get_greater(kk,ee);
				conjtrans(SL, SLtmp); // SLtmp = conjugateTranspose(SL); 
				SLtmp += SL;
				conjtrans(SG, SGtmp); // SGtmp = conjugateTranspose(SG); 
				SGtmp += SG;
				double SLnorm = negf_math::matrix_norm(SLtmp);
				double SGnorm = negf_math::matrix_norm(SGtmp);
				NEGF_FASSERT(SLnorm < constants::antiherm_check, "SL was not anti-hermitian: |SL-(-SL+)|=%e",SLnorm);
				NEGF_FASSERT(SGnorm < constants::antiherm_check, "SG was not anti-hermitian: |SG-(-SG+)|=%e",SGnorm);
			}
		}
	}
);}

void SEOpticalPhonon::communicate_As()
{//STACK_TRACE(
	logmsg->emit_noendl(LOG_INFO,"communicating As: ");
	int root = constants::mpi_master_rank;
	int my_rank = mpi->get_rank();
	
	vector< vector<int> > processes_needed; 
	vector<int> pp_E0_minus_hw_idx; 	pp_E0_minus_hw_idx  .resize(mpi->get_num_procs(), -1);
	vector<int> pp_Emax_minus_hw_idx; 	pp_Emax_minus_hw_idx.resize(mpi->get_num_procs(), -1);
	vector<int> pp_E0_idx;				pp_E0_idx           .resize(mpi->get_num_procs(), -1);
	vector<int> pp_Emax_idx;			pp_Emax_idx         .resize(mpi->get_num_procs(), -1);
	vector<int> pp_E0_plus_hw_idx;	    pp_E0_plus_hw_idx   .resize(mpi->get_num_procs(), -1);
	vector<int> pp_Emax_plus_hw_idx;	pp_Emax_plus_hw_idx .resize(mpi->get_num_procs(), -1);
	processes_needed.resize(mpi->get_num_procs());
	if (my_rank==root) 
	{
		for (int pp=0; pp < mpi->get_num_procs(); pp++) 
		{
			// process pp needs energies E0_minus_hw_idx...Emax_minus_hw_idx  (possibly not the last one)
		    //                   and     E0_plus_hw_idx...Emax_plus_hw_idx    (possibly not the first one)
			
			// get these variables
			if (pp==my_rank) {
				pp_E0_minus_hw_idx[pp]  = this->E0_minus_hw_idx;
				pp_Emax_minus_hw_idx[pp]= this->Emax_minus_hw_idx;
				pp_E0_idx[pp]           = this->E0_idx;
				pp_Emax_idx[pp]         = this->Emax_idx;
				pp_E0_plus_hw_idx[pp]   = this->E0_plus_hw_idx;
				pp_Emax_plus_hw_idx[pp] = this->Emax_plus_hw_idx;
			} else {
				int source = pp;
				int tag = 1; mpi->recv(pp_E0_minus_hw_idx[pp],  source, tag);
				    tag = 2; mpi->recv(pp_Emax_minus_hw_idx[pp],source, tag);
				    tag = 3; mpi->recv(pp_E0_idx[pp],           source, tag);
				    tag = 4; mpi->recv(pp_Emax_idx[pp],         source, tag);
				    tag = 5; mpi->recv(pp_E0_plus_hw_idx[pp],   source, tag);
				    tag = 6; mpi->recv(pp_Emax_plus_hw_idx[pp], source, tag);
			}
			
			// determine which processes compute the energies in question
			// and add them to the list
			for (int ee=pp_E0_minus_hw_idx[pp]; ee <= pp_Emax_minus_hw_idx[pp]; ee++) 
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
			for (int ee=pp_E0_plus_hw_idx[pp]; ee <= pp_Emax_plus_hw_idx[pp]; ee++) 
			{
				int pp2 = energies->get_process_computing(uint(ee));
				if (pp2==pp) continue; // this is possible, ee could be Emax_idx!
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
			logmsg->emit_noendl_all(LOG_INFO_L3,"p%d needs ",pp);
			for (uint ii=0; ii < processes_needed[pp].size(); ii++) {
				logmsg->emit_noendl_all(LOG_INFO_L3,"p%d, ",processes_needed[pp][ii]);
			}
			logmsg->emit_all(LOG_INFO_L3,"");
		}
	} else {
		int E0_idx_copy           = this->E0_idx;
		int Emax_idx_copy         = this->Emax_idx;
		int E0_minus_hw_idx_copy  = this->E0_minus_hw_idx;
		int Emax_minus_hw_idx_copy= this->Emax_minus_hw_idx;
		int E0_plus_hw_idx_copy   = this->E0_plus_hw_idx;
		int Emax_plus_hw_idx_copy = this->Emax_plus_hw_idx;
		int tag = 1; mpi->send(E0_minus_hw_idx_copy,  root, tag);
			tag = 2; mpi->send(Emax_minus_hw_idx_copy,root, tag);
			tag = 3; mpi->send(E0_idx_copy,           root, tag);
			tag = 4; mpi->send(Emax_idx_copy,         root, tag);
			tag = 5; mpi->send(E0_plus_hw_idx_copy,   root, tag);
			tag = 6; mpi->send(Emax_plus_hw_idx_copy, root, tag);
	}
	
	logmsg->emit(LOG_INFO_L2,"Broadcasting processes_needed...    ");
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
	logmsg->emit(LOG_INFO_L2,"Broadcasting energy index arrays...   ");
	mpi->broadcast(pp_E0_minus_hw_idx  , root);
	mpi->broadcast(pp_Emax_minus_hw_idx, root);
	mpi->broadcast(pp_E0_idx           , root);
	mpi->broadcast(pp_Emax_idx         , root);
	mpi->broadcast(pp_E0_plus_hw_idx   , root);
	mpi->broadcast(pp_Emax_plus_hw_idx , root);
		
	mpi->synchronize_processes();
	
	// ---------------------------------------------------------------------------------
	// create an array of zipped data (for every energy all k's) for all OWN energies
	// ---------------------------------------------------------------------------------
	logmsg->emit(LOG_INFO_L2,"Creating zipped data arrays...   ");
	//uint myNE = nE_self;
	NEGF_ASSERT(this->myNE==energies->get_my_number_of_points() && this->myNE==this->nE_self && this->myNE==this->Emax_idx-this->E0_idx+1, "inconsistent myNE.");
	double * AL_real_data[myNE];			// will store copies of the complex AL split into real and imag parts
	double * AG_real_data[myNE];
	double * AR_real_data[myNE];
	double * AL_imag_data[myNE];
	double * AG_imag_data[myNE];
	double * AR_imag_data[myNE];
	unsigned char * AL_real_char[myNE];		// will store pointer to same memory as AL_...._data, but different type
	unsigned char * AG_real_char[myNE];
	unsigned char * AR_real_char[myNE];
	unsigned char * AL_imag_char[myNE];
	unsigned char * AG_imag_char[myNE];
	unsigned char * AR_imag_char[myNE];
	unsigned long num_chars[myNE];			// will store how many chars the array of complex matrices corresponds to
	unsigned long num_chars_AR[myNE];
	unsigned long AL_real_comp_size[myNE];	// will store the size of the zipped data
	unsigned long AG_real_comp_size[myNE];
	unsigned long AR_real_comp_size[myNE];
	unsigned long AL_imag_comp_size[myNE];
	unsigned long AG_imag_comp_size[myNE];
	unsigned long AR_imag_comp_size[myNE];
	unsigned char * AL_real_compressed[myNE];	// will store pointers to the zipped data
	unsigned char * AG_real_compressed[myNE];
	unsigned char * AR_real_compressed[myNE];
	unsigned char * AL_imag_compressed[myNE];
	unsigned char * AG_imag_compressed[myNE];
	unsigned char * AR_imag_compressed[myNE];
	for (uint ee2=0; ee2<myNE; ee2++) 
	{
		uint ee = energies->get_global_index(ee2);
		
		// compress AL[ee] and AG[ee]	(anti-hermitian)
		NEGF_FASSERT(control[ee2+nE_below]==ee, "control failed (4): control[ee2+nE_below]=%d, ee=%d.",control[ee2+nE_below],ee);	
#ifdef ANTIHERM
		negf_math::do_compress_antiherm(AL[ee2+nE_below], AL_real_data[ee2], AL_imag_data[ee2], AL_real_char[ee2], AL_imag_char[ee2], num_chars[ee2], 
#else
		negf_math::do_compress(AL[ee2+nE_below], AL_real_data[ee2], AL_imag_data[ee2], AL_real_char[ee2], AL_imag_char[ee2], num_chars[ee2],
#endif
		        AL_real_comp_size[ee2], AL_imag_comp_size[ee2], AL_real_compressed[ee2], AL_imag_compressed[ee2]);

		unsigned long tmp = num_chars[ee2];
#ifdef ANTIHERM
		negf_math::do_compress_antiherm(AG[ee2+nE_below], AG_real_data[ee2], AG_imag_data[ee2], AG_real_char[ee2], AG_imag_char[ee2], num_chars[ee2], 
#else
		negf_math::do_compress(AG[ee2+nE_below], AG_real_data[ee2], AG_imag_data[ee2], AG_real_char[ee2], AG_imag_char[ee2], num_chars[ee2],
#endif
		        AG_real_comp_size[ee2], AG_imag_comp_size[ee2], AG_real_compressed[ee2], AG_imag_compressed[ee2]);
		NEGF_ASSERT(num_chars[ee2]==tmp, "something went wrong: num_chars[ee2]!=tmp.");

		if (complicated_retarded || nemo_retarded) {			
			negf_math::do_compress(AR[ee2+nE_below], AR_real_data[ee2], AR_imag_data[ee2], AR_real_char[ee2], AR_imag_char[ee2], num_chars_AR[ee2], 
					AR_real_comp_size[ee2], AR_imag_comp_size[ee2], AR_real_compressed[ee2], AR_imag_compressed[ee2]);	
		}
	}	
	mpi->synchronize_processes();
	
	// --------------------------------------------------------------------------------------
	// determine how much data is sent in total
	// we do this by going through the same algorithm as when sending the data, just w/o sending
	// necessary for the amount of buffer that needs to be allocated for a safe operation
	// --------------------------------------------------------------------------------------
	logmsg->emit(LOG_INFO_L2,"Computing total amount of data to be sent...    ");
	NEGF_FASSERT(this->Emax_idx-this->E0_idx+1==myNE,"Emax_idx=%d, E0_idx=%d, myNE=%d",this->Emax_idx,this->E0_idx,myNE);
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
			if (!i_am_sender) { // this is possible!
			} else {
				
			// PART 1: sender sends his A's to processes storing energies BELOW his energy interval
			for (uint ee=this->E0_idx; ee<=this->Emax_idx; ee++) 
			{			
				// determine receiver(s) BELOW current process - if there is none, continue to next energy
				for (int pp=0; pp<mpi->get_num_procs(); pp++) 
				{
					bool pp_receives = receiver[pp] && pp_E0_plus_hw_idx[pp]<=int(ee) && pp_Emax_plus_hw_idx[pp]>=int(ee);
					if (!pp_receives) continue;				
				
					// size of AX_xxxx_comp_size is only myNE!
					buffer_size_needed += AL_real_comp_size[ee-this->E0_idx] + AL_imag_comp_size[ee-this->E0_idx] + 2*sizeof(unsigned long)+ 8*MPI_BSEND_OVERHEAD;
					buffer_size_needed += AG_real_comp_size[ee-this->E0_idx] + AG_imag_comp_size[ee-this->E0_idx] + 2*sizeof(unsigned long)+ 8*MPI_BSEND_OVERHEAD;
					if (complicated_retarded || nemo_retarded) {
						buffer_size_needed += AR_real_comp_size[ee-this->E0_idx] + AR_imag_comp_size[ee-this->E0_idx] + 2*sizeof(unsigned long)+ 8*MPI_BSEND_OVERHEAD;
					}
				}			
			}
			
			// PART 2: sender sends his A's to processes storing energies ABOVE his energy interval
			for (uint ee=this->E0_idx; ee<=this->Emax_idx; ee++) 
			{			
				// determine receiver ABOVE current process - if there is none, continue to next energy
				for (int pp=0; pp<mpi->get_num_procs(); pp++) {
					bool pp_receives = receiver[pp] && pp_E0_minus_hw_idx[pp]<=int(ee) && pp_Emax_minus_hw_idx[pp]>=int(ee);
					if (!pp_receives) continue;
					
					// size of AX_xxxx_comp_size is only myNE!
					buffer_size_needed += AL_real_comp_size[ee-this->E0_idx] + AL_imag_comp_size[ee-this->E0_idx] + 2*sizeof(unsigned long) + 8*MPI_BSEND_OVERHEAD;
					buffer_size_needed += AG_real_comp_size[ee-this->E0_idx] + AG_imag_comp_size[ee-this->E0_idx] + 2*sizeof(unsigned long) + 8*MPI_BSEND_OVERHEAD;
					if (complicated_retarded || nemo_retarded) {
						buffer_size_needed += AR_real_comp_size[ee-this->E0_idx] + AR_imag_comp_size[ee-this->E0_idx] + 2*sizeof(unsigned long) + 8*MPI_BSEND_OVERHEAD;
					}
				}		
			}
				
			} // if(i_am_sender)
		}
	
		// mark receivers as computed. needs to be performed in ALL threads (senders, receivers and those which are neither)
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
	logmsg->emit(LOG_INFO_L2,"Allocating MPI buffer...   ");
	unsigned long buffersize_long = buffer_size_needed+1000;
	unsigned long max_array_size = 2147483648UL;
	NEGF_FASSERT(buffersize_long < max_array_size, "Buffer size (%d) does not fit into an int!!@!",buffersize_long);
	int buffersize = int(buffersize_long);
	char * buffer = new char[buffersize];
	int err = MPI_Buffer_attach(buffer, buffersize);
	NEGF_FASSERT(err==0, "MPI_Buffer_attach gave error %d",err);
	mpi->synchronize_processes();
	
	// ------------------------------
	// do MPI communication!
	// ------------------------------
	logmsg->emit(LOG_INFO_L2,"Communicate!");
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
		//mpi->synchronize_processes();
		
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
            // ----------------------------------------------------------------------------------
            // RECEIVER
            // ----------------------------------------------------------------------------------

			// PART 1: receiver receives missing A's ABOVE his energy interval
			for (uint ee=this->E0_plus_hw_idx; ee<=this->Emax_plus_hw_idx; ee++) 
			{
				int sender_id = energies->get_process_computing(ee);
				if (sender_id==my_rank) { continue;	}	// this is possible, ee could be Emax_idx!
				NEGF_FASSERT(sender[sender_id] && !receiver[sender_id],"p%d expected sender %d but something went wrong.", my_rank, sender_id);
				int tag = ee;
				
				int local_idx = ee - (this->E0_plus_hw_idx+shift) + nE_below + nE_self;
				NEGF_FASSERT(local_idx < int(this->AL.size()) && local_idx < int(this->AG.size()),"set up ALG first: local_idx=%d, AL.size()=%d.",local_idx, this->AL.size());
				NEGF_FASSERT(control[local_idx]==ee, "p%d control failed (5): control[%d]=%d, ee=%d.",mpi->get_rank(),local_idx,control[local_idx],ee);		
				vector<SEMat> & ALmat = this->AL[local_idx];
				vector<SEMat> & AGmat = this->AG[local_idx];
				for (uint kk=0; kk<Nk; kk++) {
					NEGF_FASSERT(ALmat[kk].num_rows()==NxNn && ALmat[kk].num_cols()==NxNn, "ee=%d,kk=%d: AL is only a %dx%d matrix",ee,kk,ALmat[kk].num_rows(),ALmat[kk].num_cols());
					NEGF_FASSERT(AGmat[kk].num_rows()==NxNn && AGmat[kk].num_cols()==NxNn, "ee=%d,kk=%d: AG is only a %dx%d matrix",ee,kk,AGmat[kk].num_rows(),AGmat[kk].num_cols());
				}
				
				logmsg->emit(LOG_INFO_L3,"p%d waits for ee=%d from p%d", my_rank, ee, sender_id);
#ifdef ANTIHERM
				mpi->recv_antihermitians(ALmat, sender_id, tag);
                mpi->recv_antihermitians(AGmat, sender_id, tag);
#else
				mpi->recv(ALmat, sender_id, tag);
				mpi->recv(AGmat, sender_id, tag);
#endif
				if (complicated_retarded || nemo_retarded) { 
					NEGF_FASSERT(local_idx < int(this->AR.size()),"set up AR first: local_idx=%d, AR.size()=%d.",local_idx, this->AR.size());
					vector<SEMat> & ARmat = this->AR[local_idx];
					mpi->recv(ARmat, sender_id, tag); // NOT antihermitian
				}
				logmsg->emit(LOG_INFO_L3,"p%d just got ee=%d from p%d", my_rank, ee, sender_id);
			}
			//mpi->synchronize_processes();
			//cout << "PHASE 1 DONE.";
			
			// PART 2: receiver receives missing A's BELOW his energy interval
			for (uint ee=this->E0_minus_hw_idx; ee<=this->Emax_minus_hw_idx; ee++) 
			{
				int sender_id = energies->get_process_computing(ee);
				if (sender_id==my_rank) { continue; } // this is possible, ee could be E0_idx!
				NEGF_FASSERT(sender[sender_id] && !receiver[sender_id],"p%d expected sender %d but something went wrong.", my_rank, sender_id);
				int tag = ee;
				
				int local_idx = ee - this->E0_minus_hw_idx;
				NEGF_FASSERT(local_idx < int(this->AL.size()) && local_idx < int(this->AG.size()),"set up ALG first. local_idx=%d, AL.size()=%d.",local_idx, this->AL.size());
				NEGF_FASSERT(control[local_idx]==ee, "p%d control failed (6): control[%d]=%d, ee=%d.",mpi->get_rank(),local_idx,control[local_idx],ee);	
				vector<SEMat> & ALmat = this->AL[local_idx];
				vector<SEMat> & AGmat = this->AG[local_idx];
				for (uint kk=0; kk<Nk; kk++) {
					NEGF_FASSERT(ALmat[kk].num_rows()==NxNn && ALmat[kk].num_cols()==NxNn, "ee=%d,kk=%d: AL is only a %dx%d matrix",ee,kk,ALmat[kk].num_rows(),ALmat[kk].num_cols());
					NEGF_FASSERT(AGmat[kk].num_rows()==NxNn && AGmat[kk].num_cols()==NxNn, "ee=%d,kk=%d: AG is only a %dx%d matrix",ee,kk,AGmat[kk].num_rows(),AGmat[kk].num_cols());
				}
				
				logmsg->emit(LOG_INFO_L3,"p%d waits for ee=%d from p%d", my_rank, ee, sender_id);
#ifdef ANTIHERM
				mpi->recv_antihermitians(ALmat, sender_id, tag);
                mpi->recv_antihermitians(AGmat, sender_id, tag);
#else
				mpi->recv(ALmat, sender_id, tag);
				mpi->recv(AGmat, sender_id, tag);
#endif
				
				if (complicated_retarded || nemo_retarded) {
					NEGF_FASSERT(local_idx < int(this->AR.size()), "set up AR first. local_idx=%d, AR.size()=%d.",local_idx, this->AR.size());
					vector<SEMat> & ARmat = this->AR[local_idx];
					mpi->recv(ARmat, sender_id, tag); // NOT antihermitian
				}
				logmsg->emit(LOG_INFO_L3,"p%d just got ee=%d from p%d", my_rank, ee, sender_id);
			}
			
		} else 
		{
			if (!i_am_sender) { // this is possible!
				//cout << "DONEp" << my_rank << endl;
			} else {


			// ----------------------------------------------------------------------------------
			// SENDER
			// ----------------------------------------------------------------------------------
		
			// PART 1: sender sends his A's to processes storing energies BELOW his energy interval
			for (uint ee=this->E0_idx; ee<=this->Emax_idx; ee++) 
			{			
				// determine receiver(s) BELOW current process - if there is none, continue to next energy
				for (int pp=0; pp<mpi->get_num_procs(); pp++) 
				{
					bool pp_receives = receiver[pp] && pp_E0_plus_hw_idx[pp]<=int(ee) && pp_Emax_plus_hw_idx[pp]>=int(ee);
					int receiver_id= pp;
					
					if (!pp_receives) continue;
					
					int tag = ee;
					int tag2 = tag+1;
										
					logmsg->emit(LOG_INFO_L3,"p%d sends ee=%d to p%d", my_rank, ee, receiver_id);
										
					// BUFFERED SEND!
					// mpi->send_antihermitians and mpi->recv_antihermitians send/receive the following 
					// quantities in order: real_char_size, imag_char_size, real_compressed, imag_compressed
					
					// ----------------------
					// send AL
					// ----------------------
					// send array lengths
					int real_char_size = int(AL_real_comp_size[ee-this->E0_idx]); // NOT +nE_below!
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
					
					// ----------------------
					// send AG
					// ----------------------
					// send array lengths
					real_char_size = int(AG_real_comp_size[ee-this->E0_idx]); // NOT +nE_below!
					imag_char_size = int(AG_imag_comp_size[ee-this->E0_idx]);
					err = MPI_Bsend(&real_char_size, 1, MPI_INT, receiver_id, tag, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of AG_real_char_size failed: err=%d",err);
					err = MPI_Bsend(&imag_char_size, 1, MPI_INT, receiver_id, tag2, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of AG_imag_char_size failed: err=%d",err);
		
					// send arrays
					err = MPI_Bsend(AG_real_compressed[ee-this->E0_idx], real_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of AG_real_compressed failed: err=%d",err);
					err = MPI_Bsend(AG_imag_compressed[ee-this->E0_idx], imag_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag2, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of AG_imag_compressed failed: err=%d",err);
					
					if (complicated_retarded || nemo_retarded) {						
						// ----------------------
						// send AR (mpi->send())
						// ----------------------
						// send array lengths
						real_char_size = int(AR_real_comp_size[ee-this->E0_idx]); // NOT +nE_below!
						imag_char_size = int(AR_imag_comp_size[ee-this->E0_idx]);
						err = MPI_Bsend(&real_char_size, 1, MPI_INT, receiver_id, tag, MPI_COMM_WORLD);
						NEGF_FASSERT(err==0, "MPI send of AG_real_char_size failed: err=%d",err);
						err = MPI_Bsend(&imag_char_size, 1, MPI_INT, receiver_id, tag2, MPI_COMM_WORLD);
						NEGF_FASSERT(err==0, "MPI send of AG_imag_char_size failed: err=%d",err);
			
						// send arrays
						err = MPI_Bsend(AR_real_compressed[ee-this->E0_idx], real_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag, MPI_COMM_WORLD);
						NEGF_FASSERT(err==0, "MPI send of AG_real_compressed failed: err=%d",err);
						err = MPI_Bsend(AR_imag_compressed[ee-this->E0_idx], imag_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag2, MPI_COMM_WORLD);
						NEGF_FASSERT(err==0, "MPI send of AG_imag_compressed failed: err=%d",err);
					}
					
					logmsg->emit(LOG_INFO_L3,"p%d just sent ee=%d to p%d", my_rank, ee, receiver_id);
				
				}			
			}
			
			//mpi->synchronize_processes();
			//cout << "PHASE 1 DONE.";
			
			// PART 2: sender sends his A's to processes storing energies ABOVE his energy interval
			for (uint ee=this->E0_idx; ee<=this->Emax_idx; ee++) 
			{			
				// determine receiver ABOVE current process - if there is none, continue to next energy
				for (int pp=0; pp<mpi->get_num_procs(); pp++) {
					bool pp_receives = receiver[pp] && pp_E0_minus_hw_idx[pp]<=int(ee) && pp_Emax_minus_hw_idx[pp]>=int(ee);
					if (!pp_receives) continue;
					int receiver_id= pp;
					
					logmsg->emit(LOG_INFO_L3,"p%d sends ee=%d to p%d", my_rank, ee, receiver_id);
					
					int tag = ee;
					int tag2 = tag+1;
					
					// *** BUFFERED SEND! ***
					// mpi->send_antihermitians and mpi->recv_antihermitians send/receive the following 
					// quantities in order: real_char_size, imag_char_size, real_compressed, imag_compressed
					
					// ----------------------
					// send AL
					// ----------------------
					// send array lengths
					int real_char_size = int(AL_real_comp_size[ee-this->E0_idx]); // NOT +nE_below!
					int imag_char_size = int(AL_imag_comp_size[ee-this->E0_idx]);
					NEGF_FASSERT(real_char_size>0, "(1) real_char_size was <=0: %d",real_char_size);
					NEGF_FASSERT(imag_char_size>0, "(1) imag_char_size was <=0: %d",imag_char_size);
					err = MPI_Bsend(&real_char_size, 1, MPI_INT, receiver_id, tag, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of AL_real_char_size failed: err=%d",err);
					err = MPI_Bsend(&imag_char_size, 1, MPI_INT, receiver_id, tag2, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of AL_imag_char_size failed: err=%d",err);
		
					// send arrays
					err = MPI_Bsend(AL_real_compressed[ee-this->E0_idx], real_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of AL_real_compressed failed: err=%d",err);
					err = MPI_Bsend(AL_imag_compressed[ee-this->E0_idx], imag_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag2, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of AL_imag_compressed failed: err=%d",err);
					
					// ----------------------
					// send AG
					// ----------------------
					// send array lengths
					real_char_size = int(AG_real_comp_size[ee-this->E0_idx]); // NOT +nE_below!
					imag_char_size = int(AG_imag_comp_size[ee-this->E0_idx]);
					NEGF_FASSERT(real_char_size>0, "(2) real_char_size was <=0: %d",real_char_size);
					NEGF_FASSERT(imag_char_size>0, "(2) imag_char_size was <=0: %d",imag_char_size);
					err = MPI_Bsend(&real_char_size, 1, MPI_INT, receiver_id, tag, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of AG_real_char_size failed: err=%d",err);
					err = MPI_Bsend(&imag_char_size, 1, MPI_INT, receiver_id, tag2, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of AG_imag_char_size failed: err=%d",err);
		
					// send arrays
					err = MPI_Bsend(AG_real_compressed[ee-this->E0_idx], real_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of AG_real_compressed failed: err=%d",err);
					err = MPI_Bsend(AG_imag_compressed[ee-this->E0_idx], imag_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag2, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of AG_imag_compressed failed: err=%d",err);
					
					if (complicated_retarded || nemo_retarded) {
						// ----------------------
						// send AR 
						// ----------------------
						// send array lengths
						real_char_size = int(AR_real_comp_size[ee-this->E0_idx]); // NOT +nE_below!
						imag_char_size = int(AR_imag_comp_size[ee-this->E0_idx]);
						err = MPI_Bsend(&real_char_size, 1, MPI_INT, receiver_id, tag, MPI_COMM_WORLD);
						NEGF_FASSERT(err==0, "MPI send of AG_real_char_size failed: err=%d",err);
						err = MPI_Bsend(&imag_char_size, 1, MPI_INT, receiver_id, tag2, MPI_COMM_WORLD);
						NEGF_FASSERT(err==0, "MPI send of AG_imag_char_size failed: err=%d",err);
			
						// send arrays
						err = MPI_Bsend(AR_real_compressed[ee-this->E0_idx], real_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag, MPI_COMM_WORLD);
						NEGF_FASSERT(err==0, "MPI send of AG_real_compressed failed: err=%d",err);
						err = MPI_Bsend(AR_imag_compressed[ee-this->E0_idx], imag_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag2, MPI_COMM_WORLD);
						NEGF_FASSERT(err==0, "MPI send of AG_imag_compressed failed: err=%d",err);
					}			
					
					logmsg->emit(LOG_INFO_L3,"p%d just sent ee=%d to p%d", my_rank, ee, receiver_id);
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
	logmsg->emit(LOG_INFO_L2,"Deallocating MPI buffer");
	err = MPI_Buffer_detach(&buffer, &buffersize);
	delete [] buffer;
	NEGF_FASSERT(err==0, "MPI_Buffer_attach gave error %d",err);
	
	// ----------------------------------------------------------------
	// release zipped data
	// ----------------------------------------------------------------
	logmsg->emit(LOG_INFO_L2,"Deallocating zipped arrays");
	for (uint ee2=0; ee2<this->myNE; ee2++) {
		delete [] AL_real_data[ee2];
		delete [] AL_imag_data[ee2];
		delete [] AG_real_data[ee2];
		delete [] AG_imag_data[ee2];
		
		delete [] AL_real_compressed[ee2];
		delete [] AL_imag_compressed[ee2];
		delete [] AG_real_compressed[ee2];
		delete [] AG_imag_compressed[ee2];
		
		if (complicated_retarded || nemo_retarded) {
			delete [] AR_real_data[ee2];
			delete [] AR_imag_data[ee2];
			
			delete [] AR_real_compressed[ee2];
			delete [] AR_imag_compressed[ee2];
		}
	}
	mpi->synchronize_processes();
/*);*/}


void SEOpticalPhonon::calculate_retarded()
{STACK_TRACE(
	logmsg->emit_small_header("calculating retarded optical phonon self-energy");
	NEGF_ASSERT(!complicated_retarded && !nemo_retarded, "Do not call this method when GR using complicated_retarded or nemo_retarded!");
	
	// possibility for the user to skip SR entirely (set LuisierSRpop to -1!)
	if (options->exists("LuisierSRpop") && options->get("LuisierSRpop")==-1) {
		// retarded SE were initialized to zero in the beginning
		return;
	}
	
	for (uint ee2=0; ee2<this->myNE; ee2++)
	{
		const uint ee = energies->get_global_index(ee2);
		for (uint kk=0; kk<Nk; kk++) 
		{
			const SEMat & SL = this->get_lesser(kk,ee);
			const SEMat & SG = this->get_greater(kk,ee);
			      SEMat & SR = this->get_retarded(kk,ee);
			
			sub(SG, SL, SR); // SR = SG - SL;
			SR *= 0.5;
		}
	}
	// no principal value!
);}


void SEOpticalPhonon::set_scaling(double new_scaling)
{STACK_TRACE(
	NEGF_ASSERT(new_scaling>=0.0 && new_scaling<=1.0, "bad scaling factor"); 
	logmsg->emit(LOG_INFO,"Setting scaling of optical phonon scattering to %.4g", new_scaling);
	this->scaling = new_scaling;
);}


void SEOpticalPhonon::output_debug_info()
{STACK_TRACE(
	logmsg->emit_header("Debug info about GL, GG and SigmaL, SigmaG of POP interaction");
	
	const double Ethresh = 0.0;

	// -----------------------------------------
	// compute stuff for own energies
	// -----------------------------------------
	logmsg->emit(LOG_INFO,"Computing...");
	vector<double> my_GL_neg; my_GL_neg.resize(4, 0.0);
	vector<double> my_GL_pos; my_GL_pos.resize(4, 0.0);
	vector<double> my_GG_neg; my_GG_neg.resize(4, 0.0);
	vector<double> my_GG_pos; my_GG_pos.resize(4, 0.0);
	vector<double> my_SL_neg; my_SL_neg.resize(4, 0.0);
	vector<double> my_SL_pos; my_SL_pos.resize(4, 0.0);
	vector<double> my_SG_neg; my_SG_neg.resize(4, 0.0);
	vector<double> my_SG_pos; my_SG_pos.resize(4, 0.0);
	vector<double> my_SLGG_neg; my_SLGG_neg.resize(4, 0.0);
	vector<double> my_SLGG_pos; my_SLGG_pos.resize(4, 0.0);
	vector<double> my_SGGL_neg; my_SGGL_neg.resize(4, 0.0);
	vector<double> my_SGGL_pos; my_SGGL_pos.resize(4, 0.0);
	Matc SLGG(NxNn,NxNn);
	Matc SGGL(NxNn,NxNn);
	for (uint ee2 = 0; ee2 < myNE; ee2++) 
	{
		uint ee = energies->get_global_index(ee2);
		double E = energies->get_energy_from_global_idx(ee);
		//if (ee % 5 == 0) logmsg->emit_noendl(LOG_INFO_L3, "p%d: GR(E=%d,:)...   ",mpi->get_rank(),ee);
		for (uint kk=0; kk<Nk; kk++) 
		{
			vector<double> tmp; tmp.resize(4, 0.0);
			
#ifdef USE_BANDED
			Matc GLm(NxNn,NxNn); GLm = gf->get_lesser(kk,ee);
#else
			const Matc & GLm = gf->get_lesser(kk,ee);
#endif
			this->get_CB_VB_norms(GLm, tmp); 
			if (E<Ethresh) { for(uint ii=0; ii<4; ii++) { my_GL_neg[ii] += tmp[ii]; } } else { for(uint ii=0; ii<4; ii++) { my_GL_pos[ii] += tmp[ii]; } }
			
#ifdef USE_BANDED
			Matc GGm(NxNn,NxNn); GGm = gf->get_greater(kk,ee);
#else
			const Matc & GGm = gf->get_greater(kk,ee);
#endif
			this->get_CB_VB_norms(GGm, tmp); 
			if (E<Ethresh) { for(uint ii=0; ii<4; ii++) { my_GG_neg[ii] += tmp[ii]; } } else { for(uint ii=0; ii<4; ii++) { my_GG_pos[ii] += tmp[ii]; } }
			
#ifdef USE_BANDED
			Matc SLm(NxNn,NxNn); SLm = this->get_lesser(kk,ee);
#else
			const Matc & SLm = this->get_lesser(kk,ee);
#endif
			this->get_CB_VB_norms(SLm, tmp); 
			if (E<Ethresh) { for(uint ii=0; ii<4; ii++) { my_SL_neg[ii] += tmp[ii]; } } else { for(uint ii=0; ii<4; ii++) { my_SL_pos[ii] += tmp[ii]; } }
			
#ifdef USE_BANDED
			Matc SGm(NxNn,NxNn); SGm = this->get_greater(kk,ee);
#else
			const Matc & SGm = this->get_greater(kk,ee);
#endif
			this->get_CB_VB_norms(SGm, tmp); 
			if (E<Ethresh) { for(uint ii=0; ii<4; ii++) { my_SG_neg[ii] += tmp[ii]; } } else { for(uint ii=0; ii<4; ii++) { my_SG_pos[ii] += tmp[ii]; } }
			
			mult(SLm, GGm, SLGG); // SLGG = SLm * GGm;
			this->get_CB_VB_norms(SLGG, tmp); 
			if (E<Ethresh) { for(uint ii=0; ii<4; ii++) { my_SLGG_neg[ii] += tmp[ii]; } } else { for(uint ii=0; ii<4; ii++) { my_SLGG_pos[ii] += tmp[ii]; } }
			
			mult(SGm, GLm, SGGL); // SGGL = SGm * GLm;
			this->get_CB_VB_norms(SGGL, tmp); 
			if (E<Ethresh) { for(uint ii=0; ii<4; ii++) { my_SGGL_neg[ii] += tmp[ii]; } } else { for(uint ii=0; ii<4; ii++) { my_SGGL_pos[ii] += tmp[ii]; } }
		}
	}
	mpi->synchronize_processes();
	
	// ------------------------------------------------
	// communicate to master process
	// ------------------------------------------------
	logmsg->emit(LOG_INFO,"Aggregating in master thread...");
	if (mpi->get_rank()==constants::mpi_master_rank) 
	{
		vector<double> total_GL_neg = my_GL_neg;
		vector<double> total_GL_pos = my_GL_pos;
		vector<double> total_GG_neg = my_GG_neg;
		vector<double> total_GG_pos = my_GG_pos;
		vector<double> total_SL_neg = my_SL_neg;
		vector<double> total_SL_pos = my_SL_pos;
		vector<double> total_SG_neg = my_SG_neg;
		vector<double> total_SG_pos = my_SG_pos;
		vector<double> total_SLGG_neg = my_SLGG_neg;
		vector<double> total_SLGG_pos = my_SLGG_pos;
		vector<double> total_SGGL_neg = my_SGGL_neg;
		vector<double> total_SGGL_pos = my_SGGL_pos;
					
		// collect the pieces
		vector<double> tmp; tmp.resize(4, 0.0);
		for (int pp=0; pp<mpi->get_num_procs(); pp++) 
		{
			if (pp==constants::mpi_master_rank) continue;
			
			// receive from other process
			logmsg->emit/*_all*/(LOG_INFO_L3,"Receiving from process %d...",pp);
			int tag = pp;
			uint size = 4;
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_GL_neg[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_GL_pos[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_GG_neg[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_GG_pos[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_SL_neg[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_SL_pos[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_SG_neg[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_SG_pos[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_SLGG_neg[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_SLGG_pos[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_SGGL_neg[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_SGGL_pos[ii] += tmp[ii]; }
			logmsg->emit/*_all*/(LOG_INFO_L3,"Receiving done.");
		}
		
		// screen output!
		logmsg->emit(LOG_INFO, "MATRIX NORMS (sum over E,k); E-range splits at %.3g",Ethresh);
		logmsg->emit(LOG_INFO, "quantity   E-range  CB-CB       CB-VB       VB-CB       VB-VB");
		logmsg->emit(LOG_INFO, "------------------------------------------------------------------");
		logmsg->emit(LOG_INFO, "GL         VB       %.3e   %.3e   %.3e   %.3e", total_GL_neg[0], total_GL_neg[1], total_GL_neg[2], total_GL_neg[3]);
		logmsg->emit(LOG_INFO, "GL         CB       %.3e   %.3e   %.3e   %.3e", total_GL_pos[0], total_GL_pos[1], total_GL_pos[2], total_GL_pos[3]);
		logmsg->emit(LOG_INFO, "GG         VB       %.3e   %.3e   %.3e   %.3e", total_GG_neg[0], total_GG_neg[1], total_GG_neg[2], total_GG_neg[3]);
		logmsg->emit(LOG_INFO, "GG         CB       %.3e   %.3e   %.3e   %.3e", total_GG_pos[0], total_GG_pos[1], total_GG_pos[2], total_GG_pos[3]);
		logmsg->emit(LOG_INFO, "SL         VB       %.3e   %.3e   %.3e   %.3e", total_SL_neg[0], total_SL_neg[1], total_SL_neg[2], total_SL_neg[3]);
		logmsg->emit(LOG_INFO, "SL         CB       %.3e   %.3e   %.3e   %.3e", total_SL_pos[0], total_SL_pos[1], total_SL_pos[2], total_SL_pos[3]);
		logmsg->emit(LOG_INFO, "SG         VB       %.3e   %.3e   %.3e   %.3e", total_SG_neg[0], total_SG_neg[1], total_SG_neg[2], total_SG_neg[3]);
		logmsg->emit(LOG_INFO, "SG         CB       %.3e   %.3e   %.3e   %.3e", total_SG_pos[0], total_SG_pos[1], total_SG_pos[2], total_SG_pos[3]);
		logmsg->emit(LOG_INFO, "SL*GG      VB       %.3e   %.3e   %.3e   %.3e", total_SLGG_neg[0], total_SLGG_neg[1], total_SLGG_neg[2], total_SLGG_neg[3]);
		logmsg->emit(LOG_INFO, "SL*GG      CB       %.3e   %.3e   %.3e   %.3e", total_SLGG_pos[0], total_SLGG_pos[1], total_SLGG_pos[2], total_SLGG_pos[3]);
		logmsg->emit(LOG_INFO, "SG*GL      VB       %.3e   %.3e   %.3e   %.3e", total_SGGL_neg[0], total_SGGL_neg[1], total_SGGL_neg[2], total_SGGL_neg[3]);
		logmsg->emit(LOG_INFO, "SG*GL      CB       %.3e   %.3e   %.3e   %.3e", total_SGGL_pos[0], total_SGGL_pos[1], total_SGGL_pos[2], total_SGGL_pos[3]);
	} else {
		// send to master process
		int dest = constants::mpi_master_rank;
		int tag = mpi->get_rank();
		mpi->send(my_GL_neg, dest, tag);
		mpi->send(my_GL_pos, dest, tag);
		mpi->send(my_GG_neg, dest, tag);
		mpi->send(my_GG_pos, dest, tag);
		mpi->send(my_SL_neg, dest, tag);
		mpi->send(my_SL_pos, dest, tag);
		mpi->send(my_SG_neg, dest, tag);
		mpi->send(my_SG_pos, dest, tag);
		mpi->send(my_SLGG_neg, dest, tag);
		mpi->send(my_SLGG_pos, dest, tag);
		mpi->send(my_SGGL_neg, dest, tag);
		mpi->send(my_SGGL_pos, dest, tag);
	}
	mpi->synchronize_processes();
);}


void SEOpticalPhonon::get_CB_VB_norms(const Matc & A, vector<double> & result)
{STACK_TRACE(
	NEGF_ASSERT(A.num_rows()==NxNn && A.num_cols()==NxNn && Nn==2, "expected A to be (2*Nx)^2 matrix");
	
	result.assign(4, 0.0);
	for (uint xx = 1; xx <= Nx; xx++)
	{
		for (uint yy = 1; yy <= Nx; yy++)
		{
			cplx CBCB = A((xx-1)*Nn+1, (yy-1)*Nn+1);
			cplx CBVB = A((xx-1)*Nn+1, (yy-1)*Nn+2);
			cplx VBCB = A((xx-1)*Nn+2, (yy-1)*Nn+1);
			cplx VBVB = A((xx-1)*Nn+2, (yy-1)*Nn+2);
			result[0] += std::abs(CBCB*CBCB);
			result[1] += std::abs(CBVB*CBVB);
			result[2] += std::abs(VBCB*VBCB);
			result[3] += std::abs(VBVB*VBVB);
		}
	}
	result[0] = negf_math::sqrt(result[0]);
	result[1] = negf_math::sqrt(result[1]);
	result[2] = negf_math::sqrt(result[2]);
	result[3] = negf_math::sqrt(result[3]);
);}


