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
#include "SEAcousticPhonon.h"
using namespace negf;


SEAcousticPhonon::SEAcousticPhonon(const Overlap * ov_,
					const Geometry * xspace_, 
					const Kspace * kspace_, 
					const Energies * energies_, 
					const Options * options_,
					const GreenFunctions * gf_,
					const MaterialDatabase * db):
	SelfEnergy(xspace_,kspace_,energies_,options_, constants::odSAC),
	ov(ov_),
	gf(gf_),
	scaling(1.0),
	new_version(true)
{STACK_TRACE(
	NEGF_ASSERT(ov!=NULL && xspace!=NULL && kspace!=NULL && energies!=NULL && options!=NULL && gf!=NULL, "null pointer encountered.");
	
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
	double lattice_constant;
	if (is_wurtzite) {
		logmsg->emit(LOG_INFO,"Taking wurtzite material parameters.");
		lattice_constant  = constants::convert_from_SI(units::length, 1e-10*biggest_mat->get("lattice_constant_c"));
	} else {
		lattice_constant  = constants::convert_from_SI(units::length, 1e-10*biggest_mat->get("lattice_constant"));
	}
	double mass_density   = constants::convert_from_SI(units::mass, 1.0) / constants::convert_from_SI(units::volume, 1.0)
						                                 * biggest_mat->get("mass_density");
	double speed_of_sound = constants::convert_from_SI(units::velocity, biggest_mat->get("longitudinal_acoustic_velocity"));
	double kT             = constants::convert_from_SI(units::energy, constants::SIkb * options->get("temperature"));
	double e_def_pot; 
	double h_def_pot;
	if (is_wurtzite) {
		e_def_pot = constants::convert_from_SI(units::energy, constants::SIec * fabs(biggest_mat->get("electron_acoustic_deformation_potential_D1")));
		h_def_pot = constants::convert_from_SI(units::energy, constants::SIec * fabs(biggest_mat->get("hole_acoustic_deformation_potential")));
		// bullshit values
	} else {
		e_def_pot = constants::convert_from_SI(units::energy, constants::SIec * fabs(biggest_mat->get("electron_acoustic_deformation_potential")));
		h_def_pot = constants::convert_from_SI(units::energy, constants::SIec * fabs(biggest_mat->get("hole_acoustic_deformation_potential")));
	}
	
	// set up prefactor for every band (w/ correct deformation potential)
	vector<uint> cb_bands;
	options->get_conduction_degrees_of_freedom(cb_bands); // starts at 0
	for (uint ii=0; ii<Nn; ii++) {
		bool electron_band = false;
		for (uint jj=0; jj<cb_bands.size(); jj++) {
			if (cb_bands[jj]==ii) {
				electron_band = true;
				break;
			}
		}
		double def_pot = (electron_band) ? e_def_pot : h_def_pot;
		logmsg->emit(LOG_INFO, "Deformation potential band %d: %.2g",ii+1,def_pot);
		
		this->prefactor.push_back(def_pot * def_pot * kT / (2.0 * constants::pi * mass_density * speed_of_sound*speed_of_sound
									* ((new_version) ? 1.0 : lattice_constant)));
	}
	logmsg->emit(LOG_INFO, "Mass density: %.2g[kg m-3]", mass_density / (constants::convert_from_SI(units::mass, 1.0) / constants::convert_from_SI(units::volume, 1.0)));
	logmsg->emit(LOG_INFO, "Speed of sound: %.2g[m/s]",  speed_of_sound / constants::convert_from_SI(units::velocity, 1.0));
	logmsg->emit(LOG_INFO, "Lattice constant: %.3g[A]",  lattice_constant / constants::convert_from_SI(units::length, 1.0) * 1e10);
	logmsg->emit(LOG_INFO, "Total AC phonon scattering prefactor (band 0): %.2g", this->prefactor[0]);
	
	// determine number of differences xi-xj
	double Nl_dbl = 1 + Nx*(Nx-1)/2.0;
	NEGF_ASSERT(fabs(ceil(Nl_dbl)-Nl_dbl)<1e-14 && fabs(floor(Nl_dbl)-Nl_dbl)<1e-14, "Nl must be integer");
	this->Nl = uint(Nl_dbl);
	
	// set up delta_ij
	this->delta_ij.resize(Nl,0.0);
	// xx==yy
	this->delta_ij[this->find_delta_index(1,1)] = 1.0 / lattice_constant;
	// yy<xx 
	for (uint xx=1; xx<Nx; xx++) {
		for (uint yy=1; yy<xx; yy++) {
			// find |xi-xj|
			uint gx = xspace->get_global_vertex_index(xx-1);
			uint gy = xspace->get_global_vertex_index(yy-1);
			double dx = xspace->get_distance(xspace->get_vertex(gx),xspace->get_vertex(gy));
			// note: sin(dx)/dx is even in dx, so the sign of dx does not matter
			
			this->delta_ij[this->find_delta_index(xx,yy)] = negf_math::sin(constants::pi*dx / lattice_constant) / (constants::pi*dx);
		}
	}
);}


/* Copied from SEOpticalPhonon::find_F_index (see description there) */
uint SEAcousticPhonon::find_delta_index(uint xx, uint yy) const 			// input is 1-based, result is 0-based
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


void SEAcousticPhonon::calculate()
{STACK_TRACE(
	logmsg->emit_small_header("calculating acoustic phonon self-energy");
	const OVMat & M = ov->get_internal_overlap();
	NEGF_ASSERT(M.num_cols()==NxNn && M.num_rows()==M.num_cols(), "wrong overlap matrix.");
	
	for (uint ee2 = 0; ee2 < myNE; ee2++)
	{
		uint ee = energies->get_global_index(ee2);
		
		// ------------------------------------------------------------
		// compute int k*dk M*G(E,k)*M, G=GL,GG,GL
		// ------------------------------------------------------------
		GLMat GL_Integral = GLMat_create(NxNn);
		GLMat GG_Integral = GLMat_create(NxNn);
		Matc  GR_Integral(NxNn, NxNn);
		Matc  tmp1(NxNn, NxNn);
		Matc  tmp2(NxNn, NxNn);
		for (uint kk = 0; kk < Nk; kk++)
		{
			double wk = kspace->get_point(kk).get_weight() / (constants::pi*2.0); // get_weight() includes 2*pi!!!
			
			const GLMat & MGLM = gf->get_overlap_augmented_lesser(kk,ee);
			add(MGLM, wk, GL_Integral); // GL_Integral += wk * MGLM;
			
			const GLMat & MGGM = gf->get_overlap_augmented_greater(kk,ee);
			add(MGGM, wk, GG_Integral); // GG_Integral += wk * MGGM;
			
			const Matc & GR = gf->get_retarded(kk,ee);
			add (GR, wk, tmp1); // tmp1 += wk * GR;
		}
		mult(tmp1, M, tmp2); // tmp2 = tmp1 * M;
		mult(M, tmp2, GR_Integral); // GR_Integral = M * tmp2;
		
		// ------------------------------------------------------------
		// assign integral to all self-energies (momentum-independent)
		// ... augmented with band-dependent prefactor
		// ------------------------------------------------------------
		/*for (uint nn=0; nn<Nn; nn++) {
			if (options->is_conduction_band(nn)) {
				NEGF_ASSERT(!options->is_valence_band(nn));
			} else {
				NEGF_ASSERT(options->is_valence_band(nn));
			}
		}*/
		for (uint kk = 0; kk < Nk; kk++)
		{
			SEMat & SL = this->get_lesser(kk,ee);
			SEMat & SG = this->get_greater(kk,ee);
			SEMat & SR = this->get_retarded(kk,ee);
			
			// self-energies are diagonal
			SL = SACMat_create(NxNn);
			SG = SACMat_create(NxNn);
			SR = SACMat_create(NxNn);
			
			if (!new_version)
			{
				for (uint xx=1; xx<=Nx; xx++) {
					for (uint nn=1; nn<=Nn; nn++) {
						//uint idx = (xx-1)*Nn+nn; 
						uint idx = get_mat_idx(xx,nn,Nx);
						SL(idx,idx) = (this->scaling * this->prefactor[nn-1]) * GL_Integral(idx,idx);
						SG(idx,idx) = (this->scaling * this->prefactor[nn-1]) * GG_Integral(idx,idx);
						SR(idx,idx) = (this->scaling * this->prefactor[nn-1]) * GR_Integral(idx,idx);
					}
				}
			} else 
			{
				for (uint xx=1; xx<=Nx; xx++) {
					for (uint yy=1; yy<=Nx; yy++) {
						uint ll = this->find_delta_index(xx,yy); // index in delta_ij-array corresponding to (i,j)
						for (uint mm=1; mm<=Nn; mm++) {
							for (uint nn=1; nn<=Nn; nn++) {
								if (   (options->is_conduction_band(mm-1) && options->is_valence_band(nn-1))
									|| (options->is_conduction_band(nn-1) && options->is_valence_band(mm-1)) ) {
									// no coupling in this case
									continue;
								}
								int ii = get_mat_idx(xx,mm,Nx); // matrix entry corresponding to (x,m)
								int jj = get_mat_idx(yy,nn,Nx); // matrix entry corresponding to (y,n)
	#ifdef USE_BANDED
								if (fabs(ii-jj) > SL.num_offdiags+1e-8) continue;
	#endif
								SL(ii,jj) = (this->scaling * this->delta_ij[ll] * this->prefactor[nn-1]) * GL_Integral(ii,jj);
								SG(ii,jj) = (this->scaling * this->delta_ij[ll] * this->prefactor[nn-1]) * GG_Integral(ii,jj);
								SR(ii,jj) = (this->scaling * this->delta_ij[ll] * this->prefactor[nn-1]) * GR_Integral(ii,jj);
							}
						}
					}
				}
			}
		}
	}
	mpi->synchronize_processes();
);}


void SEAcousticPhonon::set_scaling(double new_scaling)
{STACK_TRACE(
	NEGF_ASSERT(new_scaling>=0.0 && new_scaling<=1.0, "."); 
	logmsg->emit(LOG_INFO,"Setting scaling of acoustic phonon scattering to %.4g", new_scaling);
	this->scaling = new_scaling;
);}

