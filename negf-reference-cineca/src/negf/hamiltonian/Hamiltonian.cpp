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
#include "Hamiltonian.h"
using namespace negf;

// re-order storage. 
// band-indices (m,n) at position (x,y) (all 1-based) are found in A((m-1)*Nx+x, (n-1)*Nx+y)

Hamiltonian::Hamiltonian(const Geometry * xspace_, const Kspace * kspace_, const Options * options_,
						 const MaterialDatabase * db_, const char * xgridfilename_) throw (Exception *):
	xspace(xspace_),
	kspace(kspace_),
	options(options_),
	db(db_),
	xgridfilename(xgridfilename_)
{STACK_TRACE(
	logmsg->emit_header("setting up Hamiltonian class");
	NEGF_ASSERT(kspace!=NULL && options!=NULL && db!=NULL, "null pointer enctountered.");
	this->Nx   = xspace->get_num_internal_vertices();
	this->NxNn = Nn*Nx;
		
	// --------------------------------------------------------------------------
	// create information desk from which TDKP will obtain material parameters
	// --------------------------------------------------------------------------
	this->infodesk = new TdkpInfoDesk(db, options->get("temperature"));
	
	double kpmethod = options->get("kp_method");
	logmsg->emit(LOG_INFO,"kpmethod=%d. Spin degeneracy: electrons: %d, holes: %d.", int(kpmethod), 
			int(get_spin_degeneracy(kpmethod, quantities::electron_density)), int(get_spin_degeneracy(kpmethod, quantities::hole_density)) ); 
	
	// --------------------------------------------------------------------------
	// create "bulk" band structure at contact 0
	// this will be needed for calculating the fermilevel in multiband models
	// but it is created always
	// --------------------------------------------------------------------------
	logmsg->emit(LOG_INFO,"Setting bulk configuration object.");
	NEGF_ASSERT(xspace->get_dimension()==1, "contact band structure other than bulk not supported.");
#ifndef NOTDKP
	this->configure_interface(this->contact_0_config, true);
#endif	
	
	if (fabs(kpmethod) > 3.0) { // HACK at the moment
#ifndef NOTDKP
		logmsg->emit(LOG_INFO,"Setting up bulk TDKP band structure for contact fermilevel.");
		switch (xspace->get_dimension()) {
		case 1:
			try {
				this->contact_0_bandstructure = tdkp::InterfaceFactory::create_bulk_radial_slc("contact_0", this->contact_0_config, 
							*this->infodesk, xspace->get_contact(0)->get_adjacent_region()->get_material()->get_name().c_str());
			} catch (std::string s) { s.append("\nwhile setting up contact 0 bandstructure."); NEGF_EXCEPTION(s.c_str()); }
			try {
				this->contact_0_bandstructure->calculate();
			} catch (std::string s) { s.append("\nwhile calculating contact 0 bandstructure."); NEGF_EXCEPTION(s.c_str()); }
			break;
		case 2:
		case 3:
		default:
			NEGF_EXCEPTION("Contact bandstructure must be bulk at the moment.");
		}
#else
		NEGF_EXCEPTION("NEGF was compiled with TDKP disabled.");
#endif
	}
	
	// --------------------------------------------------------------------------
	// if kpmethod=0 or 1, use own effective mass implementation
	// for case 0, the dumb orthogonal delta-TB-basis is used
	// ALSO CHECK COMPUTE_SPECTRAL_XDENSITY AND COMPUTE_SPECTRAL_CURRENT WHEN SWITCHING!!!!!
	// --------------------------------------------------------------------------
	if (fabs(kpmethod - 0.0) < 1e-14) {
		logmsg->emit(LOG_INFO,"Using single-band effective mass model and orthogonal (Datta) basis.");
		this->interface = new InterfaceEffMassOrtho(xspace, db, options->get("temperature"));
		vector<double> potential; potential.resize(xspace->get_num_vertices(), 0.0);
		interface->set_potential(potential);
		// we're already done
		return;
	}
	if (fabs(kpmethod - 1.0) < 1e-14) {
		logmsg->emit(LOG_INFO,"Using single-band effective mass model and FEM basis.");
		this->interface = new InterfaceEffMass(xspace, db, options->get("temperature"));
		vector<double> potential; potential.resize(xspace->get_num_vertices(), 0.0);
		interface->set_potential(potential);
		return;
	}
	if (fabs(kpmethod - 3.0) < 1e-14) {
		logmsg->emit(LOG_INFO,"Using 2-band effective mass model and orthogonal (Datta) basis.");
		this->interface = new InterfaceEffMassOrtho2Band(xspace, db, options->get("temperature"));
		logmsg->emit(LOG_INFO,"Done creating interface.");
		vector<double> potential; potential.resize(xspace->get_num_vertices(), 0.0);
		interface->set_potential(potential);
		return;
	}
	
#ifndef NOTDKP
	// -------------------------------------------
	// configure TDKP interface
	// -------------------------------------------
	logmsg->emit(LOG_INFO,"Setting configuration object for TDKP Hamiltonian.");
	this->configure_interface(this->config, false);
	
	// -------------------------------
	// create interface!
	// -------------------------------
	logmsg->emit(LOG_INFO_L1,"Creating TDKP interface...");
	char logf[1000];
	sprintf(logf,"%sthread%d.tdkplog", fnames->get_outfiledirectory().c_str(),mpi->get_rank());
	string logfilename(logf);
	this->xgridfilename.append(".grd");
	try {
	this->interface = tdkp::InterfaceNEGFWell::factory(
							logfilename, 
							this->config, 
							*this->infodesk, 
							this->xgridfilename.c_str());
	} catch (std::string s) { NEGF_EXCEPTION(s.c_str()); }
	
	// -------------------------------
	// set potential to zero
	// -------------------------------
	logmsg->emit(LOG_INFO_L1,"Setting potential to zero...");
	vector<double> potential; potential.resize(xspace->get_num_vertices(), 0.0);
	try {
	interface->set_potential(potential);
	} catch (std::string s) { NEGF_EXCEPTION(s.c_str()); }
		
	mpi->synchronize_processes();
#else
	NEGF_EXCEPTION("NEGF was compiled with TDKP disabled.");
#endif
);}
			

#ifndef NOTDKP
void Hamiltonian::configure_interface(tdkp::InterfaceConfiguration & conf, bool is_contact)
{STACK_TRACE(
	// set kp method
	double kpmethod = options->get("kp_method");
	if (fabs(kpmethod - 8.0) < 1e-14) {
		conf.set_model(tdkp::InterfaceConfiguration::kp8x8); 
	} else if (fabs(kpmethod - 6.0) < 1e-14) {
		conf.set_model(tdkp::InterfaceConfiguration::kp6x6);
	} else if (fabs(kpmethod - 4.0) < 1e-14) {
		conf.set_model(tdkp::InterfaceConfiguration::kp4x4);
	} else if (fabs(kpmethod - 18.0) < 1e-14) {
		conf.set_model(tdkp::InterfaceConfiguration::kp8x8WZ);
	} else if (fabs(kpmethod - 16.0) < 1e-14) {
		conf.set_model(tdkp::InterfaceConfiguration::kp6x6WZ);
	} else if (fabs(kpmethod - 2.0) < 1e-14 || fabs(kpmethod - 3.0) < 1e-14 || fabs(kpmethod - 1.0) < 1e-14 || fabs(kpmethod - 0.0) < 1e-14) {
		// cases 0 and 1 are needed for bulk object for contact 0 fermilevel
		conf.set_model(tdkp::InterfaceConfiguration::EffectiveMass);
	} else {
		NEGF_EXCEPTION("kp method not understood. must be 8x8 (8), 6x6 (6), 4x4 (4), 8x8WZ (18) or 6x6WZ (16).");
	}
	
	// set length unit of real-space grid
	tdkp::InterfaceConfiguration::LengthUnits lengthunit = tdkp::InterfaceConfiguration::micrometers; 
	string lengthinfo = "um";
	conf.set_grid_unit(lengthunit);
	
	// set kmin, kmax in (1/nm!!!), number of of k points
	NEGF_ASSERT(kspace->get_dimension()==2, "at the moment only 2D transversal k-vectors are implemented.");
	const double conv = 1.0;	// assume 1/nm storage!
	conf.set_k_range(kspace->get_kmin()/conv, kspace->get_kmax()/conv);
	if (is_contact) {
		conf.set_number_of_k_points(100);
	} else {
		conf.set_number_of_k_points(kspace->get_number_of_points());
	}
	
	/** set coordinate axes
     * excerpt from tdkp::Interface.h
     * has different meaning for 1D/2D and 3D problems:
     * 	 3D: axis 0 denotes the crystal axis of the geometrical x axis,
     *       axis 1 ...                                         y axis,
     *       axis 3 ...                                         z axis,
     *
     *   2D: axis 0 denotes the cyrstal axis of the quantized direction
     *              corresponding to the x axis,
     *       axis 1 ... quantized y axis,
     *       axis 2 denotes the crystal axis of the transversal,
     *              free direction
     *
     *   1D: axis 0 denotes the crystal axis of the quantized direction
     * 	     axis 1 denotes the crystal axis of the free transversal direction
     *
     *   0D: (bulk material) bandstructure is calculated in radial 
     *       approximation. 
     *       axis 0 (x,y,z) determines the direction of the bandstructure
     * 
     *  it is clear that all axes must be orthogonal.
     */
	conf.set_axis(0, 1.0, 0.0, 0.0);
	conf.set_axis(1, 0.0, 1.0, 0.0); 
	conf.set_axis(2, 0.0, 0.0, 1.0); 
	
	// set number of bands (does this matter at all in our case???)
	conf.set_max_number_of_cb_subbands(20);
	conf.set_max_number_of_vb_subbands(20);
	
	// set amount of screen output (-1 --> nothing, 3-> a little, 6->full)
	conf.set_log_output_level(1);
	
	// irrelevant interface functions:
	//conf.set_stop_calculation_at_maximal_edge(bool stop_calculation_at_max);
	conf.use_strains(false);
	//conf.set_tdkp_input_dir(const char* input_directory);
	//conf.set_tdkp_output_dir(const char* output_directory);
	//conf.dump() const;
);}
#endif


void Hamiltonian::get(const DomainPoint & kpoint, Matc & result) const throw (Exception *)
{STACK_TRACE(
	uint Nvert = xspace->get_num_vertices();
	NEGF_FASSERT(result.num_rows()==Nvert*Nn && result.num_cols()==Nvert*Nn, 
			"size of matrix handed over must be Nvert*Nn=%dx%d! instead %dx%d",
			Nvert,Nn,result.num_rows(),result.num_cols());
	
	double kk = kpoint.get_coord_abs();
	NEGF_ASSERT(!isnan(kk), "|k| was NaN!");
	double kk_nm = 1e-9 * kk/constants::convert_from_SI(units::density_1d, 1.0);
	try {
	interface->assemble_hamiltonian(kk_nm);
	GEMatrix tmp(Nvert*Nn,Nvert*Nn);
	interface->get_hamiltonian(tmp);
	for (uint xx=1; xx<=Nvert; xx++) {
		for (uint yy=1; yy<=Nvert; yy++) {
			for (uint mm=1; mm<=Nn; mm++) {
				for (uint nn=1; nn<=Nn; nn++) {
				    NEGF_ASSERT(   !isnan(tmp((xx-1)*Nn+mm, (yy-1)*Nn+nn).real())
                                && !isnan(tmp((xx-1)*Nn+mm, (yy-1)*Nn+nn).imag()), "NaN encountered!");
					//result.flens_matrix((mm-1)*Nvert+xx, (nn-1)*Nvert+yy) = tmp((xx-1)*Nn+mm, (yy-1)*Nn+nn);
					result.flens_matrix(get_mat_idx(xx,mm,Nvert), get_mat_idx(yy,nn,Nvert)) = tmp((xx-1)*Nn+mm, (yy-1)*Nn+nn);
				}
			}
		}
	}
	} catch (std::string s) { NEGF_EXCEPTION(s.c_str()); }
	NEGF_ASSERT(result.num_rows()==Nvert*Nn && result.num_cols()==Nvert*Nn, "inconsistent matrix dimension.");
	
	/* UNIT CONVERSION:
	 * tdkp gives H_ij = \int dx phi_i(x) H(x) \phi_j(x)
	 * - the FEM shapefunctions phi are unitless (they range between 0 and 1)
	 * - H(x) is an energy
	 * - therefore H_ij has units energy*length
	 * - tdkp calulates in eV and nm
	 */
	NEGF_ASSERT(xspace->get_dimension()==1, "need different conversion factor for xspace=2D,3D");
	double conv = /* SI --> NEGF */ constants::convert_from_SI(units::energy, 1.0) * constants::convert_from_SI(units::length, 1.0)
					  * (/* tdkp --> SI */ constants::SIec * 1e-9);
	if (constants::old_orthogonal
		&& (options->get("kp_method")==0.0 || options->get("kp_method")==3.0)) {
		// orthogonal basis Hamiltonian, units energy
		conv = constants::convert_from_SI(units::energy, 1.0) *(/* tdkp --> SI */ constants::SIec);
	}
	result *= conv;
);}


// can't make it const because of work_ham
void Hamiltonian::get_internal(const DomainPoint & kpoint, Matc & result) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(result.num_rows()==NxNn && result.num_cols()==NxNn, "size of matrix handed over must be Nx*Nn!");
	uint Nvert = xspace->get_num_vertices();
	
	// get full Hamiltonian and store in temporary matrix
	Matc work_ham(Nvert*Nn, Nvert*Nn);
	this->get(kpoint, work_ham);
	
	// construct the Hamiltonian of the interior points only and store in result
	for (uint xx=1; xx<=Nx; xx++) {
		for (uint yy=1; yy<=Nx; yy++) {
			uint gx = xspace->get_global_vertex_index(xx-1) + 1;
			uint gy = xspace->get_global_vertex_index(yy-1) + 1;
//#ifdef REORDER
			for (uint mm=1; mm<=Nn; mm++) {
				for (uint nn=1; nn<=Nn; nn++) {
					result(get_mat_idx(xx,mm,Nx),get_mat_idx(yy,nn,Nx)) = work_ham(get_mat_idx(gx,mm,Nvert), get_mat_idx(gy,nn,Nvert));
				}
			}
//#else
//			result.fill_block(xx, yy, H, gx, gy);
//#endif
		}
	}
);}


void Hamiltonian::get_overlap(Matd & result) const throw (Exception *)
{STACK_TRACE(
	uint Nvert = xspace->get_num_vertices(); // NOT only internal vertices!
	
	result = Matd(Nvert,Nvert);
	try {
	interface->assemble_hamiltonian(0); // k-vector should be irrelevant
	// REORDER flag is irrelevant here because matrix is only Nvert*Nvert, not "bands-resolved"
	interface->get_overlap(result.flens_matrix);
	} catch (std::string s) { NEGF_EXCEPTION(s.c_str()); }
	NEGF_ASSERT(result.num_rows()==Nvert && result.num_cols()==Nvert, "inconsistent matrix size.");
	
	// UNIT CONVERSION: \int dx phi_i(x) \phi_j(x) has units length. tdkp calulates in nm
	NEGF_ASSERT(xspace->get_dimension()==1, "need different conversion factor for xspace=2D,3D");
	double conv = /* SI --> NEGF */ constants::convert_from_SI(units::length, 1.0) * /* tdkp --> SI */ 1e-9;
	if (constants::old_orthogonal
		&& (options->get("kp_method")==0.0 || options->get("kp_method")==3.0)) {
		// orthogonal basis Hamiltonian, units 1
		conv = 1.0;
	}
	result *= conv;
);}

/** assign new electrostatic potential to TDKP / Interface object 
 *  @param elstat_potential the electrostatic potential defined on all Poisson vertices */
void Hamiltonian::set_electrostatic_potential(const vector<double> & elstat_potential) throw (Exception *)
{STACK_TRACE(
	this->electrostatic_potential = elstat_potential;
	
	// multiply electrostatic potential with -ec to get energetic potential
	vector<double> pot = electrostatic_potential;
	const double ec = constants::convert_from_SI(units::charge, constants::SIec);
	for (uint ii = 0; ii < pot.size(); ii++) {
		pot[ii] = pot[ii] * (-ec);
		NEGF_ASSERT(!isnan(pot[ii]) && !isinf(pot[ii]), "electrostatic potential was NaN!");
	}
	NEGF_FASSERT(pot.size()==xspace->get_num_vertices(), "inconsistent potential size (%d instead of %d)",
				pot.size(),xspace->get_num_vertices());
	// assign
	try {
	interface->set_potential(pot);
	} catch (std::string s) { NEGF_EXCEPTION(s.c_str()); }
);}


void Hamiltonian::set_strain(StrainPolarization * strainpol)
{STACK_TRACE(
     interface->set_strain(strainpol);
);}


