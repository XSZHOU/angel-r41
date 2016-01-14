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
#include "TdkpInfoDesk.h"
using namespace negf;

namespace negf {

#ifndef NOTDKP
double get_spin_degeneracy(tdkp::Interface * kpresult, quantities::PhysicalQuantity e_or_h) 
{STACK_TRACE(
	if (e_or_h==quantities::electron_density && 
			(   kpresult->get_config().model!=tdkp::InterfaceConfiguration::kp8x8 
	 		 && kpresult->get_config().model!=tdkp::InterfaceConfiguration::kp8x8WZ)) {
		return 2.0;
	}
	if (e_or_h==quantities::hole_density && kpresult->get_config().model==tdkp::InterfaceConfiguration::EffectiveMass) {
		return 2.0;
	}
	return 1.0;
);}
#endif

double get_spin_degeneracy(const double & kpmethod, quantities::PhysicalQuantity e_or_h /*true-->holes*/) 
{STACK_TRACE(
	if (fabs(kpmethod - 0.0) < 1e-14) {	// single band effective mass, orthogonal (Datta) basis
		return 2.0;
	}
	if (fabs(kpmethod - 1.0) < 1e-14) {	// single band effective mass, FEM basis
		return 2.0;
	}
	if (fabs(kpmethod - 2.0) < 1e-14) {	// CB&VB band effective mass, FEM basis
		return 2.0;
	}
	if (fabs(kpmethod - 3.0) < 1e-14) {	// CB&VB band effective mass, orthogonal basis
		return 2.0;
	}
	if (e_or_h==quantities::hole_density) {	// kp multiband models, holes --> always nondegenerate
		return 1.0;
	}
	NEGF_ASSERT(e_or_h==quantities::electron_density, "electron or hole density expected.");
	if (fabs(kpmethod - 8.0) < 1e-14 || fabs(kpmethod - 18.0) < 1e-14) { // kp 8x8, electrons
		return 1.0;
	}
	if (fabs(kpmethod - 6.0) < 1e-14 || fabs(kpmethod - 16.0) < 1e-14 
		|| fabs(kpmethod - 4.0) < 1e-14) { 	// kp 4x4,6x6, electrons
		return 2.0;
	}
	NEGF_EXCEPTION("Unknown kpmethod.");
);}

}


double TdkpInfoDesk::get_property(const string & materialname, const string & propertyname) const
{INTERFACE_STACK_TRACE( 
	PropertyContainer<double> * mat = db->get_material(materialname.c_str());
	logmsg->emit(LOG_INFO_L3,"*** TDKP requests %s",propertyname.c_str());
	const double eV = constants::convert_from_SI(units::energy, constants::SIec); // 1 eV in NEGF units
	
	// -----------------------------------------
	// some hacks for R. Veprek's TDKP solver  
	// -----------------------------------------
	if (propertyname.compare("conduction_band_edge")==0) {
		return this->convert_to_tdkp_units(units::energy, 
				this->get_cbedge(mat, this->temperature, this->db) / eV);
	}
	if (propertyname.compare("bandgap")==0) {
		return this->convert_to_tdkp_units(units::energy, 
				(this->get_cbedge(mat, this->temperature, this->db) - eV*mat->get("valence_band_edge"))  /eV);
	}
	// tdkp requests "transverse masses", we assume spherical band structures (also in wurtzite case...)
	if (propertyname.compare("electron_effective_mass_transverse")==0) {
		return mat->get("electron_effective_mass");
	}
	if (propertyname.compare("hole_effective_mass_transverse")==0) {
		return mat->get("hole_effective_mass");
	}
	
	// -----------------------------
	// normal cases
	// -----------------------------
	if (   propertyname.compare("valence_band_edge")==0			// is already stored in eV
	    || propertyname.substr(0,19)=="luttinger_parameter"		// unitless
	    || propertyname.substr(0,16)=="strain_potential"		// is already stored in eV
	    || propertyname.compare("spin_orbit_splitting")==0		// is already stored in eV
	    || propertyname.compare("optical_matrix_element")==0 	// is already stored in eV
	    || propertyname.compare("electron_effective_mass")==0	// unitless
		|| propertyname.compare("hole_effective_mass")==0){		// unitless
		return mat->get(propertyname.c_str());
	}
	if (propertyname.substr(0,16)=="elastic_constant") {
		return this->convert_to_tdkp_units(units::pressure, mat->get(propertyname.c_str()) );
	}
	NEGF_FEXCEPTION("Could not convert unit for parameter %s", propertyname.c_str());
	return mat->get(propertyname.c_str()); 
);}


double TdkpInfoDesk::get_slc_param(const std::string & propertyname) const
{INTERFACE_STACK_TRACE(
	if (propertyname=="homogeneous_broadening") return 1.0e13;  // in units of angular frequency 1/s; 3e13 corresponds to 20meV
	if (propertyname=="refractive_index") return 3.5;
	if (propertyname=="geometry_output_division_factor") return 1.0; // WE DO NOT USE THIS FEATURE!!! --> SLC solution will lack this factor
	// note that this factor needs to be in tdkp units and therefore is equivalent to 1 nm^d where d is the dimensionality of the problem
	if (propertyname=="omega_num") return constants::slc_num_points;
	NEGF_FEXCEPTION("TdkpInfoDesk: Don't know how to determine property \"%s\".",propertyname.c_str());
);}


bool TdkpInfoDesk::is_set(const string & materialname, const string & propertyname) const
{INTERFACE_STACK_TRACE( 
	PropertyContainer<double> * mat = db->get_material(materialname.c_str());
	// some hacks: R. Veprek's TDKP solver requests conduction_band_edge, valence_band_edge
	if (propertyname.compare("conduction_band_edge")==0 || propertyname.compare("bandgap")==0) {
		return mat->is_set("valence_band_edge") && mat->is_set("bandgap_0K")
				 && mat->is_set("bandgap_alpha") && mat->is_set("bandgap_beta");
	}
	// this is the normal case
	return mat->is_set(propertyname); 
);}


/** tdkp expects nm - eV - s 
 *  we assume that "value" is in SI units, except for energies where it shall be in eV */
double TdkpInfoDesk::convert_to_tdkp_units(units::UnitType unit, const double & value) const
{INTERFACE_STACK_TRACE( 
	switch (unit) {
	case units::energy: 	return value; 			break;	// assume storage in eV!!!
	case units::length: 	return value * 1e-9; 	break;
	case units::time:		return value;			break;
	case units::pressure: 	return value;			break;	// only relative sizes of elastic constants needed
	default: NEGF_EXCEPTION("this unit is not implemented."); return 0.0; break;
	}
);}

/** workaround to use bowing parameter instead of interpolation of alpha, beta etc. for ternary materials */
double TdkpInfoDesk::get_cbedge(const PropertyContainer<double> * mat, const double & T/*emperature*/, 
							const MaterialDatabase * const database) /*const*/
{STACK_TRACE(
	int ternary = -1;
	const string & matname = mat->get_name();
	for (uint ii = 0; ii < Constants.ternary_names.size(); ii++) {
		// attention: material name was appended w/ molefraction at this stage
		if (matname.length()>=Constants.ternary_names[ii].length() && 
			matname.substr(0,Constants.ternary_names[ii].length()).compare(Constants.ternary_names[ii]) == 0) {
			ternary = ii;
			break;
		}
	}
	
	 // valence band edge is correct for both binary and ternary materials
	double vbedge = mat->get("valence_band_edge");
	
	// --------------------------
	// case of ternary material
	// --------------------------
	if (ternary!=-1) {
		const TernaryPropertyContainer<double> * matt = database->get_ternary_material(mat->get_name().c_str());
		double molefraction = matt->get_molefraction();
		
		const PropertyContainer<double> * mat0 = matt->get_first_pure_material();
		const PropertyContainer<double> * mat1 = matt->get_second_pure_material();
		NEGF_FASSERT(mat0->is_set("bandgap_0K"), 	"property \"bandgap_0K\" not found in material %s.", 	mat0->get_name().c_str());
		NEGF_FASSERT(mat0->is_set("bandgap_alpha"), "property \"bandgap_alpha\" not found in material %s.", mat0->get_name().c_str());
		NEGF_FASSERT(mat0->is_set("bandgap_beta"), 	"property \"bandgap_beta\" not found in material %s.", 	mat0->get_name().c_str());
		NEGF_FASSERT(mat1->is_set("bandgap_0K"), 	"property \"bandgap_0K\" not found in material %s.", 	mat1->get_name().c_str());
		NEGF_FASSERT(mat1->is_set("bandgap_alpha"), "property \"bandgap_alpha\" not found in material %s.", mat1->get_name().c_str());
		NEGF_FASSERT(mat1->is_set("bandgap_beta"), 	"property \"bandgap_beta\" not found in material %s.", 	mat1->get_name().c_str());
		
		double bandgap_0K_mat0 	= mat0->get("bandgap_0K");
		double alpha_mat0 		= mat0->get("bandgap_alpha");
		double beta_mat0 		= mat0->get("bandgap_beta");
		double bandgap_0K_mat1  = mat1->get("bandgap_0K");
		double alpha_mat1 		= mat1->get("bandgap_alpha");
		double beta_mat1 		= mat1->get("bandgap_beta");
		
		double bandgap_bowing = 0.0;
		if (matt->is_set("bandgap_bowing")) {
			bandgap_bowing = matt->get("bandgap_bowing");
			// please note that any molefraction dependence of the bowing parameter itself is treated in TernaryPropertyContainer.h
		}
		
		double bandgap_TK_mat0 =bandgap_0K_mat0 - alpha_mat0*T*T / (beta_mat0+T);
		double bandgap_TK_mat1 =bandgap_0K_mat1 - alpha_mat1*T*T / (beta_mat1+T);
		
		double bandgap_TK_alloy = (1.0-molefraction)*bandgap_TK_mat0 + molefraction*bandgap_TK_mat1
								 - molefraction*(1.0-molefraction) * bandgap_bowing;
								 
		return vbedge + bandgap_TK_alloy;
	} else {
	// --------------
	// normal case
	// --------------
		NEGF_FASSERT(mat->is_set("bandgap_0K"), "property \"bandgap_0K\" not found in material %s.", mat->get_name().c_str());
		NEGF_FASSERT(mat->is_set("bandgap_alpha"), "property \"bandgap_alpha\" not found in material %s.", mat->get_name().c_str());
		NEGF_FASSERT(mat->is_set("bandgap_beta"), "property \"bandgap_beta\" not found in material %s.", mat->get_name().c_str());
		
		double bandgap = mat->get("bandgap_0K");
		double alpha   = mat->get("bandgap_alpha");
		double beta    = mat->get("bandgap_beta");
		// assume storage in Kelvin, Tpar=300K, bandgap = Egap @ 0K !!!
		bandgap += /*alpha*300*300 / (beta+300)*/ - alpha*T*T / (beta+T);
		
		return vbedge + bandgap;
	}
);}
