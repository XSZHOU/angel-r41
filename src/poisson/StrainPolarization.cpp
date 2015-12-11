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
#include "StrainPolarization.h"
using namespace negf;


StrainPolarization::StrainPolarization(const Geometry * grid_, 
										const MaterialDatabase * db_,
										const double pol_decreaser):
	grid(grid_),
	db(db_)
{STACK_TRACE(
	NEGF_ASSERT(grid!=0 && db!=0, "encountered null pointer.");
	NEGF_ASSERT(grid->get_dimension()==1, "the whole thing works only for 1D grids.");
	logmsg->emit_small_header("Setting up strain and polarization");
	if (pol_decreaser!=1.0) {
		logmsg->emit(LOG_INFO,"Using a decreasing factor of %g", pol_decreaser);
	}
	
	// initialize arrays
	this->exx.resize(grid->get_num_regions(), 0.0);
	this->ezz.resize(grid->get_num_regions(), 0.0);
	this->pol.resize(grid->get_num_regions(), 0.0);		
	this->piezo.resize(grid->get_num_regions(), 0.0);		
	this->sheet_charge.resize(grid->get_num_vertices(), 0.0);
	this->is_wurtzite.resize(grid->get_num_regions(), false);
	this->alattice.resize(grid->get_num_regions(), 0.0);
	this->C11.resize(grid->get_num_regions(), 0.0);
	this->C12.resize(grid->get_num_regions(), 0.0);
	this->C13.resize(grid->get_num_regions(), 0.0);
	this->C33.resize(grid->get_num_regions(), 0.0);
	this->e31.resize(grid->get_num_regions(), 0.0);
	this->e33.resize(grid->get_num_regions(), 0.0);
	this->psp.resize(grid->get_num_regions(), 0.0);
	
	// ------------------------------------------------------------------------
	// find unstrained material name and lattice constants. 
	// unstrained material is material next to a contact w/ largest extension
	// keep in mind: contact regions are extended
	// ------------------------------------------------------------------------
	string unstrained_material = "";
	double extension = 0.0;
	double ai = 0.0;
	for (uint ii=0; ii<grid->get_num_contacts(); ii++) {
		// expecting exactly 1 vertex adjacent to a different region
		Vertex * v = 0;
		for (uint jj=0; jj<grid->get_contact(ii)->get_num_contact_vertices(); jj++) {
			if (grid->get_regions_near(grid->get_contact(ii)->get_contact_vertex(jj)).size() > 1) {
				NEGF_ASSERT(v==0, "there seemed to be more than 1 contact vertex adjacent to the device");
				v = grid->get_contact(ii)->get_contact_vertex(jj);
			}
		}
		
		// expecting exactly 2 regions near v
		const vector<Region *> & regs_near_v = grid->get_regions_near(v);
		NEGF_ASSERT(regs_near_v.size()==2, "expected exactly 2 regions near v.");
		
		// get material of region other than contact region (<-- identifier "contact_")
		int regidx = -1;
		if (regs_near_v[0]->get_name().substr(0,8)=="contact_") {
			regidx = 1;
		} else {
			NEGF_ASSERT(regs_near_v[1]->get_name().substr(0,8)=="contact_", "expected contact_ identifier.");
			regidx = 0;
		}
		const PropertyContainer<double> * material_near_contact = regs_near_v[regidx]->get_material();
		bool wurtzite = (material_near_contact->get_name().find("GaN")!=string::npos);
		double xmin = 1e100;
		double xmax = -1e100;
		for (uint jj=0; jj<grid->get_num_vertices(); jj++) {
			bool v_near_material_near_contact = false;
			const vector<Region *> & regs_near = grid->get_regions_near(grid->get_vertex(ii));
			for (uint kk=0; kk<regs_near.size(); kk++) {
				if (regs_near[kk]==regs_near_v[regidx]) {
					v_near_material_near_contact = true;
					break;
				}
			}
			if (v_near_material_near_contact) {
				if (grid->get_vertex(ii)->get_coordinate(0) > xmax) {
					xmax = grid->get_vertex(ii)->get_coordinate(0);
				}
				if (grid->get_vertex(ii)->get_coordinate(0) < xmin) {
					xmin = grid->get_vertex(ii)->get_coordinate(0);
				}
			}
		}
		double new_ext = xmax - xmin;
		
		// expect the materials near all contacts to be the same
		if (unstrained_material=="") {
			unstrained_material = material_near_contact->get_name();
			ai = constants::convert_from_SI(units::length, 1e-10 * 
					((wurtzite) ? material_near_contact->get("lattice_constant_c") : material_near_contact->get("lattice_constant")) );
			extension = new_ext;
		} else {
			//NEGF_ASSERT(unstrained_material == material_near_contact->get_name(), "expected same material near all contacts.");
			if (unstrained_material!=material_near_contact->get_name()) {
				// check size of material_near_contact
				if (new_ext > extension) {
					unstrained_material = material_near_contact->get_name();
					ai = constants::convert_from_SI(units::length, 1e-10 * 
							((wurtzite) ? material_near_contact->get("lattice_constant_c") : material_near_contact->get("lattice_constant")) );
					extension = new_ext;
				}
			}
		}
	}
	logmsg->emit(LOG_INFO,"Unstrained material: %s", unstrained_material.c_str());
	
	// ---------------------------------------------------------------------
	// find lattice constants, C13, C12, e31, e33, psp from material DB
	// ---------------------------------------------------------------------
	const double ec = constants::convert_from_SI(units::charge, constants::SIec);
	for (uint ii=0; ii<grid->get_num_regions(); ii++) {
		const PropertyContainer<double> * material = grid->get_region(ii)->get_material();
		this->is_wurtzite[ii] = (material->get_name().find("GaN")!=string::npos);
		if (is_wurtzite[ii]) {
			this->alattice[ii] = constants::convert_from_SI(units::length, 1e-10 * material->get("lattice_constant_c") );
			this->C11[ii] = 0.0;
			this->C12[ii] = 0.0;
			this->C13[ii] = constants::convert_from_SI(units::pressure, 		   material->get("elastic_constant_C13"));
			this->C33[ii] = constants::convert_from_SI(units::pressure, 		   material->get("elastic_constant_C33"));
			this->e31[ii] = constants::convert_from_SI(units::density_2d,   	   material->get("piezo_coefficient_e31")) * ec/constants::SIec;
			this->e33[ii] = constants::convert_from_SI(units::density_2d,   	   material->get("piezo_coefficient_e33")) * ec/constants::SIec;
			this->psp[ii] = constants::convert_from_SI(units::density_2d, 		   material->get("spontaneous_polarization")) * ec/constants::SIec;
		} else {
			this->alattice[ii] = constants::convert_from_SI(units::length, 1e-10 * material->get("lattice_constant") );
			this->C11[ii] = constants::convert_from_SI(units::pressure, 		   material->get("elastic_constant_C11"));
			this->C12[ii] = constants::convert_from_SI(units::pressure, 		   material->get("elastic_constant_C12"));
			this->C13[ii] = 0.0;
			this->C33[ii] = 0.0;
			this->e31[ii] = 0.0;
			this->e33[ii] = 0.0;
			this->psp[ii] = 0.0;
		}
	}
	
	// ----------------------------------------------------------------------------------
	// calculate exx=eyy, ezz, piezo-polarization, total polarization (all regionwise)
	// multiply total polarization with pol_decreaser, scaling the result
	// ----------------------------------------------------------------------------------
	for (uint ii=0; ii<grid->get_num_regions(); ii++) {
		this->exx[ii]   = (alattice[ii] - ai) / ai;
		logmsg->emit(LOG_INFO_L2, "region %d: ai=%g, alattice=%g, exx=%g", ii, ai, alattice[ii], exx[ii]);
		if (is_wurtzite[ii]) {
			this->ezz[ii]   = -2.0 * C13[ii]/C33[ii] * exx[ii];
			this->piezo[ii] = e31[ii] * (exx[ii]+ezz[ii]) + e33[ii] * ezz[ii];
			this->pol[ii]   = pol_decreaser * (piezo[ii] + psp[ii]);
		} else {
			this->ezz[ii]   = -2.0 * C12[ii]/C11[ii] * exx[ii];
			this->piezo[ii] = 0.0;
			this->pol[ii]   = 0.0;			
		}
	}
	
	// ----------------------------------------------------------------------------------
	// calculate sheet charge at each vertex
	// ----------------------------------------------------------------------------------
	for (uint ii=0; ii<grid->get_num_vertices(); ii++) {
		Vertex * v = grid->get_vertex(ii);
		const vector<Region *> & regs_near_v = grid->get_regions_near(grid->get_vertex(ii));
		
		if (regs_near_v.size()==1) {
			this->sheet_charge[ii] = 0.0;
			continue;
		}
		NEGF_ASSERT(regs_near_v.size()==2, "expected exactly 2 adjacent regions.");
		
		// which region is left of vertex
		const vector<Edge *> & edges_near_v = grid->get_edges_near(v);
		NEGF_ASSERT(edges_near_v.size()==2, "expected 2 edges near v.");
		Vertex * v0 = (edges_near_v[0]->get_lower_vertex()==v) ? edges_near_v[0]->get_upper_vertex() : edges_near_v[0]->get_lower_vertex();
		Vertex * v1 = (edges_near_v[1]->get_lower_vertex()==v) ? edges_near_v[1]->get_upper_vertex() : edges_near_v[1]->get_lower_vertex();
		NEGF_ASSERT(grid->get_regions_near(v0).size()==1 && grid->get_regions_near(v1).size()==1, "expected neighbouring vertices to be inside regions");
		double  x = v->get_coordinate(0);
		double x0 = v0->get_coordinate(0);
		double x1 = v1->get_coordinate(0);
		NEGF_ASSERT((x0>x && x1<x) || (x0<x && x1>x), "something went wrong.");
		Region * reg0 = grid->get_regions_near(v0)[0];
		Region * reg1 = grid->get_regions_near(v1)[0];
		NEGF_ASSERT(reg0!=reg1, "expected different regions.");
		
		// calculate sheed density w/ correct sign
		if (x0>x) {
			this->sheet_charge[ii] = + (pol[reg0->get_index()] - pol[reg1->get_index()]);
		} else {
			this->sheet_charge[ii] = + (pol[reg1->get_index()] - pol[reg0->get_index()]);
		}
	}
	
	// ----------------------------------------------------------------------------------
	// screen output
	// ----------------------------------------------------------------------------------
	// find leftmost vertex
	Vertex * v = 0;
	double x = 1e100;
	for (uint ii=0; ii < grid->get_num_vertices(); ii++) {
		if (grid->get_vertex(ii)->get_coordinate(0) < x) {
			x = grid->get_vertex(ii)->get_coordinate(0);
			v = grid->get_vertex(ii);
		}
	}
	NEGF_ASSERT(v!=0, "leftmost vertex could not be found.");
	logmsg->emit(LOG_INFO,"region                   exx=eyy       ezz     piezo       psp   tot_pol");
	logmsg->emit(LOG_INFO,"------------------------------------------------------------------------");
	// go through vertices by their coordinate
	Region * actual_region = 0;
	while (true) {
		const vector<Edge *> & edges_near_v = grid->get_edges_near(v);
		bool no_right_vertex = true;
		for (uint ii=0; ii<edges_near_v.size(); ii++) {
			Vertex * v2 = (edges_near_v[ii]->get_lower_vertex()==v) ? edges_near_v[ii]->get_upper_vertex() : edges_near_v[ii]->get_lower_vertex();
			if (v2->get_coordinate(0) > x) {
				x = v2->get_coordinate(0);
				v = v2;
				no_right_vertex = false;
				break;
			}
		}
		if (no_right_vertex) {
			logmsg->emit(LOG_INFO_L3,"vertex loop finished.");
			break;
		}
		const vector<Region *> & regs_near_v = grid->get_regions_near(v);
		if (regs_near_v.size()==1) {
			if (regs_near_v[0] != actual_region) {
				actual_region = regs_near_v[0];
				uint idx = actual_region->get_index();
				logmsg->emit(LOG_INFO,"%22s %9.3g %9.3g %9.3g %9.3g %9.3g", actual_region->get_name().c_str(), 
						exx[idx], ezz[idx], piezo[idx], psp[idx], pol[idx]);
			}
		} else {
			NEGF_ASSERT(regs_near_v.size()==2, "more than 2 adjacent regions encountered.");
			NEGF_ASSERT(this->get_sheet_charge_between(regs_near_v[0],regs_near_v[1])==this->get_sheet_charge(v), "something went wrong.");
			logmsg->emit(LOG_INFO,"           sheetcharge                   %7.3g", this->get_sheet_charge(v));
		}
	}
);}


double StrainPolarization::get_sheet_charge_between (const Region * reg1, const Region * reg2) const
{STACK_TRACE(
	// find vertex that has reg1 and reg2 adjacent
	Vertex * thevertex = 0;
	for (uint ii=0; ii<grid->get_num_vertices(); ii++) {
		const vector<Region *> & regs_near_v = grid->get_regions_near(grid->get_vertex(ii));
		if (regs_near_v.size()!=2) continue;
		
		if ((regs_near_v[0]==reg1 || regs_near_v[0]==reg2) && (regs_near_v[1]==reg1 || regs_near_v[1]==reg2)) {
			NEGF_ASSERT(thevertex==0, "found a second vertex at boundary reg1-reg2!");
			thevertex = grid->get_vertex(ii);
			// break;
		}
	}
	
	NEGF_ASSERT(thevertex!=0, "did not find any vertex where both regions are adjacent.");
	return this->sheet_charge[thevertex->get_index_global()];
);}


