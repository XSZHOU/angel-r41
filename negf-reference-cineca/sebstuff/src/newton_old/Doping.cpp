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
#include "Doping.h"
using namespace negf;

#ifndef NODFISE

Doping::Doping( const  Geometry * grid_,
				const  BoxMethod * bm_,
				const  string filename, 
				uint   options_DopingUnit, 		// 0 --> .dat-file is in cm-3,    1 --> .dat file is in m-3
				int    options_DopingSign,		// -1 --> multiply .dat file values by -1
				bool   options_UseGrainBoundaryDoping,	// 0 --> e.g. fields like ArsenicGrainBoundary are not included
				uint   options_dim):
	grid(grid_),
	bm(bm_)
{STACK_TRACE(
	NEGF_ASSERT(grid!=NULL && bm!=NULL, "encountered null pointer.");
	NEGF_ASSERT(options_DopingUnit==0 || options_DopingUnit==1, "DopingUnit option must be 0 or 1.");
	NEGF_ASSERT(options_DopingSign==1 || options_DopingSign==-1, "DopingSign option must be 1 or -1.");
	
	dependencies.clear();				// no dependencies
	this->its_type = quantities::hole_density;
	this->number_of_variables = grid->get_num_vertices();
	this->timestamp = 99999;			// makes sure update is never needed
	current_variable_values.clear();
	current_variable_values.resize(grid->get_num_vertices(), 0.0);
	
	// --------------------------------------------------------------------------
	// list of possible fieldnames
	// note: n-doping means positively ionized dopants, so same sign as holes
	//       n-doping means that the variable value is >0!
	// --------------------------------------------------------------------------
	vector<string> 							doping_fields;		// fieldnames
	vector<quantities::PhysicalQuantity> 	n_or_p;				// n or p doping
	
	doping_fields.push_back("Arsenic");     n_or_p.push_back(quantities::hole_density);
	doping_fields.push_back("Phosphorus");  n_or_p.push_back(quantities::hole_density);
	doping_fields.push_back("Antimony");    n_or_p.push_back(quantities::hole_density);
	doping_fields.push_back("Boron");       n_or_p.push_back(quantities::electron_density);
	doping_fields.push_back("Indium");      n_or_p.push_back(quantities::electron_density);
	if (options_UseGrainBoundaryDoping) {
		uint num_fields = doping_fields.size();
		for (uint ii = 0; ii < num_fields; ii++) {
			string fieldname = doping_fields[ii] + "GrainBoundary";
			doping_fields.push_back(fieldname);
			n_or_p.push_back((n_or_p[ii]==quantities::electron_density)
								? quantities::hole_density
								: quantities::electron_density);	// changes sign according to bipolar3d.dat!!! (don't know why)
//								? quantities::electron_density
//								: quantities::hole_density);
		}
	}
	// doping_fields.push_back("Doping");      n_or_p.push_back(quantities::hole_density);
	
	units::UnitType dopingunit;
	switch(options_dim) {
		case 1: dopingunit = units::density_1d; break;
		case 2: dopingunit = units::density_2d; break;
		case 3: dopingunit = units::density_3d; break;
		default: NEGF_EXCEPTION("dimension must be 1, 2 or 3.");
	}
	
	// ------------------------------------
	// read in fields
	// ------------------------------------
	InputParser parser;
	bool some_doping_was_found = false;
	for (uint ff = 0; ff < doping_fields.size(); ff++)
	{
		// check if field exists in the file; if not, continue to next field
		string fieldname = doping_fields[ff]; fieldname.append("ActiveConcentration");
		if (!parser.check_field_existence(filename.c_str(), fieldname)) {
			fieldname = doping_fields[ff]; fieldname.append("Concentration");
			if (!parser.check_field_existence(filename.c_str(), fieldname)) {
				continue;
			}
		}
		some_doping_was_found = true;
		
		// read in field
		logmsg->emit(LOG_INFO_L2, "reading doping field %s",fieldname.c_str());
		vector<double> tempvec;	 // stores values in DF-ISE numbering  NO!!! The numbering is already the same as in NEGF
		parser.read_dfise_dat(filename.c_str(), filename.c_str(), grid, tempvec, fieldname, "vertex");
		
		// convert units and sign
		vector<double> tempvec2; // stores converted values in the numbering used here
		tempvec2.resize(grid->get_num_vertices(), 0.0);
		int npsign;
		if (n_or_p[ff]==quantities::electron_density) {
			npsign = -1;
		} else if (n_or_p[ff]==quantities::hole_density) {
			npsign = 1;
		} else {
			NEGF_EXCEPTION("n_or_p must be an electron or a hole density.");
		}
		for (uint ii = 0; ii < grid->get_num_vertices(); ii++) {
			//NEGF_FASSERT(tempvec.size() > grid->get_vertex(ii)->get_index_external(), 
			//		"external vertex index (%d) exceeds field size (%d).", grid->get_vertex(ii)->get_index_external(), tempvec.size());
			//tempvec2[ii] = npsign * constants::convert_from_SI(dopingunit, tempvec[grid->get_vertex(ii)->get_index_external()]);
			tempvec2[ii] = npsign * constants::convert_from_SI(dopingunit, tempvec[grid->get_vertex(ii)->get_index_global()]);
			if (options_DopingUnit==0)	// DESSIS units, doping concentration in cm^-3
				tempvec2[ii] = tempvec2[ii] * 1e6;
			if (options_DopingSign==-1)
				tempvec2[ii] = -tempvec2[ii];
		}
		
		// add the field to the total doping field
		for (uint ii = 0; ii <  grid->get_num_vertices(); ii++) {
			this->current_variable_values[ii] += tempvec2[ii];
		}
		
		// correct doping in case oxide is around! (doping=0 in case of oxide...)
		for (uint ii = 0; ii <  grid->get_num_vertices(); ii++) {
			Vertex * vertex = grid->get_vertex(ii);
			const vector<Region*> & reg = grid->get_regions_near(vertex);
			double frac = 0.0;
			for (uint rr = 0; rr < reg.size(); rr++) {
				for (uint nn = 0; nn < Constants.oxide_names.size(); nn++) {
					if (reg[rr]->get_material_name()==Constants.oxide_names[nn]) {
						frac += bm->get_fraction_of_surroundment(vertex,reg[rr]);
						break;
					}
				}
			}
			// do not correct doping if vertex is at a contact and not enirely surrounded
			if (vertex->is_at_contact() && fabs(frac - 1.0) > 1e-5) {
				continue;
			}
			current_variable_values[ii] = (1-frac) * current_variable_values[ii];
		}
	}
	if (!some_doping_was_found) {
		logmsg->emit(LOG_WARN, "*** WARING *** no doping field was found. be sure to specify the individual dopants.");
	}
);}

#endif // NODFISE


Doping::Doping( const  Geometry * grid_,
                const  BoxMethod * bm_,
                const map< string, PropertyContainer<double> * > * cmdfile):
    grid(grid_),
    bm(bm_)
{STACK_TRACE(
    NEGF_ASSERT(grid!=NULL && bm!=NULL, "encountered null pointer.");

    // Equation stuff
    dependencies.clear();               // no dependencies
    this->its_type = quantities::hole_density;
    this->number_of_variables = grid->get_num_vertices();
    this->timestamp = 99999;            // makes sure update is never needed
    current_variable_values.clear();
    current_variable_values.resize(grid->get_num_vertices(), 0.0);

    units::UnitType dopingunit = units::density_3d;

    PropertyContainer<double> * regs = 0;
    for (map< string, PropertyContainer<double> * >::const_iterator it = (*cmdfile).begin(); it!=(*cmdfile).end(); it++) {
        if (it->first=="regions") {
            regs = it->second;
            break;
        }
    }
    NEGF_ASSERT(regs!=0, "regions section was not found.");


    vector<double> tempvec(grid->get_num_vertices(), 0.0);
    for (uint rr = 0; rr < grid->get_num_regions(); rr++) {
        // first and last regions are contacts
        int ridx = rr;
        if (rr==0) ridx = 1; // left contact
        if (rr==grid->get_num_regions()-1) ridx = grid->get_num_regions()-2; // right contact

        char buf[1000]; sprintf(buf, "region%d_doping", ridx-1); // -1 because right contact was inserted as a region
        NEGF_ASSERT(regs->is_set(buf), "could not find doping for a region.");
        double doping = constants::convert_from_SI(dopingunit, 1e6 * regs->get(buf));

        const vector<Element *> & reg_elems = grid->get_region(rr)->get_elements();
        for (uint ee=0; ee<reg_elems.size(); ee++) {
            uint v0 = reg_elems[ee]->get_vertex(0)->get_index_global();
            this->current_variable_values[v0] = doping;
            uint v1 = reg_elems[ee]->get_vertex(1)->get_index_global();
            this->current_variable_values[v1] = doping;
            // TODO:
            // problem with interface vertices... we should weigh them
        }
    }

    logmsg->emit(LOG_INFO_L1, "DOPING:");
    for (uint ii=0; ii<current_variable_values.size(); ii++) {
        logmsg->emit(LOG_INFO_L1, "  x=%5g: dop=%5g", grid->get_vertex(ii)->get_coordinate(0), current_variable_values[ii]);
    }
);}


