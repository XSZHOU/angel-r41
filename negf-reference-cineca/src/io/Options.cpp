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
#include "Options.h"
using namespace negf;

Options::Options() throw (Exception *)
{STACK_TRACE(
	logmsg->emit_header("reading options");
	logmsg->emit(LOG_INFO,"File: %s.cmd",fnames->get_filename().c_str());
	InputParser parser;
    //read_cmd_file, here problems !!!(S.Z. first try)
	this->sim = parser.read_cmd_file(fnames->get_filename());
             
	bool found = false;
    
	for (map<string, PropertyContainer<double> * >::iterator it = sim.begin(); it != sim.end(); ++it)
	{
		string::size_type loc = it->first.find(constants::options_name);
		if (loc==0 && it->first.size()==constants::options_name.length()) {
			NEGF_FASSERT(found==false, "duplicate identifier \"%s\" detected.", constants::options_name.c_str());
			this->opts = it->second;
			found = true;
		} else {
		    // OLD: clean up the rest. NEW: do nothing
			// delete it->second;
			// it->second = NULL;
		}
	}
	NEGF_FASSERT(found, "identifier \"%s\" was not found in .cmd-file.", constants::options_name.c_str());
	
	// -------------------------
	// SET DEGREES OF FREEDOM
	// at the moment this corresponds to the number of bulk kp bands
	// see also Hamiltonian.cpp
	// -------------------------	
	uint dofs = 0;
	double kpmethod = this->get("kp_method");
             std::cout<<"kpmethod Get !"<<kpmethod<<std::endl;
	if (fabs(kpmethod - 8.0) < 1e-14) {			// zincblende 8x8
		dofs = 8; 
	} else if (fabs(kpmethod - 6.0) < 1e-14) { 	// zincblende 6x6
		dofs = 6;
	} else if (fabs(kpmethod - 4.0) < 1e-14) { 	// zincblende 4x4
		dofs = 4;
	} else if (fabs(kpmethod - 18.0) < 1e-14) {	// wurtzite 8x8
		dofs = 8;
	} else if (fabs(kpmethod - 16.0) < 1e-14) {	// wurtzite 6x6
		dofs = 6;
	} else if (fabs(kpmethod - 1.0) < 1e-14) {	// single band effective mass, FEM basis
		dofs = 1;
	} else if (fabs(kpmethod - 0.0) < 1e-14) {	// single band effective mass, orthogonal basis
		dofs = 1;
	} else if (fabs(kpmethod - 2.0) < 1e-14) {	// CB&VB effective mass, FEM
		dofs = 2;
	} else if (fabs(kpmethod - 3.0) < 1e-14) {	// CB&VB effective mass, orthogonal basis
		dofs = 2;
	} else {
		NEGF_EXCEPTION("kp method not understood. must be 8x8 (8), 6x6 (6), 4x4 (4), 8x8WZ (18) or 6x6WZ (16).");
	}
	
	// Nn IS EXCLUSIVELY SET HERE!
	negf::Nn = dofs;
#ifdef FLENS
	work1 = GEMatrix(Nn,Nn);
	work2 = GEMatrix(Nn,Nn);
	work3 = GEMatrix(Nn,Nn);
#endif
	
	this->conduction_bands.clear(); this->get_conduction_degrees_of_freedom(conduction_bands);
	this->valence_bands.clear();    this->get_valence_degrees_of_freedom(valence_bands);
);}
	

bool Options::exists(const string & propertyname) const throw (Exception *)
{STACK_TRACE(
	return opts->is_set(propertyname);
);}	

double Options::get(const string & propertyname) const throw (Exception *)
{STACK_TRACE(
	NEGF_FASSERT(opts->is_set(propertyname), "option \"%s\" not found.",propertyname.c_str());
	return opts->get(propertyname);
);}

void Options::get_valence_degrees_of_freedom(vector<uint> & result) const throw (Exception *)
{STACK_TRACE(
	// see also Hamiltonian.cpp
	result.clear();
	
	double kpmethod = this->get("kp_method");
	if (fabs(kpmethod - 8.0) < 1e-14 
	 || fabs(kpmethod - 18.0) < 1e-14) {		// 8x8
		// S,X,Y,Z,S,X,Y,Z
		result.push_back(1); 
		result.push_back(2); 
		result.push_back(3); 
		result.push_back(5); 
		result.push_back(6); 
		result.push_back(7); 
		return; 
	} else if (fabs(kpmethod - 6.0) < 1e-14 
			|| fabs(kpmethod - 16.0) < 1e-14) {	// 6x6
		NEGF_EXCEPTION("unknown.");
		return;
	} else if (fabs(kpmethod - 4.0) < 1e-14) {	// 4x4
		NEGF_EXCEPTION("unknown.");
		return;
	} else if (fabs(kpmethod - 1.0) < 1e-14 
			|| fabs(kpmethod - 0.0) < 1e-14) {	// 1x1
		return;
	} else if (fabs(kpmethod - 2.0) < 1e-14
			|| fabs(kpmethod - 3.0) < 1e-14) {
		result.push_back(1); 	
		return;
	} else {
		NEGF_EXCEPTION("kp method not understood. must be 8x8 (8), 6x6 (6), 4x4 (4), 8x8WZ (18) or 6x6WZ (16).");
	}
);}


void Options::get_conduction_degrees_of_freedom(vector<uint> & result) const throw (Exception *)
{STACK_TRACE(
	// see also Hamiltonian.cpp
	result.clear();
	
	double kpmethod = this->get("kp_method");
	if (fabs(kpmethod - 8.0) < 1e-14 || fabs(kpmethod - 18.0) < 1e-14) {		// 8x8
		// S,X,Y,Z,S,X,Y,Z
		result.push_back(0); 
		result.push_back(4);
		return; 
	} else if (fabs(kpmethod - 6.0) < 1e-14 || fabs(kpmethod - 16.0) < 1e-14) {	// 6x6
		NEGF_EXCEPTION("unknown.");
		return;
	} else if (fabs(kpmethod - 4.0) < 1e-14) {									// 4x4
		NEGF_EXCEPTION("unknown.");
		return;
	} else if (fabs(kpmethod - 1.0) < 1e-14 
			|| fabs(kpmethod - 0.0) < 1e-14 
			|| fabs(kpmethod - 2.0) < 1e-14
			|| fabs(kpmethod - 3.0) < 1e-14) {	// 1x1
		result.push_back(0); 
		return;
	} else {
		NEGF_EXCEPTION("kp method not understood. must be 8x8 (8), 6x6 (6), 4x4 (4), 8x8WZ (18) or 6x6WZ (16).");
	}
);}


bool Options::is_conduction_band(const uint nn) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(nn<Nn, "invalid index nn (must be 0-based)");
	for (uint ii=0; ii<conduction_bands.size(); ii++) {
		if (conduction_bands[ii]==nn) {
			return true;
		}
	}
	return false;
);}


bool Options::is_valence_band(const uint nn) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(nn<Nn, "invalid index nn (must be 0-based)");
	for (uint ii=0; ii<valence_bands.size(); ii++) {
		if (valence_bands[ii]==nn) {
			return true;
		}
	}
	return false;
);}
