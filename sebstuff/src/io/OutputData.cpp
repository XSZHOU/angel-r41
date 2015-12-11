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
#include "OutputData.h"
using namespace negf;


OutputData::OutputData( const Geometry * const grid_, const string & resultfilename_) throw (Exception *):
	grid(grid_),
	parser(new InputParser()),
	resultfilename(resultfilename_)
{STACK_TRACE(
	NEGF_ASSERT(grid!=NULL && parser!=NULL, "Tried to assign null pointer.");
	NEGF_ASSERT(resultfilename.length() > 0, "empty filename.");
	this->dat_equations.clear();
	this->dat_unittypes.clear();
	this->plt_equations.clear();
	this->plt_unittypes.clear();
	this->plt_datanames.clear();
	this->plt_values.clear();		
);}


void OutputData::add(Equation * eqn, units::UnitType unit)  throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(eqn!=0, "Tried to assign null pointer.");
	
	// look if eqn is to be written in a .dat-file (->vertices, elements, edges) or in a .plt-file
	uint D = grid->get_dimension();
	if (   eqn->get_num_variables()!=this->grid->get_num_vertices()
		&& eqn->get_num_variables()!=this->grid->get_num_elements()
		&& eqn->get_num_variables()!=this->grid->get_num_edges()
		&& eqn->get_num_variables()!=this->grid->get_num_vertices()*D ) 
	{
		if (eqn->get_num_variables()==this->grid->get_num_edges())
			NEGF_EXCEPTION("Although DF-ISE supports location=edge, TECPLOT cannot display it.");
		logmsg->emit(LOG_INFO_L2,"Equation \"%30s\" will be output in a .plt-file.",eqn->get_name().c_str());
		this->plt_equations.push_back(eqn);
		if (eqn->get_num_variables() > 50)
			NEGF_FEXCEPTION("More than 50 datasets would be created from a single equation (%s), which does not make sense.",
							eqn->get_name().c_str());
		if (eqn->get_num_variables() > 20)
			logmsg->emit(LOG_WARN,"*** WARNING: *** equation \"%s\" has more than 20 datasets!",
							eqn->get_name().c_str());
		if (eqn->get_num_variables()==grid->get_num_contacts())
			logmsg->emit(LOG_INFO_L2,"Equation \"%30s\" is recognized as being defined for each contact.",
							eqn->get_name().c_str());
		if (eqn->get_num_variables()==grid->get_num_regions())
			logmsg->emit(LOG_INFO_L2,"Equation \"%30s\" is recognized as being defined for each region.",
							eqn->get_name().c_str());
		for (uint ii = 0; ii < eqn->get_num_variables(); ii++) {
			ostringstream name; name << eqn->get_name() << "_";
			if (eqn->get_num_variables()==grid->get_num_contacts()) {
				name << grid->get_contact(ii)->get_name();
			} else if (eqn->get_num_variables()==grid->get_num_regions()) {
				name << grid->get_region(ii)->get_name();
			} else {
				name << ii;
			}
			this->plt_datanames.push_back(name.str());		// create dataset names
			this->plt_unittypes.push_back(unit);			// create dataset unittype
			vector<double> empty_vec;
			this->plt_values.push_back(empty_vec);
		}
	} else {
		logmsg->emit(LOG_INFO_L2,"Equation \"%30s\" will be output in a .dat-file.",eqn->get_name().c_str());
		this->dat_equations.push_back(eqn);
		this->dat_unittypes.push_back(unit);
	}
);}


/* Create a new set of values to be saved from the current variable values of the plt-equations */
void OutputData::plt_snapshot() throw (Exception *)
{STACK_TRACE(
	uint counter = 0;
	for (uint ii = 0; ii < plt_equations.size(); ii++)
	{
		for (uint jj = 0; jj < plt_equations[ii]->get_num_variables(); jj++)
		{
			this->plt_values[counter].push_back( plt_equations[ii]->get_value(jj) );
			counter++;
		}
	}
	NEGF_ASSERT(counter==plt_values.size(),
			"Did the number of vars in an eqn change? Did the number of assigned eqns change?");
);}


void OutputData::set_filename(string resultfilename_) throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(resultfilename_.length() > 0, "empty filename.");
	resultfilename = resultfilename_;
);}


#ifndef NODFISE
#include "OutputDataDFISE.cpp"
#endif


void OutputData::write_plt() const throw (Exception *)
{STACK_TRACE(
	if (this->plt_equations.size()==0) {
		// there is nothing to do
		return;
	}
	NEGF_ASSERT(this->plt_datanames.size()==this->plt_unittypes.size(), 
					"unittype and dataname arrays do not have the same size");
	NEGF_ASSERT(plt_values.size()>0, "there is nothing to do.");
	NEGF_ASSERT(plt_values.size()==plt_datanames.size(), 
					"values-array has not same length as datanames array.");
	for (uint ii = 1; ii < plt_values.size(); ii++) 
		NEGF_ASSERT(plt_values[ii].size()==plt_values[0].size(), 
					"number of values in the individual datasets must be equal.");
					
	// prepare things needed by DF-ISE output routines
	uint 			num_datasets = plt_datanames.size();
	uint			num_values   = plt_values[0].size();
	vector<double> 	values;	values.resize(num_values*num_datasets, 0.0);
	for (uint ii = 0; ii < num_datasets; ii++)
		for (uint jj = 0; jj < num_values; jj++)
			values[jj*num_datasets+ii] = plt_values[ii][jj];

#ifndef NODFISE
	// write to .plt-file
	this->parser->write_dfise_plt(
			(this->resultfilename + ".plt").c_str(),
			num_values, 
			num_datasets,
			plt_datanames,
			plt_unittypes,
			values );
#endif
	this->parser->write_xy(
			(this->resultfilename + ".xy").c_str(),
			num_values, 
			num_datasets,
			plt_datanames,
			plt_unittypes,
			values );
);}

// needed in NewtonSolver
Equation * OutputData::get_dat_equation(uint ii) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(ii < this->get_num_dat_equations(), "invalid index.");
	return this->dat_equations[ii];
);}

Equation * OutputData::get_plt_equation(uint ii) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(ii < this->get_num_plt_equations(), "invalid index.");
	return this->plt_equations[ii];
);}

