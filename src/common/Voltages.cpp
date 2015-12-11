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
#include "Voltages.h"
using namespace negf;


Voltages::Voltages(const Geometry * mastergrid) throw (Exception *):
	grid(mastergrid)
{STACK_TRACE(
	logmsg->emit_header("voltage settings, which contacts are ramped by how much");
	NEGF_ASSERT(grid!=NULL, "tried to assign NULL pointer.");
	
	// GET ALL THE INFORMATION
	logmsg->emit(LOG_INFO,"Reading %s...", fnames->get_filename().c_str());
	InputParser parser;
	map<string, PropertyContainer<double>* > sim = parser.read_cmd_file(fnames->get_filename());
	
	// determine all the experiments
	map< int, PropertyContainer<double>* > experiments;
	ostringstream msg; msg << "Experiments with numbers ";
	map< string, PropertyContainer<double>* >::iterator it;
	for (it = sim.begin(); it != sim.end(); ++it)
	{
		string::size_type loc = it->first.find("experiment_");
		if (loc==0) {
			int exp_no = boost::lexical_cast<int>( it->first.substr(11,string::npos) );
			NEGF_ASSERT(experiments.find(exp_no)==experiments.end(), "duplicate experiment detected");
			experiments[exp_no] = it->second;
			msg << exp_no << ", ";
		}
	}
	msg << "have been detected. Will perform them in sorted order."; logmsg->emit(LOG_INFO, msg.str().c_str());
	
	// sort the experiments
	vector< PropertyContainer<double>* > sorted_experiments;
	sorted_experiments.resize(experiments.size(), NULL);
	map< int, PropertyContainer<double>* >::iterator it2, it3;
	for (it2 = experiments.begin(); it2 !=experiments.end(); ++it2)
	{
		// determine how many other experiments have smaller number
		uint smaller_number = 0;
		for (it3 = experiments.begin(); it3 !=experiments.end(); ++it3)
		{
			if (it3->first < it2->first) {
				smaller_number++;
			}
		}
		sorted_experiments[smaller_number] = it2->second;
	}
	for (uint ii = 0; ii < sorted_experiments.size(); ii++) {
		NEGF_ASSERT(sorted_experiments[ii]!=NULL, "did not find an experiment.");
	}
	
	// add the experiments to the step list, one after the other
	this->settings.clear();
	this->ramping_contacts.clear();
	for (uint exp = 0; exp < sorted_experiments.size(); exp++) {
		this->add_experiment(sorted_experiments[exp]);
	}
	NEGF_ASSERT(settings.size()>0, "No voltage settings were found.");
	
	// check for series resistance
	PropertyContainer<double> * options = 0;
	for (it = sim.begin(); it != sim.end(); ++it)
	{
		string::size_type loc = it->first.find(constants::options_name);
		if (loc==0 && it->first.size()==constants::options_name.length()) {
			NEGF_FASSERT(options==0, "duplicate identifier \"%s\" detected.", constants::options_name.c_str());
			options = it->second;
		}
	}
	NEGF_FASSERT(options!=0, "identifier \"%s\" was not found in .cmd-file.", constants::options_name.c_str());
	if (options->is_set("SeriesResistance")) {
		NEGF_ASSERT(grid->get_num_contacts()==2, "\"SeriesResistance\" can only be used at the moment when exactly 2 contacts exist");
		this->series_resistance = options->get("SeriesResistance");
		logmsg->emit(LOG_INFO,"Series resistance of %e Ohms detected.", series_resistance);
	} else {
		this->series_resistance = 0.0;
	}
	
	// clean up the created PropertyContainers
	for (it = sim.begin(); it != sim.end(); ++it) {
		delete it->second;
	}
	
	// standard eqn stuff
	//this->its_type = quantities::potential;
	this->number_of_variables = grid->get_num_contacts();
	//this->dependencies.clear();	// no dependencies
	//this->set_name("voltage");
	this->timestamp = 0;
	NEGF_ASSERT(settings[0].size()==this->number_of_variables, "something is wrong.");
	this->current_variable_values = settings[0];
);}


// overwritten from Equation class ---> get_value(ii) returns settings[timestamp][ii]
// ATTENTION: MAKE SURE TO CORRECTLY SET THE TIMESTAMP BEFORE ACCESSING!
double Voltages::get_value(uint contact_idx) const throw(Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->timestamp < settings.size(), "timetamp exceeds number of voltage steps.");
	NEGF_ASSERT(contact_idx < settings[this->timestamp].size(), "wrong index.");
	return settings[this->timestamp][contact_idx];
);}


void Voltages::add_experiment(PropertyContainer<double> * ex)
{STACK_TRACE(
	
	/** strategy:
	- check for each mastergrid contact if a property "<name>_voltage" exists.
	- if no, then "<name>_min", "<name>_max" and "<name>_step" must be defined
	- ramping can occur at exactly one contact for now!
	*/
		
	/* --------------------------------------------------------------------------
	 * determine which contact is ramped
	 * check for that contact if ramping info is available
	 * check also that when a contact is not ramped, ramping info is not available.
	 * ----------------------------------------------------------------------------*/
	Contact * ramping_contact = NULL;
	double ramp_min = 0.0;
	double ramp_max = 0.0;
	double ramp_step = 0.0;
	NEGF_ASSERT( grid->get_num_contacts() > 0, "No contacts!");
	for (uint cc = 0; cc < grid->get_num_contacts(); cc++)
	{
		Contact * contact = grid->get_contact(cc);
		string name_voltage = contact->get_name() + "_voltage";
		string name_min     = contact->get_name() + "_min";
		string name_max     = contact->get_name() + "_max";
		string name_step    = contact->get_name() + "_step";
		if (!ex->is_set(name_voltage)) {
			NEGF_ASSERT(ramping_contact == NULL, "Can only ramp one contact.");
			ramping_contact = contact;
			
			// check if ramping info for that contact is available; read this info
			NEGF_ASSERT(ex->is_set(name_min), "<name>_min not found.");
			NEGF_ASSERT(ex->is_set(name_max), "<name>_max not found.");
			NEGF_ASSERT(ex->is_set(name_step), "<name>_step not found.");
			ramp_min = ex->get(name_min);
			ramp_max = ex->get(name_max);
			ramp_step = ex->get(name_step);
		} else {
			// ramping info for that contact MUST NOT be available
			NEGF_ASSERT(!ex->is_set(name_min), "<name>_min found even though <name>_voltage is defined.");
			NEGF_ASSERT(!ex->is_set(name_max), "<name>_max found even though <name>_voltage is defined.");
			NEGF_ASSERT(!ex->is_set(name_step), "<name>_step found even though <name>_voltage is defined.");
		}
	}
	NEGF_ASSERT(ramping_contact!=NULL, 
			"Fixed voltages were defined for all contacts. Did not find any contact to ramp.");
	uint ramp = ramping_contact->get_index();
	
	NEGF_FASSERT(ramp_step!=0.0 && negf_math::sign(ramp_step)==negf_math::sign(ramp_max-ramp_min), 
			"invalid ramping command: ramp_step=%e, ramp_max-ramp_min=%e", ramp_step, ramp_max-ramp_min);
	
	// add the settings produced by the ramping during this experiment to the global settings list
	vector<double> setting; setting.resize(grid->get_num_contacts(), 0.0);
	
	double experiment_underrelaxation = -1.0;
	if (ex->is_set("PotentialUnderrelaxation")) {
		experiment_underrelaxation = ex->get("PotentialUnderrelaxation");
		NEGF_ASSERT(experiment_underrelaxation>=0.0-1e-10 && experiment_underrelaxation<=1.0+1e-10, "underrelaxation must be between 0 and 1.");
	}
	double experiment_late_underrelax = -1.0;
	if (ex->is_set("LateUnderrelaxation")) {
		experiment_late_underrelax = ex->get("LateUnderrelaxation");
		NEGF_ASSERT(experiment_late_underrelax>=0.0-1e-10 && experiment_late_underrelax<=1.0+1e-10, "late underrelaxation must be between 0 and 1.");
	}

	for (double volt = ramp_min; (ramp_step > 0) ? (volt <= ramp_max) : (volt >= ramp_max); volt += ramp_step)
	{
		ostringstream msg; msg << "   Step "<<settings.size()<<": ";
		for (uint cc = 0; cc < grid->get_num_contacts(); cc++)
		{
			if (cc != ramp) {
				setting[cc] = ex->get(grid->get_contact(cc)->get_name() + "_voltage");
			} else {
				setting[cc] = volt;
			}
			msg << grid->get_contact(cc)->get_name() <<"="<<setting[cc]<<", ";
		}
		this->ramping_contacts.push_back(ramping_contact);
		this->underrelaxation.push_back(experiment_underrelaxation);
		this->late_underrelax.push_back(experiment_late_underrelax);
		this->settings.push_back(setting);
		logmsg->emit(LOG_INFO_L2,msg.str().c_str());
	}	
);}


uint Voltages::get_num_steps() const throw (Exception *)
{STACK_TRACE(
	return this->settings.size();
);}


const vector<double> & Voltages::get_setting(uint step) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(settings.size()>step, "invalid step index.");
	NEGF_ASSERT((settings[step]).size()==grid->get_num_contacts(), "wrong voltage vector.");
	
	return settings[step];
);}


string Voltages::get_suffix(uint step) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(step<settings.size(), "invalid step index.")
	NEGF_ASSERT(ramping_contacts.size()==settings.size(), "something is wrong.");
	
	char buf[1000];
	sprintf(buf,"_%s%05.3fV", ramping_contacts[step]->get_name().c_str(), 
							 settings[step][ramping_contacts[step]->get_index()]);
	return buf;
);}


double Voltages::get_underrelaxation(uint step)
{STACK_TRACE(
	NEGF_ASSERT(step<settings.size(), "invalid step index.")
	return underrelaxation[step];
);}

double Voltages::get_late_underrelax(uint step)
{STACK_TRACE(
	NEGF_ASSERT(step<settings.size(), "invalid step index.")
	return late_underrelax[step];
);}

void Voltages::emit_header(uint step) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT((settings[step]).size()==grid->get_num_contacts(), "invalid settings matrix");
	NEGF_ASSERT(ramping_contacts.size()==settings.size(), "something is wrong.");
	
	string msg  = ramping_contacts[step]->get_name();
	double volt = settings[step][ramping_contacts[step]->get_index()];
	
	if (mpi->get_rank()==constants::mpi_master_rank)
	{
		logmsg->emit(LOG_INFO,"");
		logmsg->emit(LOG_INFO,"|=======================================================================================|");
		logmsg->emit(LOG_INFO,"|********************                                               ********************|");
		logmsg->emit(LOG_INFO,"|********************      %16s: %6.3fV                ********************|",msg.c_str(),volt);
		logmsg->emit(LOG_INFO,"|********************                                               ********************|");
		logmsg->emit(LOG_INFO,"|=======================================================================================|");
	}
);}


double Voltages::get_ramped_voltage(uint step) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT((settings[step]).size()==grid->get_num_contacts(), "invalid settings matrix");
	NEGF_ASSERT(ramping_contacts.size()==settings.size(), "something is wrong.");
	
	return settings[step][ramping_contacts[step]->get_index()];
);}


double Voltages::compute_value(uint line) const 
{  
	/* do nothing, but do NOT throw Exception */ 
	return 0.0;
}

