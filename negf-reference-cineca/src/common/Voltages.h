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
#ifndef VOLTAGES_H_NEGF
#define VOLTAGES_H_NEGF

#include "all.h"

#include "PropertyContainer.h"
#include "InputParser.h"
#include "Geometry.h"
#include <boost/lexical_cast.hpp>	// $(INC_PATH) must have boost included!

using namespace std;

namespace negf {

	
	/** Determines the applied voltages (in Volts) at each contact for each simulation step
	 *  and determines which contact is ramped. */
	class Voltages
	{

	public:
	
		Voltages(const Geometry * mastergrid) throw (Exception *);	//!< constructor already defines the entire settings
		~Voltages() {}

		// --------------------------------
		// access functions
		// --------------------------------

		uint get_num_steps() const throw (Exception *);		//!< return the number of simulation steps (1 step = 1 setting of contact voltages)
		const vector<double> & get_setting(unsigned int step) const throw (Exception *);		//!< return the setting of voltages for a specific step. units: Volt
		string get_suffix(unsigned int step) const throw (Exception *);	//!< return a suffix to be appended to the filename for saving the results of a specific step
		void emit_header(unsigned int step) const throw (Exception *);	//!< screen output
		double get_ramped_voltage(uint step) const throw (Exception *);	//!< get the value of the voltage which is ramped in a certain step
		
		double get_underrelaxation(uint step);	//!< get the electrostatic potential underrelaxation parameter for the current step/experiment. returns -1.0 when not set
		double get_late_underrelax(uint step);	//!< get the electrostatic potential underrelaxation parameter for the current step/experiment. returns -1.0 when not set

		double get_value(uint contact_idx) const throw(Exception *);	//!< overwritten functions of the Equation class; returns settings[timestamp][idx]
		
		// from Equation
        uint get_num_variables() const { return this->number_of_variables; }
		uint get_timestamp() const { return timestamp; }
		void set_timestamp(uint new_timestamp) { this->timestamp = new_timestamp; }


	protected:
		const Geometry * 			grid;				//!< for getting the contact names
		double 						series_resistance;	//!< needed for VoltagesTotal?
		
		vector< vector<double> > 	settings;			//!< experiments, voltage ramping within each experiment
		vector<Contact *> 			ramping_contacts;	//!< same length as settings
		vector<double> 				underrelaxation;	//!< same length as settings
		vector<double> 				late_underrelax;	//!< same length as settings
		
		void add_experiment(PropertyContainer<double> * e);	//!< helper function
		
		double compute_value(uint line) const;
		
		// from Equation
		uint timestamp;
		uint number_of_variables;
		vector<double> current_variable_values;
	};
	
} // namespace negf

#endif /*VOLTAGES_H_NEGF*/
