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
#ifndef OPTIONS_H_NEGF
#define OPTIONS_H_NEGF

#include "all.h"

#include "InputParser.h"
#include "PropertyContainer.h"

namespace negf {
	
	class Options {
		
	public:
	
		Options() throw (Exception *);	// will read options from fnames.get_filename()+".cmd"
		~Options() { delete opts; }
				
		/** return indices that correspond to valence degrees of freedom (0-based) */
		void get_valence_degrees_of_freedom(vector<uint> & result) const throw (Exception *);
		
		/** return indices that correspond to conduction degrees of freedom (complementary to valence) (0-based) */
		void get_conduction_degrees_of_freedom(vector<uint> & result) const throw (Exception *);
		
		/** check whether a band is conduction or valence */
		bool is_conduction_band(const uint nn) const throw (Exception *);
		bool is_valence_band(const uint nn) const throw (Exception *);
		
		/** check for the existence of some property */
		bool exists(const string & propertyname) const throw (Exception *);
		
		/** get some property */
		double get(const string & propertyname) const throw (Exception *);

		const map<string, PropertyContainer<double> * > * get_cmdfile() const { return &sim; }
		PropertyContainer<double> * get_container() { return opts; }
	
	
	protected:
		
		map<string, PropertyContainer<double> * > sim;
		PropertyContainer<double> * opts;
		
		vector<uint> conduction_bands;
		vector<uint> valence_bands;
	};
	
} // end of namespace

#endif /*OPTIONS_H_NEGF*/
