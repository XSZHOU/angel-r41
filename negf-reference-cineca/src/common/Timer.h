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
#ifndef TIMER_H_NEGF
#define TIMER_H_NEGF

#include <vector>
#include <string>
//#include <sys/time.h>
#include <time.h>

//#include "all.h"
#include "Logger.h"
#include "Exception.h"

using namespace std;

#ifdef WIN32
	#include <sys/timeb.h>
#endif

namespace negf {
	
	
	/** This class stores the time of "clicks", or events
	 *  Each click has 1. a name, 2. a number of seconds since 1.1.1970, 3. a formatted timestring
	 *  All system calls related to time are encapsulated by this class 
	 */
	class Timer{
		
	public:
		Timer();
		~Timer();
		
		// ------------------
		// access functions
		// ------------------

		// returns time since 1.1.1970 in seconds
		double get_starttime() { return clicks[0]; }			//!< returns time since 1.1.1970 in seconds
		
		double get_click(uint ii);								//!< get click number ii
		string get_clickname(uint ii);							//!< get the name of click ii
		double get_click(string name);							//!< get the click with a certain name
		
		double get_seconds_since_start(string name);			//!< get the seconds between the creation of the object an click "name"
		double get_seconds_since_last_click(string name);		//!< get the seconds between the click "name" and the click before
		double get_seconds_between(string name1, string name2);	//!< get the seconds between the clicks "name1" and "name2"
		
		string get_formatted_time(string clickname);			//!< get a nicely formatted string containing the time of click "name"
		
		double get_time();										//!< time in seconds since 1.1.1970 from system
		
		// -------------------
		// setup functions
		// -------------------
		void click(string name);								//!< create a new click with a certain name
		
	protected:
		uint get_index(string name);							//!< get the index in the list of clicks of a click with a certain name
		
		string get_formatted_timestring();						//!< get a nicely formatted string containing the system time
	
		vector<double> clicks;									//!< time in seconds since 1.1.1970 for each click
		vector<string> formatted_times;							//!< for each click a string containing the formatted time
		vector<string> clicknames;								//!< list of click names
			
	};	
		
}	// namespace gain	

#endif /*TIMER_H_NEGF*/
