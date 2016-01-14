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
#ifndef INTERRUPT_H_NEGF
#define INTERRUPT_H_NEGF

#include <csignal>
#include "all.h"

using namespace std;

namespace negf {

	/** This class makes the program terminatable (Ctrl-C) by creating a signal handler
	 *  Usually this is normally performed in main.cpp, but since the program is a library thats not possible 
	 *  Now a global instance is created in all.cpp; do not create this class more than once. */
	class Interrupt{
		
	public:
		Interrupt();		//!< created in all.cpp
		~Interrupt() {}
		
		static void terminate(int sig) throw (Exception *);	//<! will be called when pressing Ctrl-C
	};	
		
}	// namespace negf	

#endif /*INTERRUPT_H_NEGF*/
