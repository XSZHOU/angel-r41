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
#include "Interrupt.h"
using namespace negf;

Interrupt::Interrupt()
{STACK_TRACE(
	static bool initialized = false;
	if (!initialized) {
		signal(SIGINT, &(Interrupt::terminate));
		initialized = true;
	}
);}

void Interrupt::terminate(int sig) throw (Exception *)
{STACK_TRACE(
	NEGF_FEXCEPTION("Interrupt Signal %d detected. Aborting.", sig);
);}
