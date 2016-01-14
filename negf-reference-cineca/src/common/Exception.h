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
#ifndef _EXCEPTION_H_NEGF
#define _EXCEPTION_H_NEGF

//#include "all.h"	// would not work! don't know why exactly (Logger etc works)
#include <sstream>
#include "stdarg.h"
#include <string.h>	// may compile without it at IIS, but not at home

using namespace std;

/* macro to trace back the hierarchy of called functions when an exception is thrown
   IMPORTANT NOTE: Since the code may contain commas which would be interpreted as 
   				   different arguments to the macro, it must be variadic
   NOTE: Info about filename, line number and compile time is automatically appended in 
         the class (since a new class is constructed at every handover and the 
         constructor appends it)
   USAGE: include this macro in every function like this:
 
 	void Class::functionname(arguments)
 	{STACK_TRACE(
		...
 	);} 
*/

#define STACK_TRACE(...) \
	try { __VA_ARGS__ } catch(Exception *e) { e->append_info(__LINE__,__FILE__,__DATE__,__TIME__, __func__); throw e; }
	//try { __VA_ARGS__ } catch(Exception *e) { NEGF_EXCEPTION(*e->get_reason()); } 

// -----------------------------------------------------
// exception macros
// use FEXCEPTION or FASSERT like sprintf
// -----------------------------------------------------

#define NEGF_EXCEPTION(msg) 		 throw new Exception((msg),                  __LINE__, __FILE__, __DATE__, __TIME__, __func__);
#define NEGF_FEXCEPTION(msg, ...)    throw new Exception(2000, __LINE__, __FILE__, __DATE__, __TIME__, __func__, (msg), __VA_ARGS__);
#define NEGF_ASSERT(cond, msg)		 if(!(cond)) { throw new Exception((msg),    __LINE__, __FILE__, __DATE__, __TIME__, __func__); }
#define NEGF_FASSERT(cond, msg, ...) if(!(cond)) { throw new Exception(2000,  __LINE__, __FILE__, __DATE__, __TIME__, __func__,(msg), __VA_ARGS__); }


namespace negf {

	/** Class to handle exceptions of any kind. Use the macros, never use <assert> or the like. */
	class Exception {
		
	public:
		Exception(const char* preason,  int line, const char* pfile, const char* pdate, const char* ptime, const char* pfunc);
		Exception(unsigned int buflen, int line, const char* pfile, const char* pdate, const char* ptime, 
					const char* pfunc, const char* preason, ...);
		~Exception();
		
		const std::string & get_reason() const	{ return this->reason; }	//!< Return reason for exception (whole string)
		void append_info(int line, const char* pfile, const char* pdate, const char* ptime, const char* pfunc);
		
	protected:
		std::string reason;	//!< string storing information about the thrown Exception
		
	};
	
}
				
#endif /*_EXCEPTION_H_NEGF*/
