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
#ifndef _LOGGER_H_NEGF
#define _LOGGER_H_NEGF

#include "all.h"

using namespace std;

namespace negf {
	
	// -----------------------------------------------
	// log message level definition
	// -----------------------------------------------
	enum LoggerLevel {
		LOG_ERROR 	= 1,			//!< Highest level, error messages MUST be shown!
		LOG_WARN	= 2,			//!< Level for warnings (may be encountered a couple of times - not every warning is bad)
		LOG_INFO	= 3,			//!< Normal output
		LOG_INFO_L1 = 4,			//!< Some detailed information
		LOG_INFO_L2 = 5,			//!< Only for debugging
		LOG_INFO_L3 = 6				//!< Ridiculous amount of output, totally messy
	};
	
	/** Logger to manage program output messages.
	    Several "listeners" (output streams) can be added to / deleted from a logger. */
	class Logger{
	public:
	
		Logger(LoggerLevel level_ = LOG_INFO_L2);
		~Logger();

		void add_listener(ostream&);	//!< add an output stream (e.g. std::cout or a file that was opened with ofstream)
		void del_listener(ostream&);	//!< remove an output stream
				
		void emit				(const Exception &) const;
		void emit				(LoggerLevel level, const char* msg, ...) const;
		void emit_noendl		(LoggerLevel level, const char* msg, ...) const;	//!< emit messages without end_of_line
		void emit_all			(LoggerLevel level, const char* msg, ...) const;	//!< not only master thread
		void emit_noendl_all	(LoggerLevel level, const char* msg, ...) const;	//!< not only master thread
		void emit_small_header	(const char* msg, ...) const;	//!< emit small header
		void emit_header		(const char* msg, ...) const;	//!< emit standard header
		void emit_big_header	(const char* msg, ...) const;	//!< emit big header
		void emit_huge_header	(const char* msg, ...) const;	//!< emit huge header

		void init_progress_bar(LoggerLevel level, const char* text, int total);		
		int  set_progress_bar(int current, int total);
		void end_progress_bar(); 
		
		void set_level(LoggerLevel level_);
		LoggerLevel get_level() const { return this->log_level; }
				
	protected:
		
		void emit(const char* msg) const;				//!< basic routine that sends a given string to all outstreams in the list
		void emit_noendl(const char* msg) const;		//!< same as emit, but without end of line
		
		std::list<std::ostream*> 	outstreams;			//!< a message will be emitted to all output streams in this list
		LoggerLevel 				log_level;			//!< the current level of the logger. only messages at or above this level will be emitted.
		LoggerLevel 				progress_bar_level;
		bool 						stdout_registered;	//!< progress bar will not be displayed when std::cout is not in the list
			
	};	
		
}

#endif /*_LOGGER_H_NEGF*/
