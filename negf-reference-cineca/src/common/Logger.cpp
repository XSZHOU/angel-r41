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
#include "Logger.h"
using namespace negf;

Logger::Logger(LoggerLevel level_):
	log_level(level_),
	stdout_registered(false)
{
}

Logger::~Logger() {
}

/** set output level */
void Logger::set_level(LoggerLevel level_) 
{STACK_TRACE(
	NEGF_ASSERT(level_ >= LOG_ERROR && level_ <= LOG_INFO_L3, "invalid loglevel set");
	this->log_level = level_;	
);}

/** attach outstream as listener (could be file or terminal) */
void Logger::add_listener(std::ostream &stream) 
{
	this->outstreams.push_back(&stream);		
	if(&stream == &std::cout) {
		this->stdout_registered = true;	
	}
}
/** detach outstream as listener */
void Logger::del_listener(std::ostream &stream) 
{
	this->outstreams.remove(&stream);	
	if(&stream == &std::cout) {
		this->stdout_registered = false;	
	}
}

void Logger::emit(const char* msg) const 
{	
	for(std::list<std::ostream*>::const_iterator it = this->outstreams.begin();	it != this->outstreams.end(); it++) {
		(*(*it)) << msg << std::endl;		
	}
}

void Logger::emit_noendl(const char* msg) const 
{
	for(std::list<std::ostream*>::const_iterator it = this->outstreams.begin();	it != this->outstreams.end(); it++) {
		(*(*it)) << msg;
		(*(*it)).flush();
	}
}

/** emit exception */
void Logger::emit(const Exception& e) const 
{
	this->emit(e.get_reason().c_str());	
}

/** Emit formatted message with specific loglevel
 * @param level  desired loglevel
 * @param buflen length of expected output string (MUST BE BIGGER THAN OUTPUT STRING)
 * @param fmt    format string for sprintf
 * @param ...    argument list 
 * */
void Logger::emit(LoggerLevel level, const char* fmt, ...) const 
{
	if (mpi->get_rank()!=constants::mpi_master_rank) return;
	
	uint buflen = 1000; 
	if(level <= this->log_level) {				
	   va_list args;
	   va_start(args,fmt);
	   char* buf = new char[buflen];
	   vsprintf(buf, fmt, args);
	   this->emit(buf);
	   delete[] buf;
	}
}

/** Emit formatted message, no end_of_line */
void Logger::emit_noendl(LoggerLevel level, const char* fmt, ...) const 
{
	if (mpi->get_rank()!=constants::mpi_master_rank) return;
	
	uint buflen = 1000; 
	if(level <= this->log_level) {				
	   va_list args;
	   va_start(args,fmt);
	   char* buf = new char[buflen];
	   vsprintf(buf, fmt, args);
	   this->emit_noendl(buf);
	   delete[] buf;
	}
}

/** Same as emit, but for every MPI process */
void Logger::emit_all(LoggerLevel level, const char* fmt, ...) const 
{	
	uint buflen = 1000; 
	if(level <= this->log_level) {				
	   va_list args;
	   va_start(args,fmt);
	   char* buf = new char[buflen];
	   vsprintf(buf, fmt, args);
	   this->emit(buf);
	   delete[] buf;
	}
}

/** same as emit_noendl, but for every MPI process */
void Logger::emit_noendl_all(LoggerLevel level, const char* fmt, ...) const 
{	
	uint buflen = 1000; 
	if(level <= this->log_level) {				
	   va_list args;
	   va_start(args,fmt);
	   char* buf = new char[buflen];
	   vsprintf(buf, fmt, args);
	   this->emit_noendl(buf);
	   delete[] buf;
	}
}


void Logger::emit_small_header(const char* fmt, ...) const
{
	uint buflen = 1000;
	va_list args;
	va_start(args,fmt);
	char* buf = new char[buflen];
	vsprintf(buf, fmt, args);

	this->emit(LOG_INFO,"---------------------------------------------------");
	this->emit(LOG_INFO,"%s",buf);
	this->emit(LOG_INFO,"---------------------------------------------------");
	mpi->synchronize_processes(); // otherwise output of non-master-processes will destroy header
}

void Logger::emit_header(const char* fmt, ...) const 
{
	uint buflen = 1000; 			
	va_list args;
	va_start(args,fmt);
	char* buf = new char[buflen];
	vsprintf(buf, fmt, args);
	
	this->emit(LOG_INFO,"|--------------------------------------------------------------------|");
	// total length: 70 characters
	int whites = 68 - strlen(buf);
	int left_whites = whites/2;
	int right_whites = whites-whites/2;
	if (whites<=0) {
		this->emit(LOG_INFO,"|%s|",buf);
	} else {
		string line("|"); line.append(buf); line.append("|");
		line.insert(1, left_whites, ' ');
		line.insert(line.length()-1, right_whites, ' ');
		this->emit(LOG_INFO,"%s",line.c_str());
	}
	this->emit(LOG_INFO,"|--------------------------------------------------------------------|");
	mpi->synchronize_processes(); // otherwise output of non-master-processes will destroy header
}


void Logger::emit_big_header(const char* fmt, ...) const 
{
	uint buflen = 1000; 			
	va_list args;
	va_start(args,fmt);
	char* buf = new char[buflen];
	vsprintf(buf, fmt, args);
	
	this->emit(LOG_INFO,"|====================================================================|");
	// total length: 70 characters
	int whites = 68 - strlen(buf);
	int left_whites = whites/2;
	int right_whites = whites-whites/2;
	if (whites<=0) {
		this->emit(LOG_INFO,"|%s|",buf);
	} else {
		string line("|"); line.append(buf); line.append("|");
		line.insert(1, left_whites, ' ');
		line.insert(line.length()-1, right_whites, ' ');
		this->emit(LOG_INFO,"%s",line.c_str());
	}
	this->emit(LOG_INFO,"|====================================================================|");
	mpi->synchronize_processes(); // otherwise output of non-master-processes will destroy header
}

void Logger::emit_huge_header(const char* fmt, ...) const 
{
	uint buflen = 1000; 			
	va_list args;
	va_start(args,fmt);
	char* buf = new char[buflen];
	vsprintf(buf, fmt, args);
	
	this->emit(LOG_INFO,"|====================================================================|");
	this->emit(LOG_INFO,"|**********                                                **********|");
	// total length: 70 characters
	int whites = 66 - strlen(buf);
	int left_whites = whites/2;
	int right_whites = whites-whites/2;
	if (whites<=0) {
		this->emit(LOG_INFO,"|%s|",buf);
	} else {
		string line("|*"); line.append(buf); line.append("*|");
		line.insert(2, left_whites, ' ');
		line.insert(line.length()-2, right_whites, ' ');
		this->emit(LOG_INFO,"%s",line.c_str());
	}
	this->emit(LOG_INFO,"|**********                                                **********|");
	this->emit(LOG_INFO,"|====================================================================|");
	mpi->synchronize_processes(); // otherwise output of non-master-processes will destroy header
}


/** Initialize a progress bar of the form [XXXXX     ].
    @param level level above an including which the progress bar is displayed
    @param text some text describing what the bar stands for
    @param total a number at which 100% is achieved */
void Logger::init_progress_bar(LoggerLevel level, const char* text, int total) 
{
	if (mpi->get_rank()!=constants::mpi_master_rank) return;
	
	this->progress_bar_level = level;
	if(this->progress_bar_level <= this->log_level) {
		if(this->stdout_registered) {
			std::cout << text << " " << total << "\n";	
			for(int ii = 0; ii < 72; ii++) {
				std::cout << " ";	
			}		
			this->set_progress_bar(0, total);
		}
	}
}

/** Set the new status of a progress bar
    @param current running counter
    @param total number at which 100% is achieved */
int Logger::set_progress_bar(int current, int total) 
{
	if (mpi->get_rank()!=constants::mpi_master_rank) {
		return total - 1;
	}
	if(this->stdout_registered) {	
		int prct_x = (current * 64) / total;
		int prct   = (current * 100) / total;
		if(this->progress_bar_level <= this->log_level) {
			// delete old line
			for(int ii = 0; ii < 72; ii++) {
				std::cout << "\b";	
			}	
			// write new line
			std::cout << " [";
			for(int ii = 0; ii < prct_x; ii++) {
				std::cout << "x";	
			}
			for(int ii = prct_x; ii < 64; ii++) {
				std::cout << " ";	
			}
				
			std::cout << "] ";
			std::cout.width(3);
			std::cout << prct << "%";	
			std::cout.flush();
		}
		// calculate next
		int next = (int)floor((((float)prct + 1.0) * (float)total) / 100.0);
		return next > current ? next : (current + 1);
	} else {
		return total - 1;	
	}			
}

void Logger::end_progress_bar() 
{
	if (mpi->get_rank()!=constants::mpi_master_rank) return;
	
	if(this->progress_bar_level <= this->log_level) {
		if(this->stdout_registered) {
			this->set_progress_bar(100,100);
			std::cout << " ... done\n";	
		}
	}
}

	
