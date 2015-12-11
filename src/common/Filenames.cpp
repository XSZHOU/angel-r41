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
#include "Filenames.h"
using namespace negf;


Filenames::Filenames() throw (Exception *):
	initialized(false)
{STACK_TRACE(
);}


/** Initialize the names of the various directories. This routine must be called exactly once at the beginning of a simulation.
    As a base serves constants::default_aqua_path.
    When the environment variable AQUADIR is set, this will be used instead.
    When the environment variable AQUAOUT is set, this will be used as a base location for the output (results+logfile)
 */
void Filenames::init(const string & name)
{STACK_TRACE(
	
	NEGF_ASSERT(!this->is_initialized(), "Class has already been initialized.");
	
	// divide input filename into bare file and directory
	this->barefile = name;
	int loc = -1;
	do
	{
		this->barefile = this->barefile.substr(loc+1, barefile.size());
		loc = this->barefile.find("/", 0);
	} while (loc != (int) string::npos);
	NEGF_ASSERT(this->barefile.size()>0, "Invalid file.");
	this->filedirectory = name.substr(0,name.size()-this->barefile.size());
	if (filedirectory=="") filedirectory="./";
	this->filename = this->filedirectory + this->barefile;

	// home directory
	this->homedir = getenv("HOME");
	homedir.append((homedir[homedir.size() - 1] == '/' ? "":"/"));
	
	// look if environment variable NEGFDIR is defined
	char * tmp_negfdir = 0;
	tmp_negfdir = getenv("NEGFDIR");	
	// set up "negfdir" which determines where material defs and output are
	if (tmp_negfdir==NULL) {
		this->negfdir = constants::default_negf_dir;
		logmsg->emit(LOG_INFO,"Environment variable NEGFDIR not found. Assuming");
		logmsg->emit(LOG_INFO,"   %s",negfdir.c_str());
		logmsg->emit(LOG_INFO,"as base directory.");
	} else {
		this->negfdir = tmp_negfdir;
		negfdir = negfdir + (negfdir[negfdir.size() - 1] == '/' ? "":"/");
		logmsg->emit(LOG_INFO,"Using %s as base directory for material definitions and output.",
						negfdir.c_str());
	}
	negfdir = negfdir + (negfdir[negfdir.size() - 1] == '/' ? "":"/");
	
	// for all output files, the log file and the material parameters, negfdir is taken as root
	string outfilebase = negfdir;
	char * tmp_outfilebase = 0;
	tmp_outfilebase = getenv("NEGFOUT");
	if (tmp_outfilebase!=NULL) {
		outfilebase = tmp_outfilebase;
		outfilebase = outfilebase + (outfilebase[outfilebase.size() - 1] == '/' ? "":"/");
	}
	
	// output file
	this->outfile_dir = outfilebase + "results/" + barefile + "/";
	logmsg->emit(LOG_INFO,"Output:    %s", outfile_dir.c_str());
	this->outfile = this->outfile_dir + barefile;
	
	// logfile
	this->logfile_dir = outfilebase + "results/logs/";
	char buf[1000];
	string timestring = timer->get_formatted_time("start");
	sprintf(buf,"%snegf_%s.out", logfile_dir.c_str(),timestring.c_str());
	logmsg->emit(LOG_INFO,"Logfile:   %s",buf);
	this->logfile = buf;
	
	// path for material definitions
	this->material_dir = negfdir + "materials/";
	logmsg->emit(LOG_INFO,"Materials: %s",material_dir.c_str());
	#ifdef AMD
	setenv("NEGF_MATERIAL_CNFPATH",material_dir.c_str(),1);  // overwritten if already existing
	#endif
	#ifdef SUN
	putenv("NEGF_MATERIAL_CNFPATH",material_dir.c_str(),1);
	#endif
	#ifdef IA64
	setenv("NEGF_MATERIAL_CNFPATH",material_dir.c_str(),1);
	#endif
	#ifdef PPC
	setenv("NEGF_MATERIAL_CNFPATH",material_dir.c_str(),1);
	#endif
	
	this->initialized = true;
);}


void Filenames::set_outfile_directory_suffix(const string & suffix)
{STACK_TRACE(
	//this->outfiledirectory = this->negfdir + "results/" + this->barefile + suffix + "/";
	string new_outfiledir = this->outfile_dir;
	// strip "/"
	if (new_outfiledir[new_outfiledir.size() - 1] == '/') {
		new_outfiledir = new_outfiledir.substr(0, new_outfiledir.size() - 1); // second argument is length, not last char position
	}
	this->outfile_dir = new_outfiledir + suffix + "/";
	logmsg->emit(LOG_INFO,"New output dir: %s", outfile_dir.c_str());
	this->outfile = this->outfile_dir + this->barefile;
);}

void Filenames::set_outfile_directory(const string & new_outdir)
{STACK_TRACE(
	this->outfile_dir = new_outdir;
	// strip "/"
	//if (outfile_dir[outfile_dir.size() - 1] == '/') {
	//	outfile_dir = outfile_dir.substr(0, outfile_dir.size() - 1); // second argument is length, not last char position
	//}
	logmsg->emit(LOG_INFO,"New output dir: %s", outfile_dir.c_str());
	this->outfile = this->outfile_dir + this->barefile;
);}

const string & Filenames::get_homedirectory() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->is_initialized(), "initialize filenames first.");
	return this->homedir;
);}
const string & Filenames::get_filedirectory() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->is_initialized(), "initialize filenames first.");
	return this->filedirectory; 
);}
const string & Filenames::get_barefile() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->is_initialized(), "initialize filenames first.");
	return this->barefile;
);} 
const string & Filenames::get_filename() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->is_initialized(), "initialize filenames first.");
	return this->filename; 
);}
const string & Filenames::get_path() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->is_initialized(), "initialize filenames first.");
	return this->negfdir; 
);}
const string & Filenames::get_logfile() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->is_initialized(), "initialize filenames first.");
	return this->logfile;
);}
const string & Filenames::get_logfiledirectory() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->is_initialized(), "initialize filenames first.");
	return this->logfile_dir;
);}
const string & Filenames::get_materialdirectory()	const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->is_initialized(), "initialize filenames first.");
	return this->material_dir; 
);}
const string & Filenames::get_outfile() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->is_initialized(), "initialize filenames first.");
	return this->outfile;
);}
const string & Filenames::get_outfiledirectory() const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->is_initialized(), "initialize filenames first.");
	return this->outfile_dir;
);}
