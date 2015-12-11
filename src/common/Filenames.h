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
#ifndef FILENAMES_H_NEGF
#define FILENAMES_H_NEGF

#include "all.h"

using namespace std;

namespace negf {
	

	/** Wrapper to store filenames, paths etc.
	 *  note: directories are always terminated with "/"
	 */
	class Filenames
	{

	public:
		Filenames() throw (Exception *);	//!< created in all.cpp
		~Filenames() {}
		
		// --------------------------------
		// setup functions
		// --------------------------------
		void init(const string & name); //!< "name" was the argument handed over to the program by the user
		void set_outfile_directory_suffix(const string & suffix);	//!< append something to the name of the output directory
		void set_outfile_directory(const string & dir); //!< change entire output directory
		
		// --------------------------------
		// access functions
		// --------------------------------
		bool 		 is_initialized() 			const { return initialized; }
		
		const string & get_homedirectory() 		const throw (Exception *);	//!< returns the directory $HOME of the user
		const string & get_filedirectory()		const throw (Exception *);	//!< returns the simulation file directory
		
		const string & get_barefile() 			const throw (Exception *);	//!< returns the simulation name, e.g. "nano" when the simulation is nano.cmd, nano.grd etc.
		const string & get_filename() 			const throw (Exception *);	//!< returns the simulation name including directory, but without ending
		
		const string & get_path() 				const throw (Exception *);
		
		const string & get_logfile() 			const throw (Exception *);	//!< returns the filename of the logfile, including directory
		const string & get_logfiledirectory() 	const throw (Exception *);	//!< returns the directory where the logfile is created
		const string & get_materialdirectory()	const throw (Exception *);	//!< returns the directory where the material definitions are stored
		const string & get_outfile() 			const throw (Exception *);
		const string & get_outfiledirectory() 	const throw (Exception *);	//!< returns the output directory
		
		
	protected:
		string homedir;				//!< home directory
		string filedirectory;		//!< directory of simulation file (.grd, .dat)
		string barefile;			//!< name of simulation file (.grd, .dat) without path
		string filename;			//!< name of simulation file (.grd, .dat) including path
		string negfdir;				//!< =NEGFDIR env variable is exists, else constants::default_aqua_path
		string logfile;				//!< file where simulator messages are stored
		string logfile_dir;			//!< directory of the logfile
		string material_dir;		//!< path for material definitions
		string outfile;				//!< simulation output ("_voltage_..." gets appended)
		string outfile_dir;			//!< directory of simulation output
	
		bool initialized;			//!< becomes true when init was called
		
	};	

} // end namespace 
		
#endif /*FILENAMES_H_NEGF*/
