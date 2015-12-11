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
#ifndef OUTPUTDATA_H_NEGF
#define OUTPUTDATA_H_NEGF

#include "all.h"

#include "Geometry.h"
#include "Equation.h"
#include "InputParser.h"

using namespace std;

namespace negf {

	/** The class collects datafields (equations) for DF-ISE output */
	class OutputData
	{
		public:
			OutputData(const Geometry * const grid_, const string & resultfilename_) throw (Exception *);
			~OutputData() {}
			
			void 			 add(Equation * eqn, units::UnitType unit) throw (Exception *);
			
			uint			 get_num_dat_equations() 				const { return this->dat_equations.size(); }
			uint			 get_num_plt_equations() 				const { return this->plt_equations.size(); }
			Equation * 		 get_dat_equation(uint ii) 				const throw (Exception *);
			Equation * 		 get_plt_equation(uint ii) 				const throw (Exception *);

            void             write_plt()        const throw (Exception *);  // write .plt-file

#ifndef NODFISE
			void 			 write_dat() 		const throw (Exception *);	// write data defined on vertices/elems
#endif

			// generate a new value entry for the plt-equations
			void			 plt_snapshot() 	throw (Exception *);

			void 			 set_filename(string resultfilename_) 	throw (Exception *);
			string 			 get_filename() 						const { return this->resultfilename; }
			
			
			
		protected:			
			vector<Equation *>				dat_equations;		// equations defined on vertices/elements
			vector<units::UnitType>			dat_unittypes;
			
			vector<Equation *>				plt_equations;	// equations to write into a plt-file
			vector<units::UnitType>			plt_unittypes;
			vector<string>					plt_datanames;
			vector< vector<double> >		plt_values;	
			
			const Geometry * const 			grid;
			const InputParser * const 		parser;
			string							resultfilename;
			
			
	};
	
} // end of namespace

#endif /*OUTPUTDATA_H_NEGF*/
