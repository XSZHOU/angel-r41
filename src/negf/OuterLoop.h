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
#ifndef OUTERLOOP_H_
#define OUTERLOOP_H_

#include "all.h"

#include "InnerLoop.h"
//#include "Equation.h"
#include "PoissonProblem.h"
#include "Resonances.h"
#include "InputParser.h"

namespace negf {
	
	
	/** Reach self-consistency for a fixed electrostatic field 
	 *  The calculation of the GF is performed explicitly in this class, whereas
	 *  The SE calculation is outsourced to SelfEnergies */
	class OuterLoop {
	public:
	
		OuterLoop(InnerLoop * innerloop_, PoissonProblem * poisson_, const string & outfilename_) throw (Exception *);
		~OuterLoop() {}
	
		void 			perform() throw (Exception *);
	
		void 			update_voltages(const vector<double> & voltages_) throw (Exception *);
		void 			update_underrelaxation(double new_underrelaxation) throw (Exception *);
		void 			update_late_underrelax(double new_late_underrelax) throw (Exception *);
		
		void 			save_everything_to_file(const char * filebase) throw (Exception *);
		
		// steering
		void 			output_debug(bool yes_or_no) { this->output_after_each_iteration = yes_or_no; }
		void 			set_max_outer_iterations(uint new_number);
		void 			set_outer_convergence_crit(double new_crit);
		
	protected:
	
		void 			iterate(uint outer_iter, const double inner_err_crit, bool divergence_avoiding_only);
		bool 			converged();
			
		void 			compute_fermilevels();
		void 			ramp_scattering(uint iter);
		
		uint 			max_outer_iterations;
		double 	        rel_error;
		bool 			output_after_each_iteration;
		
		double 			potential_underrelaxation;
		double 			late_underrelax;
		vector< vector<double> > poisson_solutions;
		vector< vector<double> > mixed_potentials;
		
		InnerLoop	   * innerloop;
		PoissonProblem * poisson;
		
		// quantities which actually belong to the inner loop, but pointers are stored here anyway for convenience
		const Geometry * xspace;
		const Kspace   * kspace;
		const Options  * options;
		Energies       * energies;
		GreenFunctions * gf;
		SelfEnergies   * se;
		PostProcessing * pp;
		Resonances     * res;
		
		vector<double> 	last_electrostatic_potential;
		bool 			ramping_finished;	// will be true after ramping of scattering mechanisms was done once
		
		vector<double>  fermilevels;	// stores last calculated fermilevels
		vector<double>  voltages;		// stores last assigned voltages
		
		string          outfilename;

		bool ok;
	};
	
} // end namespace negf

#endif /*OUTERLOOP_H_*/
