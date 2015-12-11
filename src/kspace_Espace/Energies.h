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
#ifndef ENERGIES_H_
#define ENERGIES_H_

#include "all.h"

#include "Options.h"
//#include "DomainPoint.h"
//#include "DomainMaster.h"
//#include "Domain.h"

namespace negf {
	
	class Energies {
		
	public:
	
		Energies(const Options * options) throw (Exception *);		// sets up energy grid
		~Energies() {}
	
		// ------------------
		// access functions
		// ------------------
		
		uint 	get_number_of_points()    const { return this->num_energy_points; }	// overall number of energy points
		uint 	get_my_number_of_points() const { return my_energy_points.size(); } // number of energy points in this process (MPI)
		uint 	get_start_global_idx()    const { return my_start_energy_idx; }
		uint 	get_stop_global_idx()     const { return my_end_energy_idx; }
		uint 	get_number_of_points(uint process_id) const throw (Exception *);    // number of energy points for a certain MPI process
		uint 	get_start_global_idx(uint process_id) const throw (Exception *);
		uint 	get_stop_global_idx(uint process_id) const throw (Exception *);
		
		double 	get_energy_from_global_idx(uint global_Eidx) const { return energy_grid[global_Eidx]; }
		double 	get_energy_from_local_idx(uint local_Eidx)   const { return my_energy_points[local_Eidx]; }
		double  get_weight_from_global_idx(uint global_Eidx) const { return weights[global_Eidx]; }
		
		// set/get whether energy points are ignored
		bool 	is_ignored(uint global_Eidx) 				const { return ignore[global_Eidx]; }
		void 	set_ignored(uint global_Eidx, bool yes_or_no) 	  { ignore[global_Eidx] = yes_or_no; }
		const vector<bool> & get_ignore_array() 			const { return ignore; }
		void 	set_ignore_array(const vector<bool> & new_values) { ignore = new_values; }
		
		// returns the global index from the index in the process energy point list
		uint	get_global_index(uint my_idx) const throw (Exception *);
		
		// returns the index in the energy point list of the current process from the global energy grid index
		uint 	get_my_idx(uint global_idx) const throw (Exception *);
		
		// returns the process that computes a certain energy index
		uint 	get_process_computing(uint energy_point_index) const throw (Exception *);
		
		// returns the process that computes a certain energy
		uint 	get_process_computing(double energy_point_value) const throw (Exception *);
		
		const vector<double> & get_energy_grid()     const { return this->energy_grid; }
		const vector<double> & get_old_energy_grid() const { return this->old_energy_grid; }
		
		// set new energy grid! includes recomputation of weights
		void 	assign_new_energy_grid(const vector<double> & all_energies);

	protected:
	
		vector<double> energy_grid;
		vector<double> weights;
		vector<double> my_energy_points;
		vector<bool>   ignore;	// gives possibility for an energy point to be ignored in the energy integration
		
		int my_process_id;
		int num_mpi_procs;
		uint my_start_energy_idx;	// index of the first energy point computed by this thread
		uint my_end_energy_idx;		// index of the last energy point computed by this thread
		uint num_energy_points;
		
		vector<int> process_for_every_energy_point;	// stores for EVERY energy point the process that computes it
		
		vector<double> old_energy_grid;	// in case of a change in the energy grid, this variable stores the last one
	};
	
} // end of namespace

#endif /*ENERGIES_H_*/
