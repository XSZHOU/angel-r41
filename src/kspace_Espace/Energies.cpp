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
#include "Energies.h"
using namespace negf;

Energies::Energies(const Options * options) throw (Exception *):
	energy_grid()
{STACK_TRACE(
	logmsg->emit_header("setting up energies");
	
	// set up min/max energy
	double E0 = constants::convert_from_SI(units::energy, options->get("min_energy") * constants::SIec);
	double E1 = constants::convert_from_SI(units::energy, options->get("max_energy") * constants::SIec);
	
	// set up number of energy points
	double num_energy_points_dbl = options->get("num_energy_points");
	int num_energy_points_int = int(num_energy_points_dbl);
	NEGF_ASSERT(num_energy_points_int>0, "invalid number of energy points.");
	this->num_energy_points = uint(num_energy_points_int);
	
	logmsg->emit(LOG_INFO,"Emin=%.2e, Emax=%.2e, nE=%d (equally spaced in the beginning)",E0,E1,num_energy_points);
		
	// initialize "ignore"-vector
	this->ignore.resize(num_energy_points, false);
	
	// determine which energy points are assigned to the current MPI process
	this->my_process_id = mpi->get_rank();
	this->num_mpi_procs = mpi->get_num_procs();
	NEGF_ASSERT(num_mpi_procs>0 && my_process_id>=0 && my_process_id<num_mpi_procs, "something is wrong.");
	
	// this creates an equally spaced energy grid - "initial guess"
	vector<double> all_energies;
	for (uint ii=0; ii < num_energy_points; ii++) {
		all_energies.push_back(E0 + ii * (E1-E0)/(num_energy_points-1));
	}
	this->assign_new_energy_grid(all_energies);
);}


/** Be sure that ignore[] array was set BEFORE calling this routine! */
void Energies::assign_new_energy_grid(const vector<double> & all_energies)
{STACK_TRACE(
	// store current energy grid in old_energy_grid
	this->old_energy_grid = this->energy_grid;
	
	// assign new energy grid
	NEGF_ASSERT(all_energies.size()==this->num_energy_points, "inconsistent number of energy points.");
	this->energy_grid = all_energies;
	
	// set up weights
	// take into account ignore[] array
	this->weights.assign(this->num_energy_points, 0.0);
	for (uint ii=0; ii < num_energy_points; ii++) {
		uint ii_minus_1 = (ii!=0) ? ii-1 : 0;
		uint ii_plus_1  = (ii!=num_energy_points-1) ? ii+1 : num_energy_points-1;
		if (ignore[ii_minus_1] && ii_minus_1 > 0) {
			ii_minus_1--;
			NEGF_ASSERT(!ignore[ii_minus_1], "There seemed to be two consecutive energy points to be ignored - can't do that (1).");
		}
		if (ignore[ii_plus_1] && ii_plus_1 < num_energy_points-1) {
			ii_plus_1++;
			NEGF_ASSERT(!ignore[ii_plus_1], "There seemed to be two consecutive energy points to be ignored - can't do that (2).");
		}
		
		double E_ii_minus_1 = energy_grid[ii_minus_1];
		double  E_ii_plus_1 = energy_grid[ii_plus_1];
		if (ignore[ii]) {
			this->weights[ii] = 0.0;
		} else {
			this->weights[ii] = 0.5 * (E_ii_plus_1 - E_ii_minus_1);
		}
	}
	
	// average number of energy points per MPI process
	double points_per_process = double(num_energy_points) / num_mpi_procs;
	
	// determine first and last own energy index
	this->my_start_energy_idx = uint(floor(my_process_id*points_per_process)); // starts with 0
	this->my_end_energy_idx   = uint(floor((my_process_id+1)*points_per_process) - 1);
	if (this->my_end_energy_idx > num_energy_points-1) { 
		this->my_end_energy_idx = num_energy_points-1;
	}
	
	// be sure that if there was already an energy grid, the same number of points is associated
	// (because of memory allocation of Green functions etc)
	if (old_energy_grid.size()==energy_grid.size()) {
		NEGF_ASSERT(my_energy_points.size()==this->my_end_energy_idx - this->my_start_energy_idx + 1, 
				"need same number of energy points per process when reassigning energy grid!");
	}
	
	// set up own energy points
	this->my_energy_points.resize(this->my_end_energy_idx - this->my_start_energy_idx + 1);
	for (uint ii=0; ii < my_energy_points.size(); ii++) {
		 my_energy_points[ii] = this->energy_grid[this->my_start_energy_idx + ii];
	}
	
	// set up an array that stores the process computing every energy point
	this->process_for_every_energy_point.resize(num_energy_points, -1);
	for (uint pp=0; pp<uint(num_mpi_procs); pp++)
	{
		for (uint ii=uint(floor(pp*points_per_process)); ii<uint(floor((pp+1)*points_per_process)); ii++) {
			process_for_every_energy_point[ii] = pp;
		}
	}
	for (uint ii=0; ii<process_for_every_energy_point.size(); ii++) {
		NEGF_FASSERT(process_for_every_energy_point[ii]!=-1, "something went wrong: process_for_every_energy_point[%d]=-1.",ii);
	}
	
	// security checks
	if (this->my_start_energy_idx > this->my_end_energy_idx) {
		if (points_per_process<1.0) {
			NEGF_FEXCEPTION("You have too many processors (%d) for the given amount of energy points (%d)", num_mpi_procs, num_energy_points);
		} else {
			NEGF_EXCEPTION("An unknown cause yielded invalid energy point ranges.");
		}
	}
	NEGF_ASSERT(this->get_my_number_of_points()==this->my_end_energy_idx-this->my_start_energy_idx+1, "inconsistent my_number_of_energies.");
	NEGF_ASSERT(this->get_number_of_points(mpi->get_rank())==this->my_end_energy_idx-this->my_start_energy_idx+1, "inconsistent my_number_of_energies.");
	NEGF_FASSERT((int)this->get_process_computing(this->my_start_energy_idx)==my_process_id,
			"problem in the process-splitting of the energies (start): get_process_computing(%d) returns %d instead of %d",
			this->my_start_energy_idx, this->get_process_computing(this->my_start_energy_idx), my_process_id);
	NEGF_FASSERT((int)this->get_process_computing(this->my_end_energy_idx)==my_process_id, 
			"problem in the process-splitting of the energies (end): get_process_computing(%d) returns %d instead of %d",
			this->my_end_energy_idx, this->get_process_computing(this->my_end_energy_idx), my_process_id);
	// check monotonicity
	for (uint ee=0; ee<num_energy_points-1; ee++) {
		NEGF_FASSERT(this->get_energy_from_global_idx(ee) < this->get_energy_from_global_idx(ee+1),
			"nonmonotonic energy grid! E[%d]=%.3e, E[%d]=%.3e", ee,this->get_energy_from_global_idx(ee), ee+1, this->get_energy_from_global_idx(ee+1));
	}
	
	// some screen output
	if (this->my_process_id == constants::mpi_master_rank) {
		for (uint ee=0; ee<num_energy_points; ee++) {
			logmsg->emit(LOG_INFO_L2,"process %2d computes energy %3d (=%.3e), weight %.3e",
					this->get_process_computing(ee), ee,
					this->get_energy_from_global_idx(ee),
					this->get_weight_from_global_idx(ee));
		}
		logmsg->emit(LOG_INFO_L2,"New energy grid: ");
		for (uint ee=0; ee<num_energy_points; ee++) {
			logmsg->emit(LOG_INFO_L2, "    ee=%d: %.3e   ", ee, energy_grid[ee]);
		}
	}
);}

uint Energies::get_process_computing(uint energy_point_index) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(energy_point_index < this->get_number_of_points(), "invalid energy point index.");
	return process_for_every_energy_point[energy_point_index];
);}


uint Energies::get_process_computing(double energy_point_value) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(energy_point_value>=energy_grid[0], "given energy is not computed (too small).");
	uint idx=0; 
	while (energy_grid[idx] < energy_point_value) {
		idx++;
		NEGF_ASSERT(idx<num_energy_points, "given energy is not computed (too big).");
	}
	NEGF_ASSERT(fabs(energy_grid[idx] - energy_point_value) < 1e-14, "energy point does not seem to be computed.");
	return this->get_process_computing(uint(idx));
);}


uint Energies::get_global_index(uint local_idx) const throw (Exception *)
{STACK_TRACE(
	NEGF_FASSERT(local_idx < this->get_my_number_of_points(), "invalid local index %d", local_idx);
	return this->my_start_energy_idx + local_idx;
);}


uint Energies::get_my_idx(uint global_idx) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(global_idx >= this->my_start_energy_idx && global_idx <= this->my_end_energy_idx, 
			"process does not seem to compute given energy point.");
	return global_idx - this->my_start_energy_idx;
);}


/** get the number of points that a certain process computes */
uint Energies::get_number_of_points(uint process_id) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT((int)process_id < this->num_mpi_procs, "invalid process id.");
	
	uint NE = this->get_number_of_points(); // total # energy points
	double ppp = double(NE) / this->num_mpi_procs;
	uint p_start_energy_idx = uint(floor(process_id*ppp));
	uint p_end_energy_idx   = uint(floor((process_id+1)*ppp) - 1);
	if (p_end_energy_idx > NE-1) { 
		p_end_energy_idx = NE-1;
	}
	return p_end_energy_idx - p_start_energy_idx + 1;
);}


uint Energies::get_start_global_idx(uint process_id) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT((int)process_id < this->num_mpi_procs, "invalid process id.");
	int NE = this->get_number_of_points(); // total # energy points
	double ppp = double(NE) / this->num_mpi_procs;
	return uint(floor(process_id*ppp));
);}


uint Energies::get_stop_global_idx(uint process_id) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT((int)process_id < this->num_mpi_procs, "invalid process id.");
	uint NE = this->get_number_of_points(); // total # energy points
	double ppp = double(NE) / this->num_mpi_procs;
	uint result = uint(floor((process_id+1)*ppp) - 1);
	if (result > NE-1) { 
		result = NE-1;
	}
	return result;
);}

