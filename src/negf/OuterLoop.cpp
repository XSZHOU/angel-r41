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
#include "OuterLoop.h"
using namespace negf;


OuterLoop::OuterLoop(InnerLoop * innerloop_, PoissonProblem * poisson_, const string & outfilename_) throw (Exception *):
	max_outer_iterations       (constants::max_outer_iters),
	rel_error                  (constants::convert_from_SI(units::potential, constants::outer_errcrit)),
	output_after_each_iteration(false),
	potential_underrelaxation  (0.0),
	late_underrelax            (0.0),
	innerloop                  (innerloop_),
	poisson                    (poisson_),
	xspace                     (innerloop->get_xspace()),
	kspace                     (innerloop->get_kspace()),
	options                    (innerloop->get_options()),
	energies                   (innerloop->get_energies()),
	gf                         (innerloop->get_green_functions()),
	se                         (innerloop->get_self_energies()),
	pp                         (innerloop->get_post_processing()),
	ramping_finished           (false),
	outfilename                (outfilename_)
{STACK_TRACE(
	logmsg->emit_header("setting up outer loop Poisson-NEGF");
	NEGF_ASSERT(innerloop!=NULL && poisson!=NULL, "null pointer encountered.");
	NEGF_ASSERT(xspace!=NULL && kspace!=NULL && options!=NULL && energies!=NULL, "null pointer encountered (2)");
	NEGF_ASSERT(gf!=NULL && se!=NULL && pp!=NULL, "null pointer encountered (3)");
		
	this->res = new Resonances(innerloop);
	
	if (options->exists("PotentialUnderrelaxation")) {
		this->potential_underrelaxation = options->get("PotentialUnderrelaxation");
		NEGF_ASSERT(this->potential_underrelaxation>=0.0 && this->potential_underrelaxation<1.0, "invalid potential underrelaxation.");
		logmsg->emit(LOG_INFO,"Electrostatic potential will be underrelaxed by %.1g",this->potential_underrelaxation);
	}
	this->late_underrelax = this->potential_underrelaxation;
	
	this->last_electrostatic_potential.resize(poisson->get_grid()->get_num_vertices(), 0.0);
	mpi->synchronize_processes();
);}


void OuterLoop::set_max_outer_iterations(uint new_number) 
{
	NEGF_ASSERT(new_number>0 && new_number<1000, "invalid maximum number of outer iterations.");
	logmsg->emit(LOG_INFO,"Maximum number of outer iterations is now %d.", new_number);
	this->max_outer_iterations = new_number;
}


void OuterLoop::set_outer_convergence_crit(double new_crit) 
{
	NEGF_ASSERT(new_crit>0.0 && new_crit<1.0, "invalid outer convergence criterion.");
	logmsg->emit(LOG_INFO,"New outer convergence criterion is %g.", new_crit);
	this->rel_error = new_crit;
}


void OuterLoop::perform() throw (Exception *)
{STACK_TRACE(
	int iter;
	for (iter = 1; iter <= (int)this->max_outer_iterations; iter++)
	{
		logmsg->emit_huge_header("Outer iteration %d...", iter);
		
		// -----------------------------------------------------------
		// maybe perform ramping in first iteration at first voltage
		// -----------------------------------------------------------
		if (iter==1 && !this->ramping_finished) {
			this->ramp_scattering(iter);
		}
		
		// ------------------------------------------------------------------------------------------------------
		// perform iteration (calculation of GR,GL,GG and then the new SigmaR, SigmaL, SigmaG) on every process
		// ------------------------------------------------------------------------------------------------------
		bool divergence_avoiding_only = (iter==1) ? true : false;
		this->iterate(iter, options->get("inner_errcrit"), divergence_avoiding_only);
		mpi->synchronize_processes();		
		
		// ----------------------------------------
		// master thread tests convergence
		// ----------------------------------------
		int root = constants::mpi_master_rank;
		int conv;
		if(mpi->get_rank()==root) {
			conv = (this->converged()) ? 1 : 0;	
			if (iter==1) {
				logmsg->emit(LOG_INFO,"First outer iteration - no convergence.");
				conv = 0;
			}
		}
		// communicate to all processes if convergence is reached
		mpi->broadcast(conv, root);	// performed by all processes, conv is send for mpi_master, recv for rest
		if (conv==1) {
			break;
		}	
		if (this->output_after_each_iteration) 
		{
			this->pp->compute_scattering_current(SEtype_contact);
			if (this->se->has_self_energy(SEtype_optical_phonon)) {
				this->pp->compute_scattering_current(SEtype_optical_phonon);
			}
			if (this->se->has_self_energy(SEtype_spont_photon)) {
				this->pp->compute_scattering_current(SEtype_spont_photon);
			}
			if (mpi->get_rank()==root) {
				char buf[1000];
				sprintf(buf,"%s_step%d",outfilename.c_str(),iter);
				this->save_everything_to_file(buf);
			}
		}
	}
	if (iter==(int)max_outer_iterations+1) 
	{
		// compute spectral density and density (the latter being stored in the master MPI process)
		this->pp->compute_local_dos();
		this->pp->compute_spectral_edensity();
		this->pp->compute_edensity();
		this->pp->compute_spectral_hdensity();
		this->pp->compute_hdensity();
		// compute luminescence
		if (se->has_self_energy(SEtype_spont_photon)) {
			pp->compute_scattering_current(SEtype_spont_photon);
			pp->compute_luminescence();
		}	
		// save to file
		if (mpi->get_rank()==constants::mpi_master_rank) {
			char buf[1000];
			sprintf(buf,"%s_NOTCONVERGED",outfilename.c_str());
			this->save_everything_to_file(buf);

			pp->get_contact_current()->snapshot(888.888);
			pp->get_contact_current()->save_to_file(buf);
			if (se->has_self_energy(SEtype_spont_photon)) {
				pp->get_luminescence2()->snapshot(888.888);
				pp->get_luminescence2()->write_power_to_file(buf);
			}
		}
		
		mpi->synchronize_processes(); // important! otherwise exiting is done before master rank can write files
		NEGF_EXCEPTION("\n\nMaximum number of outer iterations reached without convergence.");
	} else {
		logmsg->emit(LOG_INFO, "Outer Convergence reached in %d iterations", iter);
	}
);}


void OuterLoop::ramp_scattering(uint iter)
{STACK_TRACE(
	// determine if there are self-energies present to ramp, and adjust their strength
	double decrease = (options->exists("ScatteringDecreaseFactor"))	? options->get("ScatteringDecreaseFactor") : constants::scattering_decrease;
	NEGF_ASSERT( decrease > 0.0 && decrease <= 1.0, "scattering decrease factor must be within (0,1].");
	bool ramp = false;
	if (se->has_self_energy(SEtype_buettiker)) {
		se->get_buettiker_selfenergy()->set_energy_parameter(se->get_buettiker_selfenergy()->get_energy_parameter() * decrease);
		ramp = true;
	}
	if (se->has_self_energy(SEtype_golizadeh_momentum)) {
		se->get_golizadeh_m_selfenergy()->set_energy_parameter(se->get_golizadeh_m_selfenergy()->get_energy_parameter() * decrease);
		ramp = true;
	}
	if (se->has_self_energy(SEtype_golizadeh_phase)) {
		se->get_golizadeh_p_selfenergy()->set_energy_parameter(se->get_golizadeh_p_selfenergy()->get_energy_parameter() * decrease);
		ramp = true;
	}
	if (se->has_self_energy(SEtype_optical_phonon)) {
		se->get_optical_phonon_selfenergy()->set_scaling(decrease);
		ramp = true;
	}
	if (se->has_self_energy(SEtype_acoustic_phonon)) {
		se->get_acoustic_phonon_selfenergy()->set_scaling(decrease);
		ramp = true;
	}
	if (se->has_self_energy(SEtype_spont_photon)) {
		se->get_spontaneous_photon_selfenergy()->set_scaling(decrease);
		ramp = true;
	}
	if (se->has_self_energy(SEtype_ion_imp)) {
		se->get_ionized_impurities_selfenergy()->set_scaling(decrease);
		ramp = true;
	}
	
	if (ramp) 
	{
		// for each damping factor perform a single outer iteration, meaning the entire inner loop followed by a Poisson update	
		
		double ramp_factor = (options->exists("ScatteringRampFactor")) ? options->get("ScatteringRampFactor") : constants::scattering_ramp;
		
		// very very first iteration: solve ballistic problem
		logmsg->emit_huge_header("Outer iteration with zero scattering");
		
		// get final scattering parameters
		double buett = (se->has_self_energy(SEtype_buettiker))          ? se->get_buettiker_selfenergy()->get_energy_parameter()   : 0.0;
		double golim = (se->has_self_energy(SEtype_golizadeh_momentum)) ? se->get_golizadeh_m_selfenergy()->get_energy_parameter() : 0.0;
		double golip = (se->has_self_energy(SEtype_golizadeh_phase))    ? se->get_golizadeh_p_selfenergy()->get_energy_parameter() : 0.0;
		
		// do an inner iteration w/ zero scattering
		if (se->has_self_energy(SEtype_buettiker))          se->get_buettiker_selfenergy()->set_energy_parameter(0.0);
		if (se->has_self_energy(SEtype_golizadeh_momentum)) se->get_golizadeh_m_selfenergy()->set_energy_parameter(0.0);
		if (se->has_self_energy(SEtype_golizadeh_phase))    se->get_golizadeh_p_selfenergy()->set_energy_parameter(0.0);
		if (se->has_self_energy(SEtype_optical_phonon))     se->get_optical_phonon_selfenergy()->set_scaling(0.0);
		if (se->has_self_energy(SEtype_acoustic_phonon))    se->get_acoustic_phonon_selfenergy()->set_scaling(0.0);
		if (se->has_self_energy(SEtype_spont_photon))       se->get_spontaneous_photon_selfenergy()->set_scaling(0.0);
		if (se->has_self_energy(SEtype_ion_imp))       		se->get_ionized_impurities_selfenergy()->set_scaling(0.0);
		this->iterate(iter, constants::inner_ramp_errcrit, true);
		
		// set scattering parameter back to initial value
		if (se->has_self_energy(SEtype_buettiker))          se->get_buettiker_selfenergy()->set_energy_parameter(buett);
		if (se->has_self_energy(SEtype_golizadeh_momentum)) se->get_golizadeh_m_selfenergy()->set_energy_parameter(golim);
		if (se->has_self_energy(SEtype_golizadeh_phase))    se->get_golizadeh_p_selfenergy()->set_energy_parameter(golip);
		if (se->has_self_energy(SEtype_optical_phonon))     se->get_optical_phonon_selfenergy()->set_scaling(decrease);
		if (se->has_self_energy(SEtype_acoustic_phonon))    se->get_acoustic_phonon_selfenergy()->set_scaling(decrease);
		if (se->has_self_energy(SEtype_spont_photon))       se->get_spontaneous_photon_selfenergy()->set_scaling(decrease);
		if (se->has_self_energy(SEtype_ion_imp))       		se->get_ionized_impurities_selfenergy()->set_scaling(decrease);
		
		while (true) 
		{
			double multiplier = min(ramp_factor, 1.0/decrease);
			NEGF_ASSERT(multiplier>= 1.0, "scattering decrease factor must be >=1.");
			decrease = decrease * multiplier;
			logmsg->emit_huge_header("Outer iteration with scattering reduced by %.3g",decrease);
			if (se->has_self_energy(SEtype_buettiker)) {
				se->get_buettiker_selfenergy()->set_energy_parameter(se->get_buettiker_selfenergy()->get_energy_parameter() * multiplier);
			}
			if (se->has_self_energy(SEtype_golizadeh_momentum)) {
				se->get_golizadeh_m_selfenergy()->set_energy_parameter(se->get_golizadeh_m_selfenergy()->get_energy_parameter() * multiplier);
			}
			if (se->has_self_energy(SEtype_golizadeh_phase)) {
				se->get_golizadeh_p_selfenergy()->set_energy_parameter(se->get_golizadeh_p_selfenergy()->get_energy_parameter() * multiplier);
			}
			if (se->has_self_energy(SEtype_optical_phonon)) {
				se->get_optical_phonon_selfenergy()->set_scaling(decrease);
			}
			if (se->has_self_energy(SEtype_acoustic_phonon)) {
				se->get_acoustic_phonon_selfenergy()->set_scaling(decrease);
			}
			if (se->has_self_energy(SEtype_spont_photon)) {
				se->get_spontaneous_photon_selfenergy()->set_scaling(decrease);
			}
			if (se->has_self_energy(SEtype_ion_imp)) {
				se->get_ionized_impurities_selfenergy()->set_scaling(decrease);
			}
			this->iterate(iter, constants::inner_ramp_errcrit, false);
			if (decrease - 1.0 > -1e-5) break;
		}
	}
	this->ramping_finished = true;
);}



void OuterLoop::iterate(uint outer_iter, const double inner_err_crit, bool divergence_avoiding_only)
{STACK_TRACE(
	// determine new appropriate energy grid
	res->determine_new_energy_grid(divergence_avoiding_only, this->fermilevels);
		
	// interpolate Green functions onto new grid
	gf->interpolate_from_old_energies();
	
	// if (outer_iter>1) this->compute_fermilevels();
	
	// do inner loop Keldysh-Dyson until convergence is reached
	bool inner_loop_converged = this->innerloop->perform(outer_iter, inner_err_crit);
	if (!inner_loop_converged) 
	{
		// output situation at the end
		if (se->has_self_energy(SEtype_spont_photon)) {
			pp->compute_scattering_current(SEtype_spont_photon);
			pp->compute_luminescence();
		}
		if (mpi->get_rank()==constants::mpi_master_rank) {
			char buf[1000];
			sprintf(buf,"%s_NOTCONVERGED",outfilename.c_str());
			this->save_everything_to_file(buf);

			pp->get_contact_current()->snapshot(888.888);
			pp->get_contact_current()->save_to_file(buf);
			if (se->has_self_energy(SEtype_spont_photon)) {
				pp->get_luminescence2()->snapshot(888.888);
				pp->get_luminescence2()->write_power_to_file(buf);
			}
		}

		// abort
		mpi->synchronize_processes();
		NEGF_EXCEPTION("Inner loop did not converge - terminating.");
	}
	
	// compute spectral density and density (the latter being stored in the master MPI process)
	this->pp->compute_local_dos();
	this->pp->compute_spectral_edensity();
	this->pp->compute_edensity();
	this->pp->compute_spectral_hdensity();
	this->pp->compute_hdensity();
	
	logmsg->emit_header("Solving Poisson equation");
	int root = constants::mpi_master_rank;
	if (mpi->get_rank()==root) 
	{		
		// store the old electrostatic potential
		this->last_electrostatic_potential = this->poisson->get_electrostatic_potential();
		
		// assign new density to poisson
		// to do so, we must enlarge the vector for the non-internal vertices
		// the value at these new vertices does not matter since the poisson eqn is replaced with Dirichlet or Neumann contitions there
		vector<double> edens = this->pp->get_edensity();
		vector<double> hdens = this->pp->get_hdensity();
		NEGF_ASSERT(edens.size()==this->xspace->get_num_internal_vertices(), "unexpected density vector size.");
		NEGF_ASSERT(hdens.size()==this->xspace->get_num_internal_vertices(), "unexpected density vector size.");
		vector<double> edens2; edens2.resize(xspace->get_num_vertices(), 0.0);
		vector<double> hdens2; hdens2.resize(xspace->get_num_vertices(), 0.0);
		for (uint ii=0; ii < xspace->get_num_internal_vertices(); ii++) {
			edens2[xspace->get_global_vertex_index(ii)] = edens[ii];
			hdens2[xspace->get_global_vertex_index(ii)] = hdens[ii];
		}
		// fill density of vertices inside contacts
		for (uint cc = 0; cc < xspace->get_num_contacts(); cc++) {
			const vector<Vertex *> contact_verts = xspace->get_contact(cc)->get_contact_vertices();
			double contact_edens = 0.0;
			double contact_hdens = 0.0;
			for (uint ii = 0; ii < contact_verts.size(); ii++) 
			{
				const vector<Edge *> edges_near_ii = xspace->get_edges_near(contact_verts[ii]);
				for (uint jj=0; jj<edges_near_ii.size(); jj++) {
					Vertex * v = (edges_near_ii[jj]->get_lower_vertex()==contact_verts[ii]) 
									? edges_near_ii[jj]->get_upper_vertex() : edges_near_ii[jj]->get_lower_vertex();
					if (!v->is_at_contact()) {
						uint idx = v->get_index_global();
						if (contact_edens!=0.0) {
							NEGF_ASSERT(contact_edens==edens2[idx], "density must be the same at all vertices within a contact");
						} else {
							contact_edens = edens2[idx];
							contact_hdens = hdens2[idx];
						}
					}
				}
			}
			logmsg->emit(LOG_INFO,"edensity at contact %d: %.3e", cc, contact_edens);
			logmsg->emit(LOG_INFO,"hdensity at contact %d: %.3e", cc, contact_hdens);
			NEGF_ASSERT(contact_edens>0.0 || contact_hdens>0.0, "densities are totally zero.");
			for (uint ii = 0; ii < contact_verts.size(); ii++) {
				uint idx = contact_verts[ii]->get_index_global();
				edens2[idx] = contact_edens;
				hdens2[idx] = contact_hdens;
			}
		}
		
		// check for irregular stuff
		const Matd & sedens = this->pp->get_entire_spectral_edensity();
		double nE_max = 0.0;
		for (uint ii=1; ii<=sedens.num_rows(); ii++) {
			nE_max = max(nE_max, sedens(ii,1));
		}
		logmsg->emit(LOG_INFO,"nE_max at internal vertex 1: %g",nE_max);
		if (nE_max > 2.0) {
			//logmsg->emit(LOG_INFO,"****** nE_max TOO BIG - DISCARDING!");
			logmsg->emit(LOG_INFO,"****** BAD nE!!! *****");
		}
		
		this->poisson->assign_new_edensity(edens2);
		this->poisson->assign_new_hdensity(hdens2);
		
		// solve poisson
		//NewtonSolver * newton = this->poisson->get_newton_solver();
		//newton->calculate_jacobian();
		//newton->calculate_rhs();
		//newton->calculate_update();
		//newton->update_newton_vars(this->poisson->get_poisson_equation()->get_timestamp() + 1);
		this->poisson->solve_one_step();

		// UNDERRELAXATION
		double underrelax = (outer_iter<=3) ? this->potential_underrelaxation : this->late_underrelax;
		if (underrelax>0.0 && outer_iter==1) {
			logmsg->emit(LOG_INFO,"First outer iteration. No underrelaxation is performed.");
		}
		if (underrelax>0.0 && outer_iter>1) 
		{
			bool pulay = false;
			if (outer_iter > 2 && options->exists("PulayMixing") && options->get("PulayMixing")==1) pulay = true;
			if (pulay)
			{
				logmsg->emit(LOG_INFO,"Underrelaxing potential (Pulay mixing with factor %.1g(new), %.1g(old))",
					1.0-underrelax, underrelax);
				uint m = this->poisson_solutions.size();
				NEGF_ASSERT(m>0, "need at leat 1 poisson solution.");
				
				// ----------------------------------------------------------------------------------------------
				// find c_l such that | sum_{l=0}^{m-1} c_l (poisson_solutions[l] - mixed_potential[l]) | = min.
				// ----------------------------------------------------------------------------------------------
				Matd Diff(m+1, m+1);
				for (uint ii=1; ii<=m; ii++) {
					for (uint jj=1; jj<=m; jj++) {
						double tmp = 0.0;
						for (uint xx=0; xx < xspace->get_num_vertices(); xx++) {
							tmp +=   (poisson_solutions[ii-1][xx] - mixed_potentials[ii-1][xx])
							       * (poisson_solutions[jj-1][xx] - mixed_potentials[jj-1][xx]);
						}
						Diff(ii,jj) = 2.0 * tmp;
					}
				}
				for (uint ii=1; ii<=m; ii++) {
					Diff(m,ii) = 1.0;
					Diff(ii,m) = 1.0;
				}
				Vecd rhs(m+1);
				rhs(m+1) = 1.0;
				
				// solve Diff*coeff = rhs
				solve_linear_problem(Diff,rhs);	// rhs is overwritten with solution
								
				vector<double> coeff; coeff.resize(m, 0.0);
				for (uint ii=0; ii<m; ii++) {
					coeff[ii] = rhs(ii+1);
				}
				
				// ----------------------------------------------------------------------------------------------
				// mixed_potentials[m] = sum_{l=0}^{m-1} c_l ((1.0-underrelax)*poisson_solutions[l] + undderrelax*mixed_potential[l])
				// ----------------------------------------------------------------------------------------------
				// please note that the new poisson solution does NOT influence the new mixed potential!
				this->poisson_solutions.push_back(this->poisson->get_electrostatic_potential());
				vector<double> mix; mix.resize(xspace->get_num_vertices(), 0.0);
				for (uint ii=0; ii < xspace->get_num_vertices(); ii++) {
					for (uint mm=0; mm < poisson_solutions.size(); mm++) {
						mix[ii] += coeff[mm] * ((1.0-underrelax) * poisson_solutions[mm][ii]
						                       +     underrelax  * mixed_potentials[mm][ii]);
					}
				}
				this->mixed_potentials.push_back(mix);
			} else
			{
				logmsg->emit(LOG_INFO,"Underrelaxing potential (Kerker mixing: phi = %.1g*phi_new+%.1g*phi_old)",
					1.0-underrelax, underrelax);
				const vector<double> & phi_new = this->poisson->get_electrostatic_potential();
				const vector<double> & phi_old = this->last_electrostatic_potential;
				vector<double> mixed_potential;
				mixed_potential.resize(xspace->get_num_vertices(), 0.0);
				for (uint xx=0; xx < xspace->get_num_vertices(); xx++) {
					mixed_potential[xx] = (1.0-underrelax) * phi_new[xx]
										+ underrelax       * phi_old[xx];
				}
				//this->poisson->get_poisson_equation()->set_values(mixed_potential, this->poisson->get_poisson_equation()->get_timestamp());
				this->poisson->set_electrostatic_potential(mixed_potential);

				this->poisson_solutions.push_back(this->poisson->get_electrostatic_potential());
				this->mixed_potentials.push_back(mixed_potential);
			}
		}
		
		// broadcast the new potential to all other MPI processes
		vector<double> pot = this->poisson->get_electrostatic_potential();
		logmsg->emit(LOG_INFO,"Value of new el.stat. potential at first vertex: %g, at last vertex: %g.", pot[0], pot[pot.size()-1]);
		mpi->broadcast(pot, root);
		
		// assign new electrostatic potential to own Hamiltonian object
		this->innerloop->get_hamiltonian()->set_electrostatic_potential(pot);

		if (options->exists("InjectingStatesCutoff")) {
			vector<double>  Ec = poisson->get_cbedge()/*->get_values()*/;
			vector<double>  Ev = poisson->get_vbedge()/*->get_values()*/;
			const double ec = constants::convert_from_SI(units::charge, constants::SIec);
			NEGF_ASSERT(Ec.size()==pot.size() && Ev.size()==pot.size(), "inconsistent sizes of Ec, Ev, phi");
			for (uint ii=0; ii<pot.size(); ii++) {
				Ec[ii] = Ec[ii] - ec*pot[ii];
				Ev[ii] = Ev[ii] - ec*pot[ii];
			}
			this->innerloop->get_self_energies()->get_contact_selfenergy()->assign_bandedges(Ec, Ev);
		}
	} else {
		// receive the freshly computed electrostatic potential from the master thread
		vector<double> pot;
		pot.resize(this->poisson->get_grid()->get_num_vertices(), 0.0);	// necessary! vector needs to have correct length
		mpi->broadcast(pot, root);
		
		// assign this potential to the own Hamiltonian object
		this->innerloop->get_hamiltonian()->set_electrostatic_potential(pot);

		if (options->exists("InjectingStatesCutoff")) {
			vector<double>  Ec = poisson->get_cbedge()/*->get_values()*/;
			vector<double>  Ev = poisson->get_vbedge()/*->get_values()*/;
			const double ec = constants::convert_from_SI(units::charge, constants::SIec);
			NEGF_ASSERT(Ec.size()==pot.size() && Ev.size()==pot.size(), "inconsistent sizes of Ec, Ev, phi");
			for (uint ii=0; ii<pot.size(); ii++) {
				Ec[ii] = Ec[ii] - ec*pot[ii];
				Ev[ii] = Ev[ii] - ec*pot[ii];
			}
			this->innerloop->get_self_energies()->get_contact_selfenergy()->assign_bandedges(Ec, Ev);
		}
	}
);}


bool OuterLoop::converged()
{STACK_TRACE(
	if(mpi->get_rank()!=constants::mpi_master_rank) {
		NEGF_EXCEPTION("OuterLoop::converged() should only be called from master thread");
	}
	
	uint Nvert = this->last_electrostatic_potential.size();
	
	// look at change in electrostatic potential
	vector<double> change;
	change.resize(Nvert, 0.0);
	const vector<double> & new_elstat_pot = this->poisson->get_electrostatic_potential();
	NEGF_ASSERT(new_elstat_pot.size()==Nvert, "something is wrong.");
	for (uint ii=0; ii < Nvert; ii++) {
		change[ii] = new_elstat_pot[ii] - this->last_electrostatic_potential[ii];
	}
	double vector_norm = 0.0;
	for (uint ii=0; ii < Nvert; ii++) {
		vector_norm += change[ii]*change[ii];
	}
	vector_norm = negf_math::sqrt(vector_norm);
	vector_norm = vector_norm / Nvert;
	if (vector_norm < this->rel_error) {
		logmsg->emit(LOG_INFO,"Outer loop converged! |phi_new - phi_old|=%.3e < %.3e", vector_norm, this->rel_error);
		return true;
	} else {
		logmsg->emit(LOG_INFO,"Outer loop did not converge: |phi_new - phi_old|=%.3e >= %.3e", vector_norm, this->rel_error);
		return false;
	}
);}


void OuterLoop::update_voltages(const vector<double> & voltages_) throw (Exception *)
{STACK_TRACE(
	logmsg->emit_header("updating voltages");
	NEGF_ASSERT(voltages_.size()==xspace->get_num_contacts(), "expected 1 voltage per contact.");
	this->voltages = voltages_;
	
	this->compute_fermilevels();
	
	// start over with Poisson mixing schemes
	this->poisson_solutions.clear();
	this->mixed_potentials.clear();
);}


void OuterLoop::update_underrelaxation(double new_underrelaxation) throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(new_underrelaxation>=0.0-1e-10 && new_underrelaxation<=1.0+1e-10, "new underrelaxation must be between 0 and 1.");
	this->potential_underrelaxation = new_underrelaxation;
	logmsg->emit(LOG_INFO,"new electrostatic potential underrelaxation parameter: %.3g", this->potential_underrelaxation);
);}

void OuterLoop::update_late_underrelax(double new_late_underrelax) throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(new_late_underrelax>=0.0-1e-10 && new_late_underrelax<=1.0+1e-10, "new late underrelaxation must be between 0 and 1.");
	this->late_underrelax = new_late_underrelax;
	logmsg->emit(LOG_INFO,"new electrostatic potential LATE underrelaxation: %.3g", this->late_underrelax);
);}

void OuterLoop::compute_fermilevels()
{STACK_TRACE(
	logmsg->emit_header("updating fermilevels");
	const Geometry * xgrid = this->poisson->get_grid();
	
	// get doping of contact 0 vertex 0
	Vertex * v = xgrid->get_contact(0)->get_contact_vertex(0);
	double doping = this->poisson->get_doping()[v->get_index_global()]/*->get_value(v->get_index_global())*/;
	
	// get fermilevel of contact 0 vertex 0
	// will be the same for all voltage settings
	this->fermilevels.clear(); fermilevels.resize(voltages.size());
	this->innerloop->get_contact_0_fermilevel()->calculate_contact0_dos_from_device_GR();
	this->innerloop->get_contact_0_fermilevel()->compute_contact_0_fermilevel(doping);
	fermilevels[0] = this->innerloop->get_contact_0_fermilevel()->get_contact_0_fermilevel();
	
	// all other fermilevels are then determined by the voltage differences
	for (uint ii=1; ii < voltages.size(); ii++) {
		fermilevels[ii] = fermilevels[0] + (voltages[ii] - voltages[0]);
	}
	
	if (this->options->exists("NonzeroBoundaryField")
		&& this->options->get("NonzeroBoundaryField")==1
		&& fermilevels.size()==2) {
		Vertex * v2 = xgrid->get_contact(1)->get_contact_vertex(0);
		double dist = xgrid->get_distance(v,v2);
		double E_field = -(fermilevels[1] - fermilevels[0]) / dist;
		logmsg->emit(LOG_INFO,"SETTING E-FIELD FOR NEUMANN CONDITION Of POISSON EQN TO %e",E_field);
		//this->poisson->get_poisson_equation()->set_neumann_electric_field(E_field);
		NEGF_EXCEPTION("Need to change this.");
	}
	
	// assign!
	for (uint cc=0; cc < xgrid->get_num_contacts(); cc++) {
		xgrid->get_contact(cc)->set_bndcond(quantities::fermilevel, bndconds::BC_Dirichlet);
		xgrid->get_contact(cc)->set_bc_num_values(quantities::fermilevel, 1);
		xgrid->get_contact(cc)->set_bnd_value(quantities::fermilevel, fermilevels[cc]);
	}
);}


/** called by master process only */
void OuterLoop::save_everything_to_file(const char * filebase) throw (Exception *)
{STACK_TRACE(
#ifndef NODFISE
    string old_output_filename = this->outfilename;
    this->poisson->get_output_data()->set_filename(filebase);
	this->poisson->get_output_data()->write_dat();
    this->poisson->get_output_data()->set_filename(old_output_filename);
#endif
	
	const vector<double> & phi   = poisson->get_electrostatic_potential();
	const vector<double> & edens = poisson->get_edensity()/*->get_values()*/;
	const vector<double> & hdens = poisson->get_hdensity()/*->get_values()*/;
	vector<double>         Ec    = poisson->get_cbedge()/*->get_values()*/;
	vector<double>         Ev    = poisson->get_vbedge()/*->get_values()*/;
	const double ec = constants::convert_from_SI(units::charge, constants::SIec);
	for (uint ii=0; ii < phi.size(); ii++) {
		Ec[ii] += -ec*phi[ii];
		Ev[ii] += -ec*phi[ii];
	}
	
	//const vector<double> & curr   = poisson->get_electrostatic_potential();
	const vector<double> & ecurr = pp->get_ecurrent(); 
	const vector<double> & hcurr = pp->get_hcurrent(); 
	vector<double> curr(phi.size(), 0.0);
	for (uint ii=0; ii<ecurr.size(); ii++) {
		if (ii>curr.size()-1) break;
		curr[ii] = ecurr[ii] + hcurr[ii]; // +? -?
	}
	
	
	logmsg->emit_noendl(LOG_INFO,"potential=");
	for (uint ii=0; ii < phi.size(); ii++) {
		logmsg->emit_noendl(LOG_INFO,"%e   ",phi[ii]);
	}
	logmsg->emit(LOG_INFO,"");
	logmsg->emit(LOG_INFO,"potential at first inner vertex: %.6g, at last inner vertex: %.6g",
				phi[xspace->get_global_vertex_index(0)], phi[xspace->get_global_vertex_index(xspace->get_num_internal_vertices()-1)]);
	
	InputParser parser;
	char buf[1000];
	sprintf(buf,"%s_phi_n_p_J_Ec_Ev",filebase);
	parser.write_phi_n_p_Ec_Ev_J(buf, xspace, phi, edens, hdens, curr, Ec, Ev);
	sprintf(buf,"%s_JEe",filebase);
	parser.write_current_matrix(buf, pp->get_entire_spectral_ecurrent(), xspace, energies->get_energy_grid());
	sprintf(buf,"%s_JEe2",filebase);
	parser.write_xE_matrix(buf, pp->get_entire_spectral_ecurrent2(), xspace, energies->get_energy_grid());
	sprintf(buf,"%s_JEh",filebase);
	parser.write_current_matrix(buf, pp->get_entire_spectral_hcurrent(), xspace, energies->get_energy_grid());
	sprintf(buf,"%s_JEh2",filebase);
	parser.write_xE_matrix(buf, pp->get_entire_spectral_hcurrent2(), xspace, energies->get_energy_grid());
	// make sure pp->compute_scattering_current(SEtype_contact) was called by all processes
	sprintf(buf,"%s_Jecont",filebase);
	parser.write_xE_matrix(buf, pp->get_entire_scattering_ecurrent(SEtype_contact), xspace, energies->get_energy_grid());
	sprintf(buf,"%s_Jhcont",filebase);
	parser.write_xE_matrix(buf, pp->get_entire_scattering_hcurrent(SEtype_contact), xspace, energies->get_energy_grid());
	if (se->has_self_energy(SEtype_optical_phonon)) {
		// make sure pp->compute_scattering_current(SEtype_optical_phonon) was called by all processes
		sprintf(buf,"%s_Jephon",filebase);
		parser.write_xE_matrix(buf, pp->get_entire_scattering_ecurrent(SEtype_optical_phonon), xspace, energies->get_energy_grid());
		sprintf(buf,"%s_Jhphon",filebase);
		parser.write_xE_matrix(buf, pp->get_entire_scattering_hcurrent(SEtype_optical_phonon), xspace, energies->get_energy_grid());
	}
	if (se->has_self_energy(SEtype_spont_photon)) {
		// make sure pp->compute_scattering_current(SEtype_optical_photon) was called by all processes
		sprintf(buf,"%s_Jephot",filebase);
		parser.write_xE_matrix(buf, pp->get_entire_scattering_ecurrent(SEtype_spont_photon), xspace, energies->get_energy_grid());
		sprintf(buf,"%s_Jhphot",filebase);
		parser.write_xE_matrix(buf, pp->get_entire_scattering_hcurrent(SEtype_spont_photon), xspace, energies->get_energy_grid());
		
		// make sure pp->compute_luminescence() was called by all processes - otherwise an old version is saved
		sprintf(buf,"%s_RphotLAK",filebase);
		pp->get_luminescence()->write_recombination_to_file(buf);
		sprintf(buf,"%s_spectrumLAK",filebase);
		pp->get_luminescence()->write_spectrum_to_file(buf);
		
		sprintf(buf,"%s_RphotGAL",filebase);
		pp->get_luminescence2()->write_recombination_to_file(buf);
		sprintf(buf,"%s_spectrumGAL",filebase);
		pp->get_luminescence2()->write_spectrum_to_file(buf);
	}
	sprintf(buf,"%s_nE",filebase);
	parser.write_xE_matrix(buf, pp->get_entire_spectral_edensity(), xspace, energies->get_energy_grid());
	sprintf(buf,"%s_pE",filebase);
	parser.write_xE_matrix(buf, pp->get_entire_spectral_hdensity(), xspace, energies->get_energy_grid());
	sprintf(buf,"%s_LDOS",filebase);
	parser.write_xE_matrix(buf, pp->get_entire_local_dos(), xspace, energies->get_energy_grid());
	sprintf(buf,"%s_LDOS_k0",filebase);
	parser.write_xE_matrix(buf, pp->get_entire_local_dos_k0(), xspace, energies->get_energy_grid());
	if (fabs(Nn-2.0) < 1e-10) {
		sprintf(buf,"%s_LDOS_VB",filebase);
		parser.write_xE_matrix(buf, pp->get_entire_local_dos_VB(), xspace, energies->get_energy_grid());
	}
	
	if (options->exists("Transmission") && options->get("Transmission")==1) {
		sprintf(buf,"%s_TE",filebase);
		pp->get_transmission()->write_to_file(buf);
	}

    if (options->exists("QuasiFermilevels") && options->get("QuasiFermilevels")==1) {
        sprintf(buf,"%s_QFL",filebase);
        pp->get_quasi_fermilevel()->write_to_file(buf);
    }
);}

