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
#include "PoissonProblem.h"
using namespace negf;

PoissonProblem::PoissonProblem(const Geometry * grid_, const MaterialDatabase * db_, const Options * options) throw (Exception *):
	grid(grid_),
	db(db_)
{STACK_TRACE(
	logmsg->emit_header("setting up Poisson problem");
	NEGF_ASSERT(grid!=NULL && db!=NULL, "null pointer encountered.");	
	
	// create Box Method
	this->bm = new BoxMethod(grid);
	
	// ---------------------------------------
	// create equation objects
	// ---------------------------------------
	
	// create temperature, eff. mass, eff. DOS, bandedges
	vector<double> vals; vals.clear();	vals.push_back(options->get("temperature"));
	this->temperature = new ArbitraryData(quantities::temperature, vals);
	this->emass       = new EffectiveMass(grid, bm, db, quantities::electron_density, "vertex");
	this->hmass       = new EffectiveMass(grid, bm, db, quantities::hole_density, "vertex");
	this->Nc          = new EffectiveDOS(temperature, emass, 3/*dos dim*/,	quantities::electron_density);
	this->Nv          = new EffectiveDOS(temperature, hmass, 3/*dos dim*/,	quantities::hole_density);
	this->cbedge      = new Bandedge(grid, bm, db, quantities::electron_density, "vertex", temperature->get_value(0));
	this->vbedge      = new Bandedge(grid, bm, db, quantities::hole_density, "vertex", temperature->get_value(0));
	
	// create dielectric constant
	this->epsilon = new EpsilonRegionwise(grid, db);
	
	// create doping
	// assume doping to be stored in SI units (m^-3)
	string doping_datfile = fnames->get_filename();

	// OLD - read in from DF-ISE
	// this->doping = new Doping(grid, bm, doping_datfile.c_str(), 1, 1, 0, 3);
	//		options_DopingUnit, options_DopingSign, options_UseGrainBoundaryDoping, options_dim);

	// NEW - regionwise constant doping from .cmd-file
	this->doping = new Doping(grid, bm, options->get_cmdfile());


	//this->doping = new ConstantDensity(quantities::hole_density, grid, 0.0);
	double doping_norm = 0.0;
	for (uint ii=0; ii<doping->get_num_variables(); ii++) {
		doping_norm += doping->get_value(ii)*doping->get_value(ii);
	}
	doping_norm = negf_math::sqrt(doping_norm);
	logmsg->emit(LOG_INFO,"Norm of doping: %e",doping_norm);
	
	// create strain and polarization
	vector<double> static_rhs;
	if (options->exists("StrainPolarization") && options->get("StrainPolarization")==1) {
		double pol_decreaser = 1.0;
		if (options->exists("PolarizationDecreaser")) {
			pol_decreaser = options->get("PolarizationDecreaser");
		}
		this->strainpol = new StrainPolarization(grid, db, pol_decreaser);
		
		// prepare sheet densities entering FEM-Poisson equation
		static_rhs.resize(grid->get_num_vertices(), 0.0);
		for (uint ii=0; ii<grid->get_num_vertices(); ii++) {
			static_rhs[ii] = strainpol->get_sheet_charge(grid->get_vertex(ii));
		}
	} else {
		this->strainpol = 0;
		//this->poisson  = new PoissonBM(grid, bm);
	}
	
	logmsg->emit_small_header("Setting up Poisson equation and density equations");
	string gridfile = fnames->get_filename() + ".grd";	// required because tdkp has its own Geometry object
	//string polarization_file = "";
	//this->poisson  = new PoissonFEM(grid, *db, gridfile.c_str(), /*polarization_file.c_str()*/static_rhs);
	this->poisson  = new PoissonFEM1D(grid, *db, gridfile.c_str(), /*polarization_file.c_str()*/static_rhs);
	
	// create edensity, hdensity equations
	vals.assign(grid->get_num_vertices(), 0.0);
	//this->hdensity = new ConstantDensity(quantities::hole_density, grid, 0.0);
	//this->edensity = new ArbitraryData(quantities::electron_density, vals);
	//this->hdensity = new ArbitraryData(quantities::hole_density, vals);
	this->edensity = new SemiclassicalDensity(this->poisson, cbedge, Nc, temperature);
	this->hdensity = new SemiclassicalDensity(this->poisson, vbedge, Nv, temperature);
	this->edensity->set_values(vals, 0);
	this->hdensity->set_values(vals, 0);
	
	// allocate densities to Poisson eqn
	this->poisson->allocate_owncharge(edensity, hdensity);
	this->poisson->allocate_doping(doping);
	this->poisson->allocate_epsilon(epsilon);
		
	this->epsilon    ->set_name("epsilon");
	this->doping     ->set_name("doping");
	this->edensity   ->set_name("edensity");
	this->hdensity   ->set_name("hdensity");
	this->poisson    ->set_name("potential");
	this->temperature->set_name("temperature");
	this->emass      ->set_name("emass");
	this->hmass      ->set_name("hmass");
	this->Nc         ->set_name("electron_effective_dos");
	this->Nv         ->set_name("hole_effective_dos");
	this->cbedge     ->set_name("cbedge");
	this->vbedge     ->set_name("vbedge");
	
	// ----------------------------------------------------------
	// assign appropriate boundary conditions to contacts!
	// ----------------------------------------------------------
	NEGF_ASSERT(grid->get_num_contacts()>0, "must have at least one contact!");
	/*
	logmsg->emit(LOG_INFO,"Assigning Dirichlet (=0) condition for potential to contact %s, Neumann elsewhere.",
					grid->get_contact(0)->get_name().c_str());
	grid->get_contact(0)->set_bndcond(quantities::potential, bndconds::BC_Dirichlet);
	grid->get_contact(0)->set_bc_num_values(quantities::potential, grid->get_contact(0)->get_num_contact_vertices());
	for (uint ii=0; ii < grid->get_contact(0)->get_num_contact_vertices(); ii++) {
		grid->get_contact(0)->set_bnd_value(quantities::potential, ii, 0.0);
	}
	for (uint ii=1; ii < grid->get_num_contacts(); ii++) {
		grid->get_contact(ii)->set_bndcond(quantities::potential, bndconds::BC_Neumann);
	}*/
	
	logmsg->emit(LOG_INFO,"Assigning Neumann condition for potential at every contact.");
	for (uint ii=0; ii < grid->get_num_contacts(); ii++) {
		grid->get_contact(ii)->set_bndcond(quantities::potential, bndconds::BC_Neumann);
	}
	
	// ---------------------------
	// initial guess
	// ---------------------------
	logmsg->emit(LOG_INFO,"Initializing potential values to zero.");
	vals.clear(); vals.resize(grid->get_num_vertices(), 0.0);
	
	// note: fermilevel at contact 0 is sought with potential=0
	const double ec = constants::convert_from_SI(units::charge, constants::SIec);
	const double kT = constants::convert_from_SI(units::energy, constants::SIkb * temperature->get_value(0));
	
	/*
	uint v_potential_zero = grid->get_contact(0)->get_contact_vertex(0)->get_index_global();
	int sign_zero = negf_math::sign(doping->get_value(v_potential_zero));
	double zero_level = sign_zero * kT/ec * negf_math::fermihalf_inverse(fabs(doping->get_value(v_potential_zero)) / Nc->get_value(v_potential_zero));
	if (mpi->get_rank()==constants::mpi_master_rank) cout << "Zero level at " << zero_level << endl;
	for (uint ii=0; ii < grid->get_num_vertices(); ii++) {
		vals[ii] = negf_math::sign(doping->get_value(ii)) * kT/ec * negf_math::fermihalf_inverse(fabs(doping->get_value(ii)) / Nc->get_value(ii)) - zero_level; 
		if (mpi->get_rank()==constants::mpi_master_rank) cout << vals[ii] << "   ";
	}*/
	
	/*
	for (uint ii=0; ii < grid->get_num_vertices(); ii++) {
		vals[ii] = cbedge->get_value(ii) + kT/ec * negf_math::log(doping->get_value(ii) / Nc->get_value(ii)); 
	}
	*/
	
	double doping1 = doping->get_value(grid->get_contact(0)->get_contact_vertex(0)->get_index_global());
	double doping2 = doping->get_value(grid->get_contact(1)->get_contact_vertex(0)->get_index_global());
	if (doping1*doping2 < 0.0) {
		logmsg->emit(LOG_INFO,"p-n situation encountered. setting potential to a linear slope:");
		
		// p-n situation
		double Ncv1 = (doping1 < 0.0) ? Nv->get_value(grid->get_contact(0)->get_contact_vertex(0)->get_index_global())
									  : Nc->get_value(grid->get_contact(0)->get_contact_vertex(0)->get_index_global());
		double Ncv2 = (doping2 < 0.0) ? Nv->get_value(grid->get_contact(1)->get_contact_vertex(0)->get_index_global())
									  : Nc->get_value(grid->get_contact(1)->get_contact_vertex(0)->get_index_global());
		
		double shift1 = kT/ec * negf_math::fermihalf_inverse(fabs(doping1) / Ncv1);
		double shift2 = kT/ec * negf_math::fermihalf_inverse(fabs(doping2) / Ncv2);
		double cbedge_eV = TdkpInfoDesk::get_cbedge(grid->get_contact(0)->get_adjacent_region()->get_material(), temperature->get_value(0), db);
		double vbedge_eV = grid->get_contact(0)->get_adjacent_region()->get_material()->get("valence_band_edge");
		double Egap = constants::convert_from_SI(units::energy, constants::SIec * (cbedge_eV - vbedge_eV));
		double delta_phi = Egap + shift1 - shift2; // signs?
		if (doping1>0.0) {
			// contact 0 is n, contact 1 is p
			delta_phi *= -1.0;
		} else {
			// contact 0 is p, contact 1 is n
		}
		
		/*cout << "doping1=" << doping1 << ", doping2=" << doping2 << endl;
		cout << "Ncv1="    << Ncv1    << ", Ncv2="    << Ncv2    << endl;
		cout << "shift1="  << shift1  << ", shift2="  << shift2  << endl;
		cout << "Ec_eV="   << cbedge_eV << ", Ev_eV=" << vbedge_eV << endl;
		cout << "Egap="    << Egap    << ", dphi="    << delta_phi  << endl;*/
		
		double dtot = grid->get_distance(grid->get_contact(1)->get_contact_vertex(0), grid->get_contact(0));
		for (uint ii=0; ii < grid->get_num_vertices(); ii++) {
			vals[ii] = delta_phi * grid->get_distance(grid->get_vertex(ii), grid->get_contact(0)) / dtot;
		}
		
		// so far, phi=0 at contact 0, phi=delta_phi at contact 1
		// however, we need to shift the potential in case doping1>0 (WHY???)
		if (doping1>0.0) {
			for (uint ii=0; ii < grid->get_num_vertices(); ii++) {
				vals[ii] += Egap/2;
			}
		}
		
		if (mpi->get_rank()==constants::mpi_master_rank) {
			for (uint ii=0; ii < grid->get_num_vertices(); ii++) {
				cout << vals[ii] << "   ";
			}
		}
	}
	
	this->poisson->set_values(vals, poisson->get_timestamp());
	
	// ---------------------------
	// create newton solver
	// ---------------------------
	vector<Equation *> handler;
	handler.push_back(temperature);
	handler.push_back(emass);
	handler.push_back(hmass);
	handler.push_back(cbedge);
	handler.push_back(vbedge);
	handler.push_back(Nc);
	handler.push_back(Nv);
	handler.push_back(epsilon);
	handler.push_back(doping);
	handler.push_back(edensity);
	handler.push_back(hdensity);
	handler.push_back(poisson);
	this->newton = new NewtonSolver(handler, 5,   1e-8,   1e-8,    false,    "poisson", true, false);
	// arguments:						max_iters, rel_err, abs_err, bank_rose, name,  verbose, debug
	
	// --------------------------
	// prepare output data object
	// --------------------------
	this->outputdata = new OutputData(grid, fnames->get_outfile());
	this->outputdata->add(poisson, units::potential);
	this->outputdata->add(edensity, units::density_3d);
	this->outputdata->add(hdensity, units::density_3d);
	this->outputdata->add(doping, units::density_3d);
	this->outputdata->add(cbedge, units::potential); // want output in eV
	this->outputdata->add(vbedge, units::potential); // want output in eV
);}


vector<Equation *> PoissonProblem::get_all_equations() const
{STACK_TRACE(
	vector<Equation *> handler;
	handler.push_back(temperature);
	handler.push_back(emass);
	handler.push_back(hmass);
	handler.push_back(cbedge);
	handler.push_back(vbedge);
	handler.push_back(epsilon);
	handler.push_back(doping);
	handler.push_back(edensity);
	handler.push_back(hdensity);
	handler.push_back(poisson);
	return handler;
);}

void PoissonProblem::assign_new_edensity(const vector<double> & new_edensity)
{STACK_TRACE(
	NEGF_FASSERT(new_edensity.size()==grid->get_num_vertices(), "inconsistent vector size (%d instead of %d)",
				new_edensity.size(),grid->get_num_vertices());
	this->edensity->set_values(new_edensity, edensity->get_timestamp());
);}

void PoissonProblem::assign_new_hdensity(const vector<double> & new_hdensity)
{STACK_TRACE(
	NEGF_FASSERT(new_hdensity.size()==grid->get_num_vertices(), "inconsistent vector size (%d instead of %d)",
				new_hdensity.size(),grid->get_num_vertices());
	this->hdensity->set_values(new_hdensity, hdensity->get_timestamp());
);}


vector<double> PoissonProblem::get_electrostatic_potential() const
{STACK_TRACE(
	vector<double> vals = this->poisson->get_values();
	// in the beginning potential was not yet computed
	// --> take zeros
	if (vals.size()!=this->poisson->get_num_variables()) {
		logmsg->emit(LOG_WARN,"Warning: potential did not seem to be computed.");
		vals.resize(this->poisson->get_num_variables());
	}
	return vals;
);}


