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
#include "Kspace.h"
using namespace negf;

Kspace::Kspace(const Geometry * xspace, const Options * options_, const MaterialDatabase * db) throw (Exception *):
	options(options_),
	my_spacing(equal), // default: equal spacing in k-space (~parabolic spacing in E-space)
	my_rule(trapez)    // default
{STACK_TRACE(
	logmsg->emit_header("setting up k-space");
	
	this->dimension = 2; // hard coded
	this->kmin = 0.0;
	this->kmax = constants::convert_from_SI(units::density_1d, 1e9*options->get("maximal_k_value")); // 1/nm assumed!!!
	
	// check options for different discretization / integration scheme
	if (options->exists("KSpaceDiscretization")) {
	    double num = options->get("KSpaceDiscretization");
	         if (fabs(num-0.0) < 1e-13) my_spacing = equal;
	    else if (fabs(num-1.0) < 1e-13) my_spacing = square_root;
	    else NEGF_EXCEPTION("unknown KSpaceDiscretization value. 0=equal, 1=square_root");
	}
    if (options->exists("KSpaceIntegration")) {
        double num = options->get("KSpaceIntegration");
             if (fabs(num-0.0) < 1e-13) my_rule = trapez;
        else if (fabs(num-1.0) < 1e-13) my_rule = romberg_simpson;
        else if (fabs(num-2.0) < 1e-13) my_rule = three_point;
        else NEGF_EXCEPTION("unknown KSpaceIntegration value. 0=trapezoidal, 1=Romberg-Simpson, 2=Three-Point");
    }

	this->Nk = uint(options->get("num_k_points"));
	if (my_rule==romberg_simpson && Nk%2==0) Nk++;	// need an odd number for romberg

	// find the bulk material, which is characterized by having the biggest volume
	vector<double> material_volumes;
	material_volumes.resize(db->get_num_materials(), 0.0);
	for (uint ii=0; ii<xspace->get_num_regions(); ii++) {
		const vector<Element *> & region_elems = xspace->get_region(ii)->get_elements();
		double region_volume = 0.0;
		for (uint jj=0; jj < region_elems.size(); jj++) {
			region_volume += region_elems[jj]->get_volume();
		}
		
		uint material_id = xspace->get_region(ii)->get_material()->get_id();
		NEGF_ASSERT(material_id < material_volumes.size(), "invalid material id.");
		material_volumes[material_id] += region_volume;
	}
	uint biggest_material_id = 0;
	for (uint ii=1; ii < material_volumes.size(); ii++) {
		if (material_volumes[ii] > material_volumes[biggest_material_id]) {
			biggest_material_id = ii;
		}
	}
	const PropertyContainer<double> * biggest_mat = db->get_material(biggest_material_id);
	logmsg->emit(LOG_INFO,"Material \"%s\" is bulk.",biggest_mat->get_name().c_str());
	double me = constants::convert_from_SI(units::mass, constants::SIm0 * biggest_mat->get("electron_effective_mass"));
	double mh = constants::convert_from_SI(units::mass, constants::SIm0 * biggest_mat->get("hole_effective_mass"));
	double m = 0.0;
	if (options->get("kp_method") > 1.0) {
		m = max(me,mh);
	} else {
		m = me;
	}
	logmsg->emit(LOG_INFO,"Testing with m=%.3gm0.",m/constants::convert_from_SI(units::mass, constants::SIm0));
	
	double dE = constants::convert_from_SI(units::energy, constants::SIec * 0.005); // wanted energy resolution
	const double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
	logmsg->emit(LOG_INFO,"Suggested kmax: %e[nm-1].",sqrt(2 * m * this->Nk * dE) / hbar * constants::convert_from_SI(units::length, 1.0)*1e-9);
	//this->kmax = sqrt(2 * m * this->Nk * dE) / hbar;
		
	DomainPoint * p = 0;
	for (uint kk=0; kk<Nk; kk++)
	{
		double k = 0.0;
		double w = 0.0;
		
		double delta_k, k0, k2, dk, Nk_m_1; // helper variables
		switch (my_spacing)
		{
		case equal:
		    delta_k = (Nk>1) ? kmax/(Nk-1) : 0.0;
			k = kk * delta_k; // equally spaced
			
			// -----------------------
			// compute weight!
			// -----------------------
			switch (my_rule)
			{
			case trapez:
				if (kk==0 || kk==Nk-1) {
					w = 1.0/2.0;
				} else {
					w = 1.0;
				}
				w = 2.0*constants::pi * k * delta_k * w;
				break;
			case romberg_simpson:
				if (kk==0 || kk==Nk-1) {
					w = 1.0/3.0;
				} else if(kk%2==0) {
					w = 2.0/3.0;
				} else {
					w = 4.0/3.0;
				}
				w = 2.0*constants::pi * k * delta_k * w;
				break;
			case three_point:
				if (kk==0 || kk==Nk-1) {
					w = 17.0/48.0;
				} else if(kk==1 || kk==Nk-2) {
					w = 59.0/48.0;
				} else if(kk==2 || kk==Nk-3) {
					w = 43.0/48.0;
				} else if(kk==3 || kk==Nk-4) {
					w = 49.0/48.0;
				} else {
					w = 1.0;
				}
				w = 2.0*constants::pi * k * delta_k * w;
				break;
			default:
			    NEGF_EXCEPTION("unknown integration rule.");
			}
			break;
		case square_root:
		    Nk_m_1 = (Nk>1) ? Nk-1 : 0.0;
			k = sqrt(kk/Nk_m_1) * kmax; // -->E(k)~k^2 will be equally spaced

			k0 = (kk!=0)    ? sqrt((kk-1.0)/Nk_m_1) * kmax : k;
			k2 = (kk!=Nk-1) ? sqrt((kk+1.0)/Nk_m_1) * kmax : k;
			dk = 0.5*(k2-k0);
			w = 2.0*constants::pi * k * dk;
			break;
		default:
		    NEGF_EXCEPTION("unknown spacing.");
		}
		
		// -------------------------------
		// create new point, add to vector
		// -------------------------------
		p = new DomainPoint(k, 0.0);
		p->set_index(kk);
		p->set_weight(w);
		
		this->k_grid.push_back(p);
	}
	this->test_accuracy(m);
);}

void Kspace::test_accuracy(const double & m)
{STACK_TRACE(
	if (mpi->get_rank()==constants::mpi_master_rank)
	{
		// test 2D k-integration accuracy
		// analytical result: F0(E) = m*kT/(2*pi*hbar^2) * log(1+exp((EF-E)/kT));
		// k-integration: sum_k  weight_k * 1.0 ./ (1+exp((E + Ek - EF)/kT)), where Ek = hbar^2 k^2 / 2m
				
		double kT   = constants::convert_from_SI(units::energy, constants::SIkb * options->get("temperature"));
		double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
		uint   nE   = 100;
		double Emin = constants::convert_from_SI(units::energy, -0.2 * constants::SIec);
		double Emax = constants::convert_from_SI(units::energy, 0.5 * constants::SIec);
		double EF   = 0.0;
				
		bool wrong_result = false;
		for (uint ee=0; ee < nE; ee++)
		{
			double E = Emin + (Emax - Emin) / (nE-1) * ee;
			
			// analytical result
			double F0 = m*kT/(2*constants::pi*hbar*hbar) * negf_math::log(1.0 + negf_math::exp((EF-E)/kT));
			
			// numerical result
			double numerical_result = 0.0;
			for (uint kk=0; kk < this->get_number_of_points(); kk++) 
			{
				double k = this->get_point(kk).get_coord_abs();
				double Ek = hbar*hbar * k*k / (2.0 * m);
				
				double wk = this->get_point(kk).get_weight();
				//double k2 = (kk == this->get_number_of_points()-1) ? k : this->get_point(kk+1).get_coord_abs();
				//double k0 = (kk == 0) ? k : this->get_point(kk-1).get_coord_abs();
				//double wk = 2.0*constants::pi * (k2 - k0)/2 * k;
				
				numerical_result += wk * 1.0/(4.0*constants::pi*constants::pi) * 
					1.0 / ( 1.0 + negf_math::exp((E + Ek - EF)/kT) );
			}
			
			// compare!
			if (fabs((numerical_result - F0) / F0) > 0.01)
			{
				wrong_result = true;
				logmsg->emit(LOG_WARN, "E=%.3e[eV]: numerical and analytical results do not agree: F0=%.4e, numerically %.4e, rel_err=%.2g",
					E, F0, numerical_result, fabs((numerical_result - F0) / F0));
			}
		}
		//NEGF_ASSERT(!wrong_result, "Aborting.");
	}
	mpi->synchronize_processes();
);}
