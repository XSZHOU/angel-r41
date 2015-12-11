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
#include "SemiclassicalDensity.h"
using namespace negf;

SemiclassicalDensity::SemiclassicalDensity(
		Equation * potential_, 
		Equation * bandedge_,
		Equation * effdos_,
		Equation * temperature_):
	potential(potential_),
	bandedge(bandedge_),
	effdos(effdos_),
	temperature(temperature_)
{STACK_TRACE(
	// security checks
	NEGF_ASSERT(potential!=NULL && bandedge!=NULL && effdos!=NULL && temperature!=NULL, "null pointer encountered.");
	NEGF_ASSERT(potential->get_num_variables()==effdos->get_num_variables()
				&& potential->get_num_variables()==bandedge->get_num_variables()
				&& temperature->get_num_variables()==1, "strange number of vars of a depcy.");
	NEGF_ASSERT(potential->get_type()==quantities::potential && bandedge->get_type()==quantities::energy
				&& (effdos->get_type()==quantities::electron_density || effdos->get_type()==quantities::hole_density)
				&& temperature->get_type()==quantities::temperature, "a depcy had the wrong type.");
	
	// phi_ref = 0.0!
	
	this->dependencies.clear();
	dependencies.push_back(potential);
	dependencies.push_back(bandedge);
	dependencies.push_back(effdos);
	dependencies.push_back(temperature);
	
	// standard eqn stuff
	this->its_type = effdos->get_type();
	this->number_of_variables = potential->get_num_variables();
	this->timestamp = 0;
	this->current_variable_values.clear();
);}


double SemiclassicalDensity::compute_value(uint line) const
{STACK_TRACE(	
	NEGF_EXCEPTION("SemiclassicalDensity::compute_value() is not possible! use set_values() instead.");
);}


void SemiclassicalDensity::direct_derivatives(uint line, const Equation * eqn, uint & nonzeros, 
										uint indices [], double coeffs []) const
{STACK_TRACE(	
	if (eqn==this->temperature) {
		nonzeros = 0;
		return;
	}
	if (eqn==this->bandedge) {
		nonzeros = 0;
		return;
	}
	if (eqn==this->effdos) {
		nonzeros = 0;
		return;
	}
	if (eqn==this->potential)
	{
		nonzeros = 1;
		indices[0] = line;
		if (coeffs==NULL) return;
		
		/** calculate 
		 *  EFn = Ec - e*phi - kT * F_{+0.5}^-1(n/Nc)   or
		 *  EFp = Ev - e*phi + kT * F_{+0.5}^-1(p/Nv) 
		 */
		//int    sign = (effdos->get_type()==quantities::electron_density) ? 1 : -1;
		double ec   = constants::convert_from_SI(units::charge, constants::SIec);
		double kT   = constants::convert_from_SI(units::energy, constants::SIkb * temperature->get_value(0));
		//double phi  = potential->get_value(line);
		double dens = this->get_value(line);
		double Ncv  = effdos->get_value(line);
		//double Ex   = bandedge->get_value(line);
		//double EF   = Ex - ec*phi- - sign * kT * negf_math::fermi_int(0.5, dens / Ncv);
			
		/** calculate
		 *  dn/dphi ~ ec/kT * Ncv * F_{-0.5}(nu_x)   where Ncv = Nc or Nv and 
		 *  nu_n = (EFn - Ec + e*phi) / kT,
		 *  nu_p = (Ev - e*phi - EFp) / kT = -(EFp - Ev + e*phi) / kT
		 * 
		 *  nu_x is calculated from nu_x = F_{+0.5}^{-1}(dens/Ncv)
		 */
		double nu_x = negf_math::fermihalf_inverse(dens/Ncv);
		int sgn = (this->get_type()==quantities::electron_density) ? +1 : -1;
		//int sgn = +1;
		coeffs[0] = sgn * ec / kT * Ncv * negf_math::fermi_int(-0.5, nu_x);
		
		return;
	}
	NEGF_EXCEPTION("Unknown equation.");
);}

