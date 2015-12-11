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
#include "EffectiveDOS.h"
using namespace negf;

EffectiveDOS::EffectiveDOS(Equation * temperature_, Equation * mass_, uint dos_dim_,
						quantities::PhysicalQuantity electron_or_hole):
	dos_dim(dos_dim_),
	temperature(temperature_),
	mass(mass_)
{STACK_TRACE(
	NEGF_ASSERT(temperature!=NULL && mass!=NULL, "null pointer encountered.");
	NEGF_ASSERT(temperature->get_num_variables()==mass->get_num_variables()
				|| temperature->get_num_variables()==1,
				"invalid number of variables.");
	NEGF_ASSERT(dos_dim==1 || dos_dim==2 || dos_dim==3, "Wrong DOS dimensionality.");
	
	// assign standard Equation stuff
	switch (electron_or_hole) {
	case quantities::electron_density: this->its_type = quantities::electron_density; break;
	case quantities::hole_density:     this->its_type = quantities::hole_density;     break;
	default: NEGF_EXCEPTION("wrong type."); break;
	}
	this->number_of_variables = mass->get_num_variables();
	this->dependencies.push_back(temperature);
	this->dependencies.push_back(mass);
	
	this->timestamp = min(temperature->get_timestamp(), mass->get_timestamp());
	this->current_variable_values.clear();
	this->compute_values(this->timestamp);
);}

double EffectiveDOS::compute_value(uint line) const
{STACK_TRACE(
	double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
	double kT = (temperature->get_num_variables()==1) 
				? constants::convert_from_SI(units::energy, constants::SIkb * temperature->get_value(0))
				: constants::convert_from_SI(units::energy, constants::SIkb * temperature->get_value(line));
	double m = mass->get_value(line);
	
	double factor = m*kT / (2*constants::pi*hbar*hbar);
	
	switch (dos_dim) {
		case 3:	return 2 * negf_math::pow(factor, 1.5); break;
//		case 2:	return factor; 			 		  break;
		case 2:	return 2 * factor; 			 	  break;
		case 1:	return std::sqrt(factor);		  break;
		default: NEGF_EXCEPTION("Strange dimensionality."); return 0.0;
	}
);}

void EffectiveDOS::direct_derivatives(uint line, const Equation * eqn, uint & nonzeros, 
										uint indices [], double coeffs[]) const
{STACK_TRACE(
	// get nonzeros, indices
	if (eqn==temperature)
	{
		nonzeros = 1;
		if (temperature->get_num_variables()==1) {
			indices[0] = 0;
		} else {
			indices[0] = line;
		}
	} else if (eqn==mass)
	{
		nonzeros = 1;
		indices[0] = line;
	} else NEGF_EXCEPTION("unknown equation.");
	
	// get coeffs
	if (coeffs==NULL)
		return;
		
	double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
	double kT = (temperature->get_num_variables()==1) 
				? constants::convert_from_SI(units::energy, constants::SIkb * temperature->get_value(0))
				: constants::convert_from_SI(units::energy, constants::SIkb * temperature->get_value(line));
	double m = mass->get_value(line);
	
	double factor = m*kT / (2*constants::pi*hbar*hbar);
	
	if (eqn==temperature)
	{
		switch (dos_dim) {
			case 3:	coeffs[0] = 2 * 1.5*negf_math::pow(factor, 0.5) * m / (2*constants::pi*hbar*hbar); 	break;
//			case 2:	coeffs[0] = m / (2*constants::pi*hbar*hbar); 			 		 			 		break;
			case 2:	coeffs[0] = m / (constants::pi*hbar*hbar); 			 		 			 			break;
			case 1:	coeffs[0] = 0.5*negf_math::pow(factor, -0.5) * m / (2*constants::pi*hbar*hbar);	 	break;
			default: NEGF_EXCEPTION("Strange dimensionality.");
		}
		coeffs[0] = coeffs[0] * constants::convert_from_SI(units::energy, constants::SIkb * 1.0);
		return;
	}
	if (eqn==mass)
	{
		switch (dos_dim) {
			case 3:	coeffs[0] = 2 * 1.5*negf_math::pow(factor, 0.5) * kT / (2*constants::pi*hbar*hbar); 	break;
//			case 2:	coeffs[0] = kT / (2*constants::pi*hbar*hbar); 			 		 			 			break;
			case 2:	coeffs[0] = kT / (constants::pi*hbar*hbar); 			 		 			 			break;
			case 1:	coeffs[0] = 0.5*negf_math::pow(factor, -0.5) * kT / (2*constants::pi*hbar*hbar);	 	break;
			default: NEGF_EXCEPTION("Strange dimensionality.");
		}
		return;
	}
	NEGF_EXCEPTION("unknown equation.");
);}
