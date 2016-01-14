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
#ifndef POISSONPROBLEM_H_NEGF
#define POISSONPROBLEM_H_NEGF

#include "all.h"

#include "Geometry.h"
#include "BoxMethod.h"
#include "MaterialDatabase.h"
#include "ImplicitEquation.h"
#include "Poisson.h"
#include "PoissonBM.h"
#include "PoissonFEM.h"
#include "PoissonFEM1D.h"
#include "ArbitraryData.h"
#include "ConstantDensity.h"
#include "Doping.h"
#include "EpsilonRegionwise.h"
#include "EffectiveMass.h"
#include "EffectiveDOS.h"
#include "Bandedge.h"
#include "SemiclassicalDensity.h"
#include "StrainPolarization.h"
#include "NewtonSolver.h"
#include "OutputData.h"
#include "Options.h"
#include "TdkpInfoDesk.h"

using namespace std;

namespace negf {

	/** Wrapper class to create, store & access everything related to the Poisson problem */
	class PoissonProblem 
	{

		public:
			
			PoissonProblem(const Geometry * grid_, const MaterialDatabase * db_,
							const Options * options) throw (Exception *);
			~PoissonProblem() {}

			const Geometry  * const get_grid() 				const { return grid; }
			const BoxMethod * const get_boxmethod() 		const { return bm; }
			Poisson * 				get_poisson_equation()  const { return poisson; }
			const Equation *		get_dielectric()  		const { return epsilon; }
			Equation * 				get_edensity() 			const { return edensity; }
			Equation * 				get_hdensity() 			const { return hdensity; }
			Equation * 				get_cbedge() 			const { return cbedge; }
			Equation * 				get_vbedge() 			const { return vbedge; }
			Equation * 				get_effective_edos()	const { return Nc; }
			Equation * 				get_effective_hdos()	const { return Nv; }
			const Equation * 		get_doping() 			const { return doping; }
			OutputData * 			get_output_data() 		const { return outputdata; }
			
			vector<Equation *> 		get_all_equations() const;
			vector<double> 			get_electrostatic_potential() const;
			NewtonSolver * 			get_newton_solver() 	const { return newton; }
			
			void 					assign_new_edensity(const vector<double> & new_edensity);
			void 					assign_new_hdensity(const vector<double> & new_hdensity);
			
		protected:
			
			const Geometry * grid;
			const BoxMethod * bm;
			const MaterialDatabase * db;
			
			Equation * temperature;
			Equation * emass;
			Equation * hmass;
			Equation * Nc;
			Equation * Nv;
			Equation * cbedge;
			Equation * vbedge;
			Equation * epsilon;
			Equation * edensity;
			Equation * hdensity;
			Equation * doping;
			StrainPolarization * strainpol;
			Poisson  * poisson;
			
			NewtonSolver * newton;
			OutputData * outputdata;
			
	};

}	// end namespace negf

#endif /*POISSONPROBLEM_H_NEGF*/
