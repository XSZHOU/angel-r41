/*
Copyright (c) 2010 Sebastian Steiger, Integrated Systems Laboratory, ETH Zurich.
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
#include "MaterialDatabase.h"
#include "PoissonFEM1D.h"
#include "StrainPolarization.h"
//#include "OutputData.h"
#include "Options.h"
#include "CSRMatrix.h"
#include "LinearSolverUmfpack.h"

using namespace std;

namespace negf {

    /** Wrapper class to create, store & access everything related to the Poisson problem */
    class PoissonProblem
    {

        public:

            PoissonProblem(const Geometry * grid_, const MaterialDatabase * db_,
                            const Options * options) throw (Exception *);
            ~PoissonProblem() {}
			
			// preferably call initial_guess(...) immediately after construction
			void initial_guess() { this->initial_guess(false, 0.0, 0.0); }
			void initial_guess(bool new_method, double EF_0, double EF_1); // more sophisticated initial guess using contact fermilevels

            const Geometry *        get_grid()              const { return grid; }
            PoissonFEM1D *          get_poisson_equation()  const { return poisson; }
            const vector<double> &  get_doping()            const { return doping; }
            const vector<double> &  get_cbedge()            const { return cbedge; }
            const vector<double> &  get_vbedge()            const { return vbedge; }
            const vector<double> &  get_Nc()                const { return Nc; }
            const vector<double> &  get_Nv()                const { return Nv; }
            const vector<double> &  get_dielectric()        const { return epsilon; }
            const vector<double> &  get_edensity()          const { return edensity; }
            const vector<double> &  get_hdensity()          const { return hdensity; }

            StrainPolarization *    get_strain_polarization() const { return strainpol; }

            const vector<double> &  get_electrostatic_potential() const { return elstat_pot; }

            void                    set_electrostatic_potential(const vector<double> & new_potential);

            void                    assign_new_edensity(const vector<double> & new_edensity);
            void                    assign_new_hdensity(const vector<double> & new_hdensity);

            static double get_cbedge(const PropertyContainer<double> * mat, const double & T/*emperature*/,
                                     const MaterialDatabase * const database);

            //void solve();
            void solve_one_step();

        protected:

            const Geometry * grid;
            const MaterialDatabase * db;

            vector<double> emass;
            vector<double> hmass;
            vector<double> cbedge;
            vector<double> vbedge;
            vector<double> Nc;
            vector<double> Nv;
            double temperature;

            vector<double> edensity;
            vector<double> hdensity;
            vector<double> doping;
            vector<double> epsilon;

            vector<double> elstat_pot;

            StrainPolarization * strainpol;
            PoissonFEM1D  * poisson;

            // Newton iteration
            CSRMatrix<double> * jacobian;
            vector<double> rhs;
            vector<double> solution;
            LinearSolverUmfpack * solver;
			
			static const double max_update_V = 0.3; // limit was determined from ballistic pn.cmd simulation, V=0
			

    };

}   // end namespace negf

#endif /*POISSONPROBLEM_H_NEGF*/
