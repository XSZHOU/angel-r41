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
#ifndef POISSONFEM1D_H_NEGF
#define POISSONFEM1D_H_NEGF

#include "all.h"

#include "Geometry.h"
#include "MaterialDatabase.h"

using namespace std;

namespace negf {

    /** 1D Nonlinear Poisson Equation using FEM discretization */
    class PoissonFEM1D {
        public:
            PoissonFEM1D(
                const Geometry * grid_,
                const vector<double> & static_rhs_          //!< contains int Ni(x) rho(x), rho are the polarization sheet charges
            );
            ~PoissonFEM1D() {}

            // ----------------------------------
            // setup
            // ----------------------------------

            void set_elstat_pot(vector<double> * potential);
            void set_edensity(vector<double> * edensity);
            void set_hdensity(vector<double> * hdensity);
            void set_doping(vector<double> * doping);
            void set_epsilon(vector<double> * epsilon);

            // for predictor-derivative
            void set_Nc(vector<double> * Nc);
            void set_Nv(vector<double> * Nv);
            void set_kT(double kT);

            void set_piezo_decreaser(double new_decreaser_value);

            // ----------------------------------
            // access
            // ----------------------------------

            vector<int> & get_icol() { return icol; }
            vector<int> & get_prow() { return prow; }

            double  get_newton_function(uint idx) const; // F_i

            void    get_newton_derivative(uint idx, uint & nonzeros,
                                      uint indices[], double coeff[]) const; // dF_i / dphi_j

        protected:

            const Geometry * grid;
            uint Nvert;

            vector<double> stiffness_matrix;    //!< int epsilon grad Ni * grad Nj, stored in CSR format
            vector<double> mass_matrix;         //!< int Ni * Nj, stored in CSR format
            vector<int> icol;                   //!< CSR pattern for mass and stiffness matrix
            vector<int> prow;                   //!< CSR pattern for mass and stiffness matrix
            vector<double> static_rhs;          //!< static rhs contribution (int Ni * rho_static)

            double piezo_decreaser;             //!< factor that multiplies the static rhs, to decrease sheet densities

            double neumann_field;

            vector<double> * elstat_pot;
            vector<double> * edensity;
            vector<double> * hdensity;
            vector<double> * doping;
            vector<double> * epsilon;
            vector<double> * Nc;
            vector<double> * Nv;
            double           kT;

    };

}   // end of namespace

#endif /*POISSONFEM1D_H_NEGF*/
