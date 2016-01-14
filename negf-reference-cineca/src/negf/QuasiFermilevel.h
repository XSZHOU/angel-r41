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
#ifndef QUASIFERMILEVEL_H_
#define QUASIFERMILEVEL_H_


#include "Geometry.h"
#include "Kspace.h"
#include "Energies.h"
#include "NEGFObject.h"
#include "SelfEnergy.h"
#include "GreenFunctions.h"

namespace negf {

    /** Quasi-Fermilevel (QFL) calculator. <BR>
     *  The class computes the electron/hole QFL out of the numerical density of states (from GR) and
     *  a given electron/hole density. <BR>
     *  Suggested mode of operations:
     *  1. set_densities(..)
     *  2. calculate()
     *  3. get_*_qfl()
     *  */
    class QuasiFermilevel {
    public:

        QuasiFermilevel(const Geometry * xspace_,
                    const Kspace * kspace_,
                    const Energies * energies_,
                    const GreenFunctions * gf_) throw (Exception *);

        ~QuasiFermilevel();

        /** set the densities, defined on all grid points of xspace */
        void set_densities(const vector<double> & edens_, const vector<double> & hdens_);

        /** calculate the quasi-fermilevels corresponding to the densities */
        void calculate() throw (Exception *);

        /** get the electron quasifermilevel on the interior x-points */
        vector<double> get_electron_qfl() const throw (Exception *);

        /** get the hole quasifermilevel on the interior x-points */
        vector<double> get_hole_qfl() const throw (Exception *);

        /** write to file */
        void write_to_file(char * filename) const throw (Exception *);

    protected:

        // helpers
        void calculate_dos();
        double calculate_density(double EF, uint ii, bool e_or_h, bool derivative);

        const Options        *  options;
        const Geometry       *  xspace;
        const Kspace         *  kspace;
        const Energies       *  energies;
        const GreenFunctions *  gf;
        uint                    Nx;
        uint                    NxNn;
        uint                    Nk;
        uint                    NE;
        uint                    myNE;

        // densities (assigned externally), size = num_internal_vertices
        vector<double> edens;
        vector<double> hdens;

        // quasi-fermilevels (only master thread), size = total num. vertices
        vector<double> electron_qfl;
        vector<double> hole_qfl;

        // DOS(ee,xx) - used by master thread
        vector< vector<double> > electron_dos;
        vector< vector<double> > hole_dos;
    };

} // end of namespace


#endif
