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
#include "PoissonProblem.h"
#include <iomanip>
using namespace negf;

PoissonProblem::PoissonProblem(const Geometry * grid_, const MaterialDatabase * db_, const Options * options) throw (Exception *):
    grid(grid_),
    db(db_)
{STACK_TRACE(
    logmsg->emit_header("setting up Poisson problem");
    NEGF_ASSERT(grid!=NULL && db!=NULL, "null pointer encountered.");

    // create temperature, eff. mass, eff. DOS, bandedges
    this->temperature = options->get("temperature");
    double kT = constants::convert_from_SI(units::energy, constants::SIkb * this->temperature);

    const double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
    const double m0   = constants::convert_from_SI(units::mass, constants::SIm0);

    uint num_verts = grid->get_num_vertices();

    // ----------------------------------------------------
    // construct vertex-based me, mh, Nc, Nv, Ec, Ev
    // ----------------------------------------------------
    this->emass.resize(num_verts, 0.0);
    this->hmass.resize(num_verts, 0.0);
    this->cbedge.resize(num_verts, 0.0);
    this->vbedge.resize(num_verts, 0.0);
    this->Nc.resize(num_verts, 0.0);
    this->Nv.resize(num_verts, 0.0);

    for (uint ii=0; ii<num_verts; ii++)
    {
        vector<Element *> elems_near_ii = grid->get_elems_near(grid->get_vertex(ii));
        double tot_measure = 0.0;
        for (uint jj=0; jj<elems_near_ii.size(); jj++) {
            NEGF_ASSERT(elems_near_ii[jj]->get_type()==element_type::interval, "expected interval.");
            double el_meas = elems_near_ii[jj]->get_edge(0)->get_length();
            tot_measure += el_meas;

            const PropertyContainer<double> * mat = elems_near_ii[jj]->get_region()->get_material();

            emass[ii] += el_meas * m0 * mat->get("electron_effective_mass");
            hmass[ii] += el_meas * m0 * mat->get("hole_effective_mass");
            cbedge[ii] += el_meas * this->get_cbedge(mat, this->temperature, this->db);
            vbedge[ii] += el_meas * mat->get("valence_band_edge");
        }
        NEGF_ASSERT(tot_measure > 0.0, "something went wrong.");
        emass[ii] /= tot_measure;
        hmass[ii] /= tot_measure;
        cbedge[ii] /= tot_measure;
        vbedge[ii] /= tot_measure;

        double Nc_factor = emass[ii]*kT / (2*constants::pi*hbar*hbar);
        double Nv_factor = hmass[ii]*kT / (2*constants::pi*hbar*hbar);
        Nc[ii] = 2 * negf_math::pow(Nc_factor, 1.5);
        Nv[ii] = 2 * negf_math::pow(Nv_factor, 1.5);
    }

    // ----------------------------------------------------
    // create dielectric constant (element-based)
    // ----------------------------------------------------
    this->epsilon.resize(grid->get_num_elements(), 0.0);
    for (uint ii=0; ii<grid->get_num_elements(); ii++) {
        epsilon[ii] = constants::convert_from_SI(units::dielectric,
                grid->get_element(ii)->get_region()->get_material()->get("static_dielectric_constant")
                        * constants::SIeps0 );
    }

    // ----------------------------------------------------
    // NEW - regionwise constant doping from .cmd-file
    // ----------------------------------------------------
    this->doping.resize(num_verts, 0.0);

    const map< string, PropertyContainer<double> * > * cmdfile = options->get_cmdfile();
    PropertyContainer<double> * regs = 0;
    for (map< string, PropertyContainer<double> * >::const_iterator it = (*cmdfile).begin(); it!=(*cmdfile).end(); it++) {
        if (it->first=="regions") { regs = it->second; break;
        }
    }
    NEGF_ASSERT(regs!=0, "regions section was not found.");

    for (uint ii=0; ii<num_verts; ii++)
    {
        vector<Element *> elems_near_ii = grid->get_elems_near(grid->get_vertex(ii));
        double tot_measure = 0.0;
        for (uint jj=0; jj<elems_near_ii.size(); jj++) {
            NEGF_ASSERT(elems_near_ii[jj]->get_type()==element_type::interval, "expected interval.");
            double el_meas = elems_near_ii[jj]->get_edge(0)->get_length();
            tot_measure += el_meas;

            uint rr = elems_near_ii[jj]->get_region()->get_index();
            uint ridx = rr;
            if (rr==0) ridx = 1; // left contact
            if (rr==grid->get_num_regions()-1) ridx = grid->get_num_regions()-2; // right contact

            char buf[1000]; sprintf(buf, "region%d_doping", ridx-1); // -1 because right contact was inserted as a region
            NEGF_ASSERT(regs->is_set(buf), "could not find doping for a region.");

            doping[ii] += el_meas * constants::convert_from_SI(units::density_3d, 1e6 * regs->get(buf));
        }
        NEGF_ASSERT(tot_measure > 0.0, "something went wrong.");
        doping[ii] /= tot_measure;
    }
    logmsg->emit(LOG_INFO_L1, "DOPING:");
    for (uint ii=0; ii<num_verts; ii++) {
        //logmsg->emit(LOG_INFO_L1, "  x=%5g: dop=%5g", grid->get_vertex(ii)->get_coordinate(0), doping[ii]);
        logmsg->emit_noendl(LOG_INFO_L1, "  %5g", grid->get_vertex(ii)->get_coordinate(0), doping[ii]);
    }

    // ----------------------------------------------------
    // edensity and hdensity are set to 0
    // ----------------------------------------------------
    this->edensity.resize(num_verts, 0.0);
    this->hdensity.resize(num_verts, 0.0);

    // --------------------------------
    // create strain and polarization
    // --------------------------------
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
    }

    // ----------------------------------------------------------
    // assign appropriate boundary conditions to contacts!
    // ----------------------------------------------------------
    NEGF_ASSERT(grid->get_num_contacts()>0, "must have at least one contact!");
    logmsg->emit(LOG_INFO,"Assigning Neumann condition for potential at every contact.");
    for (uint ii=0; ii < grid->get_num_contacts(); ii++) {
        grid->get_contact(ii)->set_bndcond(quantities::potential, bndconds::BC_Neumann);
    }

    // ---------------------------
    // initial guess for potential
    // ---------------------------
    logmsg->emit(LOG_INFO,"Initializing potential values to zero.");
    elstat_pot.assign(num_verts, 0.0);

	// see initial_guess() - called later

    // --------------------------------
    // create Poisson equation
    // --------------------------------
    logmsg->emit_small_header("Setting up Poisson equation and density equations");
    this->poisson = new PoissonFEM1D(grid, static_rhs);

    poisson->set_elstat_pot(&elstat_pot);
    poisson->set_edensity(&edensity);
    poisson->set_hdensity(&hdensity);
    poisson->set_doping(&doping);
    poisson->set_epsilon(&epsilon);
    poisson->set_Nc(&Nc);
    poisson->set_Nv(&Nv);
    poisson->set_kT(kT);

    // -----------------------------------------------------
    // set up matrix and linear solver for Newton iteration
    // -----------------------------------------------------
    this->jacobian = new CSRMatrix<double>(num_verts, nonsymmetric_matrix);
    // ATTENTION: prow, icol sparsity is wrong for contact vertices - OK because not used there in PoissonFEM1D.
    // --> we should really take get_newton_derivative sparsity.
    uint nonzeros=0; uint indices[100];
    for (uint ii=0; ii<num_verts; ii++) {
        poisson->get_newton_derivative(ii, nonzeros, indices, NULL);
        NEGF_ASSERT(nonzeros>0, "must have at least 1 nonzero per line.");

        for (uint jj=0; jj<nonzeros; jj++) {
            jacobian->announce(ii, indices[jj]);
        }
    }
    jacobian->set_structure(); // sets arrays prow, icol, nonzeros

    this->rhs.resize(num_verts, 0.0);
    this->solution.resize(num_verts, 0.0);

    this->solver = new LinearSolverUmfpack(jacobian, &(rhs[0]), &(solution[0]));

    // --------------------------
    // prepare output data object
    // --------------------------
//    this->outputdata = new OutputData(grid, fnames->get_outfile());
);}


void PoissonProblem::initial_guess(bool new_method, double EF_0, double EF_1)
{STACK_TRACE(
    // ---------------------------
    // initial guess for potential
    // ---------------------------
    const double kT = constants::convert_from_SI(units::energy, constants::SIkb * this->temperature);
    const double ec   = constants::convert_from_SI(units::charge, constants::SIec);

    // note: fermilevel at contact 0 is sought with potential=0

    uint cvidx0 = grid->get_contact(0)->get_contact_vertex(0)->get_index_global();
    uint cvidx1 = grid->get_contact(1)->get_contact_vertex(0)->get_index_global();
	logmsg->emit(LOG_INFO,"Contact 0 is near x=%g, Contact 1 near x=%g.", 
		grid->get_vertex(cvidx0)->get_coordinate(0), grid->get_vertex(cvidx1)->get_coordinate(0) );

    double doping0 = doping[cvidx0];
    double doping1 = doping[cvidx1];
	double Ncv_0 = (doping0 < 0.0) ? Nv[cvidx0] : Nc[cvidx0];
	double Ncv_1 = (doping1 < 0.0) ? Nv[cvidx1] : Nc[cvidx1];
			
	double Nc_0 = Nc[cvidx0];
	double Nc_1 = Nc[cvidx1];
	double Nv_0 = Nv[cvidx0];
	double Nv_1 = Nv[cvidx1];
	double Ec_0 = cbedge[cvidx0];
	double Ev_0 = vbedge[cvidx0];
 	double Ec_1 = cbedge[cvidx1];
	double Ev_1 = vbedge[cvidx1];
	
	
	if (new_method)
	{
		logmsg->emit(LOG_INFO, "Initial guess for potential (new method using EF0=%.3f, EF1=%.3fV)...", EF_0, EF_1); 
		// 1. determine phi(0) and phi(N) from semiclassical relationship
		// 2. interpolate phi(i) linearly in between
		// 3. subtract phi(0) from the result, because the class ContactFermilevel assumes zero electrostatic potential at contact 0...
		
		// note that for n-doping, phi(0)~0 is natural. however, for p-doping that essentially introduces a constant shift value of the potential! 
		
		// doping>0:       doping = Nc * F_{0.5}(nu_n), nu_n = (EF - Ec + e*phi) / kT
		//            -->  e*phi  = Ec - EF + kT*F_{0.5}^{-1}(doping/Nc)
		//
		// doping<0:      -doping = Nv * F_{0.5}(nu_p), nu_p = (Ev - e*phi - EF) / kT
		//            -->  e*phi  = Ev - EF - kT*F_{0.5}^{-1}(-doping/Nv)
				
		double phi0 = (doping0>0.0)  
							?  (Ec_0 - EF_0 + kT*negf_math::fermihalf_inverse( doping0 / Nc_0) ) / ec
							:  (Ev_0 - EF_0 - kT*negf_math::fermihalf_inverse(-doping0 / Nv_0) ) / ec  ;
		double phi1 = (doping1>0.0)  
							?  (Ec_1 - EF_1 + kT*negf_math::fermihalf_inverse( doping1 / Nc_1) ) / ec
							:  (Ev_1 - EF_1 - kT*negf_math::fermihalf_inverse(-doping1 / Nv_1) ) / ec  ;
		double dphi = phi1 - phi0;
						
	    double dtot = grid->get_distance(grid->get_contact(1)->get_contact_vertex(0), grid->get_contact(0));
	    for (uint ii=0; ii < grid->get_num_vertices(); ii++) {
		   double frac = grid->get_distance(grid->get_vertex(ii), grid->get_contact(0)) / dtot;
 	       elstat_pot[ii] = dphi * frac;
	    }
		
	} else {	
		 if (doping0*doping1 < 0.0) {
	        logmsg->emit(LOG_INFO,"p-n situation encountered. setting potential to a linear slope:");
	
	        // p-n situation
	        double Egap_0 = Ec_0 - Ev_0;
	        //double Egap_1 = Ec_1 - Ev_1;

	        double nu_0 = kT/ec * negf_math::fermihalf_inverse(fabs(doping0) / Ncv_0);
	        double nu_1 = kT/ec * negf_math::fermihalf_inverse(fabs(doping1) / Ncv_1);
		
	        double delta_phi = Egap_0 + nu_0 - nu_1; // signs?
	        if (doping0>0.0) {
    	        // contact 0 is n, contact 1 is p
        	    delta_phi *= -1.0;
        	} else {
         	   // contact 0 is p, contact 1 is n
        	}

	        double dtot = grid->get_distance(grid->get_contact(1)->get_contact_vertex(0), grid->get_contact(0));
	        for (uint ii=0; ii < grid->get_num_vertices(); ii++) {
 	           elstat_pot[ii] = delta_phi * grid->get_distance(grid->get_vertex(ii), grid->get_contact(0)) / dtot;
	        }

	        // so far, phi=0 at contact 0, phi=delta_phi at contact 1
	        // however, we need to shift the potential in case doping0>0 (WHY???)
	        if (doping0>0.0) {
	            for (uint ii=0; ii < grid->get_num_vertices(); ii++) {
	                elstat_pot[ii] += Egap_0/2;
	            }
	        }
	    }
	}
	logmsg->emit(LOG_INFO,"inital guess of el.stat.potential:");
	for (uint ii=0; ii<grid->get_num_vertices(); ii++) {
		logmsg->emit_noendl(LOG_INFO,"%g   ", elstat_pot[ii]);
	} 
	logmsg->emit(LOG_INFO,"");
);}



void PoissonProblem::solve_one_step()
{STACK_TRACE(
    bool verbose = false;
    if (verbose) logmsg->emit(LOG_INFO_L1,"performing a single Newton step...");
    uint num_verts = grid->get_num_vertices();

    double max_update = constants::convert_from_SI(units::energy, constants::SIec * this->max_update_V)
                        / constants::convert_from_SI(units::charge, constants::SIec);

    uint nonzeros = 0;
    uint indices[100];
    double coeffs[100];
    for (uint ii=0; ii<num_verts; ii++) {
        // calculate Jacobian
        poisson->get_newton_derivative(ii, nonzeros, indices, coeffs);
        for (uint jj=0; jj<nonzeros; jj++) {
            jacobian->set(ii, indices[jj], coeffs[jj]);
        }

        // calculate RHS
        rhs[ii] = poisson->get_newton_function(ii);
    }

    // solve J*x=rhs
    solver->solve();


    // phi(n+1) = phi(n) - x, x=J^-1*F(phi(n))
    for (uint ii=0; ii<num_verts; ii++) {
        double update = (fabs(solution[ii])<max_update) ? solution[ii] : negf_math::sign(solution[ii]) * max_update;
        elstat_pot[ii] -= update;
    }

);}


void PoissonProblem::assign_new_edensity(const vector<double> & new_edensity)
{STACK_TRACE(
    NEGF_FASSERT(new_edensity.size()==grid->get_num_vertices(), "inconsistent vector size (%d instead of %d)",
                new_edensity.size(),grid->get_num_vertices());
    this->edensity = new_edensity;
);}

void PoissonProblem::assign_new_hdensity(const vector<double> & new_hdensity)
{STACK_TRACE(
    NEGF_FASSERT(new_hdensity.size()==grid->get_num_vertices(), "inconsistent vector size (%d instead of %d)",
                new_hdensity.size(),grid->get_num_vertices());
    this->hdensity = new_hdensity;
);}


void PoissonProblem::set_electrostatic_potential(const vector<double> & new_potential)
{STACK_TRACE(
    NEGF_ASSERT(new_potential.size()==elstat_pot.size(), "inconsistent vector size.");
    this->elstat_pot = new_potential;
);}


/** Workaround to use bowing parameter instead of interpolation of alpha, beta etc. for ternary materials. <BR>
    function is static --> can also be used in other classes */
double PoissonProblem::get_cbedge(const PropertyContainer<double> * mat, const double & T/*emperature*/,
                            const MaterialDatabase * const database)
{STACK_TRACE(
    int ternary = -1;
    const string & matname = mat->get_name();
    for (uint ii = 0; ii < Constants.ternary_names.size(); ii++) {
        // attention: material name was appended w/ molefraction at this stage
        if (matname.length()>=Constants.ternary_names[ii].length() &&
            matname.substr(0,Constants.ternary_names[ii].length()).compare(Constants.ternary_names[ii]) == 0) {
            ternary = ii;
            break;
        }
    }

     // valence band edge is correct for both binary and ternary materials
    double vbedge = mat->get("valence_band_edge");

    // --------------------------
    // case of ternary material
    // --------------------------
    if (ternary!=-1) {
        const TernaryPropertyContainer<double> * matt = database->get_ternary_material(mat->get_name().c_str());
        double molefraction = matt->get_molefraction();

        const PropertyContainer<double> * mat0 = matt->get_first_pure_material();
        const PropertyContainer<double> * mat1 = matt->get_second_pure_material();
        NEGF_FASSERT(mat0->is_set("bandgap_0K"),    "property \"bandgap_0K\" not found in material %s.",    mat0->get_name().c_str());
        NEGF_FASSERT(mat0->is_set("bandgap_alpha"), "property \"bandgap_alpha\" not found in material %s.", mat0->get_name().c_str());
        NEGF_FASSERT(mat0->is_set("bandgap_beta"),  "property \"bandgap_beta\" not found in material %s.",  mat0->get_name().c_str());
        NEGF_FASSERT(mat1->is_set("bandgap_0K"),    "property \"bandgap_0K\" not found in material %s.",    mat1->get_name().c_str());
        NEGF_FASSERT(mat1->is_set("bandgap_alpha"), "property \"bandgap_alpha\" not found in material %s.", mat1->get_name().c_str());
        NEGF_FASSERT(mat1->is_set("bandgap_beta"),  "property \"bandgap_beta\" not found in material %s.",  mat1->get_name().c_str());

        double bandgap_0K_mat0  = mat0->get("bandgap_0K");
        double alpha_mat0       = mat0->get("bandgap_alpha");
        double beta_mat0        = mat0->get("bandgap_beta");
        double bandgap_0K_mat1  = mat1->get("bandgap_0K");
        double alpha_mat1       = mat1->get("bandgap_alpha");
        double beta_mat1        = mat1->get("bandgap_beta");

        double bandgap_bowing = 0.0;
        if (matt->is_set("bandgap_bowing")) {
            bandgap_bowing = matt->get("bandgap_bowing");
            // please note that any molefraction dependence of the bowing parameter itself is treated in TernaryPropertyContainer.h
        }

        double bandgap_TK_mat0 =bandgap_0K_mat0 - alpha_mat0*T*T / (beta_mat0+T);
        double bandgap_TK_mat1 =bandgap_0K_mat1 - alpha_mat1*T*T / (beta_mat1+T);

        double bandgap_TK_alloy = (1.0-molefraction)*bandgap_TK_mat0 + molefraction*bandgap_TK_mat1
                                 - molefraction*(1.0-molefraction) * bandgap_bowing;

        return vbedge + bandgap_TK_alloy;
    } else {
    // --------------
    // normal case
    // --------------
        NEGF_FASSERT(mat->is_set("bandgap_0K"), "property \"bandgap_0K\" not found in material %s.", mat->get_name().c_str());
        NEGF_FASSERT(mat->is_set("bandgap_alpha"), "property \"bandgap_alpha\" not found in material %s.", mat->get_name().c_str());
        NEGF_FASSERT(mat->is_set("bandgap_beta"), "property \"bandgap_beta\" not found in material %s.", mat->get_name().c_str());

        double bandgap = mat->get("bandgap_0K");
        double alpha   = mat->get("bandgap_alpha");
        double beta    = mat->get("bandgap_beta");
        // assume storage in Kelvin, Tpar=300K, bandgap = Egap @ 0K !!!
        bandgap += /*alpha*300*300 / (beta+300)*/ - alpha*T*T / (beta+T);

        return vbedge + bandgap;
    }
);}

