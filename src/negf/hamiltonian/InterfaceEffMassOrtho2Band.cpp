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
#include "InterfaceEffMassOrtho2Band.h"
using namespace negf;

InterfaceEffMassOrtho2Band::InterfaceEffMassOrtho2Band(const Geometry * xspace_,
						 		const MaterialDatabase * db,
						 		const double temperature) throw (Exception *):
	xspace(xspace_),
	strainpol(NULL)
{STACK_TRACE(
	NEGF_ASSERT(xspace!=NULL, "null pointer encountered.");
	NEGF_ASSERT(xspace->get_dimension()==1, "class is only designed for 1D.");
	NEGF_ASSERT(Nn==2, "expected Nn=2.");
	uint Nvert = xspace->get_num_vertices();
	
	const double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
	
	// set up weights dx[ii] for each vertex which are used in new discretization
	this->dx.assign(Nvert, 0.0);
	for (uint xx=2; xx<=Nvert-1; xx++) {
		const vector<Element *> elems_near_x = xspace->get_elems_near(xspace->get_vertex(xx-1));
		NEGF_ASSERT(elems_near_x.size()==2, "expected 2 edges near x.");
		this->dx[xx-1] = 0.5 * (elems_near_x[0]->get_edge(0)->get_length() + elems_near_x[1]->get_edge(0)->get_length());
	}
	this->dx[      0] = /* 2.0* */xspace->get_elems_near(xspace->get_vertex(      0))[0]->get_edge(0)->get_length();
	this->dx[Nvert-1] = /* 2.0* */xspace->get_elems_near(xspace->get_vertex(Nvert-1))[0]->get_edge(0)->get_length();
		
	this->overlap          = DGEMatrix(Nvert,Nvert);
	this->core_hamiltonian = GEMatrix(2*Nvert,2*Nvert);
	
	// -----------------------------------------
	// rows of interior vertices
	// -----------------------------------------
	for (uint xx=2; xx<=Nvert-1; xx++) {		
		Vertex * vx = xspace->get_vertex(xx-1);
		const vector<Element *> elems_near_x = xspace->get_elems_near(vx);
		NEGF_ASSERT(elems_near_x.size()==2, "expected 2 edges near x.");
		
		if (constants::old_orthogonal) {
			this->overlap(xx,xx) = 1.0;	
		} else {
			this->overlap(xx,xx) = dx[xx-1];
		}
		
		for (uint ii = 0; ii < elems_near_x.size(); ii++) 
		{
			NEGF_ASSERT(elems_near_x[ii]->get_num_edges()==1, "expected exactly 1 edge.");
			Edge * edge = elems_near_x[ii]->get_edge(0);
			Vertex * vy = (edge->get_lower_vertex()==vx) ? edge->get_upper_vertex() : edge->get_lower_vertex();
			uint     yy = vy->get_index_global() + 1;
			double  lxy = edge->get_length();
			
			// band edges: linear combination of surrounding element values
			const PropertyContainer<double> * mat = elems_near_x[ii]->get_region()->get_material();
			double Ec = TdkpInfoDesk::get_cbedge(mat, temperature, db);
			double Ev = constants::convert_from_SI(units::energy, mat->get("valence_band_edge") * constants::SIec);
			double total_length = 0.0;
			for (uint jj = 0; jj < elems_near_x.size(); jj++) {
				total_length += elems_near_x[jj]->get_edge(0)->get_length();
			}
			
			// kinetic part of Hamiltonian
			double me = constants::convert_from_SI(units::mass, mat->get("electron_effective_mass") * constants::SIm0);
			double mh = constants::convert_from_SI(units::mass, mat->get("hole_effective_mass") * constants::SIm0);
			double efac = -hbar*hbar / (2.0 * me);
			double hfac = +hbar*hbar / (2.0 * mh);
			
			if (constants::old_orthogonal) {
				// band edges
				this->core_hamiltonian((xx-1)*Nn+1,(xx-1)*Nn+1) += Ec * lxy/total_length;
				this->core_hamiltonian((xx-1)*Nn+2,(xx-1)*Nn+2) += Ev * lxy/total_length;
				
				// symmetric pseudo-laplace, but gives shit when the discretization changes
				this->core_hamiltonian((xx-1)*Nn+1,(xx-1)*Nn+1) -= efac / (lxy*lxy);
				this->core_hamiltonian((xx-1)*Nn+1,(yy-1)*Nn+1) =  efac / (lxy*lxy);
				this->core_hamiltonian((xx-1)*Nn+2,(xx-1)*Nn+2) -= hfac / (lxy*lxy);
				this->core_hamiltonian((xx-1)*Nn+2,(yy-1)*Nn+2) =  hfac / (lxy*lxy);
			} else {
				// band edges
				this->core_hamiltonian((xx-1)*Nn+1,(xx-1)*Nn+1) += Ec * lxy/total_length * dx[xx-1];
				this->core_hamiltonian((xx-1)*Nn+2,(xx-1)*Nn+2) += Ev * lxy/total_length * dx[xx-1];
				
				// nonsymmetric, but consistent discretization of laplace operator 
				// nonsymmetricity gets lifted because of left-multiplication w/ diagonal overlap matrix
				this->core_hamiltonian((xx-1)*Nn+1,(xx-1)*Nn+1) -= efac / lxy;	// *dx[xx-1]/dx[xx-1]
				this->core_hamiltonian((xx-1)*Nn+1,(yy-1)*Nn+1) =  efac / lxy;
				this->core_hamiltonian((xx-1)*Nn+2,(xx-1)*Nn+2) -= hfac / lxy;
				this->core_hamiltonian((xx-1)*Nn+2,(yy-1)*Nn+2) =  hfac / lxy;
			}
		}
	}
	// --------------------------------------------------------------
	// handling of xx=1,Nvert rows (though I think it is not used)
	// --------------------------------------------------------------
	if (constants::old_orthogonal) {
		this->overlap(    1,    1) = 1.0;	
		this->overlap(Nvert,Nvert) = 1.0;	
	} else {
		this->overlap(    1,    1) = dx[      0];	
		this->overlap(Nvert,Nvert) = dx[Nvert-1];
	}
	// offdiagonal entry is given by Hermiticity requirement
	this->core_hamiltonian((    1-1)*Nn+1,(    2-1)*Nn+1) = this->core_hamiltonian((    2-1)*Nn+1,(    1-1)*Nn+1);
	this->core_hamiltonian((    1-1)*Nn+2,(    2-1)*Nn+2) = this->core_hamiltonian((    2-1)*Nn+2,(    1-1)*Nn+2);
	this->core_hamiltonian((Nvert-1)*Nn+1,(Nvert-2)*Nn+1) = this->core_hamiltonian((Nvert-2)*Nn+1,(Nvert-1)*Nn+1);
	this->core_hamiltonian((Nvert-1)*Nn+2,(Nvert-2)*Nn+2) = this->core_hamiltonian((Nvert-2)*Nn+2,(Nvert-1)*Nn+2);
	// diagonal entry: material and spacing are certainly the same as xx=2,Nvert-1
	this->core_hamiltonian((    1-1)*Nn+1,(    1-1)*Nn+1) = this->core_hamiltonian((    2-1)*Nn+1,(    2-1)*Nn+1);
	this->core_hamiltonian((    1-1)*Nn+2,(    1-1)*Nn+2) = this->core_hamiltonian((    2-1)*Nn+2,(    2-1)*Nn+2);
	this->core_hamiltonian((Nvert-1)*Nn+1,(Nvert-1)*Nn+1) = this->core_hamiltonian((Nvert-2)*Nn+1,(Nvert-2)*Nn+1);
	this->core_hamiltonian((Nvert-1)*Nn+2,(Nvert-1)*Nn+2) = this->core_hamiltonian((Nvert-2)*Nn+2,(Nvert-2)*Nn+2);
	
	
	// checks
	for (uint ii=1; ii<=Nvert*Nn; ii++) {
		for (uint jj=1; jj<ii; jj++) {
			NEGF_FASSERT(abs(this->core_hamiltonian(ii,jj)-this->core_hamiltonian(jj,ii))<1e-13, "H should be symmetric; instead H(%d,%d)=%.15e, H(%d,%d)=%.15e",
					ii, jj, abs(this->core_hamiltonian(ii,jj)), jj, ii, abs(this->core_hamiltonian(jj,ii)));
		}
	}
	DGEMatrix tmp1;
	tmp1 = transpose(this->overlap) - this->overlap;
	NEGF_ASSERT(negf_math::matrix_norm(tmp1) < 1e-14, "M should be symmetric.");
	GEMatrix tmp2;
	tmp2 = transpose(this->core_hamiltonian) - this->core_hamiltonian;	// real...
	NEGF_ASSERT(negf_math::matrix_norm(tmp2) < 1e-14, "H should be symmetric.");
	
	this->potential_augmented_hamiltonian       = this->core_hamiltonian;
	this->potential_and_k_augmented_hamiltonian = this->core_hamiltonian;
	
);}


// add phi to diagonal!
void InterfaceEffMassOrtho2Band::set_potential(const vector<double>& node_potential)
{STACK_TRACE(
	uint Nvert = xspace->get_num_vertices();
	NEGF_ASSERT(node_potential.size()==Nvert, "inconsistent potential vector.");

	this->potential_augmented_hamiltonian = this->core_hamiltonian;
	for (uint xx=1; xx<=Nvert; xx++) {
		for (uint nn=1; nn<=Nn; nn++) {
			if (constants::old_orthogonal) {
				this->potential_augmented_hamiltonian((xx-1)*Nn+nn,(xx-1)*Nn+nn) += node_potential[xx-1];
			} else {
				this->potential_augmented_hamiltonian((xx-1)*Nn+nn,(xx-1)*Nn+nn) += node_potential[xx-1] * dx[xx-1];
			}
		}
	}
);}


// H(k) = H(0) + hbar^2 k^2 / (2m) * M
void InterfaceEffMassOrtho2Band::assemble_hamiltonian(const double& k_transversal_nm)
{STACK_TRACE(
	double k_transversal = constants::convert_from_SI(units::density_1d, 1e9*k_transversal_nm);
	uint Nvert = xspace->get_num_vertices();
	const double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
	
	this->potential_and_k_augmented_hamiltonian = this->potential_augmented_hamiltonian;
	for (uint xx=1; xx<=Nvert; xx++) 
	{
		Vertex * vx = xspace->get_vertex(xx-1);
		const vector<Element *> elems_near_x = xspace->get_elems_near(vx);
		NEGF_ASSERT(elems_near_x.size()==1 || elems_near_x.size()==2, "unexpected number of edges near x.");
		for (uint ii = 0; ii < elems_near_x.size(); ii++) 
		{
			NEGF_ASSERT(elems_near_x[ii]->get_num_edges()==1, "expected exactly 1 edge.");
			Edge * edge = elems_near_x[ii]->get_edge(0);
			double lxy =  edge->get_length();
			double total_length = 0.0;
			for (uint jj = 0; jj < elems_near_x.size(); jj++) {
				total_length += elems_near_x[jj]->get_edge(0)->get_length();
			}
			
			// kinetic part of Hamiltonian
			const PropertyContainer<double> * mat = elems_near_x[ii]->get_region()->get_material();
			double me   = constants::convert_from_SI(units::mass, mat->get("electron_effective_mass") * constants::SIm0);
			double mh   = constants::convert_from_SI(units::mass, mat->get("hole_effective_mass") * constants::SIm0);
			double face = hbar*hbar * k_transversal*k_transversal / (2.0 * me);
			double fach = hbar*hbar * k_transversal*k_transversal / (2.0 * mh);
			
			if (constants::old_orthogonal) {
				this->potential_and_k_augmented_hamiltonian((xx-1)*Nn+1,(xx-1)*Nn+1) += face * lxy/total_length;
				this->potential_and_k_augmented_hamiltonian((xx-1)*Nn+2,(xx-1)*Nn+2) -= fach * lxy/total_length;
			} else {
				this->potential_and_k_augmented_hamiltonian((xx-1)*Nn+1,(xx-1)*Nn+1) += face * lxy/total_length * dx[xx-1];
				this->potential_and_k_augmented_hamiltonian((xx-1)*Nn+2,(xx-1)*Nn+2) -= fach * lxy/total_length * dx[xx-1];
			}
		}
	}
);}


void InterfaceEffMassOrtho2Band::set_strain(StrainPolarization * strainpol_)
{STACK_TRACE(
    NEGF_ASSERT(this->strainpol==NULL, "at the moment strain can be assigned only once.");
    this->strainpol = strainpol_;
    if (strainpol==NULL) return; // there is no strain

    uint Nvert = xspace->get_num_vertices();
    vector<double> delta_Ec(Nvert, 0.0);
    vector<double> delta_Ev(Nvert, 0.0);
    for (uint xx=1; xx<=Nvert; xx++)
    {
        Vertex * vx = xspace->get_vertex(xx-1);
        const vector<Element *> elems_near_x = xspace->get_elems_near(vx);
        NEGF_ASSERT(elems_near_x.size()==1 || elems_near_x.size()==2, "unexpected number of edges near x.");
        double total_length = 0.0;
        for (uint jj = 0; jj < elems_near_x.size(); jj++) {
            total_length += elems_near_x[jj]->get_edge(0)->get_length();
        }

        for (uint ii = 0; ii < elems_near_x.size(); ii++)
        {
            NEGF_ASSERT(elems_near_x[ii]->get_num_edges()==1, "expected exactly 1 edge.");
            Region * reg = elems_near_x[ii]->get_region();
            double lxy = elems_near_x[ii]->get_edge(0)->get_length();

            // calculate band edge shift for this region
            const PropertyContainer<double> * mat = reg->get_material();
            double exx = strainpol->get_exx_eyy(reg);
            double eyy = exx;
            double ezz = strainpol->get_ezz(reg);
            double ac  = constants::convert_from_SI(units::energy, constants::SIec * mat->get("strain_potential_ac"));
            double av  = constants::convert_from_SI(units::energy, constants::SIec * mat->get("strain_potential_av"));

            // the material file assumes total_shift = ac-av!
            double Ec_shift =  ac * (exx+eyy+ezz);
            double Ev_shift = -av * (exx+eyy+ezz); // HH/LH splitting is neglected.

            // add to core hamiltonian
            if (constants::old_orthogonal) {
                this->core_hamiltonian((xx-1)*Nn+1,(xx-1)*Nn+1) += Ec_shift * lxy/total_length;
                this->core_hamiltonian((xx-1)*Nn+2,(xx-1)*Nn+2) += Ec_shift * lxy/total_length;
            } else {
                this->core_hamiltonian((xx-1)*Nn+1,(xx-1)*Nn+1) += Ev_shift * lxy/total_length * dx[xx-1];
                this->core_hamiltonian((xx-1)*Nn+2,(xx-1)*Nn+2) += Ev_shift * lxy/total_length * dx[xx-1];
            }

            delta_Ec[xx-1] += Ec_shift * lxy/total_length;
            delta_Ev[xx-1] += Ev_shift * lxy/total_length;
        }
    }

    logmsg->emit(LOG_INFO,"strain-induced CB shift:");
    for (unsigned int ii=0; ii<delta_Ec.size(); ii++) {
        logmsg->emit_noendl(LOG_INFO, "%g   ", delta_Ec[ii]);
    }
    logmsg->emit(LOG_INFO,"\nstrain-induced VB shift:");
    for (unsigned int ii=0; ii<delta_Ev.size(); ii++) {
        logmsg->emit_noendl(LOG_INFO, "%g   ", delta_Ev[ii]);
    }
    logmsg->emit(LOG_INFO,"");
    mpi->synchronize_processes();
);}

