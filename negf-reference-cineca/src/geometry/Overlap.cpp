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
#include "Overlap.h"
using namespace negf;

Overlap::Overlap(Hamiltonian * ham_, const Geometry * xspace_) throw (Exception *):
	ham   (ham_),
	xspace(xspace_),
	Nx    (xspace->get_num_internal_vertices()),
	NxNn  (Nx*Nn),
	Nvert (xspace->get_num_vertices()),
	global_overlap_banded  (BMatc(Nvert*Nn, Nvert*Nn, constants::odOV)),
	internal_overlap_banded(BMatc(NxNn, NxNn, constants::odOV))
{STACK_TRACE(
	// default constructors for global_overlap, internal_overlap, overlap, small_overlap --> 1x1-matrices
	logmsg->emit_header("setting up overlap matrices");
	NEGF_ASSERT(ham!=NULL && xspace!=NULL, "null pointer encountered.");
	NEGF_FASSERT(Nx>0 && Nn>0 && Nvert>0, "something is wrong with Nx=%d, Nn=%d or Nvert=%d",Nx,Nn,Nvert);
	
	// obtain overlap matrix at all vertices from Hamiltonian class
	ham->get_overlap(this->overlap);
	this->global_overlap   = Matc(Nvert*Nn, Nvert*Nn);	// blown-up to num_degrees_of_freedom per site
	this->internal_overlap = Matc(NxNn, NxNn);			// internal vertices only
	this->small_overlap    = Matd(Nx,Nx);				// internal vertices, w/o other degrees of freedom
	this->global_overlap_banded   = BMatc(Nvert*Nn, Nvert*Nn, constants::odOV);
	this->internal_overlap_banded = BMatc(NxNn, NxNn, constants::odOV);
	
	Matc val_ij(Nn,Nn);
	for (uint xx=1; xx<=Nvert; xx++) 
	{
		for (uint yy=1; yy<=Nvert; yy++) 
		{
			// block matrix is diagonal:
			// M_{im,jn} = <im�jn> = <i�j><m�n> = \int dx phi_i(x) phi_j(x) delta_mn
			// delta function arises because Bloch basis vectors m,n are orthonormal
			for (uint nn=1; nn<=Nn; nn++) {
				val_ij(nn,nn) = overlap(xx,yy);
			}
			
			// fill global overlap matrix
			for (uint mm=1; mm<=Nn; mm++) {
				for (uint nn=1; nn<=Nn; nn++) {
					//global_overlap((mm-1)*Nvert+xx,(nn-1)*Nvert+yy) = val_ij(mm,nn);
					global_overlap(get_mat_idx(xx,mm,Nvert),get_mat_idx(yy,nn,Nvert)) = val_ij(mm,nn);
				}
			}
			
			// fill internal overlap matrix if appropriate
			if (xspace->get_vertex(xx-1)->get_index_internal()>=0 && xspace->get_vertex(yy-1)->get_index_internal()>=0) 
			{
				uint ii = xspace->get_internal_vertex_index(xx-1) + 1;
				uint jj = xspace->get_internal_vertex_index(yy-1) + 1;
				
				internal_overlap.fill_block(ii,jj,val_ij,Nx);
				
				this->small_overlap(ii,jj) = overlap(xx,yy);
			}
		}
	}
	 
	this->global_overlap_banded   = this->global_overlap;	// discards far-offdiagonal terms
	this->internal_overlap_banded = this->internal_overlap;	// discards far-offdiagonal terms
);}


const OVMat & Overlap::get_overlap() const
{STACK_TRACE(
#ifdef USE_BANDED
	return this->global_overlap_banded;
#else 
	return this->global_overlap;
#endif
);}

const OVMat & Overlap::get_internal_overlap() const
{STACK_TRACE(
#ifdef USE_BANDED
	return this->internal_overlap_banded;
#else 
	return this->internal_overlap;
#endif
);}
