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
#include "NEGFObject.h"
using namespace negf;

BandedNEGFObject::BandedNEGFObject(const Geometry * xspace_, const Kspace * kspace_, 
						const Energies * energies_, const Options * opts_, uint offdiags):
	xspace(xspace_),
	kspace(kspace_),
	energies(energies_),
	opts(opts_),
	Nx   (xspace==NULL   ? 0 : xspace->get_num_internal_vertices()),
	NxNn (Nx*Nn),
	Nk   (kspace==NULL   ? 0 : kspace->get_number_of_points()),
	NE   (energies==NULL ? 0 : energies->get_number_of_points()),
	myNE (energies==NULL ? 0 : energies->get_my_number_of_points())
{STACK_TRACE(
	NEGF_ASSERT(xspace!=NULL && kspace!=NULL && energies!=NULL && opts!=NULL, "null pointer encountered.");
	
	// allocate memory
	logmsg->emit(LOG_INFO,"Allocating memory for a banded (2*%d offdiags) NEGF object with Nx=%d, Nn=%d, Nk=%d",offdiags,Nx,Nn,Nk);
	BMatc empty_matrix = BMatc(NxNn, NxNn, offdiags);
	vector<BMatc> empty_matrices; empty_matrices.resize(Nk, empty_matrix);
	this->matrices.resize(myNE, empty_matrices);
	NEGF_ASSERT(matrices.capacity() >= myNE, "Too little space for NEGF object.");
);}

BandedNEGFObject::~BandedNEGFObject() { }

/** STORAGE:
 * matrices[ee2][kk] stores NEGF object belonging to k(kk), E(ee), ee2 = LOCAL energy index = energies->get_my_idx(ee)
 * kk = 0...nk-1, ee = 0...nE-1
 */
BMatc & BandedNEGFObject::get(uint kidx, uint global_Eidx)
{STACK_TRACE(
	NEGF_ASSERT(global_Eidx < energies->get_number_of_points() && kidx < kspace->get_number_of_points(), "invalid index");
	NEGF_FASSERT((int)energies->get_process_computing(global_Eidx)==mpi->get_rank(), 
			"energy %.3e (point %d) seems to belong to process %d (this one: %d).",
			energies->get_energy_from_global_idx(global_Eidx), global_Eidx, energies->get_process_computing(global_Eidx), mpi->get_rank());
	uint ee2 = energies->get_my_idx(global_Eidx);
	//NEGF_ASSERT(kidx*this->myNE + ee2 < matrices.size(), "wrong index.");
	//return this->matrices[kidx*this->myNE + ee2];
	NEGF_ASSERT(ee2 < matrices.size() && kidx < matrices[ee2].size(), "wrong index.");
	return this->matrices[ee2][kidx];
);}

vector<BMatc> & BandedNEGFObject::get(uint global_Eidx)
{STACK_TRACE(
	NEGF_ASSERT(global_Eidx < energies->get_number_of_points(), "invalid index");
	NEGF_FASSERT((int)energies->get_process_computing(global_Eidx)==mpi->get_rank(), 
			"energy %.3e (point %d) seems to belong to process %d (this one: %d).",
			energies->get_energy_from_global_idx(global_Eidx), global_Eidx, energies->get_process_computing(global_Eidx), mpi->get_rank());
	uint ee2 = energies->get_my_idx(global_Eidx);
	NEGF_ASSERT(ee2 < matrices.size(), "wrong index.");
	return this->matrices[ee2];
);}

	
FullNEGFObject::FullNEGFObject(const Geometry * xspace_, const Kspace * kspace_, 
						const Energies * energies_, const Options * opts_):
	xspace(xspace_),
	kspace(kspace_),
	energies(energies_),
	opts(opts_),
	Nx   (xspace==NULL   ? 0 : xspace->get_num_internal_vertices()),
	NxNn (Nx*Nn),
	Nk   (kspace==NULL   ? 0 : kspace->get_number_of_points()),
	NE   (energies==NULL ? 0 : energies->get_number_of_points()),
	myNE (energies==NULL ? 0 : energies->get_my_number_of_points())
{STACK_TRACE(
	NEGF_ASSERT(xspace!=NULL && kspace!=NULL && energies!=NULL && opts!=NULL, "null pointer encountered.");
	
	// allocate memory
	logmsg->emit(LOG_INFO,"Allocating memory for a NEGF object with Nx=%d, Nn=%d, Nk=%d",Nx,Nn,Nk);
	Matc empty_matrix(NxNn, NxNn);
	vector<Matc> empty_matrices; empty_matrices.resize(Nk, empty_matrix);
	this->matrices.resize(myNE, empty_matrices);
	NEGF_ASSERT(matrices.capacity() >= myNE, "Too little space for NEGF object.");
);}


FullNEGFObject::~FullNEGFObject() { }

Matc & FullNEGFObject::get(uint kidx, uint global_Eidx)
{STACK_TRACE(
	NEGF_ASSERT(global_Eidx < energies->get_number_of_points() && kidx < kspace->get_number_of_points(), "invalid index");
	NEGF_FASSERT((int)energies->get_process_computing(global_Eidx)==mpi->get_rank(), 
			"energy %.3e (point %d) seems to belong to process %d (this one: %d).",
			energies->get_energy_from_global_idx(global_Eidx), global_Eidx, energies->get_process_computing(global_Eidx), mpi->get_rank());
	uint ee2 = energies->get_my_idx(global_Eidx);
	//NEGF_ASSERT(kidx*this->myNE + ee2 < matrices.size(), "wrong index.");
	//return this->matrices[kidx*this->myNE + ee2];
	NEGF_ASSERT(ee2 < matrices.size() && kidx < matrices[ee2].size(), "wrong index.");
	return this->matrices[ee2][kidx];
);}

vector<Matc> & FullNEGFObject::get(uint global_Eidx)
{STACK_TRACE(
	NEGF_ASSERT(global_Eidx < energies->get_number_of_points(), "invalid index");
	NEGF_FASSERT((int)energies->get_process_computing(global_Eidx)==mpi->get_rank(), 
			"energy %.3e (point %d) seems to belong to process %d (this one: %d).",
			energies->get_energy_from_global_idx(global_Eidx), global_Eidx, energies->get_process_computing(global_Eidx), mpi->get_rank());
	uint ee2 = energies->get_my_idx(global_Eidx);
	NEGF_ASSERT(ee2 < matrices.size(), "wrong index.");
	return this->matrices[ee2];
);}
	
