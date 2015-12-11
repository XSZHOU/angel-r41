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
#include "SelfEnergy.h"
using namespace negf;

SelfEnergy::SelfEnergy(const Geometry * xspace_, 
					   const Kspace * kspace_, 
					   const Energies * energies_,
					   const Options * options_, uint num_offdiags):
	options(options_),
	xspace(xspace_),
	kspace(kspace_),
	energies(energies_),
	Nx   (xspace==NULL   ? 0 : xspace->get_num_internal_vertices()),
	NxNn (Nx*Nn),
	Nk   (kspace==NULL   ? 0 : kspace->get_number_of_points()),
	NE   (energies==NULL ? 0 : energies->get_number_of_points()),
	myNE (energies==NULL ? 0 : energies->get_my_number_of_points()),
	Nvert(xspace==NULL   ? 0 : xspace->get_num_vertices())
{STACK_TRACE(	
	NEGF_ASSERT(options!=NULL && xspace!=NULL && kspace!=NULL && energies!=NULL, "null pointer encountered.");
	
	// NEGFObject constructor creates the matrices for the energies of the calling process 
#ifdef USE_BANDED
	this->SigmaR = new BandedNEGFObject(xspace, kspace, energies, options, num_offdiags);
	this->SigmaL = new BandedNEGFObject(xspace, kspace, energies, options, num_offdiags);
	this->SigmaG = new BandedNEGFObject(xspace, kspace, energies, options, num_offdiags);
#else
	this->SigmaR = new FullNEGFObject(xspace, kspace, energies, options);
	this->SigmaL = new FullNEGFObject(xspace, kspace, energies, options);
	this->SigmaG = new FullNEGFObject(xspace, kspace, energies, options);
#endif
);}

SelfEnergy::~SelfEnergy()
{STACK_TRACE(	
	delete SigmaR;
	delete SigmaL;
	delete SigmaG;
);}
