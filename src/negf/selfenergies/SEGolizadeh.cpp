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
#include "SEGolizadeh.h"
using namespace negf;


SEGolizadeh::SEGolizadeh(const Hamiltonian * ham_,
						const Overlap * ov_,
						const Geometry * xspace_, 
						const Kspace * kspace_, 
						const Energies * energies_, 
						const Options * options_,
						const GreenFunctions * gf_,
						const double & energy_parameter_,
						const bool momentum_relaxation_):
	SelfEnergy(xspace_,kspace_,energies_,options_, constants::odSGol),
	ham(ham_),
	ov(ov_),
	gf(gf_),
	energy_parameter(energy_parameter_),
	momentum_relaxation(momentum_relaxation_),
	security_checking(false)
{STACK_TRACE(
	NEGF_ASSERT(ham!=NULL && ov!=NULL && xspace!=NULL && kspace!=NULL && energies!=NULL && options!=NULL && gf!=NULL,
					"null pointer encountered.");
	NEGF_ASSERT(energy_parameter>=0.0, "negative energhy parameter encountered.");
);}


void SEGolizadeh::set_energy_parameter(const double & new_parameter)
{STACK_TRACE(
	NEGF_ASSERT(new_parameter>=0.0, "negative energy parameter encountered.");
	logmsg->emit(LOG_INFO,"Setting new Golizadeh energy parameter to %.3e",new_parameter);
	this->energy_parameter = new_parameter;
);}
	
	
void SEGolizadeh::calculate()
{STACK_TRACE(
	logmsg->emit_small_header("calculating Golizadeh self-energies");
	this->calculate_retarded();
	this->calculate_lesser_greater();
);}


void SEGolizadeh::calculate_retarded()
{STACK_TRACE(	
	const OVMat & M = ov->get_internal_overlap(); // defined on internal indices only, (Nx*Nn)^2
	NEGF_ASSERT(M.num_rows()==NxNn && M.num_cols()==NxNn, "inconsistent overlap matrix size.");
	
	// helper matrices
#ifdef USE_BANDED
	GLMat GR_diag = GLMat_create(NxNn);
#else
	Matc GR_diag(NxNn,NxNn);
#endif
	Matc GRM(NxNn,NxNn);
	
	for (uint ee2 = 0; ee2 < myNE; ee2++)
	{
		uint ee = energies->get_global_index(ee2);
		for (uint kk = 0; kk < Nk; kk++)
		{
#ifdef USE_BANDED
			// crop GR so that its sparsity is consistent with  GL,GG
			GLMat GR = GLMat_create(NxNn); GR = gf->get_retarded(kk,ee);
#else
			const Matc & GR = gf->get_retarded(kk,ee);
#endif
			SEMat & SR = this->get_retarded(kk,ee);
			
			if (this->momentum_relaxation) {
				// re-initialization of GR_diag not necessary since every nonzero is overwritten
				for (uint xx=1; xx<=Nx; xx++) {
					if (options->exists("GoliQWOnly") && options->get("GoliQWOnly")==1.0) {
						// only include Golizadeh broadening in QW
						NEGF_ASSERT(options->exists("GoliQWLeft") && options->exists("GoliQWRight"), "need GoliQWLeft and GoliQWRight options.");
						double xleft = options->get("GoliQWLeft");
						double xright = options->get("GoliQWRight");
						
						double xcoord = xspace->get_vertex(xspace->get_global_vertex_index(xx-1))->get_coordinate(0);
						if (xcoord<xleft || xcoord>xright) {
							continue;
						}
					}
					GR_diag.fill_block(xx, xx, GR, xx, xx, Nx);
					for (uint mm=1; mm<=Nn; mm++) {
						NEGF_ASSERT(abs(GR_diag(get_mat_idx(xx,mm,Nx),get_mat_idx(xx,mm,Nx)) - GR(get_mat_idx(xx,mm,Nx),get_mat_idx(xx,mm,Nx))) < 1e-14, "something went wrong!");
						//for (uint nn=1; nn<=Nn; nn++) {
							//GR_diag(get_mat_idx(xx,mm,Nx),get_mat_idx(xx,nn,Nx)) = GR(get_mat_idx(xx,mm,Nx),get_mat_idx(xx,nn,Nx));
						//}
					}
				}
				mult(GR_diag, M, GRM);  // tmp = GR_diag * M;
			} else {
				mult(GR, M, GRM); 		// tmp = GR * M;
			}
			mult(M, GRM, SR); 			// SR = M * tmp;
			SR *= this->energy_parameter;	
		}
	}	
	mpi->synchronize_processes();
);}

void SEGolizadeh::calculate_lesser_greater()
{STACK_TRACE(	
	const OVMat & M = ov->get_internal_overlap(); // defined on internal indices only, (Nx*Nn)^2
	NEGF_ASSERT(M.num_rows()==NxNn && M.num_cols()==NxNn, "inconsistent overlap matrix size.");
	
	// helper matrices
	GLMat GL_diag = GLMat_create(NxNn);
	Matc  GLM(NxNn,NxNn);
	SEMat SAmat = SGolMat_create(NxNn);
	
	for (uint ee2 = 0; ee2 < myNE; ee2++)
	{
		uint ee = energies->get_global_index(ee2);
		for (uint kk = 0; kk < Nk; kk++)
		{
			// --------------------------------------------
			// calculate lesser self-energy
			// --------------------------------------------
			const GLMat & GL =   gf->get_lesser(kk,ee);
			      SEMat & SL = this->get_lesser(kk,ee);
			
			if (this->momentum_relaxation) {
				// re-initialization of GL_diag not necessary since every nonzero is overwritten
				for (uint xx=1; xx<=Nx; xx++) {
					if (options->exists("GoliQWOnly") && options->get("GoliQWOnly")==1.0) {
						// only include Golizadeh broadening in QW
						NEGF_ASSERT(options->exists("GoliQWLeft") && options->exists("GoliQWRight"), "need GoliQWLeft and GoliQWRight options.");
						double xleft = options->get("GoliQWLeft");
						double xright = options->get("GoliQWRight");
						
						double xcoord = xspace->get_vertex(xspace->get_global_vertex_index(xx-1))->get_coordinate(0);
						if (xcoord<xleft || xcoord>xright) {
							continue;
						}
					}
					GL_diag.fill_block(xx, xx, GL, xx, xx, Nx);
				}
				mult(GL_diag, M, GLM);  // GLM = GL_diag * M;
			} else {
				mult(GL, M, GLM); 		// GLM = GL * M;
			}
			mult(M, GLM, SL); 			// SL = M * GLM;
			SL *= this->energy_parameter;	
			
			// --------------------------------------------
			// calculate greater self-energy
			// SG = SR - SA + SL = SR - SR* + SL
			// --------------------------------------------
			      SEMat & SG = this->get_greater(kk,ee);
			const SEMat & SR = this->get_retarded(kk,ee);
			conjtrans(SR, SAmat); 		// SAmat = conjugateTranspose(SR);
			sub(SR, SAmat, SG);   		// SG = SR - SAmat;
			SG += SL;
		}
	}
	mpi->synchronize_processes();
);}



