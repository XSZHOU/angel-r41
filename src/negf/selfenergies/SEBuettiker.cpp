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
#include "SEBuettiker.h"
using namespace negf;


SEBuettiker::SEBuettiker(const Hamiltonian * ham_,
						const Overlap * ov_,
						const Geometry * xspace_, 
						const Kspace * kspace_, 
						const Energies * energies_, 
						const Options * options_,
						const GreenFunctions * gf_,
						const SelfEnergy * se_contact_,
						const double & energy_parameter_):
	SelfEnergy(xspace_,kspace_,energies_,options_, constants::odSGol),
	ham(ham_),
	ov(ov_),
	gf(gf_),
	se_contact(se_contact_),
	energy_parameter(energy_parameter_),
	kT(constants::convert_from_SI(units::energy, constants::SIkb * options->get("temperature"))),
	how_many_neighbours(2)	// 2nd-NEAREST NEIGHBOUR!
{STACK_TRACE(
	NEGF_ASSERT(ham!=NULL && ov!=NULL && xspace!=NULL && kspace!=NULL && energies!=NULL && options!=NULL && gf!=NULL,
					"null pointer encountered.");
	
	this->i_am_master = (mpi->get_rank()==constants::mpi_master_rank) ? true : false;
	
	// ---------------------------------------------
	// determine the sites with fixed fermilevels
	// ---------------------------------------------
	this->fixed_fermilevels.resize(Nx, NULL);
	for (uint xx=0; xx < Nx; xx++)
	{
		Vertex * v = xspace->get_vertex(xspace->get_global_vertex_index(xx));
		// get all adjacent vertices/edges
		const vector<Edge*> & edges_near_v = xspace->get_edges_near(v);
		for (uint ii=0; ii < edges_near_v.size(); ii++) {
			Vertex * v2 = edges_near_v[ii]->get_other_vertex(v);
			if (v2->is_at_contact()) {
				this->fixed_fermilevels[xx] = v2->get_contact();
				break;
			}
		}
	}
	
	// --------------------------------------------------------------
	// initial_guess for fermilevels of Buettiker probes is zero
	// --------------------------------------------------------------
	// this->interpolate_fermilevels(); don't do it because contact fermilevels not yet set
	this->site_fermilevels.resize(Nx, 0.0);
	this->first_calculation = true;
	
	// ------------------------------------------------------------------
	// compute retarded self-energy (which is independent of everything)
	// ------------------------------------------------------------------
	this->calculate_retarded();
	
	// ------------------------------------------------------------------
	// compute lesser / greater self-energy for initial fermilevels
	// does not make sense, because Green functions were not yet computed
	// ------------------------------------------------------------------
	//this->calculate_lesser_greater();
);}


void SEBuettiker::set_energy_parameter(const double & new_parameter)
{STACK_TRACE(
	NEGF_ASSERT(new_parameter>=0.0, "negative energy parameter encountered.");
	logmsg->emit(LOG_INFO,"Setting new Buettiker energy parameter to %.3e",new_parameter);
	this->energy_parameter = new_parameter;
);}


void SEBuettiker::calculate()
{STACK_TRACE(
	// retarded GF does not have to be recomputed, was already done in constructor 
		
	// first find self-consistent Ferilevels s.th. current is constant within device
	// ... for fixed retarded Green functions!
	if (this->first_calculation) {
		this->interpolate_fermilevels();
	}
	this->compute_fermilevels();
	// then compute the lesser/greater self-energies for fixed fermilevels
	this->calculate_lesser_greater();
	
	this->first_calculation = false;
);}


void SEBuettiker::interpolate_fermilevels()
{STACK_TRACE(	
	this->site_fermilevels.resize(xspace->get_num_internal_vertices(), 0.0);
	for (uint xx=0; xx < xspace->get_num_internal_vertices(); xx++)
	{
		// vertices near contacts
		if (this->fixed_fermilevels[xx]!=NULL) {
			this->site_fermilevels[xx] = this->fixed_fermilevels[xx]->get_bc_value(quantities::fermilevel);
			continue;
		}
		
		// other vertices
		Vertex * v = xspace->get_vertex(xspace->get_global_vertex_index(xx));
		double total_distance = 0.0;
		for (uint ii=0; ii < xspace->get_num_contacts(); ii++) {
			double dist = xspace->get_distance(v, xspace->get_contact(ii));
			total_distance += dist;
		}
		NEGF_ASSERT(total_distance > 0.0, "something went wrong (no contacts???");
		for (uint ii=0; ii < xspace->get_num_contacts(); ii++) {
			double dist = xspace->get_distance(v, xspace->get_contact(ii));
			double contact_fermilevel = xspace->get_contact(ii)->get_bc_value(quantities::fermilevel);	
			this->site_fermilevels[xx] += dist / total_distance * contact_fermilevel;
			//this->site_fermilevels[xx] = 0.9;
		}
	}
);}


void SEBuettiker::calculate_retarded()
{STACK_TRACE(	
	const OVMat & M = ov->get_internal_overlap(); // defined on internal indices only, (Nx*Nn)^2
	NEGF_ASSERT(M.num_rows()==NxNn && M.num_cols()==NxNn, "inconsistent overlap matrix size.");
	
	for (uint ee2 = 0; ee2 < myNE; ee2++)
	{
		uint ee = energies->get_global_index(ee2);
		for (uint kk = 0; kk < Nk; kk++)
		{
			SEMat & SR = this->get_retarded(kk,ee);
			mult(M, -constants::imag_unit * energy_parameter, SR); // SR = -constants::imag_unit * energy_parameter * M;
		}
	}
	mpi->synchronize_processes();
);}


void SEBuettiker::calculate_lesser_greater()
{STACK_TRACE(
	logmsg->emit_header("calculating lesser and greater Buettiker self-energies");
	
	const OVMat & M = ov->get_internal_overlap(); // defined on internal indices only, (Nx*Nn)^2
	NEGF_ASSERT(M.num_rows()==NxNn && M.num_cols()==NxNn, "inconsistent overlap matrix size.");
		
	SEMat SL = SGolMat_create(NxNn);
	SEMat SG = SGolMat_create(NxNn);
	
	for (uint ee2 = 0; ee2 < myNE; ee2++)
	{
		uint ee = energies->get_global_index(ee2);
		logmsg->emit_noendl_all(LOG_INFO_L2, "p%d: SR_buett(E=%d,k=:)...  ",mpi->get_rank(),ee);
		
		// ------------------------------------------------------------------------------------
		// lesser self-energy:  SL_ij = 2i*eta*0.5(fi+fj)*M_ij
		// greater self-energy: SG_ij = (SR - SA + SL)_ij = eta*(-i - i + 2i* 0.5(fi+fj))*M_ij
		//                            = -2i*eta*(1 - 0.5(fi+fj)) M_ij
		// ... independent of k-vector
		// ------------------------------------------------------------------------------------
		const double E = energies->get_energy_from_global_idx(ee);
		SL = M;
		SG = M;
		for (uint xx=1; xx<=Nx; xx++) 
		{
			const double mu_x = site_fermilevels[xx-1];
			double fermi_x = 1.0 / (1.0 + negf_math::exp((E - mu_x) / this->kT));
			for (uint yy=1; yy<=Nx; yy++) 
			{
				if (fabs(double(xx)-double(yy)) > this->how_many_neighbours + 0.0001) continue;
				
				const double mu_y = site_fermilevels[yy-1];
				double fermi_y = 1.0 / (1.0 + negf_math::exp((E - mu_y) / this->kT));
				
				SL.multiply_block(xx,yy,  2.0*constants::imag_unit  * this->energy_parameter * 0.5*(fermi_x+fermi_y)        , Nx);
				SG.multiply_block(xx,yy, -2.0*constants::imag_unit * this->energy_parameter * (1.0 - 0.5*(fermi_x+fermi_y)) , Nx);
			}
		}
		
		for (uint kk = 0; kk < Nk; kk++) {
			SEMat & SL_kk = this->get_lesser(kk,ee);
			SEMat & SG_kk = this->get_greater(kk,ee);
			
			SL_kk = SL;
			SG_kk = SG;
		}
	}
	mpi->synchronize_processes();
);}


void SEBuettiker::compute_fermilevels()
{STACK_TRACE(
	logmsg->emit_header("computing Buettiker Fermilevels s.th. current is conserved");
	uint max_iters = 20;
	uint iter = 0;
	bool converged = false;
	for (iter = 1; iter <= max_iters; iter++) {
		logmsg->emit(LOG_INFO,"Newton step %d...",iter);
		converged = this->newton_step();
		if (converged) break;
	}
	if (converged) {
		logmsg->emit(LOG_INFO,"Newton iteration for Buettiker fermilevels converged after %d iterations.",iter);
	} else {
		NEGF_FEXCEPTION("Newton iteration for Buettiker fermilevels did not converge after %d iterations.",max_iters);
	}
);}


bool SEBuettiker::newton_step()
{STACK_TRACE(
	const double max_update = 10.0*this->kT;
	// convergence: when every fermilevel does not change more than 1 ueV
	const double converged_update = constants::convert_from_SI(units::energy, 1e-6*constants::SIec);
	
	// -----------------------------
	// NEWTON ITERATION
	// -----------------------------
	// 1. COMPUTE F(x)
	logmsg->emit(LOG_INFO,"  computing F(x)...");
	Vecd F(Nx);
	this->compute_newton_function(F); // the master thread has the entire result (all energies)
	// 1a. DEBUGGING - OUTPUT CURRENT W/ CURRENT FERMILEVELS
	this->compute_buettiker_current();
	// 2. COMPUTE J_F(x) - not sparse
	logmsg->emit(LOG_INFO,"  computing Jacobian...");
	Matd JF(Nx,Nx);
	this->compute_newton_derivative(JF);
	mpi->synchronize_processes();
	// STEPS 3-6 ARE PERFORMED BY THE MASTER THREAD ONLY
	int converged = 1;
	if (i_am_master) 
	{
		cout << "Jacobian(1,:)=\n"; for(uint ii=1; ii<=Nx; ii++) { cout << JF(1,ii) << "   "; } cout << endl;
		cout << "Jacobian(2,:)=\n"; for(uint ii=1; ii<=Nx; ii++) { cout << JF(2,ii) << "   "; } cout << endl;
		cout << "RHS=\n";           for(uint ii=1; ii<=Nx; ii++) { cout <<    F(ii) << "   "; } cout  << endl;
		
		// STEP 3: invert JF
		logmsg->emit(LOG_INFO,"    inverting...");
		invert(JF);
		
		// STEP 4: compute JF^-1*F
		logmsg->emit(LOG_INFO,"    computing update...");
		Vecd update(Nx);
		mult(JF, F, update); // update = JF*F;
		
		// STEP 5: limit the update to some maximum value
		// ... and check for convergence
		for (uint xx=1; xx<Nx; xx++) {
			if (abs(update(xx))>converged_update) {
				converged = 0;
			}
			double fac = abs(update(xx)) / max_update;
			if (fac>1.0) {
				cout << "changed fermilevel update at vertex " << xx << " from " << update(xx) << " to ";
				update(xx) = update(xx) / fac;
				cout << update(xx) << endl;
			}
		}
		
		// STEP 6: update the chemical potentials
		NEGF_ASSERT(site_fermilevels.size()==Nx, "wrong site_fermilevels size.");
		cout << "old fermilevels: "; for (uint xx=1; xx<=Nx; xx++) cout << this->site_fermilevels[xx-1] << "   ";
		cout << "\nupdate: ";		 for (uint xx=1; xx<=Nx; xx++) cout << -update(xx) << "   ";
		for (uint xx=1; xx <= Nx; xx++) {
			this->site_fermilevels[xx-1] -= update(xx);
		}
		cout << "\nnew fermilevels: "; for (uint xx=1; xx<=Nx; xx++) cout << this->site_fermilevels[xx-1] << "   ";
		
		// communicate the new fermilevels to all other MPI threads
		int root = mpi->get_rank();
		NEGF_ASSERT(root==constants::mpi_master_rank, "strange master thread.");
		mpi->broadcast(this->site_fermilevels, root);
		
		// communicate whether convergence was reached
		mpi->broadcast(converged, root);
	} else {
		// receive the new fermilevels from the master thread
		this->site_fermilevels.assign(Nx,0.0);
		int root = constants::mpi_master_rank;
		mpi->broadcast(this->site_fermilevels, root);
		
		// receive whether the Newton iteration is converged
		mpi->broadcast(converged, root);
	}
	mpi->synchronize_processes();
	if (converged==0) {
		return false;
	} else if (converged==1) {
		return true;
	} else { NEGF_EXCEPTION("Expected converged=0 or 1."); return false; }
);}

/** compute the following:
 *  Fnewton(i) = sum_E sum_n sum_k (H*GR*M*F*GA - GR*M*F*GA*H)_ii,   i=line
 */
void SEBuettiker::compute_newton_function(Vecd & result) const
{STACK_TRACE(
	double kspace_factor = 1.0 / (negf_math::pow(2.0*constants::pi, kspace->get_dimension()));
	const double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
	
	const OVMat & M = ov->get_internal_overlap(); // defined on internal indices only, (Nx*Nn)^2
	NEGF_ASSERT(M.num_rows()==NxNn && M.num_cols()==NxNn, "inconsistent overlap matrix size.");
	
	if (i_am_master) {
		for (uint ii=0; ii<Nx; ii++) {
			cout << "FL[" << ii << "]=" << this->site_fermilevels[ii] << "    ";
		} cout << endl;
	}
	
	// --------------------------------------------------------------
	// compute contribution to newton function from own energies
	// --------------------------------------------------------------
	result = Vecd(Nx);  // will store sum_E sum_k sum_n (H*GR*M*F*GA - GR*M*F*GA*H)_line,line, sum_E just over thread-local energies 
	
	// initialize helper matrices
	Matc Hsmall(NxNn,NxNn);
	Matc HminusEM(NxNn,NxNn);
	Matc F(NxNn,NxNn);
	Matc SL(NxNn,NxNn);
	Matc GRSL(NxNn,NxNn);
	Matc GRSLGA(NxNn,NxNn);
	Matc HGRSLGA(NxNn,NxNn);
	Matc GRSLGAH(NxNn,NxNn);
	Matc GRSLGA_ij(Nn,Nn);
	
	double t12 = 0.0;
	double t23 = 0.0;
	double t34 = 0.0;
	
	// the outermost loop is k because Hamiltonian assembly takes some time
	for (uint kk = 0; kk < Nk; kk++)
	{
		// get the weight of the k=point
		double wk = kspace->get_point(kk).get_weight() * kspace_factor;
		
		// get Hamiltonian of the interior points only
		ham->get_internal(kspace->get_point(kk), Hsmall);
		NEGF_ASSERT(Hsmall.num_rows()==Nx*Nn, "something went wrong.");
		
		for (uint ee2 = 0; ee2 < myNE; ee2++)
		{
			Timer * t = new Timer();
			t->click("1");
			
			uint ee = energies->get_global_index(ee2);
			
			// calculate the weight of the current energy point
			const double E = energies->get_energy_from_global_idx(ee);
			/*double E2 = (ee < energies->get_number_of_points()-1) 
						? energies->get_energy_from_global_idx(ee+1)
						: E;
			double E0 = (ee > 0) 
						? energies->get_energy_from_global_idx(ee-1)
						: E;
			double wE = 0.5 * (E2 - E0);*/
			double wE = energies->get_weight_from_global_idx(ee) / (2.0*constants::pi);
			
			// NON_ORTHOGONAL BASIS: H --> H - EM
			mult(M, -E, HminusEM); // HminusEM = (-E)*M;
			HminusEM += Hsmall; // Hsmall+HminusEM would fail!
			
			// compute the matrix MF(E) (with longitudinal energies only entering the Fermi function)
			// MF_xy(E) = 0.5*(fermi_x+fermi_y) * M_xy
			// MF does not have to be re-initialized because all nonzero entries are overwritten
			for (uint yy=1; yy<=Nx; yy++) {
				const double mu_y = site_fermilevels[yy-1];
				const double fermi_y = 1.0 / (1.0 + negf_math::exp((E - mu_y) / this->kT));
				
				for (uint xx=1; xx<=Nx; xx++) {
					NEGF_ASSERT(xspace->get_dimension()==1, "nearest-neighbour stuff works only in 1D.");
					if (fabs(double(xx)-double(yy)) > this->how_many_neighbours + 0.0001) continue;
					
					const double mu_x = site_fermilevels[xx-1];
					const double fermi_x = 1.0 / (1.0 + negf_math::exp((E - mu_x) / this->kT));
					
					// could also do it with get_block... which is faster???
					/*for (uint nx=1; nx<=Nn; nx++) {
						for (uint ny=1; ny<=Nn; ny++) {
							SL((xx-1)*Nn+nx, (yy-1)*Nn+ny) = 0.5*(fermi_x+fermi_y) * M((xx-1)*Nn+nx, (yy-1)*Nn+ny);
						}
					}*/
					Matc aM_xy(Nn,Nn);
					M.get_block(xx,yy, aM_xy, Nx);
					aM_xy *= (0.5*(fermi_x+fermi_y));
					SL.fill_block(xx, yy, aM_xy, Nx);
		
					// HACK!!!
					/*if (xx==1 || yy==1 || xx==Nx || yy==Nx) {
						SL_xy *= 0.0;
					}*/
				}
			}
			SL *= (2.0 * constants::imag_unit * this->energy_parameter);
			
			/*
			// check that MF is symmetric
			Matc MF_tmp(NxNn,NxNn);
			conjtrans(SL, MF_tmp);
			MF_tmp -= SL;
			double MFnorm = negf_math::matrix_norm(MF_tmp);
			NEGF_FASSERT(MFnorm < 1e-10, "MF was not symmetric! norm of residual=%e",MFnorm);
			*/
			
			// add contact self-energy to 2i*eta*MF
			const SEMat & SL_other = se_contact->get_lesser(kk,ee);	// always exists thanks to SelfEnergies::initial_guess()
			if (kk==0) cout << "E=" << E <<": |SL_Buett|=" << negf_math::matrix_norm(SL) << ", |SL_cont|=" << negf_math::matrix_norm(SL_other) << endl;
			SL += SL_other;
			
			// DO NOT ADD ANYTHING IF THE CONTACT SELF_ENERGY IS ZERO
			if (negf_math::matrix_norm(SL_other) < 1e-10) continue;
			
			// get retarded and advanced GF
			const Matc & GR = gf->get_retarded(kk,ee);
			const Matc & GA = gf->get_advanced(kk,ee);
			
			t->click("2");
			mult(GR, SL, GRSL);
			mult(GRSL, GA, GRSLGA);
			mult(HminusEM, GRSLGA, HGRSLGA);
			mult(GRSLGA, HminusEM, GRSLGAH);
			 
			t->click("3");			
			for (uint line=1; line<=Nx; line++) 
			{
				// ---------------------------------------------------------------------------------------
				// if the vertex is adjacent to a contact, the Newton function is given by mu - mu_c = 0
				// ---------------------------------------------------------------------------------------
				if (this->fixed_fermilevels[line-1] != NULL) {
					result(line) = this->site_fermilevels[line-1] - this->fixed_fermilevels[line-1]->get_bc_value(quantities::fermilevel);
					// is performed nE*nk times, but this does not matter
					continue;
				}
				
				// get the entries (...)_ii
				Matc result1(Nn,Nn); HGRSLGA.get_block(line, line, result1, Nx);
				Matc result2(Nn,Nn); GRSLGAH.get_block(line, line, result2, Nx);
				
				// take the trace over the degrees of freedom of result1-result2
				cplx trace = 0.0;
				for (uint nn=1; nn<=Nn; nn++) {
					trace += result1(nn,nn) - result2(nn,nn);
				}
				trace = 1.0 / hbar * trace;
				NEGF_FASSERT(fabs(trace.imag()) < 100.0*constants::imag_err, "encountered complex trace: (%.4e,%.4e)",
						trace.real(), trace.imag());
				double trace_dbl = trace.real();
				
				// add to the newton function with the appropriate k-weight and E-weight
				result(line) += wE * wk * trace_dbl;
			}
			t->click("4");
			
			t12 += t->get_seconds_between("1", "2");
			t23 += t->get_seconds_between("2", "3");
			t34 += t->get_seconds_between("3", "4");
			delete t;
		}
	}
	logmsg->emit_all(LOG_INFO,"Time used (p%d): 1-2:%.3g[s], 2-3:%.3g[s], 3-4:%.3g[s]",mpi->get_rank(),t12,t23,t34);
	logmsg->emit(LOG_INFO,"");
	mpi->synchronize_processes();
	
	// -------------------------------------------------------------------
	// communicate result to master thread which adds up everything
	// it is also the master thread alone which does the Newton iteration
	// -------------------------------------------------------------------
	if (i_am_master) {
		// collect the pieces
		for (int pp=0; pp<mpi->get_num_procs(); pp++) 
		{
			if (pp==constants::mpi_master_rank) continue;
			vector<double> other_thread_result;
			int tag = pp;
			uint siz = Nx;
			mpi->recv(other_thread_result, siz, pp, tag);
			NEGF_ASSERT(other_thread_result.size()==Nx, "something went wrong.");
			// add to own result
			for (uint xx = 0; xx < Nx; xx++) {
				// ... EXCEPT if the vertex is near a contact!
				if (this->fixed_fermilevels[xx] != NULL) continue;
				
				result(xx+1) += other_thread_result[xx];
			}
		}
	} else {		
		// send to master thread
		int dest = constants::mpi_master_rank;
		int tag = mpi->get_rank();
		vector<double> result_vector; 
		result_vector.resize(Nx, 0.0);
		for (uint xx = 0; xx < Nx; xx++) {
			result_vector[xx] = result(xx+1);
		}
		mpi->send(result_vector, dest, tag);
	}
	mpi->synchronize_processes();
	
	return;
);}


/** compute the Jacobian  dFnewton(i)/dmu_j
 *    = sum_E sum_n sum_k f'(E,mu_j) ((H*GR*M)_ij*GA_ji - (GR*M)_ij*(GA*H)_ji)_ii
 */
void SEBuettiker::compute_newton_derivative(Matd & result) const
{STACK_TRACE(	
	double kspace_factor = 1.0 / (negf_math::pow(2.0*constants::pi, kspace->get_dimension()));
	const double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
	
	result = Matd(Nx,Nx); // initialize to zero
		
	// --------------------------------------------------------------
	// compute contribution to newton function from own energies
	// --------------------------------------------------------------
	
	const OVMat & M = ov->get_internal_overlap(); // defined on internal indices only, (Nx*Nn)^2
	NEGF_ASSERT(M.num_rows()==NxNn && M.num_cols()==NxNn, "inconsistent overlap matrix size.");
	
	// initialize helper matrices
	Matc Hsmall(NxNn,NxNn);
	Matc HminusEM(NxNn,NxNn);
	Matc HGR(NxNn,NxNn);
	Matc GAH(NxNn,NxNn);
	Matc GRM(NxNn,NxNn);
	Matc MGA(NxNn,NxNn);
	Matc MGAH(NxNn,NxNn);
	Matc HGRM(NxNn,NxNn);
	Matc result1(Nn,Nn);
	Matc result2(Nn,Nn);
	Matc result3(Nn,Nn);
	Matc result4(Nn,Nn);
	Matc tmp(Nn,Nn);
	
	double t12 = 0.0;
	double t23 = 0.0;
	
	int next = 0; int count = 0;
	logmsg->init_progress_bar(LOG_INFO_L1, "   progress of master thread, Nk=", Nk);
	for (uint kk = 0; kk < Nk; kk++)
	{
		// get the weight of the k=point
		double wk = kspace->get_point(kk).get_weight() * kspace_factor;
		
		// get Hamiltonian of the interior points only
		ham->get_internal(kspace->get_point(kk), Hsmall);
		NEGF_ASSERT(Hsmall.num_rows()==Nx*Nn, "something went wrong.");
		
		for (uint ee2 = 0; ee2 < myNE; ee2++)
		{
			Timer * t = new Timer();
			t->click("1");
			uint ee = energies->get_global_index(ee2);
			
			// calculate the weight of the current energy point
			const double E = energies->get_energy_from_global_idx(ee);
			/*double E2 = (ee < energies->get_number_of_points()-1) 
						? energies->get_energy_from_global_idx(ee+1)
						: E;
			double E0 = (ee > 0) 
						? energies->get_energy_from_global_idx(ee-1)
						: E;
			double wE = 0.5 * (E2 - E0);*/
			double wE = energies->get_weight_from_global_idx(ee) / (2.0*constants::pi);
			
			// NON_ORTHOGONAL BASIS: H --> H - EM
			mult(M, -E, HminusEM); // HminusEM = (-E)*M;
			HminusEM += Hsmall;
						
			// IF THE CONTACT SELF_ENERGY IS ZERO NOTHING WAS ADDED TO THE NEWTON FUNCTION
			// HENCE DO NOT ADD ANY CONTRIBUTION TO THE JACOBIAN FROM THIS ENERGY
			const SEMat & SL_other = se_contact->get_lesser(kk,ee);	// always exists thanks to SelfEnergies::initial_guess()
			if (negf_math::matrix_norm(SL_other) < 1e-10) {
				for (uint xx=1; xx<=Nx; xx++) 
				{
					// ---------------------------------------------------------------------------------------
					// if the vertex is adjacent to a contact, the Newton function is given by mu - mu_c = 0
					// 
					// in the end everything is added to the Jacobian of the master thread. but these entries
					// need to be 1.0 on the diagonal in the master thread. hence the master thread shall 
					// fill in 1.0 himself and any other thread shall stay with zeros
					// ---------------------------------------------------------------------------------------
					if (this->fixed_fermilevels[xx-1] != NULL) {
						if (i_am_master) {
							result(xx,xx) = 1.0; // performed nE*nk times, but this does not matter
						} else {
							// entire line is zero
						}
						continue;
					}
				}
				continue;
			}
			
			// get retarded and advanced GF
			const Matc & GR = gf->get_retarded(kk,ee);
			const Matc & GA = gf->get_advanced(kk,ee);
			
			// compute all neccessary NnNx*NnNx - matrices
			/*Matc M2(NxNn,NxNn);
			M2 = M;
			Matc M2fill(Nn,Nn);
			for (uint yy=1; yy<=Nx; yy++) {
				M2.fill_block( 1,yy,M2fill);
				M2.fill_block(Nx,yy,M2fill);
			}
			for (uint xx=2; xx<=Nx-1; xx++) {
				M2.fill_block(xx, 1,M2fill);
				M2.fill_block(xx,Nx,M2fill);
			}
			mult(GR,M2,GRM);
			mult(M2,GA,MGA);*/
			
			mult(GR,M,GRM); // GRM  = GR*M;
			mult(M,GA,MGA); // MGA  = M*GA;
			
			mult(HminusEM,GR,HGR);   	// HGR  = HminusEM*GR;
			mult(GA,HminusEM,GAH);   	// GAH  = GA*HminusEM;
			mult(HminusEM,GRM,HGRM); 	// HGRM = HminusEM*GRM;
			mult(MGA, HminusEM, MGAH);  // MGAH = MGA*HminusEM;
			
			t->click("2");
			for (uint xx=1; xx<=Nx; xx++) 
			{
				// ---------------------------------------------------------------------------------------
				// if the vertex is adjacent to a contact, the Newton function is given by mu - mu_c = 0
				// 
				// in the end everything is added to the Jacobian of the master thread. but these entries
				// need to be 1.0 on the diagonal in the master thread. hence the master thread shall 
				// fill in 1.0 himself and any other thread shall stay with zeros
				// ---------------------------------------------------------------------------------------
				if (this->fixed_fermilevels[xx-1] != NULL) {
					if (i_am_master) {
						result(xx,xx) = 1.0; // performed nE*nk times, but this does not matter
					} else {
						// entire line is zero
					}
					continue;
				}
				// ---------------------------------------------
				// otherwise, find line xx of the Jacobian
				// ---------------------------------------------
				for (uint yy=1; yy<=Nx; yy++) 
				{
					// compute result1 = (H*GR)_xy * (M*GA)_yx
					Matc HGR_xy(Nn,Nn); HGR.get_block(xx,yy, HGR_xy, Nx);
					Matc MGA_yx(Nn,Nn); MGA.get_block(yy,xx, MGA_yx, Nx);
					mult(HGR_xy, MGA_yx, result1); // result1 = HGR_xy * MGA_yx;
					
					// compute result2 = (H*GR*M)_xy * GA_yx
					Matc HGRM_xy(Nn,Nn); HGRM.get_block(xx,yy, HGRM_xy, Nx);
					Matc GA_yx(Nn,Nn);   GA  .get_block(yy,xx, GA_yx  , Nx);
					mult(HGRM_xy, GA_yx, result2); // result2 = HGRM_xy * GA_yx;
					
					// compute result3 = GR_xy * (M*GA*H)_yx
					Matc GR_xy(Nn,Nn);   GR  .get_block(xx,yy, GR_xy  , Nx);
					Matc MGAH_yx(Nn,Nn); MGAH.get_block(xx,yy, MGAH_yx, Nx);
					mult(GR_xy, MGAH_yx, result3); // result3 = GR_xy * MGAH_yx;
					
					// compute result4 = (GR*M)_xy*(GA*H)_yx
					Matc GRM_xy(Nn,Nn); GRM.get_block(xx,yy, GRM_xy, Nx);
					Matc GAH_yx(Nn,Nn); GAH.get_block(xx,yy, GAH_yx, Nx);
					mult(GRM_xy, GAH_yx, result4); // result4 = GRM_xy * GAH_yx;
					
					// take the trace of result1-result2
					cplx trace = 0.0;
					for (uint nn=1; nn <= Nn; nn++) {
						trace += result1(nn,nn) + result2(nn,nn) - result3(nn,nn) - result4(nn,nn);
					}
					trace = constants::imag_unit * this->energy_parameter / hbar * trace; // inclusion of imag_unit is especially important!!!
					NEGF_FASSERT(fabs(trace.imag()) < 100.0*constants::imag_err, "encountered complex trace: (%.4e,%.4e)",
						trace.real(), trace.imag());
					double trace_dbl = trace.real();
					
					// find the derivative of the fermi function w.r.t. the chemical potential:
					double nu = (E - site_fermilevels[yy-1]) / this->kT;
					double fermi_der = 1.0/(2.0+negf_math::exp(nu)+negf_math::exp(-nu)) * 1.0/this->kT;
					
					// multiply by the derivative of the fermi function, the k-weight, and the E-weight
					// and add to result[jj-1]
					result(xx,yy) += wE * wk * fermi_der * trace_dbl;
				}
			}
			t->click("3");
			
			t12 += t->get_seconds_between("1", "2");
			t23 += t->get_seconds_between("2", "3");
			delete t;
		}
		if(next == count++) next = logmsg->set_progress_bar(count, Nk);
	}	
	logmsg->end_progress_bar();
	logmsg->emit_all(LOG_INFO,"Time used by process %d: 1-2:%.3g[s], 2-3: %.3g[s]",mpi->get_rank(),t12,t23);
	logmsg->emit(LOG_INFO,"");
	mpi->synchronize_processes();
	
	// -------------------------------------------------------------------
	// communicate result to master thread which adds up everything
	// it is also the master thread alone which does the Newton iteration
	// -------------------------------------------------------------------
	if (i_am_master) {
		// collect the pieces
		for (int pp=0; pp<mpi->get_num_procs(); pp++) 
		{
			if (pp==constants::mpi_master_rank) continue;
			Matd other_thread_result(Nx,Nx);
			int tag = pp;
			mpi->recv(other_thread_result, pp, tag);
			
			// add to own result
			// entries for vertices near contacts were set to zero for threads other than the master thread!
			result += other_thread_result;
		}
	} else {		
		// send matrix to master thread
		int dest = constants::mpi_master_rank;
		int tag = mpi->get_rank();
		mpi->send(result, dest, tag);
	}
	mpi->synchronize_processes();
);}


void SEBuettiker::compute_buettiker_current()
{STACK_TRACE(
	logmsg->emit_header("computing Buettiker current");
	double kspace_factor = 1.0 / (negf_math::pow(2.0*constants::pi, kspace->get_dimension()));
	const double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
	
	const OVMat & M = ov->get_internal_overlap(); // defined on internal indices only, (Nx*Nn)^2
	NEGF_ASSERT(M.num_rows()==NxNn && M.num_cols()==NxNn, "inconsistent overlap matrix size.");
		
	// --------------------------------------------------------------
	// compute contribution to newton function from own energies
	// --------------------------------------------------------------
	Vecd result(Nx-1);  
	// will store sum_E sum_k sum_n T_i,i+1 (GR*SL*GA)_i+1,i - (GR*SL*GA)_i+1,i T_i,i+1, sum_E just over thread-local energies 
	// Here T=H-EM
	
	// compute SL
	this->calculate_lesser_greater();
		
	// initialize helper matrices
	Matc H(Nvert*Nn,Nvert*Nn);
	Matc Hsmall(NxNn,NxNn);
	Matc HminusEM(NxNn,NxNn);
	Matc SL(NxNn,NxNn);
	Matc GLtmp(NxNn,NxNn);
	Matc GL(NxNn,NxNn);
	Matc tmp1(Nn,Nn);
	Matc tmp2(Nn,Nn);
	
	// the outermost loop is k because Hamiltonian assembly takes some time
	for (uint kk = 0; kk < Nk; kk++)
	{
		// get the weight of the k=point
		double wk = kspace->get_point(kk).get_weight() * kspace_factor;
		
		// get Hamiltonian of the interior points only
		ham->get_internal(kspace->get_point(kk), Hsmall);
		NEGF_ASSERT(Hsmall.num_rows()==Nx*Nn, "something went wrong.");
				
		for (uint ee2 = 0; ee2 < myNE; ee2++)
		{			
			uint ee = energies->get_global_index(ee2);
			
			// calculate the weight of the current energy point
			const double E = energies->get_energy_from_global_idx(ee);
			/*double E2 = (ee < energies->get_number_of_points()-1) 
						? energies->get_energy_from_global_idx(ee+1)
						: E;
			double E0 = (ee > 0) 
						? energies->get_energy_from_global_idx(ee-1)
						: E;
			double wE = 0.5 * (E2 - E0) / (2.0*constants::pi);*/
			double wE = energies->get_weight_from_global_idx(ee) / (2.0*constants::pi);
			
			// NON_ORTHOGONAL BASIS: H --> H - EM
			mult(M, -E, HminusEM); // HminusEM = (-E)*M;
			HminusEM += Hsmall; // Hsmall+HminusEM would fail!
			
			// get retarded and advanced GF
			const Matc & GR = gf->get_retarded(kk,ee);
			const Matc & GA = gf->get_advanced(kk,ee);
			
			
						
			// IF THE CONTACT SELF_ENERGY IS ZERO DO NOT ADD SHIT
			const SEMat & SL_other = se_contact->get_lesser(kk,ee);	// always exists thanks to SelfEnergies::initial_guess()
			
			SL = this->get_lesser(kk,ee);
			if (kk==0 && ee%10==0) cout << "E=" << E <<": |SL_Buett|=" << negf_math::matrix_norm(SL) << ", |SL_cont|=" << negf_math::matrix_norm(se_contact->get_lesser(kk,ee));
			SL += SL_other;
			
			//if (ee==100 && kk==0) cout << "SL(100,0)=" << endl << SL << endl;
			
			// calculate GL
			mult(SL, GA, GLtmp); // GLtmp = SL*GA;
			mult(GR, GLtmp, GL); // GL = GR*GLtmp;
			if (kk==0 && ee%10==0) cout << ", |SL_tot|=" << negf_math::matrix_norm(SL) << ", |GL|=" << negf_math::matrix_norm(GL) <<  endl;
			
			if (negf_math::matrix_norm(SL_other) < 1e-10) continue;
			
			// calculate current
			for (uint xx=1; xx<=Nx-1; xx++) {
				multiply_blocks(HminusEM, xx, xx+1, GL, xx+1, xx, tmp1, Nx);
				multiply_blocks(HminusEM, xx+1, xx, GL, xx, xx+1, tmp2, Nx);
				
				cplx trace = 0.0;
				for (uint nn=1; nn<=Nn; nn++) {
					trace += tmp1(nn,nn) - tmp2(nn,nn);
				}
				NEGF_FASSERT(fabs(trace.imag()) < 100.0*constants::imag_err, "encountered complex trace: (%.4e,%.4e)",
						trace.real(), trace.imag());
				double trace_dbl = trace.real();
				result(xx) += wE * wk * 1.0/hbar * trace_dbl;
			}
		}
	}
	logmsg->emit(LOG_INFO,"");
	mpi->synchronize_processes();
	
	// -------------------------------------------------------------------
	// communicate result to master thread which adds up everything
	// it is also the master thread alone which does the Newton iteration
	// -------------------------------------------------------------------
	if (i_am_master) {
		// collect the pieces
		for (int pp=0; pp<mpi->get_num_procs(); pp++) 
		{
			if (pp==constants::mpi_master_rank) continue;
			vector<double> other_thread_result;
			int tag = pp;
			uint siz = Nx-1;
			mpi->recv(other_thread_result, siz, pp, tag);
			NEGF_ASSERT(other_thread_result.size()==Nx-1, "something went wrong.");
			// add to own result
			for (uint xx = 0; xx < Nx-1; xx++) {				
				result(xx+1) += other_thread_result[xx];
			}
		}
		
		// OUTPUT!
		cout << "Master current:" << endl;
		for (uint ii=1; ii<=Nx-1; ii++) {
			cout << result(ii) << "    ";
		} cout << endl;
	} else {		
		// send to master thread
		int dest = constants::mpi_master_rank;
		int tag = mpi->get_rank();
		vector<double> result_vector; 
		result_vector.resize(Nx-1, 0.0);
		for (uint xx = 0; xx < Nx-1; xx++) {
			result_vector[xx] = result(xx+1);
		}
		mpi->send(result_vector, dest, tag);
	}
	mpi->synchronize_processes();
	
	return;
);}
