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
#include "SEContacts.h"
using namespace negf;

#define INCOHERENT_CONTACTS

SEContacts::SEContacts(const Hamiltonian * ham_,
						const Overlap * ov_,
						const Geometry * xspace_, 
						const Kspace * kspace_, 
						const Energies * energies_, 
						const Options * options_):
	SelfEnergy(xspace_,kspace_,energies_,options_, constants::odSC),
	ham(ham_),
	ov(ov_),
	calculated(false),
	security_checking(false)
{STACK_TRACE(
	NEGF_ASSERT(ham!=NULL && ov!=NULL && xspace!=NULL && kspace!=NULL && energies!=NULL && options!=NULL,
					"null pointer encountered.");
	this->temperature = options->get("temperature");
	//this->contact_fermilevels.clear();
	
	// --------------------------------------------------------------------------
	// determine for each contact the vertices directly adjacent to the device
	// --------------------------------------------------------------------------
	this->interface_vertices.resize(xspace->get_num_contacts()); // vector of empty vectors
	for (uint cc = 0; cc < xspace->get_num_contacts(); cc++) 
	{
		const vector<Vertex *> & contact_verts = xspace->get_contact(cc)->get_contact_vertices();
		for (uint ii=0; ii < contact_verts.size(); ii++) 
		{
			Vertex * v = contact_verts[ii];
			// look if there is a device vertex near the current contact vertex
			// device vertices distinguish themselves in that (they have an internal index >=0) they are_at_contact()!
			const vector<Edge*> &  edges_near_ii = xspace->get_edges_near(v);
			bool device_vertex_found = false;
			for (uint jj=0; jj<edges_near_ii.size(); jj++)
			{
				Vertex * v2 = (edges_near_ii[jj]->get_lower_vertex()==v) ? edges_near_ii[jj]->get_upper_vertex() : edges_near_ii[jj]->get_lower_vertex();
				if (!v2->is_at_contact()) {
					device_vertex_found = true;
					break;
				}
			}
			if (device_vertex_found) {
				this->interface_vertices[cc].push_back(v);
			}
		}
		NEGF_ASSERT(this->interface_vertices[cc].size()>0, "no contact vertices at the interface to the device found!");
		if (xspace->get_dimension()==1) {
			NEGF_ASSERT(this->interface_vertices[cc].size()==1, "exactly 1 interface vertex expected for a 1D x-grid.");
		}
	}
	
	// --------------------------------------------------------------------------------------------
	// determine for each interface vertex the vertex connected to it lying perpendicular to the 
	// contact-device interface inside the contact
	// --------------------------------------------------------------------------------------------
	this->second_row_vertices.resize(xspace->get_num_contacts()); // vector of empty vectors
	for (uint cc = 0; cc < xspace->get_num_contacts(); cc++) 
	{
		this->second_row_vertices[cc].resize(this->interface_vertices[cc].size(), NULL);
		const vector<Vertex *> & contact_verts = xspace->get_contact(cc)->get_contact_vertices();
		for (uint ii=0; ii < this->interface_vertices[cc].size(); ii++) 
		{
			// the sought vertex has the following properties:
			// 1. it is connected with an edge to interface_vertices[cc][ii]
			// 2. it is within the contact
			// 3. it is itself not in the list interface_vertices[cc]
			// In 1D this should determine exactly one vertex
			// In 2D/3D the following property should determine a unique vertex:
			// 4. the connecting edge is perpendicular to all interface edges around the interface vertex
			Vertex * v = contact_verts[ii];
			const vector<Edge*> &  edges_near_ii = xspace->get_edges_near(v);
			for (uint jj=0; jj<edges_near_ii.size(); jj++)
			{
				Vertex * v2 = (edges_near_ii[jj]->get_lower_vertex()==v) ? edges_near_ii[jj]->get_upper_vertex() : edges_near_ii[jj]->get_lower_vertex();
				// criterion 2
				if (!v2->is_at_contact()) {
					continue;
				}
				// criterion 3
				bool interface_vertex = false;
				for (uint kk=0; kk < this->interface_vertices[cc].size(); kk++) {
					if (this->interface_vertices[cc][kk]==v2) {
						interface_vertex = true;
						break;
					}
				}
				if (interface_vertex) {
					continue;
				}
				if (xspace->get_dimension()==1) {
					// we have found the sought vertex
					NEGF_ASSERT(this->second_row_vertices[cc][ii]==NULL, "duplicate second row vertex in 1D.");
					this->second_row_vertices[cc][ii] = v2;
				} else {
					// criterion 4
					bool edge_is_perp_to_all_interface_edges = true;
					double e1[3];
					e1[0] = v2->get_coordinate(0) - v->get_coordinate(0);
					e1[1] = v2->get_coordinate(1) - v->get_coordinate(1);
					if (xspace->get_dimension()==3) {
						e1[2] = v2->get_coordinate(2) - v->get_coordinate(2);
					} else {
						e1[2] = 0.0;
					}
					for (uint kk=0; kk < edges_near_ii.size(); kk++) 
					{
						if (kk==jj) continue;
						bool interface_edge = false;
						Vertex * v3 = (edges_near_ii[kk]->get_lower_vertex()==v) ? edges_near_ii[kk]->get_upper_vertex() : edges_near_ii[kk]->get_lower_vertex();
						for (uint ll=0; ll<this->interface_vertices[cc].size(); ll++) {
							if (this->interface_vertices[cc][ll]==v3) {
								interface_edge = true;
								break;
							}
						}
						if (!interface_edge) continue;
						
						double e2[3];
						e2[0] = v3->get_coordinate(0) - v->get_coordinate(0);
						e2[1] = v3->get_coordinate(1) - v->get_coordinate(1);
						if (xspace->get_dimension()==3) {
							e2[2] = v3->get_coordinate(2) - v->get_coordinate(2);
						} else {
							e2[2] = 0.0;
						}
						double inner_product = negf_math::vector_scalar_product_3d(e1, e2);
						if (fabs(inner_product) > edges_near_ii[jj]->get_length()*edges_near_ii[kk]->get_length()*1e-6) {
							edge_is_perp_to_all_interface_edges = false;
							break;
						}
					}
					if (edge_is_perp_to_all_interface_edges) {
						// we have found the sought vertex
						NEGF_ASSERT(this->second_row_vertices[cc][ii]==NULL, "duplicate second row vertex in 2D/3D.");
						this->second_row_vertices[cc][ii] = v2;
					}
				}
			}
			NEGF_ASSERT(this->second_row_vertices[cc][ii]!=NULL, "second row vertex not found!");
		}
	}
	
	// --------------------------------------------------------------------------------------------
	// determine dx (>0)
	// --------------------------------------------------------------------------------------------
	this->dx.resize(xspace->get_num_contacts(), 0.0);
	for (uint cc = 0; cc < xspace->get_num_contacts(); cc++) 
	{
		Vertex * v1 = interface_vertices[cc][0];
		Vertex * v2 = second_row_vertices[cc][0];
		const vector<Edge *> & edges_near_v1 = xspace->get_edges_near(v1);
		double dist_v1_v2 = 0.0;
		double dist_v1_notv2 = 0.0;
		NEGF_ASSERT(edges_near_v1.size()==2, "expected exactly 2 edges near interface vertex.");
		for (uint ii=0; ii < edges_near_v1.size(); ii++) {
			Vertex * v3 = (edges_near_v1[ii]->get_lower_vertex()==v1) ? edges_near_v1[ii]->get_upper_vertex() : edges_near_v1[ii]->get_lower_vertex();
			if (v3==v2) {
				dist_v1_v2 = xspace->get_distance(v1,v3);
			} else {
				dist_v1_notv2 = xspace->get_distance(v1,v3);
			}
		}
		NEGF_ASSERT(dist_v1_v2 > 0.0, "zero distance encountered?!");
		NEGF_ASSERT(fabs(dist_v1_v2 - dist_v1_notv2) < 1e-12, "distance of interface vertex to adjacent vertices must be equal!");
		dx[cc] = dist_v1_v2;
		/*for (uint ii=1; ii < interface_vertices[cc].size(); ii++) {
			double dx2 = xspace->get_distance(interface_vertices[cc][ii], second_row_vertices[cc][ii]);
			NEGF_ASSERT(fabs(dx[cc]-dx2) < constants::convert_from_SI(units::length, 1e-11), "inconsistent dx!");
		}*/
	}
	
	// --------------------------------------------------------------------------------------------
	// determine for each interface vertex the vertex connected to it lying perpendicular to the 
	// contact-device interface inside the DEVICE (--> there we will have self-energy)
	// --------------------------------------------------------------------------------------------
	this->device_vertices.resize(xspace->get_num_contacts()); // vector of empty vectors
	for (uint cc = 0; cc < xspace->get_num_contacts(); cc++) 
	{
		this->device_vertices[cc].resize(this->interface_vertices[cc].size(), NULL);
		const vector<Vertex *> & contact_verts = xspace->get_contact(cc)->get_contact_vertices();
		for (uint ii=0; ii < this->interface_vertices[cc].size(); ii++) 
		{
			// the sought vertex has the following properties:
			// 1. it is connected with an edge to interface_vertices[cc][ii]
			// 2. it is NOT within the contact
			// In 1D this should determine exactly one vertex
			// In 2D/3D the following property should determine a unique vertex:
			// 3. the connecting edge is perpendicular to all interface edges around the interface vertex
			Vertex * v = contact_verts[ii];
			const vector<Edge*> &  edges_near_ii = xspace->get_edges_near(v);
			for (uint jj=0; jj<edges_near_ii.size(); jj++)
			{
				Vertex * v2 = (edges_near_ii[jj]->get_lower_vertex()==v) ? edges_near_ii[jj]->get_upper_vertex() : edges_near_ii[jj]->get_lower_vertex();
				// criterion 2
				if (v2->is_at_contact()) {
					continue;
				}
				if (xspace->get_dimension()==1) {
					// we have found the sought vertex
					NEGF_ASSERT(this->device_vertices[cc][ii]==NULL, "duplicate device vertex in 1D.");
					this->device_vertices[cc][ii] = v2;
				} else {
					// criterion 3
					bool edge_is_perp_to_all_interface_edges = true;
					double e1[3];
					e1[0] = v2->get_coordinate(0) - v->get_coordinate(0);
					e1[1] = v2->get_coordinate(1) - v->get_coordinate(1);
					if (xspace->get_dimension()==3) {
						e1[2] = v2->get_coordinate(2) - v->get_coordinate(2);
					} else {
						e1[2] = 0.0;
					}
					for (uint kk=0; kk < edges_near_ii.size(); kk++) 
					{
						if (kk==jj) continue;
						bool interface_edge = false;
						Vertex * v3 = (edges_near_ii[kk]->get_lower_vertex()==v) ? edges_near_ii[kk]->get_upper_vertex() : edges_near_ii[kk]->get_lower_vertex();
						for (uint ll=0; ll<this->interface_vertices[cc].size(); ll++) {
							if (this->interface_vertices[cc][ll]==v3) {
								interface_edge = true;
								break;
							}
						}
						if (!interface_edge) continue;
						
						double e2[3];
						e2[0] = v3->get_coordinate(0) - v->get_coordinate(0);
						e2[1] = v3->get_coordinate(1) - v->get_coordinate(1);
						if (xspace->get_dimension()==3) {
							e2[2] = v3->get_coordinate(2) - v->get_coordinate(2);
						} else {
							e2[2] = 0.0;
						}
						double inner_product = negf_math::vector_scalar_product_3d(e1, e2);
						if (fabs(inner_product) > edges_near_ii[jj]->get_length()*edges_near_ii[kk]->get_length()*1e-6) {
							edge_is_perp_to_all_interface_edges = false;
							break;
						}
					}
					if (edge_is_perp_to_all_interface_edges) {
						// we have found the sought vertex
						NEGF_ASSERT(this->device_vertices[cc][ii]==NULL, "duplicate device vertex in 2D/3D.");
						this->device_vertices[cc][ii] = v2;
					}
				}
			}
			NEGF_ASSERT(this->device_vertices[cc][ii]!=NULL, "device vertex not found!");
		}
	}	
	
	
	// initialize frey_broadening array if necessary
	// stores for each own energy ee2, k-point kk, band index nn and contact cc (in this order) the complex optical potential
	if (options->exists("FreyModel") && options->get("FreyModel")>=1) 
	{
		vector<cplx>                     tmp1; tmp1.resize(2, 0.0);
		vector< vector<cplx> >           tmp2; tmp2.resize(Nn, tmp1);
		vector< vector< vector<cplx> > > tmp3; tmp3.resize(Nk, tmp2);
		this->frey_broadening.resize(myNE, tmp3);
	}
);}


void SEContacts::assign_bandedges(const vector<double> & cb_, const vector<double> & vb_)
{STACK_TRACE(
	logmsg->emit(LOG_INFO,"Assigning band edges to contact self-energies for filtering gap states later on");
	NEGF_ASSERT(cb_.size()==xspace->get_num_vertices() && vb_.size()==xspace->get_num_vertices(), "inconsistent size");
	// convert all vertices --> internal vertices
	this->cb.resize(xspace->get_num_internal_vertices(), 0.0);
	this->vb.resize(xspace->get_num_internal_vertices(), 0.0);
	for (uint ii=0; ii<xspace->get_num_internal_vertices(); ii++) {
		this->cb[ii] = cb_[xspace->get_global_vertex_index(ii)];
		this->vb[ii] = vb_[xspace->get_global_vertex_index(ii)];
	}
);}


void SEContacts::calculate()
{STACK_TRACE(
	logmsg->emit_small_header("calculating contact self-energies");
	//NEGF_ASSERT(this->contact_fermilevels.size()==xspace->get_num_contacts(), "set up contact fermilevels first!");
	if (fabs(options->get("kp_method") - 0.0) < 1e-14 && false) {
		this->calculate_easy_retarded(); // possible for one-band effective mass model
	} else {
		this->calculate_retarded();
	}
	this->calculate_lesser_greater();
);}


/** see Datta tutorial */
void SEContacts::calculate_easy_retarded()
{STACK_TRACE(
	NEGF_ASSERT(Nn==1, "only possible for single-band models!");
	
	for (uint ee2 = 0; ee2 < myNE; ee2++)
	{
		uint ee = energies->get_global_index(ee2);
		logmsg->emit_noendl_all(LOG_INFO_L2, "p%d: SR_cont(E=%d,k=:)...  ",mpi->get_rank(),ee);
		for (uint kk = 0; kk < Nk; kk++)
		{	
			Matc H(Nvert,Nvert);
			ham->get(kspace->get_point(kk), H); // 0 or kk? kk!
			NEGF_ASSERT(H.num_rows()==Nvert && H.num_cols()==Nvert, "inconsistent Hamiltonian size.");
			
			SEMat & SR = this->get_retarded(kk,ee);
			SR = SCMat_create(Nx);
			for (uint cc=0; cc < xspace->get_num_contacts(); cc++) 
			{
				NEGF_ASSERT(this->device_vertices[cc].size()==1, "inconsistent number of interface vertices.");
				NEGF_ASSERT(this->interface_vertices[cc].size()==1, "inconsistent number of interface vertices.");
				
				uint xx = this->device_vertices[cc][0]->get_index_internal();
				// computed entry will be SR(xx+1,xx+1)
				
				// -----------------------------------------------------------------------------------------
				// assumption: parabolic band structure in the lead of the form
				// E = Ec+Ek+U + 2t(1-cos(k1*a)),
				// arising from an EMA Schroedinger equation of the form -t*Psi_0 + (Ec+Ek+U+2t)*Psi_1 - t*Psi_2 = E * Psi_1
				// hence k1 = 1/a*arccos(1-(E-Ec-Ek-U)/2t)
				// (here U is the Hartree potential, Ek=hbar^2 k^2 / 2m)
				// -----------------------------------------------------------------------------------------
				uint ii = this->device_vertices   [cc][0]->get_index_global() + 1;
				uint jj = this->interface_vertices[cc][0]->get_index_global() + 1;
				//cout << "vi=" << ii-1 << ", vj=" << jj-1 << endl;
				cplx t1 = -H(ii,jj);
				
				cplx t2;
				const vector<Edge *> edges_near_ii = xspace->get_edges_near(xspace->get_vertex(ii-1));
				NEGF_ASSERT(edges_near_ii.size()==2, "something went wrong.");
				vector<Vertex *> verts_near_ii;
				for (uint edge_ii = 0; edge_ii < edges_near_ii.size(); edge_ii++) {
					verts_near_ii.push_back((edges_near_ii[edge_ii]->get_lower_vertex()==this->device_vertices[cc][0])
											? edges_near_ii[edge_ii]->get_upper_vertex()
											: edges_near_ii[edge_ii]->get_lower_vertex());
				}
				uint jj2;
				if (verts_near_ii[0]==this->interface_vertices[cc][0]) {
					jj2 = verts_near_ii[1]->get_index_global() + 1;
					t2 = -H(ii,jj2);
				} else {
					NEGF_ASSERT(verts_near_ii[1]==this->interface_vertices[cc][0], "something is inconsistent.");
					jj2 = verts_near_ii[0]->get_index_global() + 1;
					t2 = -H(ii,jj2);
				}
				NEGF_FASSERT(t1.real()>0.0 && t2.real()>0.0, 
						"coupling constants should be positive; instead, t1 = -H(%d,%d) = (%.3e, %.3e) and t1 = -H(%d,%d) = (%.3e, %.3e)",
						ii, jj, t1.real(), t1.imag(), ii, jj2, t2.real(), t2.imag());
				NEGF_ASSERT(fabs(t1.imag())<constants::imag_err && fabs(t2.imag())<constants::imag_err, "imaginary t?");
				double a = dx[cc];
				double E = energies->get_energy_from_global_idx(ee);
				cplx Ec_plus_U_plus_Ek = H(ii,ii) - t1 - t2;
				NEGF_ASSERT(fabs(Ec_plus_U_plus_Ek.imag()) < constants::imag_err, "imaginary E+U?");
				cplx z = 1.0 - (E-Ec_plus_U_plus_Ek)/(t1+t2);
				cplx k1 = 1.0/a * negf_math::acos(z);
				//cout << "a=" << a << "[nm], |k|=" << kspace->get_point(kk).get_coord_abs() << "[1/nm], |k1|=" << abs(k1) << "[1/nm], t=" << t1.real() << ", Ec+U+Ek=" << Ec_plus_U_plus_Ek.real() << endl; 
				
				SR(xx+1,xx+1) = -t1 * std::exp(constants::imag_unit*k1*a);
			}
			//if (kk==0) logmsg->emit(LOG_INFO,"|SR(E=%d,k=0)|=%.3e",ee,negf_math::matrix_norm(SR));
		}
	}
	mpi->synchronize_processes();
);}


/** SR (retarded self-energy) is hard */
void SEContacts::calculate_retarded()
{STACK_TRACE(	
	const OVMat & M = ov->get_overlap(); // defined on all indices (not just internal), (Nvert*Nn)^2
	NEGF_ASSERT(M.num_rows()==Nvert*Nn && M.num_cols()==Nvert*Nn, "inconsistent overlap matrix size.");

	double Ec_min =  (cb.size()>0) ?  1e100 : -1e100;
	double Ev_max =  (vb.size()>0) ? -1e100 :  1e100;
	if (options->exists("InjectingStatesCutoff")) {
		//NEGF_ASSERT(cb.size()>0, "No CB was assigned"); // what about the beginning of the simulation
		double delta = options->get("InjectingStatesCutoff");
		for (uint jj=0; jj<cb.size(); jj++) {
			Ec_min = min(Ec_min, cb[jj]);
		}
		for (uint jj=0; jj<vb.size(); jj++) {
			Ev_max = max(Ev_max, vb[jj]);
		}
		logmsg->emit(LOG_INFO,"CB states cutoff at %e, VB states cutoff at %e", Ec_min-delta, Ev_max+delta);
	}

	for (uint ee2 = 0; ee2 < myNE; ee2++)
	{	
                
		uint ee = energies->get_global_index(ee2);
		logmsg->emit_noendl_all(LOG_INFO_L2, "p%d: SR_cont(E=%d,k=:)...  ",mpi->get_rank(),ee);
		for (uint kk = 0; kk < Nk; kk++)
		{	
			// get Hamiltonian including energy from transversal k-vector and electrostatic potential!
			Matc H(Nvert*Nn,Nvert*Nn);
			ham->get(kspace->get_point(kk), H);	
			NEGF_ASSERT(H.num_rows()==Nvert*Nn && H.num_cols()==Nvert*Nn, "inconsistent Hamiltonian size.");
			
			SEMat & SR = this->get_retarded(kk,ee);
			SR = SCMat_create(NxNn);
                               // std::cout<<"test0... S.Z."<<std::endl;
			for (uint cc=0; cc < xspace->get_num_contacts(); cc++) 
			{
				// ---------------------------------------------------------------------------------------------------------
				// construct the matrix Hc of size 2NcNn*2NcNn, Nc = number of vertices at the interface contact-device
				// Hc(      1:NcNn,      1:NcNn) will be denoted H00
				// Hc(NcNn+1:2NcNn,      1:NcNn) will be denoted H10
				// Hc(      1:NcNn,NcNn+1:2NcNn) will be denoted H01
				// Hc(NcNn+1:2NcNn,NcNn+1:2NcNn) will be denoted H11
				// ---------------------------------------------------------------------------------------------------------
				uint Nc = this->interface_vertices[cc].size();
				
				Matc Hc(2*Nc*Nn, 2*Nc*Nn);
				Matc Mc(2*Nc*Nn, 2*Nc*Nn);
				for (uint ii=1; ii<=Nc; ii++) 
				{
					uint gii1 = this->interface_vertices[cc][ii-1]->get_index_global() + 1;
					uint gii2 = this->second_row_vertices[cc][ii-1]->get_index_global() + 1;
					
					// gii1, gii2 are the "block matrix indices" in the matrix of the Hamiltonian of the whole system
					// corresponding to the coupling of the interface vertex ii and its adjacent second-row vertex

					// ordering is always (xx-1)*Nn+nn in Hc, Mc!
					for (uint mm=1; mm<=Nn; mm++) {
						for (uint nn=1; nn<=Nn; nn++) {
							Hc(      (ii-1)*Nn+mm,       (ii-1)*Nn+nn) = H(get_mat_idx(gii1,mm,Nvert), get_mat_idx(gii1,nn,Nvert));
							Hc(      (ii-1)*Nn+mm, Nc*Nn+(ii-1)*Nn+nn) = H(get_mat_idx(gii1,mm,Nvert), get_mat_idx(gii2,nn,Nvert));
							Hc(Nc*Nn+(ii-1)*Nn+mm,       (ii-1)*Nn+nn) = H(get_mat_idx(gii2,mm,Nvert), get_mat_idx(gii1,nn,Nvert));
							Hc(Nc*Nn+(ii-1)*Nn+mm, Nc*Nn+(ii-1)*Nn+nn) = H(get_mat_idx(gii2,mm,Nvert), get_mat_idx(gii2,nn,Nvert));
							
							Mc(      (ii-1)*Nn+mm,       (ii-1)*Nn+nn) = M(get_mat_idx(gii1,mm,Nvert), get_mat_idx(gii1,nn,Nvert));
							Mc(      (ii-1)*Nn+mm, Nc*Nn+(ii-1)*Nn+nn) = M(get_mat_idx(gii1,mm,Nvert), get_mat_idx(gii2,nn,Nvert));
							Mc(Nc*Nn+(ii-1)*Nn+mm,       (ii-1)*Nn+nn) = M(get_mat_idx(gii2,mm,Nvert), get_mat_idx(gii1,nn,Nvert));
							Mc(Nc*Nn+(ii-1)*Nn+mm, Nc*Nn+(ii-1)*Nn+nn) = M(get_mat_idx(gii2,mm,Nvert), get_mat_idx(gii2,nn,Nvert));
						}
					}
				}	
				Matc Hc_00(Nc*Nn,Nc*Nn); Hc.get_submatrix(      1,   Nc*Nn,       1,   Nc*Nn, Hc_00);	
				Matc Hc_01(Nc*Nn,Nc*Nn); Hc.get_submatrix(      1,   Nc*Nn, Nc*Nn+1, 2*Nc*Nn, Hc_01); // UPPER diagonal
				Matc Hc_10(Nc*Nn,Nc*Nn); Hc.get_submatrix(Nc*Nn+1, 2*Nc*Nn,       1,   Nc*Nn, Hc_10); // LOWER diagonal
				Matc Hc_11(Nc*Nn,Nc*Nn); Hc.get_submatrix(Nc*Nn+1, 2*Nc*Nn, Nc*Nn+1, 2*Nc*Nn, Hc_11);
				
				Matc Mc_00(Nc*Nn,Nc*Nn); Mc.get_submatrix(      1,   Nc*Nn,       1,   Nc*Nn, Mc_00);	
				Matc Mc_01(Nc*Nn,Nc*Nn); Mc.get_submatrix(      1,   Nc*Nn, Nc*Nn+1, 2*Nc*Nn, Mc_01);
				Matc Mc_10(Nc*Nn,Nc*Nn); Mc.get_submatrix(Nc*Nn+1, 2*Nc*Nn,       1,   Nc*Nn, Mc_10);
				
				// ----------------------------------------------
				// construct D  = E*Mc_00-Hc_00
				//           Tu = E*Mc_01-Hc_01
				//           Tl = E*Mc_10-Hc_10
				// ----------------------------------------------
				cplx E = energies->get_energy_from_global_idx(ee);
				
				// BROADENING / SCATTERING IN CONTACTS
				// Frey model overrules optical potential
				bool frey_model        = (options->exists("FreyModel") && options->get("FreyModel")>=1);
				bool optical_potential = (options->exists("IncoherentContacts") && options->get("IncoherentContacts")==1);
				
				Matc D(Nc*Nn, Nc*Nn);
				if (!frey_model && optical_potential) {
					double eta_eV = options->exists("IncoherentContactBroadening") ? options->get("IncoherentContactBroadening") : constants::contact_eta;
					//if (kk==0) logmsg->emit(LOG_INFO,"eta=%g", eta_eV);
					E += constants::imag_unit * constants::convert_from_SI(units::energy, constants::SIec * eta_eV);
				}		
				if (frey_model) {
					// just a little bit broadening to remove divergences
					E += constants::imag_unit * constants::convert_from_SI(units::energy, constants::SIec * 1e-4);
				}
				mult(Mc_00,E,D);    // D  = E*Mc_00;
				
				if (frey_model) {
					// get band-, E-, k-, contact-dependent broadening (real+imag), 
					// multiply this diagonal matrix with diagonal part of Mc_00 and add to D
					// ordering is (xx-1)*Nn+nn
					Matc eta_frey(Nc*Nn,Nc*Nn); 
					for (uint nn=0; nn<Nn; nn++) {
						NEGF_ASSERT(frey_broadening             .size() == myNE, "wrong array size");
						NEGF_ASSERT(frey_broadening[ee2]        .size() ==   Nk, "wrong array size");
						NEGF_ASSERT(frey_broadening[ee2][kk]    .size() ==   Nn, "wrong array size");
						NEGF_ASSERT(frey_broadening[ee2][kk][nn].size() ==    2, "wrong array size");
						cplx eta  = frey_broadening[ee2][kk][nn][cc]; // NOT multiplied by imag_unit!
						for (uint xx=1; xx<=Nc; xx++) {
							uint idx = (xx-1)*Nn+nn + 1;	// xx and idx are 1-based, ee2/kk/nn are 0-based
							//eta_frey(idx,idx) = eta * Mc_00(idx,idx);
							eta_frey(idx,idx) = eta; // length is already contained
						}
					}
					D += eta_frey; // eta must be negative
				}

				D -= Hc_00;
				Matc Tu(Nc*Nn, Nc*Nn);
				mult(Mc_01, E, Tu); // Tu  = E*Mc_01;
				Tu -= Hc_01;
				Matc Tl(Nc*Nn, Nc*Nn);
				mult(Mc_10, E, Tl); // Tl  = E*Mc_10;
				Tl -= Hc_10;
				
				// ----------------------------------------------
				// construct A = -Tu^-1*D
				// ----------------------------------------------
				// find tmp=Tu^-1
				Matc tmp(Nc*Nn, Nc*Nn);
				tmp = Tu;
                                //std::cout<<"test1... S.Z."<<std::endl;
				invert(tmp);
				
				// multiply, negative
				Matc A(Nc*Nn, Nc*Nn);
				mult(tmp, D, A); // A = tmp * D
				A *= -1;
				
				// ----------------------------------------------
				// solve eigenvalue problem A*psi=lambda*psi
				// ----------------------------------------------
				uint n = A.num_rows();
				cplx AA[n*n];				// SOOOOOOO INEFFICIENT
				for (uint ii=1; ii <= n; ii++) {
					for (uint jj=1; jj <= n; jj++) {
						AA[(jj-1)*n+(ii-1)] = A(ii,jj);
					}
				}
				cplx lambda[n];			// eigenvalues
				cplx vr[n*n];
				negf_math::zgeev(n, AA, lambda, vr); 
				vector< vector<cplx> > psi;
				psi.resize(2*n);
				for (uint ii=0; ii<n; ii++) {
					psi[2*ii]  .resize(n);				// psi[2*ii]   stores eigenvector ii      --> +k
					psi[2*ii+1].resize(n);				// psi[2*ii+1] also stores eigenvector ii --> -k
					for (uint jj=0; jj<n; jj++) {
						psi[2*ii]  [jj] = vr[ii*n+jj];	// vr(jj,ii)
						psi[2*ii+1][jj] = vr[ii*n+jj];	// vr(jj,ii)
					}
				}
				
				// ------------------------------------------------------------------------
				// get the k-vectors (in transport direction) from the eigenvalues lambda
				// using analytic continuation of cosine
				// there are always 2 solutions: k and -k!
				// ------------------------------------------------------------------------
				vector<cplx> k;
				k.resize(2*n, 0.0); 
				for (uint ii=0; ii < n; ii++) 
				{
					cplx z = lambda[ii]/2.0;
					
					// note: the branch cut of acos is at |Re(z)|>1, Im(z)=0
					// points near this branch cut have Re(acos(z))~0 or 2pi and Im(acos(z))~+-pi/2
					// hence, when we add +acos and -acos to the list, we have made the simulation insensitive					
					k[ii*2]   = +1.0/dx[cc] * negf_math::acos(z);	// dx>0
					k[ii*2+1] = -1.0/dx[cc] * negf_math::acos(z);
				}
				
				// ------------------------------------------------------------------------
				// get dE/dk (in transport direction) 
				// (H_00-E(k)M_00 - 2cos(k*dx) T_u) psi = 0, T_u = EM_01 - H_01
				//  --> (-dE/dk M_00 + 2*dx*sin(k*dx) T_u - 2cos(k*dx) dE/dk M_01) psi = 0
				//  -->   dE/dk (M_00 + 2cos(k*dx) M_01) = 2*dx*sin(k*dx) * T_u
				//  -->   dE/dk = (M_00 + 2cos(k*dx) M_01)^-1 * (2*dx*sin(k*dx) * T_u)
				//              = dEdk_helper2                * dEdk_helper1
				// ------------------------------------------------------------------------
				
				vector<cplx> dEdk; dEdk.resize(k.size(), 0.0);
				Matc dEdk_helper1(Nc*Nn, Nc*Nn);
				Matc dEdk_helper2(Nc*Nn, Nc*Nn);
				Matc dEdk_helper3(Nc*Nn, Nc*Nn);
				for (uint ii=0; ii < k.size(); ii++) 
				{
					mult(   Tu, +2.0*dx[cc]*sin(k[ii]*dx[cc]), dEdk_helper1); // dEdk_helper1 = 2*dx[cc]*sin(k[ii]*dx[cc]) * Tu;
					
					mult(Mc_01, +2.0*cos(k[ii]*dx[cc]), dEdk_helper2);        // dEdk_helper2 = 2*cos(k[ii]*dx[cc]) * Mc_01;
					dEdk_helper2 += Mc_00;
					invert(dEdk_helper2);
					
					mult(dEdk_helper2, dEdk_helper1, dEdk_helper3); // dEdk_helper3 = dEdk_helper2 * dEdk_helper1;					
										
					// now we have dEdk_helper3 * psi = dE/dk * psi
					// we find the scalar dE/dk by calculating the LHS and comparing it to the RHS whenever psi_i!=0
					vector<cplx> dEdk_psi;
					dEdk_psi.resize(Nc*Nn, 0.0);
					for (uint jj=0; jj<Nc*Nn; jj++) {
						for (uint ll=0; ll<Nc*Nn; ll++) {
							dEdk_psi[jj] += dEdk_helper3(jj+1,ll+1) * psi[ii][ll];
						}
					}
					uint num_candidates = 0;
					//if (kk==0 && cc==0 && ee%20==0) cout << "   candidates for dE/dk: ";
					for (uint jj=0; jj<Nc*Nn; jj++) {
						if (abs(psi[ii][jj]) > 1e-13) {
							dEdk[ii] = dEdk_psi[jj] / psi[ii][jj];
							//if (kk==0 && cc==0 && ee%20==0) cout << dEdk[ii] << "   ";
							//if (kk==0 && ee%20==0) cout << "E=" << E.real() << ", kt=0, cc=" << cc << ", n=" << ii << ", k=" << k[ii] << ", dEdk=" << dEdk[ii] << endl;
							num_candidates++;
						}
					}
					//if (kk==0 && cc==0 && ee%20==0) cout << endl;
					NEGF_ASSERT(num_candidates>0, "is the eigenvector entirely 0???");
					NEGF_ASSERT(num_candidates==1, "there seemed to be more than one possibility for dE/dk!");
				}
				
                                //std::cout<<"test2... S.Z."<<std::endl;
				// ----------------------------------------------
				// filter out unwanted k-vectors
				// ----------------------------------------------
				vector<cplx> k_filtered;
				vector< vector<cplx> > psi_filtered;
				const double eps = 1e-13;
				for (uint ii=0; ii < k.size(); ii++) {
					const cplx & num = /*k[ii]*/ dEdk[ii];
					bool filter_criterion = false;
					if (true) {
						//filter_criterion =  /* (fabs(num.real())<eps &&      num.imag() <-eps)  // state decays to the right
						//			 	  ||*/ (     num.real() >eps && fabs(num.imag())<eps); // state propagates to the right
						
						// allow states which travel to the right (velocity > 0) and which decay, but not too fast
						const double max_cplx = constants::convert_from_SI(units::density_1d, 1.0/(constants::contact_min_xdecay * 1e-9));
						const cplx & num2 = k[ii];
						filter_criterion =  (num.real() > eps) && (num2.imag()>-eps && num2.imag()<max_cplx);
					}
					if (false) {
						filter_criterion =   (fabs(num.real())<eps &&      num.imag() >eps)  // state decays to the left
									 	  || (     num.real() <eps && fabs(num.imag())<eps); // state propagates to the left
					}
					if (options->exists("IncoherentContacts") && options->get("IncoherentContacts")==1) {
						//filter_criterion = (num.real()>-eps && num.imag()>-eps);
						
						// I think there is the problem that dEdk has become meaningless since k is complex now.
						// Hence we must orient ourselves with k again.
						const cplx & num2 = k[ii];	
						
						// upper limit for imaginary part of k-vector: 1/XXX [1/nm]
						// --> do not include states decaying over XXX nm or less 
						const double max_cplx = constants::convert_from_SI(units::density_1d, 1.0/(constants::contact_min_xdecay * 1e-9));
						
						if (Nn==1) {
							filter_criterion = (num2.real()>-eps && num2.imag()>-eps);
						}
						if (Nn==2) {
							NEGF_ASSERT(Nc==1, "Nc=1 expected.");
							double cbval = abs(psi[ii][0]);
							double vbval = abs(psi[ii][1]);
							NEGF_ASSERT(fabs(cbval+vbval - 1.0) < 1e-10, "expected unity norm.");
							if (cbval > 0.9) {
								//filter_criterion = (num2.real()>-eps && num2.imag()>-eps);
								filter_criterion = (num2.real()>-eps && num2.imag()>-eps && num2.imag()<max_cplx);
							} else {
								//filter_criterion = (num2.imag()>-eps);
								//filter_criterion = (num2.imag()>-eps && num2.imag()<max_cplx);
								filter_criterion = (num2.real()<eps && num2.imag()>-eps && num2.imag()<max_cplx); // Re(k) negative, Im(k)>=0 and not too much
								//filter_criterion = (num2.real()<eps && num2.imag()>-max_cplx && num2.imag()<eps);   // Re(k) negative, Im(k)<=0 and not too much:  DOES NOT WORK
								//if (filter_criterion) { cout << "VB-good: k=" << k[ii] << ", dEdk=" << dEdk[ii] << endl; }
							}
						}
						NEGF_ASSERT(Nn<=2, "Nn>2 contact broadening NYI.");
					}
					if (options->exists("IncludeImaginaryContactStates") && options->get("IncludeImaginaryContactStates")==1) {
						NEGF_ASSERT(xspace->get_contact(0)->get_contact_vertex(0)->get_coordinate(0) < xspace->get_contact(1)->get_contact_vertex(0)->get_coordinate(0),
								"contact 0 must be on the left, contact 1 on the right!");
						const cplx & num2 = k[ii] /*dEdk[ii]*/;
						if (true) {
							if (cc==0) { filter_criterion = filter_criterion || (fabs(num2.real())<eps && num2.imag()<-eps); }
							if (cc==1) { filter_criterion = filter_criterion || (fabs(num2.real())<eps && num2.imag()>eps); }
						}
						if (false) {
							if (cc==0) { filter_criterion = filter_criterion || (fabs(num2.real())<eps && num2.imag()>eps); }
							if (cc==1) { filter_criterion = filter_criterion || (fabs(num2.real())<eps && num2.imag()<-eps); }
						}
					}
					if (options->exists("InjectingStatesCutoff")) {
						//NEGF_ASSERT(cb.size()>0, "No CB was assigned");
						double delta = options->get("InjectingStatesCutoff");
						double EEE = energies->get_energy_from_global_idx(ee);
						if (Nn==1) { // only conduction band
							filter_criterion = filter_criterion && (EEE > Ec_min - delta);
						}
						if (Nn==2) {
							double cbval = abs(psi[ii][0]);
							double vbval = abs(psi[ii][1]);
							NEGF_ASSERT(fabs(cbval+vbval - 1.0) < 1e-10, "expected unity norm.");
							if (cbval > 0.9) {
								filter_criterion = filter_criterion && (EEE > Ec_min - delta);
							} else {
								filter_criterion = filter_criterion && (EEE < Ev_max + delta);
							}
						}
						NEGF_ASSERT(Nn<=2, "Nn>2 criterion NYI.");
					}
					if (filter_criterion) {
						k_filtered.push_back(k[ii]);
						psi_filtered.push_back(psi[ii]);
						//if (kk==0 && ee%20==0) cout << "E=" << E.real() << ", kt=0, cc=" << cc << ", n=" << ii << " filtered: k=" << k[ii] << ", dEdk=" << dEdk[ii] << endl;
					} else {
						// throw away
					}
				}
				uint Nkx = k_filtered.size(); 
				logmsg->emit_all(LOG_INFO_L2,"%d out of %d k-vectors remain after filtering.",Nkx, k.size());
				if (Nkx==0) {
					continue;
				}
				
                                //std::cout<<"test3... S.Z."<<std::endl;
				// -----------------------------------------------
				// construct the matrix P(-1)
				// -----------------------------------------------
				Matc P(Nkx,Nkx);
				for (uint ii=0; ii<Nkx; ii++) {
					P(ii+1,ii+1) = exp(-constants::imag_unit*k_filtered[ii]*(-dx[cc])); // SIGN?
				}
				
				// -----------------------------------------------------------
				// construct the matrix Phi: columns of Phi are eigenvectors
				// -----------------------------------------------------------
				Matc Phi(Nc*Nn,Nkx);
				for (uint ii=0; ii<Nc*Nn; ii++) {
					for (uint ll=0; ll<Nkx; ll++) {
						Phi(ii+1,ll+1) = psi_filtered[ll][ii];
					}
				}
			
				Matc PhiT(Nkx,Nc*Nn);
				conjtrans(Phi, PhiT); // PhiT = conjugateTranspose(Phi);
				if (security_checking) { // check invertibility of PhiT*Phi
					Matc check(Nkx,Nkx);
					mult(Phi, PhiT, check); // check = Phi*PhiT;
					invert(check); // would throw error if not possible
				}
				
				// ---------------------------------------------------------
				// check for band-orthogonality of eigenvectors!
				// ---------------------------------------------------------
				if (true && Nn==2) {
					NEGF_ASSERT(Nc==1, "Nc=1 expected.");
					for (uint ii=1; ii<=Nkx; ii++) {
						double cbval = abs(Phi(1,ii));
						double vbval = abs(Phi(2,ii));
						if (cbval > 0.9) {
							NEGF_FASSERT(vbval<1e-15, "E=%e, kk=%d: Detected mixing contact eigenstate with |cbval|=%e, |vbval|=%e.",E.real(), kk, cbval,vbval);
						} else {
							NEGF_FASSERT(cbval<1e-15, "E=%e, kk=%d: Detected mixing contact eigenstate with |cbval|=%e, |vbval|=%e.",E.real(), kk, cbval,vbval);
						}
						NEGF_ASSERT(fabs(cbval+vbval - 1.0) < 1e-10, "expected unity norm.");
					}
				}
				
				// ---------------------------------------------------------------------------
				// construct the matrix g_tilde^-1 = PhiT*(D*Phi+Tu*Phi*P)
				//                                 = PhiT*(D*Phi+Tu*tmp3)
				//                                 = PhiT*(D*Phi+tmp4)
				//                                 = PhiT*tmp5
				// ---------------------------------------------------------------------------
				
				Matc tmp3(Nkx,Nkx);
				mult(Phi, P, tmp3); // tmp3 = Phi*P;
				
				Matc tmp4(Nkx,Nkx);
				mult(Tu, tmp3, tmp4); // tmp4 = Tu*tmp3;
				
				Matc tmp5(Nkx,Nkx);
				mult(D, Phi, tmp5); // tmp5 = D*Phi;
				tmp5 += tmp4; 
				
				Matc g_tilde(Nkx,Nkx);
				mult(PhiT, tmp5, g_tilde); // g_tilde = PhiT * tmp5;
				
				// -----------------
				// invert
				// -----------------
				invert(g_tilde);
				
				// -----------------------------------------------
				// FINALLY construct the self-energy
				// -----------------------------------------------
				// construct SRcc = Tl * Phi * g_tilde * PhiT * Tu
				//                = Tl * Phi * tmp7           * Tu
				//                = Tl * g00                  * Tu
				//                = Tl * tmp8
				
				//                = Tl * Phi * g_tilde * tmp6
				//                = Tl * Phi * tmp7
								
				Matc tmp7(Nc*Nn,Nc*Nn);	
				mult(g_tilde, PhiT, tmp7); // tmp7 = g_tilde * PhiT;
				
				Matc g00(Nc*Nn,Nc*Nn);	// new
				mult(Phi, tmp7, g00); // g00 = Phi * tmp7;
				
				Matc tmp8(Nc*Nn,Nc*Nn);	
				mult(g00, Tu, tmp8); // tmp8 = g00 * Tu;
				
				Matc SRcc(Nc*Nn,Nc*Nn);	
				mult(Tl, tmp8, SRcc); // SRcc = Tl * tmp8;
								
				// assign the guys to the right place in the total retarded contact self-energy
				for (uint ii=1; ii<=Nc; ii++) 
				{
					uint gii3 = this->device_vertices[cc][ii-1]->get_index_internal() + 1;
					for (uint jj=1; jj<=Nc; jj++) 
					{
						uint gii4 = this->device_vertices[cc][jj-1]->get_index_internal() + 1;
						
						// ordering is always (xx-1)*Nn+nn in Hc, Mc!
						for (uint mm=1; mm<=Nn; mm++) {
							for (uint nn=1; nn<=Nn; nn++) {
#ifdef USE_BANDED
								// m!=n will probably be outside band
								if (fabs(get_mat_idx(gii3,mm,Nx)-get_mat_idx(gii4,nn,Nx)) > SR.num_offdiags+1e-8) continue;
#endif
								SR(get_mat_idx(gii3,mm,Nx),get_mat_idx(gii4,nn,Nx)) = SRcc((ii-1)*Nn+mm, (jj-1)*Nn+nn);
							}
						}
						
						// SR.fill_block(gii3, gii4, SRcc, ii, jj);
					}
				}
			}
		}
	}
	mpi->synchronize_processes();
);}


/** SL and SG are much easier than SR
 *  Make sure this is called DIRECTLY AFTER SR was calculated */
void SEContacts::calculate_lesser_greater()
{STACK_TRACE(
	for (uint ee2 = 0; ee2 < myNE; ee2++)
	{	
		uint ee = energies->get_global_index(ee2);
		logmsg->emit_noendl_all(LOG_INFO_L2, "p%d: SL_cont,SG_cont(E=%d,k=:)...  ",mpi->get_rank(),ee);
		for (uint kk = 0; kk < Nk; kk++)
		{	
			
			SEMat & SR = this->get_retarded(kk,ee);
			
			SEMat & SL = this->get_lesser(kk,ee);
			SEMat & SG = this->get_greater(kk,ee);
			
			for (uint cc=0; cc < xspace->get_num_contacts(); cc++) 
			{
				double f = this->get_fermi(cc, ee); // 1/(1+exp(...))

	            // if there is only 1 k-point, assume ballistic calculation and parabolic bands
	            // --> use Fermi-integral of order 0 instead of simple Fermi-Dirac factor
				if (Nk==1) f = this->get_F0(cc, ee);

				vector<uint> verts;
				this->get_internal_idx_of_vertices_near_contact(cc, verts);
				NEGF_ASSERT(verts.size()>0, "need device vertices.");
				for (uint ii = 0; ii < verts.size(); ii++) 
				{
					// FLENS indices start with 1
					// Block matrix on diagonal for vertex i starts at A((i-1)*Nn+1, (i-1)*Nn+1) and stops at A(i*Nn, i*Nn)
					uint xx = verts[ii] + 1;
					Matc SRmSAcc(Nn,Nn); SR.get_block(xx,xx, SRmSAcc, Nx);
					
					Matc SAcc(Nn,Nn);
					conjtrans(SRmSAcc, SAcc); // SAcc = conjugateTranspose(SRcc);
					SRmSAcc -= SAcc;
					
					Matc SLcc(Nn,Nn);
					mult(SRmSAcc, -f, SLcc);
					
					Matc SGcc(Nn,Nn);
					mult(SRmSAcc, 1.0-f, SGcc);
					
					// ATTENTION:SPIN DEGENERACY NOT HERE!!! OTHERWISE GG=GR-GA+GL WILL BE VIOLATED...
					SL.fill_block(xx,xx, SLcc, Nx); // SLcc =      -f * (SRcc - SAcc);
					SG.fill_block(xx,xx, SGcc, Nx); // SGcc = (1.0-f) * (SRcc - SAcc);
					
					// security check
					if (security_checking) {
						for (uint ii1 = 1; ii1 <= Nn; ii1++) {
							for (uint ii2 = 1; ii2 <= Nn; ii2++) {
								NEGF_ASSERT(SLcc(ii1,ii2).real() < constants::imag_err, "lesser self-energy should be pure imaginary.");
								NEGF_FASSERT(SLcc(ii1,ii2).imag() == SLcc(ii2,ii1).imag(), "lesser self-energy should be symmetric; instead %e <> %e",
										SLcc(ii1,ii2).imag(),SLcc(ii2,ii1).imag());
								NEGF_ASSERT(SGcc(ii1,ii2).real() < constants::imag_err, "greater self-energy should be pure imaginary.");
								NEGF_FASSERT(SGcc(ii1,ii2).imag() == SGcc(ii2,ii1).imag(), "greater self-energy should be symmetric; instead %e <> %e",
										SGcc(ii1,ii2).imag(),SGcc(ii2,ii1).imag());
							}
						}
					}
				}
			}
			
			// test for anti-Hermiticity
			if (security_checking) {
				SEMat tmp = SCMat_create(NxNn); conjtrans(SL, tmp); tmp += SL;
				double tmp_norm = negf_math::matrix_norm(tmp);
				NEGF_FASSERT(tmp_norm < 1e-10, "SEContacts: tmp failed (delta=%e) : kk=%d, ee=%d=%.4g (ee2=%d)", tmp_norm, kk, ee, energies->get_energy_from_global_idx(ee), ee2);
			}
		}
	}
	mpi->synchronize_processes();	
	this->calculated = true;
);}


double SEContacts::get_fermi(uint contact_idx, uint global_Eidx) const
{STACK_TRACE(
	NEGF_ASSERT(contact_idx<xspace->get_num_contacts(), "invalid index");
	double kT   = constants::convert_from_SI(units::energy, constants::SIkb * this->temperature);
	double E    = energies->get_energy_from_global_idx(global_Eidx);
	double mu   = xspace->get_contact(contact_idx)->get_bc_value(quantities::fermilevel);

    double fac = 1.0 / (1.0 + negf_math::exp((E - mu) / kT));
	return fac;
);}


double SEContacts::get_F0(uint contact_idx, uint global_Eidx) const
{STACK_TRACE(
    NEGF_ASSERT(Nk==1, "only Nk=1 should use this.");
    NEGF_ASSERT(Nn==1, "Nk=1 possible only for effective mass band model");
    NEGF_ASSERT(contact_idx<xspace->get_num_contacts(), "invalid index");
    double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
    double kT   = constants::convert_from_SI(units::energy, constants::SIkb * this->temperature);
    double E    = energies->get_energy_from_global_idx(global_Eidx);
    double mu   = xspace->get_contact(contact_idx)->get_bc_value(quantities::fermilevel);

    // use Fermi integral of order 0 instead of Fermi-Dirac factor
    // multiplication w/ spin degeneracy factor is ***not*** done here!
    double m = constants::convert_from_SI(units::mass, constants::SIm0 *
               0.067); // AAAAAAH NO NO NO NO TODO

    double fac =   m * kT / (constants::pi * hbar * hbar)
                 * negf_math::log(1+negf_math::exp((mu - E) / kT));
    return fac;

);}




void SEContacts::get_internal_idx_of_vertices_near_contact(uint contact_idx, vector<uint> & verts) const
{STACK_TRACE(
	verts.clear();
	if (xspace->get_dimension()==1) {
		// we expect a single vertex to exist which is connected to a contact vertex and not have internal index -1
		Vertex * v = 0;
		const vector<Vertex *> & cont_verts = xspace->get_contact(contact_idx)->get_contact_vertices();
		for (uint ii=0; ii < cont_verts.size(); ii++) {
			Vertex * v1 = cont_verts[ii];
			const vector<Edge *> edges_near_ii = xspace->get_edges_near(v1);
			for (uint jj=0; jj < edges_near_ii.size(); jj++) {
				Vertex * v2 = (edges_near_ii[jj]->get_lower_vertex()==v1) ? edges_near_ii[jj]->get_upper_vertex() : edges_near_ii[jj]->get_lower_vertex();
				if (v2->get_index_internal()!=-1) {
					NEGF_FASSERT(!v2->is_at_contact(), "vertex %d (x=%.3e) has internal index %d but is at contact.", 
							v2->get_index_global(), v2->get_coordinate(0), v2->get_index_internal());
					v = v2;
					break;
				}
			}
			if (v!=0) break;
		}
		uint xx = v->get_index_internal();
		verts.push_back(xx);
		return;
	}
	if (xspace->get_dimension()==2 || xspace->get_dimension()==3) {
		NEGF_EXCEPTION("NYI.");
	}
	NEGF_EXCEPTION("Strange dimensionality.");
);}


void SEContacts::assign_broadening(vector< vector< vector< vector<cplx> > > > & frey_broadening_)
{STACK_TRACE(
	// check correctness of array
	NEGF_ASSERT(frey_broadening_.size()==myNE, "wrong array size.");
	for (uint ee2=0; ee2<myNE; ee2++) {
		NEGF_ASSERT(frey_broadening_[ee2].size()==Nk, "wrong array size.");
		for (uint kk=0; kk<Nk; kk++) {
			NEGF_ASSERT(frey_broadening_[ee2][kk].size()==Nn, "wrong array size.");
			for (uint nn=0; nn<Nn; nn++) {
				NEGF_ASSERT(frey_broadening_[ee2][kk][nn].size()==2/*num_contacts*/, "wrong array size.");				
			}
		}
	}
	
	this->frey_broadening = frey_broadening_;
);}

