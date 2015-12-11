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
#include "SEPhotonSpontaneous.h"
using namespace negf;

#define INCLUDE_UPCOUPLING

SEPhotonSpontaneous::SEPhotonSpontaneous(const Overlap * ov_,
					const Geometry * xspace_, 
					const Kspace * kspace_, 
					const Energies * energies_, 
					const Options * options_,
					const GreenFunctions * gf_,
					const MaterialDatabase * db) throw (Exception *):
	SelfEnergy(xspace_,kspace_,energies_,options_, constants::odSPhot),
	ov(ov_),
	gf(gf_),
	scaling(1.0),
	complicated_retarded((options->exists("LuisierSRphot") && options->get("LuisierSRphot")==1) ? true : false),
	security_checking(false)
{STACK_TRACE(
	NEGF_ASSERT(ov!=NULL && xspace!=NULL && kspace!=NULL && energies!=NULL && options!=NULL && gf!=NULL, "null pointer encountered.");
	
	// ---------------------------------------------------------------------------------------------------------
	// set the minimum and maximum emission energies, chosen from the lowest band gap present in the structure
	// user has the possibility to specify the window in terms of the distance above/below this quantity
	// ---------------------------------------------------------------------------------------------------------
	this->hwmin = 1e100;
	this->hwmax = 0.0;
	double kane_EP = 0.0;
	double delta_hwmax_Egap = constants::convert_from_SI(units::energy, constants::SIec * constants::delta_hwmax_Egap);
	double delta_hwmin_Egap = constants::convert_from_SI(units::energy, constants::SIec * constants::delta_hwmin_Egap);
	if (options->exists("DeltaEgapHWmin")) {
		delta_hwmin_Egap = constants::convert_from_SI(units::energy, constants::SIec * options->get("DeltaEgapHWmin"));
	}
	if (options->exists("DeltaEgapHWmax")) {
		delta_hwmax_Egap = constants::convert_from_SI(units::energy, constants::SIec * options->get("DeltaEgapHWmax"));
	}
	
	double Egap_min = 1e100;
	for (uint ii=0; ii<xspace->get_num_regions(); ii++) 
	{
		const PropertyContainer<double> * mat = xspace->get_region(ii)->get_material();		
		double Ec = constants::convert_from_SI(units::energy, constants::SIec * TdkpInfoDesk::get_cbedge(mat, options->get("temperature"), db));
		double Ev = constants::convert_from_SI(units::energy, constants::SIec * mat->get("valence_band_edge"));
		double Egap = Ec - Ev;
		
		if (Egap < Egap_min) {
			Egap_min = Egap;
			hwmax = Egap_min + delta_hwmax_Egap;
			hwmin = Egap_min - delta_hwmin_Egap;
			kane_EP = constants::convert_from_SI(units::energy, constants::SIec * mat->get("optical_matrix_element"));
		}
	}
	
	// optical matrix element (Kane parameter P) is chosen from the LOWest bandgap material
	logmsg->emit(LOG_INFO, "Spontaneous emission happens in the window [%.3geV...%.3geV] (min. Eg=%.3geV).\n Kane parameter P=%.3geV", hwmin,hwmax,Egap_min,kane_EP);
	
	// ---------------------------------------------------------------------------------
	// set up the factor with which Mtilde*G is multiplied in the energy integration
	// ---------------------------------------------------------------------------------
	const double ec   = constants::convert_from_SI(units::charge, constants::SIec);
	const double hbar = constants::convert_from_SI(units::hbar, constants::SIhbar);
	const double eps0 = constants::convert_from_SI(units::dielectric, constants::SIeps0);
	const double c0   = constants::convert_from_SI(units::velocity, constants::SI_speed_of_light);
	const double m0   = constants::convert_from_SI(units::mass, constants::SIm0);
	// Chuang, p.370: m0/6*EP = Mb^2, Mb=p_cv; we need p_cv/m_0 = Mb/m0 = sqrt(m0/6*EP)/m0 = sqrt(EP/(6*m0))
	double my_P = negf_math::sqrt(kane_EP/(6.0*m0));
	this->prefactor   = ec*ec*my_P*my_P/(8.0*constants::pi*constants::pi * eps0 * hbar * c0*c0*c0);
	// prefactor is unitless
	
	/* ---------------------------------------------------------------------------------
	kpmethod=0/1: single-band effective mass
	kpmethod=2: 2-band effective mass
	kpmethod=3: 2-band effective mass, orthogonal basis
	kpmethod=4: zincblende kp 4x4
	kpmethod=6: zincblende kp 6x6
	kpmethod=8: zincblende kp 8x8
	kpmethod=16: wurtzite kp 6x6
	kpmethod=18: wurtzite kp 8x8 
	--------------------------------------------------------------------------------- */
	double kpmethod = options->get("kp_method");
	if (fabs(kpmethod-0.0)<1e-10 || fabs(kpmethod-1.0)<1e-10) {
		NEGF_EXCEPTION("Please turn off photon self-energy for single-band simulations");
	}
	
	// --------------------------------------------------------------------
	// Set up the Nn*Nn boolean matrix nonzero_bandcoupling 
	// nonzero_bandcoupling[m][n] tells us if the entry \Pi_mn is nonzero 
	// --------------------------------------------------------------------
#ifdef INCLUDE_UPCOUPLING
	logmsg->emit(LOG_INFO,"Upcoupling is INCLUDED in the photon emission.");
#else
	logmsg->emit(LOG_INFO,"Upcoupling is NOT INCLUDED in the photon emission.");
#endif
	vector< vector<bool> > nonzero_bandcoupling;
	nonzero_bandcoupling.resize(Nn);
	for (uint mm=0; mm<Nn; mm++) {
		nonzero_bandcoupling[mm].resize(Nn, false);
	}
	if (fabs(kpmethod-2.0)<1e-10 || fabs(kpmethod-3.0)<1e-10) {
#ifdef INCLUDE_UPCOUPLING
		nonzero_bandcoupling[0][1] = true;
#endif
		nonzero_bandcoupling[1][0] = true;
	} else if (fabs(kpmethod-4.0)<1e-10) {
		nonzero_bandcoupling[0][1] = true; nonzero_bandcoupling[0][2] = true; nonzero_bandcoupling[0][3] = true; nonzero_bandcoupling[0][4] = true;
		nonzero_bandcoupling[1][0] = true; nonzero_bandcoupling[2][0] = true; nonzero_bandcoupling[3][0] = true; nonzero_bandcoupling[4][0] = true;
	} else if (fabs(kpmethod-6.0)<1e-10) {
		nonzero_bandcoupling[0][1] = true; nonzero_bandcoupling[0][2] = true; nonzero_bandcoupling[0][3] = true; 
		nonzero_bandcoupling[0][4] = true; nonzero_bandcoupling[0][5] = true; nonzero_bandcoupling[0][6] = true;
		nonzero_bandcoupling[1][0] = true; nonzero_bandcoupling[2][0] = true; nonzero_bandcoupling[3][0] = true; 
		nonzero_bandcoupling[4][0] = true; nonzero_bandcoupling[5][0] = true; nonzero_bandcoupling[6][0] = true;
	} else if (fabs(kpmethod-8.0)<1e-10) {
		nonzero_bandcoupling[0][1] = true; nonzero_bandcoupling[0][2] = true; nonzero_bandcoupling[0][3] = true; 
		nonzero_bandcoupling[0][5] = true; nonzero_bandcoupling[0][6] = true; nonzero_bandcoupling[0][7] = true;
		nonzero_bandcoupling[1][0] = true; nonzero_bandcoupling[2][0] = true; nonzero_bandcoupling[3][0] = true; 
		nonzero_bandcoupling[5][0] = true; nonzero_bandcoupling[6][0] = true; nonzero_bandcoupling[7][0] = true;
	} else { NEGF_EXCEPTION("kp method nyi!"); }
#ifdef INCLUDE_UPCOUPLING
	// check that this is symmetric
	for (uint ii=0; ii<Nn; ii++) {
		for (uint jj=0; jj<ii; jj++) {
			NEGF_ASSERT(nonzero_bandcoupling[ii][jj]==nonzero_bandcoupling[jj][ii], "expected symmetric nonzero bandcoupling.");
		}
	}
#endif
	
	// ---------------------------------------------------------------------
	// Set up the Nn*Nn*Nn*Nn-tensor Gamma
	// \Gamma_mnop = \sum_lambda int_0^2\pi d\phi int_0^\pi \sin\theta d\theta \Pi_mo(\theta,\phi,\lambda) \Pi_np(\theta,\phi,\lambda)
	// ---------------------------------------------------------------------
	vector< vector< vector< vector<double> > > > Gamma;
	Gamma.resize(Nn);
	for (uint mm=0; mm<Nn; mm++) {
		Gamma[mm].resize(Nn);
		for (uint nn=0; nn<Nn; nn++) {
			Gamma[mm][nn].resize(Nn);
			for (uint oo=0; oo<Nn; oo++) {
				Gamma[mm][nn][oo].resize(Nn, 0.0);
			}
		}
	}
	if (fabs(kpmethod-2.0)<1e-10 || fabs(kpmethod-3.0)<1e-10) {
		// 2-band model: 
		//double tmp = 8.0*constants::pi; 		// 2*4pi, 2--> polarization
		//double tmp = constants::pi/2;   		// don't know anymore why I did this
		double tmp = 8.0*constants::pi / 3.0;	// 2*4pi, but only 1/3 because of direction
#ifdef INCLUDE_UPCOUPLING
		// only 4 nonzero entries (out of 16) - assumption: band 0 is CB, band 1 is VB
		Gamma[0][0][1][1] = tmp;
		Gamma[0][1][1][0] = tmp;
		Gamma[1][0][0][1] = tmp;
#endif
		Gamma[1][1][0][0] = tmp;
	} else if(fabs(kpmethod-4.0)<1e-10) {
		NEGF_EXCEPTION("KP 4x4 not yet implemented!");
	} else if(fabs(kpmethod-6.0)<1e-10) {
		// S = Sup + Sdown, NOT *0.5...
		double tmp = 8.0*constants::pi/3.0;
		Gamma[0][0][1][1] = tmp; Gamma[0][0][1][4] = tmp; Gamma[0][0][2][2] = tmp; Gamma[0][0][2][5] = tmp;  Gamma[0][0][3][3] = tmp; Gamma[0][0][3][6] = tmp; 
		Gamma[0][0][4][1] = tmp; Gamma[0][0][4][4] = tmp; Gamma[0][0][5][2] = tmp; Gamma[0][0][5][5] = tmp;  Gamma[0][0][6][3] = tmp; Gamma[0][0][6][6] = tmp; 
		
		Gamma[0][1][1][0] = tmp; Gamma[0][1][4][0] = tmp;
		Gamma[0][2][2][0] = tmp; Gamma[0][2][5][0] = tmp;
		Gamma[0][3][3][0] = tmp; Gamma[0][3][6][0] = tmp;
		Gamma[0][4][1][0] = tmp; Gamma[0][4][4][0] = tmp;
		Gamma[0][5][2][0] = tmp; Gamma[0][5][5][0] = tmp;
		Gamma[0][6][3][0] = tmp; Gamma[0][6][6][0] = tmp;
		
		Gamma[1][0][0][1] = tmp; Gamma[1][0][0][4] = tmp; Gamma[1][1][0][0] = tmp; Gamma[1][4][0][0] = tmp;
		Gamma[2][0][0][2] = tmp; Gamma[2][0][0][5] = tmp; Gamma[2][2][0][0] = tmp; Gamma[2][5][0][0] = tmp;
		Gamma[3][0][0][3] = tmp; Gamma[3][0][0][6] = tmp; Gamma[3][3][0][0] = tmp; Gamma[3][6][0][0] = tmp;
		Gamma[4][0][0][1] = tmp; Gamma[4][0][0][4] = tmp; Gamma[4][1][0][0] = tmp; Gamma[4][4][0][0] = tmp;
		Gamma[5][0][0][2] = tmp; Gamma[5][0][0][5] = tmp; Gamma[5][2][0][0] = tmp; Gamma[5][5][0][0] = tmp;
		Gamma[6][0][0][3] = tmp; Gamma[6][0][0][6] = tmp; Gamma[6][3][0][0] = tmp; Gamma[6][6][0][0] = tmp;
	} else if(fabs(kpmethod-8.0)<1e-10) {
		double tmp = 8.0*constants::pi/3.0;
		Gamma[0][0][1][1] = tmp; Gamma[0][0][2][2] = tmp; Gamma[0][0][3][3] = tmp;
		Gamma[0][1][1][0] = tmp;
		Gamma[0][2][2][0] = tmp;
		Gamma[0][3][3][0] = tmp;
		Gamma[0][4][1][5] = tmp; Gamma[0][4][2][6] = tmp; Gamma[0][4][3][7] = tmp;
		Gamma[0][5][1][4] = tmp;
		Gamma[0][6][2][4] = tmp;
		Gamma[0][7][3][4] = tmp;
		
		Gamma[1][0][0][1] = tmp; Gamma[1][1][0][0] = tmp; Gamma[1][4][0][5] = tmp; Gamma[1][5][0][4] = tmp;
		Gamma[2][0][0][2] = tmp; Gamma[2][2][0][0] = tmp; Gamma[2][4][0][6] = tmp; Gamma[2][6][0][4] = tmp;
		Gamma[3][0][0][3] = tmp; Gamma[3][3][0][0] = tmp; Gamma[3][4][0][7] = tmp; Gamma[3][7][0][4] = tmp;
		
		Gamma[4][0][5][1] = tmp; Gamma[4][0][6][2] = tmp; Gamma[4][0][7][3] = tmp;
		Gamma[4][1][5][0] = tmp;
		Gamma[4][2][6][0] = tmp;
		Gamma[4][3][7][0] = tmp;
		Gamma[4][4][5][5] = tmp; Gamma[4][4][6][6] = tmp; Gamma[4][4][7][7] = tmp;
		Gamma[4][5][5][4] = tmp;
		Gamma[4][6][6][4] = tmp;
		Gamma[4][7][7][4] = tmp;
		
		Gamma[5][0][4][1] = tmp; Gamma[5][1][4][0] = tmp; Gamma[5][4][4][5] = tmp; Gamma[5][5][4][4] = tmp;
		Gamma[6][0][4][2] = tmp; Gamma[6][2][4][0] = tmp; Gamma[6][4][4][6] = tmp; Gamma[6][6][4][4] = tmp;
		Gamma[7][0][4][3] = tmp; Gamma[7][3][4][0] = tmp; Gamma[7][4][4][7] = tmp; Gamma[7][7][4][4] = tmp;
	}
	// check Gamma[m][n][o][p] = Gamma[p][o][n][m] = Gamma[o][p][m][n] (Nn is small -> no big thing)
	logmsg->emit(LOG_INFO,"Checking structure of Gamma...");
	for (uint mm=0; mm<Nn; mm++) {
		for (uint nn=0; nn<Nn; nn++) {	
			for (uint oo=0; oo<Nn; oo++) {
				for (uint pp=0; pp<Nn; pp++) {
#ifdef INCLUDE_UPCOUPLING
					NEGF_ASSERT(fabs(Gamma[mm][nn][oo][pp]-Gamma[pp][oo][nn][mm]) < 1e-14, "expected Gamma[m][n][o][p] = Gamma[p][o][n][m].");
					NEGF_ASSERT(fabs(Gamma[mm][nn][oo][pp]-Gamma[oo][pp][mm][nn]) < 1e-14, "expected Gamma[m][n][o][p] = Gamma[o][p][m][n].");
#endif
					NEGF_ASSERT(fabs(Gamma[mm][nn][oo][pp]-Gamma[nn][mm][pp][oo]) < 1e-14, "expected Gamma[m][n][o][p] = Gamma[n][m][p][o]."); // for anti-herm.
				}
			}
		}
	}
	
	const Matd & M = ov->get_vertex_internal_overlap();
	NEGF_ASSERT(M.num_rows()==Nx && M.num_cols()==Nx, "Expected Nx * Nx matrix");
	
	// ---------------------------------------------------------------------
	// Set up Mt(a,b,c,d) = M_ik * M_lj * \Gamma_mnop
	// ... where each argument takes Nx*Nn values and a=(i,m), b=(j,n), c=(k,o), d=(l,p)
	// Storage: this quantity is (fortunately) very sparse, we use a
	// vector< vector< CSRMatrix<double> *> > such that
	// Mt(a,b,c,d) is stored in Mt[a][b]->get(c,d)
	// ---------------------------------------------------------------------
	logmsg->emit(LOG_INFO,"Computing M[a][b][c][d]...");
	this->Mt.resize(NxNn);
	// aa,bb,cc,dd (general matrix indices) are 0-based; ii,jj,kk,ll (position indices) and mm,nn,oo,pp (band indices) are 1-based!
	for (uint ii = 1; ii <= Nx; ii++) {
		for (uint mm = 1; mm <= Nn; mm++) {
			uint aa = get_mat_idx(ii,mm,Nx) - 1;
			NEGF_ASSERT(aa < Mt.size(), "aa<Mt.size() failed.");
			this->Mt[aa].resize(NxNn, NULL);
			
			for (uint jj = 1; jj <= Nx; jj++) {
				for (uint nn = 1; nn <= Nn; nn++) {
					uint bb = get_mat_idx(jj,nn,Nx) - 1;
					NEGF_ASSERT(bb < Mt[aa].size(), "bb<Mt[aa].size() failed.");
					this->Mt[aa][bb] = new CSRMatrix<double>(NxNn, negf::nonsymmetric_matrix); // false -> nonsymmetric
					
					// find the sparse structure of the matrix Mt[aa][bb]
					for (uint kk = 1; kk <= Nx; kk++) {
						if (abs(M(ii,kk)) < 1e-10) continue;
						for (uint oo = 1; oo <= Nn; oo++) {
							NEGF_ASSERT(mm-1 < nonzero_bandcoupling.size() && oo-1 < nonzero_bandcoupling[mm-1].size(), " invalid nonzero_bandcoupling.");
							if (!nonzero_bandcoupling[mm-1][oo-1]) continue;
							uint cc = get_mat_idx(kk,oo,Nx) - 1;
							
							for (uint ll = 1; ll <= Nx; ll++) {
								if (abs(M(ll,jj)) < 1e-10) continue;
								for (uint pp = 1; pp <= Nn; pp++) {
									//if (!nonzero_bandcoupling[pp-1][nn-1]) continue;
									if (!nonzero_bandcoupling[nn-1][pp-1]) continue;
									uint dd = get_mat_idx(ll,pp,Nx) - 1;
									
									this->Mt[aa][bb]->announce(cc,dd); // announce a nonzero entry
									// note that cc and dd start at zero and that's also what announce(.,.) expects
								}
							}
						}
					}
					this->Mt[aa][bb]->set_structure(); 	// finalize the matrix structure. afterwards: set, add, get
					
					// fill the matrix
					for (uint kk = 1; kk <= Nx; kk++) {
						if (abs(M(ii,kk)) < 1e-10) continue;
						for (uint oo = 1; oo <= Nn; oo++) {
							if (!nonzero_bandcoupling[mm-1][oo-1]) continue;
							uint cc = get_mat_idx(kk,oo,Nx) - 1;
							
							for (uint ll = 1 ; ll <= Nx; ll++) {
								if (abs(M(ll,jj)) < 1e-10) continue;
								for (uint pp = 1; pp <= Nn; pp++) {
#ifdef INCLUDE_UPCOUPLING
									if (!nonzero_bandcoupling[pp-1][nn-1]) continue;
#endif
									if (!nonzero_bandcoupling[nn-1][pp-1]) continue;
									uint dd = get_mat_idx(ll,pp,Nx) - 1;
									
									this->Mt[aa][bb]->set(cc, dd, M(ii,kk) * M(ll,jj) * Gamma[mm-1][nn-1][oo-1][pp-1]); 
									// note that cc and dd start at zero and that's also what set(.,.) expects
									// Gamma[][][][] also expects 0-based indices, but M(,) has 1-based indices
								}
							}
						}
					}
				}
			}
		}
	}
	mpi->synchronize_processes();
	
	// check Mt[a][b][c][d] = Mt[b][a][d][c] - necessary for anti-hermiticity of lesser and greater self-energies
	logmsg->emit(LOG_INFO,"Checking structure of Mt...");
	for (uint aa=0; aa<NxNn; aa++) 
	{
		for (uint bb=0; bb<NxNn; bb++) 
		{
			CSRMatrix<double> * Mt_ab = Mt[aa][bb];
			CSRMatrix<double> * Mt_ba = Mt[bb][aa];
			NEGF_ASSERT(Mt_ab->get_size()==Mt_ba->get_size(), "Mt_ab and Mt_ba should have the same dimension!");
			NEGF_ASSERT(Mt_ab->get_num_nonzeros()==Mt_ba->get_num_nonzeros(), "Mt_ab and Mt_ba should have the same #nonzeros!");
					
			int fidx = Mt_ab->get_index_base(); // the Fortran problem... linear solvers need index base 1! so fidx=1
			NEGF_ASSERT(fidx==Mt_ba->get_index_base(), "you're stupid.");
			
			for (uint row=0; row<Mt_ab->get_size(); row++) 		// row index, 0-based!
			{
				int    pos = Mt_ab->get_prow()[row];			// position in nonzeros array, fidx-based
				int    col = Mt_ab->get_icol()[pos-fidx];		// column index				 , fidx-based
				double val = Mt_ab->get_nonzeros()[pos-fidx];	// value
							
				// check if M[b][a][d][c] is announced as being nonzero
				int    pos2= Mt_ba->get_prow()[col-fidx];	
				int    pos3= Mt_ba->get_prow()[col-fidx+1];		// position of next row in nonzeros array, fidx-based
				bool announced = false;
				for (int pp=pos2-fidx; pp<pos3-fidx; pp++) {	// pp is 0-based, like row
					if (Mt_ba->get_icol()[pp]-fidx==int(row)) {
						announced = true;
						break;
					}
				}
				NEGF_FASSERT(announced, "entry M[a][b][c][d] was announced as nonzero but M[b][a][d][c] was not! a=%d,b=%d,c=%d,d=%d",aa,bb,row,col-fidx);
				
				// check for same value
				double val2 = Mt_ba->get(col-fidx,row);			// will NOT throw error if not announced; note that get(.,.) expects 0-based indices
				NEGF_FASSERT(fabs(val - val2) < 1e-14, "expected same value at M[a][b][c][d] and M[b][a][d][c], but %e<>%e!",val,val2);
			}	
		}
	}
	
	// --------------------------------------------------------
	// allocate memory for arrays used in zlib decompression
	// --------------------------------------------------------
	SEMat dummymat = SPhotMat_create(NxNn);
	vector<SEMat> dummyvec; dummyvec.resize(Nk, dummymat);
	uint num_doubles = negf_math::check_and_get_num_doubles(dummyvec, false); 
	// false --> not antihermitian matrices --> for antihermitian matrices there is certainly enough memory
	int doublesize = sizeof(double);
	//int charsize   = sizeof(char);
	unsigned long num_charsss = num_doubles * doublesize;
	this->decompress_array_real = new unsigned char[num_charsss];
	this->decompress_array_imag = new unsigned char[num_charsss];
	// assumption: compressed arrays will be smaller or almost equal to uncompressed arrays
	this->compressed_array_real = new unsigned char[num_charsss+10000];
	this->compressed_array_imag = new unsigned char[num_charsss+10000];

	// ------------------------------------------------------------------
	// set up which processes need which (will also be called later on)
	// ------------------------------------------------------------------
	this->determine_mpi_stuff();
);}
					
					
SEPhotonSpontaneous::~SEPhotonSpontaneous()
{STACK_TRACE(
	for (uint aa=0; aa<NxNn; aa++) {
		for (uint bb=0; bb<NxNn; bb++) {
			delete this->Mt[aa][bb];
		}
	}
	delete [] this->compressed_array_real;
	delete [] this->compressed_array_imag;
	delete [] this->decompress_array_real;
	delete [] this->decompress_array_imag;
);}


void SEPhotonSpontaneous::calculate() throw (Exception *)
{STACK_TRACE(
	this->determine_mpi_stuff();
	this->calculate_lesser_greater();
	if (!complicated_retarded) {	// otherwise the retarded self-energy was computed in calculate_lesser_greater
		this->calculate_retarded();
	}
);}


void SEPhotonSpontaneous::calculate_retarded()
{STACK_TRACE(
	logmsg->emit_small_header("calculating retarded spontaneous photon self-energy");
	NEGF_ASSERT(!complicated_retarded, "do not call this method when using complicated_retarded.");
	for (uint ee2=0; ee2<this->myNE; ee2++)
	{
		const uint ee = energies->get_global_index(ee2);
		for (uint kk=0; kk<Nk; kk++) 
		{
			const SEMat & SL = this->get_lesser(kk,ee);
			const SEMat & SG = this->get_greater(kk,ee);
			      SEMat & SR = this->get_retarded(kk,ee);
			
			sub(SG, SL, SR); // SR = SG - SL;
			SR *= 0.5;
		}
	}
	// no principal value!
);}


void SEPhotonSpontaneous::calculate_lesser_greater()
{STACK_TRACE(
	if (!complicated_retarded) {
		logmsg->emit_small_header("calculating lesser and greater spontaneous photon emission self-energy");
	} else {
		logmsg->emit_small_header("calculating SL,SG and Luisier-SR of spontaneous photon emission self-energy");
	}
	
	if (this->scaling==0.0) {	// shortcut...
		logmsg->emit(LOG_INFO,"Skipping detailed computation because scaling=0");
		for (uint ee2=0; ee2<this->myNE; ee2++) {
			uint ee=energies->get_global_index(ee2);
			for (uint kk=0; kk<Nk; kk++) {
				SEMat & SL = this->get_lesser(kk,ee);
				SEMat & SG = this->get_greater(kk,ee);
				SL = SEMat_create(NxNn);
				SG = SEMat_create(NxNn);
				if (complicated_retarded) {
					SEMat & SR = this->get_retarded(kk,ee);
					SR = SPhotMat_create(NxNn);
				}
			}
		}
		return;
	}
		
	for (uint aa=0; aa < this->AL.size(); aa++) {
		for (uint kk=0; kk<Nk; kk++) {
			(this->AL[aa][kk])(1,1) = 888.888;
			(this->AG[aa][kk])(1,1) = 888.888;
			if (complicated_retarded) {
				(this->AR[aa][kk])(1,1) = 888.888;
			}
		}
	}
	
	// --------------------------------------------
	// calculate own energies
	// --------------------------------------------
	logmsg->emit(LOG_INFO,"Own energies...");
	SEMat ALtmp = SPhotMat_create(NxNn);
	SEMat AGtmp = SPhotMat_create(NxNn);
	SEMat ARtmp = SPhotMat_create(NxNn);
	for (uint kk=0; kk<Nk; kk++)
	{
		for (uint ee2=0; ee2<this->myNE; ee2++)
		{
			uint ee=energies->get_global_index(ee2);
						
			const GLMat & GL = gf->get_lesser(kk,ee);
			const GLMat & GG = gf->get_greater(kk,ee);
			
			NEGF_ASSERT(AL.size() > ee-this->E0_idx, "invalid AL size");
			SEMat & ALmat = this->AL[ee-this->E0_idx][kk];
			SEMat & AGmat = this->AG[ee-this->E0_idx][kk];
			
			// re-init to zero
			ALmat = SPhotMat_create(NxNn);
			AGmat = SPhotMat_create(NxNn);
			
			// AL_ab[E][k] = \sum_cd M_abcd GL_cd
			// AG_ab[E][k] = \sum_cd M_abcd GG_cd
			for (int aa=1; aa<=int(NxNn); aa++) 
			{
				for (int bb=1; bb<=int(NxNn); bb++) 
				{
#ifdef USE_BANDED
					if (fabs(aa-bb) > ALmat.num_offdiags+1e-8) continue;
					if (fabs(aa-bb) > GG.num_offdiags+1e-8) continue; // make sure GG ist really stored banded - otherwise contributions would be lost with the preceding clause
#endif
					const CSRMatrix<double> * Mtilde = Mt[aa-1][bb-1];
					uint fidx = Mtilde->get_index_base(); // the Fortran problem (linear solvers need index base 1). fidx=1 by default
					
					double * nz = Mtilde->get_nonzeros();
					int *  icol = Mtilde->get_icol();  // array of length nonzeros, gives column index of element ii at nonzeros(ii) 
					int *  prow = Mtilde->get_prow();  // array of length size + 1, gives position of row ii in nonzeros 
					for (uint cc=0; cc<NxNn; cc++) {				 // cc is 0-based
						for (int ii=prow[cc]; ii<prow[cc+1]; ii++) { // ii is fidx-based
							int dd = icol[ii-fidx];			 		 // dd is fidx-based
							ALmat(aa,bb) += nz[ii-fidx] * GL(cc+1, dd-fidx+1);
							//AGmat(aa,bb) += nz[ii-fidx] * GG(cc+1, dd-fidx+1);	// old
#ifdef USE_BANDED
							if (fabs((cc+1)-(dd-fidx+1)) > AGmat.num_offdiags+1e-8) continue;
#endif
							AGmat(cc+1, dd-fidx+1) += nz[ii-fidx] * GG(aa, bb);		// new (see report): SG(c,d) = \sum_{a,b} Mtilde(c,d,a,b) GG(a,b)
						}
					}
				}
			}
			ALmat *= this->prefactor*this->scaling;
			AGmat *= this->prefactor*this->scaling;
			
			// security check: check anti-hermiticity
			if (security_checking) {
				conjtrans(ALmat, ALtmp); // ALtmp = conjugateTranspose(ALmat); 
				ALtmp += ALmat;
				conjtrans(AGmat, AGtmp); // AGtmp = conjugateTranspose(AGmat); 
				AGtmp += AGmat;
				double ALnrm = negf_math::matrix_norm(ALtmp);
				double AGnrm = negf_math::matrix_norm(AGtmp);
				NEGF_FASSERT(ALnrm < constants::antiherm_check, "AL was not anti-hermitian: |AL-(-AL+)|=%e, ee=%d=%.4g, kk=%d", ALnrm, ee, energies->get_energy_from_global_idx(ee), kk);
				NEGF_FASSERT(AGnrm < constants::antiherm_check, "AG was not anti-hermitian: |AG-(-AG+)|=%e, ee=%d=%.4g, kk=%d", AGnrm, ee, energies->get_energy_from_global_idx(ee), kk);
			}
			// compute norms
			this->AL_norm[ee-this->E0_idx][kk] = negf_math::matrix_norm(ALmat);
			this->AG_norm[ee-this->E0_idx][kk] = negf_math::matrix_norm(AGmat);
			
			if (complicated_retarded) 
			{				
				const Matc & GR = gf->get_retarded(kk,ee);
				SEMat & ARmat = this->AR[ee-this->E0_idx][kk];
				ARmat = SPhotMat_create(NxNn);
				for (int aa=1; aa<=int(NxNn); aa++) 
				{
					for (int bb=1; bb<=int(NxNn); bb++) 
					{
#ifdef USE_BANDED
						if (fabs(aa-bb) > ARmat.num_offdiags+1e-8) continue;
#endif
						const CSRMatrix<double> * Mtilde = Mt[aa-1][bb-1];
						uint fidx = Mtilde->get_index_base(); // the Fortran problem (linear solvers need index base 1). fidx=1 by default
						
						double * nz = Mtilde->get_nonzeros();
						int *  icol = Mtilde->get_icol();  // array of length nonzeros, gives column index of element ii at nonzeros(ii) 
						int *  prow = Mtilde->get_prow();  // array of length size + 1, gives position of row ii in nonzeros 
						for (uint cc=0; cc<NxNn; cc++) {				 // cc is 0-based
							for (int ii=prow[cc]; ii<prow[cc+1]; ii++) { // ii is fidx-based
								int dd = icol[ii-fidx];			 		 // dd is fidx-based
								ARmat(aa,bb) += nz[ii-fidx] * GR(cc+1, dd-fidx+1);
							}
						}
					}
				}
				ARmat *= this->prefactor*this->scaling;
				this->AR_norm[ee-this->E0_idx][kk] = negf_math::matrix_norm(ARmat);
			}
		}
	}
	/*double max_AL_norm = 0.0;
	double max_AG_norm = 0.0;
	double max_AR_norm = 0.0;
	for (uint kk=0; kk<Nk; kk++) {
		for (uint ee2=0; ee2<this->myNE; ee2++) {
			max_AL_norm = max(max_AL_norm, AL_norm[ee2][kk]);
			max_AG_norm = max(max_AG_norm, AG_norm[ee2][kk]);
			if (complicated_retarded) {
				max_AR_norm = max(max_AR_norm, AR_norm[ee2][kk]);
			}
		}
	}
	cout << "p" << mpi->get_rank() << ": max_AL=" << max_AL_norm << ", max_AG=" << max_AG_norm << ", max_AR=" << max_AR_norm << endl;
	*/
	logmsg->emit_noendl_all(LOG_INFO_L2,"p%d  ",mpi->get_rank());
	mpi->synchronize_processes();
	
	// -----------------------------------------------------------------
	// communicate AL's and add to self-energies at the same time!!!!!
	// -----------------------------------------------------------------
	this->communicate_As();
	logmsg->emit_noendl_all(LOG_INFO,"p%d  ",mpi->get_rank());
);}

 
void SEPhotonSpontaneous::determine_mpi_stuff()
{STACK_TRACE(
	logmsg->emit(LOG_INFO,"Determining MPI communication stuff for photons");
	const uint nE = energies->get_number_of_points();
	
	// MPI: set up the array of matrices GL, GG which is needed to calculate the own self-energies
	
	this->E0_idx = energies->get_start_global_idx(mpi->get_rank());
	double E0    = energies->get_energy_from_global_idx(E0_idx);
	this->E1_idx = energies->get_stop_global_idx(mpi->get_rank());
	double E1    = energies->get_energy_from_global_idx(E1_idx);
	
	// -----------------------------------------------------------------------------
	// compute the following quantities:
	// 1. E0_minus_hwmax_idx - the energy just below E0-hwmin, at least zero
	// 2. E1_minus_hwmin_idx - the energy just above E1-hwmax, at most E0
	// 3. E0_plus_hwmin_idx  - the energy just below E0+hwmin, at least E1
	// 4. E1_plus_hwmax_idx  - the energy just above E1+hwmax, at most the highest energy
	// 5. nE_below = E1_minus_hwmin_idx-E0_minus_hwmax_idx
	// 6. nE_above = E1_plus_hwmax_idx-E0_plus_hwmin_idx
	// -----------------------------------------------------------------------------
	
	// determine lowest needed global energy index BELOW interval
	this->E0_minus_hwmax_idx = E0_idx;
	while (E0_minus_hwmax_idx > 0 
		   && energies->get_energy_from_global_idx(E0_minus_hwmax_idx) > E0-hwmax) {
		E0_minus_hwmax_idx--;
	}
	// determine highest needed global energy index BELOW interval
	this->E1_minus_hwmin_idx = E0_idx;
	while (E1_minus_hwmin_idx > 0 
		   && energies->get_energy_from_global_idx(E1_minus_hwmin_idx) > E1-hwmin) {
		E1_minus_hwmin_idx--;
	}
	if (E1_minus_hwmin_idx<E0_idx) E1_minus_hwmin_idx++;
	// determine lowest needed global energy index ABOVE interval
	this->E0_plus_hwmin_idx = E1_idx;
	while (E0_plus_hwmin_idx+1 < nE 
		   && energies->get_energy_from_global_idx(E0_plus_hwmin_idx) < E0+hwmin) {
		E0_plus_hwmin_idx++;
	}
	if (E0_plus_hwmin_idx > E1_idx) E0_plus_hwmin_idx--;
	// determine highest needed global energy index ABOVE interval
	this->E1_plus_hwmax_idx = E1_idx;
	while (E1_plus_hwmax_idx+1 < nE 
		   && energies->get_energy_from_global_idx(E1_plus_hwmax_idx) < E1+hwmax) {
		E1_plus_hwmax_idx++;
	}
	logmsg->emit_all(LOG_INFO_L3,"p%d (energy indices %d...%d = %.3geV...%.3geV) needs indices %d(%.3geV)...%d(%.3geV) and %d(%.3geV)...%d(%.3geV)",
			mpi->get_rank(), E0_idx, E1_idx, E0, E1,
			E0_minus_hwmax_idx,energies->get_energy_from_global_idx(E0_minus_hwmax_idx),
			E1_minus_hwmin_idx,energies->get_energy_from_global_idx(E1_minus_hwmin_idx),
			E0_plus_hwmin_idx, energies->get_energy_from_global_idx(E0_plus_hwmin_idx),
			E1_plus_hwmax_idx, energies->get_energy_from_global_idx(E1_plus_hwmax_idx));
	
	this->nE_below = E1_minus_hwmin_idx-E0_minus_hwmax_idx+1;
	this->nE_above = E1_plus_hwmax_idx-E0_plus_hwmin_idx+1;
	this->nE_self  = E1_idx - E0_idx + 1;
	logmsg->emit_all(LOG_INFO_L2,"p%d: nE_self=%d, nE_above=%d, nE_below=%d",mpi->get_rank(), nE_self, nE_above, nE_below);
	
	// we want to get rid of the points E0_idx and hwmax_idx in the auxiliary arrays because we don't need them
	if (E1_minus_hwmin_idx==E0_idx) {
		nE_below--;
	}
	if (E0_plus_hwmin_idx==E1_idx) {
		nE_above--;
	}
	mpi->synchronize_processes();
		
	// allocate space; AL[...][kj] will store the lesser and greater matrices at point (...,kj)
	// remember: assign(..) = { clear(...); resize(...); }
	SEMat tmp = SPhotMat_create(NxNn);
	vector<SEMat> Atmp; 
	Atmp.resize(Nk, tmp);
	this->AL.assign(nE_self, Atmp);
	this->AG.assign(nE_self, Atmp);
	vector<double> Anormtmp; Anormtmp.resize(Nk,0.0);
	this->AL_norm.assign(nE_self, Anormtmp); 
	this->AG_norm.assign(nE_self, Anormtmp);
	if (complicated_retarded) {
		this->AR.assign(nE_self, Atmp);
		this->AR_norm.assign(nE_self, Anormtmp); 
	}
	
	// clean up heap if it was previously used
	for (uint ii=0; ii< this->AL_real_data.size(); ii++) {
		if (this->AL_real_data[ii]!=NULL) { delete [] this->AL_real_data[ii]; this->AL_real_data[ii] = NULL; }
		if (this->AG_real_data[ii]!=NULL) { delete [] this->AG_real_data[ii]; this->AG_real_data[ii] = NULL; }
		if (this->AL_imag_data[ii]!=NULL) { delete [] this->AL_imag_data[ii]; this->AL_imag_data[ii] = NULL; }
		if (this->AG_imag_data[ii]!=NULL) { delete [] this->AG_imag_data[ii]; this->AG_imag_data[ii] = NULL; }
		if (this->AL_real_compressed[ii]!=NULL) { delete [] this->AL_real_compressed[ii]; this->AL_real_compressed[ii] = NULL; }
		if (this->AG_real_compressed[ii]!=NULL) { delete [] this->AG_real_compressed[ii]; this->AG_real_compressed[ii] = NULL; }
		if (this->AL_imag_compressed[ii]!=NULL) { delete [] this->AL_imag_compressed[ii]; this->AL_imag_compressed[ii] = NULL; }
		if (this->AG_imag_compressed[ii]!=NULL) { delete [] this->AG_imag_compressed[ii]; this->AG_imag_compressed[ii] = NULL; }
		if (complicated_retarded) {
			if (this->AR_real_data[ii]!=NULL) { delete [] this->AR_real_data[ii]; this->AR_real_data[ii] = NULL; }
			if (this->AR_imag_data[ii]!=NULL) { delete [] this->AR_imag_data[ii]; this->AR_imag_data[ii] = NULL; }
			if (this->AR_real_compressed[ii]!=NULL) { delete [] this->AR_real_compressed[ii]; this->AR_real_compressed[ii] = NULL; }
			if (this->AR_imag_compressed[ii]!=NULL) { delete [] this->AR_imag_compressed[ii]; this->AR_imag_compressed[ii] = NULL; }
		}
	}
	// resize compressed data vectors for later use 
	this->AL_real_data.assign(nE_self, NULL);	// AL_real_data[ee2] will store copies of the complex AL[ee2] (for all k) split into real and imag parts
	this->AG_real_data.assign(nE_self, NULL);
	this->AL_imag_data.assign(nE_self, NULL);
	this->AG_imag_data.assign(nE_self, NULL);
	this->AL_real_char.assign(nE_self, NULL);	// will store pointer to same memory as AL_...._data, but different type
	this->AG_real_char.assign(nE_self, NULL);
	this->AL_imag_char.assign(nE_self, NULL);
	this->AG_imag_char.assign(nE_self, NULL);
	this->num_chars        .assign(nE_self, 0);	// will store how many chars the array of complex matrices corresponds to
	this->AL_real_comp_size.assign(nE_self, 0);	// will store the size of the zipped data
	this->AG_real_comp_size.assign(nE_self, 0);
	this->AL_imag_comp_size.assign(nE_self, 0);
	this->AG_imag_comp_size.assign(nE_self, 0);
	this->AL_real_compressed.assign(nE_self, NULL);	// will store pointers to the zipped data
	this->AG_real_compressed.assign(nE_self, NULL);
	this->AL_imag_compressed.assign(nE_self, NULL);
	this->AG_imag_compressed.assign(nE_self, NULL);
	if (complicated_retarded) {
		this->AR_real_data.assign(nE_self, NULL);
		this->AR_imag_data.assign(nE_self, NULL);
		this->AR_real_char.assign(nE_self, NULL);
		this->AR_imag_char.assign(nE_self, NULL);
		this->num_chars_AR     .assign(nE_self, 0); // AR is not antihermitian...
		this->AR_real_comp_size.assign(nE_self, 0);	
		this->AR_imag_comp_size.assign(nE_self, 0);
		this->AR_real_compressed.assign(nE_self, NULL);
		this->AR_imag_compressed.assign(nE_self, NULL);
	}
	
	mpi->synchronize_processes();
);}


void SEPhotonSpontaneous::determine_needed_processes(
		vector< vector<int> > & processes_needed,
		vector<int> & pp_E0_minus_hwmax_idx,
		vector<int> & pp_E1_minus_hwmin_idx,
		vector<int> & pp_E0_idx,
		vector<int> & pp_E1_idx,
		vector<int> & pp_E0_plus_hwmin_idx,
		vector<int> & pp_E1_plus_hwmax_idx) const
{STACK_TRACE(
	logmsg->emit_noendl(LOG_INFO,"Setting up processes_needed...   ");
	int root = constants::mpi_master_rank;
	int my_rank = mpi->get_rank();
	
	processes_needed.clear(); processes_needed.resize(mpi->get_num_procs());
	pp_E0_minus_hwmax_idx.assign(mpi->get_num_procs(), -1);
	pp_E1_minus_hwmin_idx.assign(mpi->get_num_procs(), -1);
	pp_E0_idx           .assign(mpi->get_num_procs(), -1);
	pp_E1_idx           .assign(mpi->get_num_procs(), -1);
	pp_E0_plus_hwmin_idx .assign(mpi->get_num_procs(), -1);
	pp_E1_plus_hwmax_idx .assign(mpi->get_num_procs(), -1);
	if (my_rank==root) 
	{
		for (int pp=0; pp < mpi->get_num_procs(); pp++) 
		{
			// process pp needs energies E0_minus_hwmax_idx...E0_idx-1 and hwmax_idx+1...E1_plus_hwmax_idx
			
			// get these variables
			if (pp==my_rank) {
				pp_E0_minus_hwmax_idx[pp]= this->E0_minus_hwmax_idx;
				pp_E1_minus_hwmin_idx[pp]= this->E1_minus_hwmin_idx;
				pp_E0_idx[pp]            = this->E0_idx;
				pp_E1_idx[pp]            = this->E1_idx;
				pp_E0_plus_hwmin_idx[pp] = this->E0_plus_hwmin_idx;
				pp_E1_plus_hwmax_idx[pp] = this->E1_plus_hwmax_idx;
			} else {
				int source = pp;
				int tag = 1; mpi->recv(pp_E0_minus_hwmax_idx[pp], source, tag);
				    tag = 2; mpi->recv(pp_E1_minus_hwmin_idx[pp], source, tag);
				    tag = 3; mpi->recv(pp_E0_idx[pp]           , source, tag);
				    tag = 4; mpi->recv(pp_E1_idx[pp]           , source, tag);
				    tag = 5; mpi->recv(pp_E0_plus_hwmin_idx[pp] , source, tag);
				    tag = 6; mpi->recv(pp_E1_plus_hwmax_idx[pp] , source, tag);
			}
			
			// determine which processes compute the energies in question
			// and add them to the list
			for (int ee=pp_E0_minus_hwmax_idx[pp]; ee <= pp_E1_minus_hwmin_idx[pp]; ee++) 
			{
				int pp2 = energies->get_process_computing(uint(ee));
				if (pp2==pp) continue; // this is possible, ee could be E0_idx!
				bool pp2_already_in_list = false;
				for (uint ii=0; ii<processes_needed[pp].size(); ii++) {
					if (processes_needed[pp][ii]==pp2) {
						pp2_already_in_list = true;
						break;
					}
				}
				if (!pp2_already_in_list) {
					processes_needed[pp].push_back(pp2);
				}
			}
			for (int ee=pp_E0_plus_hwmin_idx[pp]; ee <= pp_E1_plus_hwmax_idx[pp]; ee++) 
			{
				int pp2 = energies->get_process_computing(uint(ee));
				if (pp2==pp) continue; // this is possible, ee could be E1_idx!
				bool pp2_already_in_list = false;
				for (uint ii=0; ii<processes_needed[pp].size(); ii++) {
					if (processes_needed[pp][ii]==pp2) {
						pp2_already_in_list = true;
						break;
					}
				}
				if (!pp2_already_in_list) {
					processes_needed[pp].push_back(pp2);
				}
			}
		}
		// screen output
		for (int pp=0; pp < mpi->get_num_procs(); pp++) {
			logmsg->emit_noendl_all(LOG_INFO_L3,"p%d needs ",pp);
			for (uint ii=0; ii < processes_needed[pp].size(); ii++) {
				logmsg->emit_noendl_all(LOG_INFO_L3,"p%d, ",processes_needed[pp][ii]);
			}
			logmsg->emit_all(LOG_INFO_L3,"");
		}
	} else {
		int E0_idx_copy           = this->E0_idx;
		int E1_idx_copy           = this->E1_idx;
		int E0_minus_hwmax_idx_copy= this->E0_minus_hwmax_idx;
		int E1_minus_hwmin_idx_copy= this->E1_minus_hwmin_idx;
		int E0_plus_hwmin_idx_copy = this->E0_plus_hwmin_idx;
		int E1_plus_hwmax_idx_copy = this->E1_plus_hwmax_idx;
		int tag = 1; mpi->send(E0_minus_hwmax_idx_copy, root, tag);
			tag = 2; mpi->send(E1_minus_hwmin_idx_copy, root, tag);
			tag = 3; mpi->send(E0_idx_copy           , root, tag);
			tag = 4; mpi->send(E1_idx_copy           , root, tag);
			tag = 5; mpi->send(E0_plus_hwmin_idx_copy , root, tag);
			tag = 6; mpi->send(E1_plus_hwmax_idx_copy , root, tag);
	}
	
	// -----------------------------------------------------
	// broadcast processed_needed and all the rest
	// -----------------------------------------------------
	logmsg->emit_noendl(LOG_INFO,"Broadcasting processes_needed...   ");
	for (int pp=0; pp<mpi->get_num_procs(); pp++) {
		int num_processes_needed = 0;
		if (my_rank==root) {
			num_processes_needed = processes_needed[pp].size();
		}
		mpi->broadcast(num_processes_needed, root);
		if (my_rank!=root) {
			processes_needed[pp].resize(num_processes_needed);
		}
		mpi->broadcast(processes_needed[pp], root);
	}
	logmsg->emit_noendl(LOG_INFO,"Broadcasting energy index arrays...   ");
	mpi->broadcast(pp_E0_minus_hwmax_idx, root);
	mpi->broadcast(pp_E1_minus_hwmin_idx, root);
	mpi->broadcast(pp_E0_idx           , root);
	mpi->broadcast(pp_E1_idx           , root);
	mpi->broadcast(pp_E0_plus_hwmin_idx , root);
	mpi->broadcast(pp_E1_plus_hwmax_idx , root);
	
	mpi->synchronize_processes();
);}


void SEPhotonSpontaneous::communicate_As()
{STACK_TRACE(
	logmsg->emit(LOG_INFO,"communicating As:");
	int root = constants::mpi_master_rank;
	int my_rank = mpi->get_rank();
	
	// ---------------------------
	// clear own self-energies
	// ---------------------------
	for (uint kk=0; kk<Nk; kk++) {
		for (uint ee2=0; ee2<this->myNE; ee2++) {
			uint ee=energies->get_global_index(ee2);
			SEMat & SL = this->get_lesser(kk,ee);
			SEMat & SG = this->get_greater(kk,ee);
			
			SL = SPhotMat_create(NxNn);
			SG = SPhotMat_create(NxNn);
			
			if (complicated_retarded) {
				SEMat & SR = this->get_retarded(kk,ee);
				SR = SPhotMat_create(NxNn);
			}
		}
	}
	
	// --------------------------------------------------------------------------------------------------------
	// set up an array which will store for every process all other processes it relies on (EXcluding itself)
	// computed only by master process
	// --------------------------------------------------------------------------------------------------------
	vector< vector<int> > processes_needed;
	vector<int> pp_E0_minus_hwmax_idx;
	vector<int> pp_E1_minus_hwmin_idx;
	vector<int> pp_E0_idx;
	vector<int> pp_E1_idx;
	vector<int> pp_E0_plus_hwmin_idx;
	vector<int> pp_E1_plus_hwmax_idx;
	this->determine_needed_processes(processes_needed, pp_E0_minus_hwmax_idx, pp_E1_minus_hwmin_idx, pp_E0_idx, pp_E1_idx, pp_E0_plus_hwmin_idx, pp_E1_plus_hwmax_idx);
		
	// ----------------------------------------------------------------
	// create an array of zipped data (for every energy all k's)
	// ----------------------------------------------------------------
	logmsg->emit_noendl(LOG_INFO,"Creating zipped data arrays...   ");
	NEGF_ASSERT(AL.size()==energies->get_my_number_of_points() && AL.size()==this->nE_self && this->myNE==this->nE_self, "inconsistent myNE.");
	NEGF_ASSERT(this->AL_real_data.size()==this->myNE, "resize AL_real_data first!");
	for (uint ee2=0; ee2<this->myNE; ee2++) 
	{
		NEGF_ASSERT(AL_real_data[ee2]==NULL, "clean up heap memory in AL_real_data first!");
		
		// compress AL[ee] and AG[ee]	(anti-hermitian)
		negf_math::do_compress_antiherm(AL[ee2], AL_real_data[ee2], AL_imag_data[ee2], AL_real_char[ee2], AL_imag_char[ee2], num_chars[ee2], 
					AL_real_comp_size[ee2], AL_imag_comp_size[ee2], AL_real_compressed[ee2], AL_imag_compressed[ee2]);
		unsigned long tmp = num_chars[ee2];
		negf_math::do_compress_antiherm(AG[ee2], AG_real_data[ee2], AG_imag_data[ee2], AG_real_char[ee2], AG_imag_char[ee2], num_chars[ee2], 
					AG_real_comp_size[ee2], AG_imag_comp_size[ee2], AG_real_compressed[ee2], AG_imag_compressed[ee2]);
		NEGF_ASSERT(num_chars[ee2]==tmp, "something went wrong.");
		
		if (complicated_retarded) {
			// compress AR[ee] (NOT anti-hermitian --> different num_chars!)
			negf_math::do_compress(AR[ee2], AR_real_data[ee2], AR_imag_data[ee2], AR_real_char[ee2], AR_imag_char[ee2], num_chars_AR[ee2], 
					AR_real_comp_size[ee2], AR_imag_comp_size[ee2], AR_real_compressed[ee2], AR_imag_compressed[ee2]);
		}
	}
	mpi->synchronize_processes();
	
	// --------------------------------------------------------------------------------------
	// determine how much data is sent in total
	// we do this by going through the same algorithm as when sending the data, just w/o sending
	// necessary for the amount of buffer that needs to be allocated for a safe operation
	// --------------------------------------------------------------------------------------
	logmsg->emit_noendl(LOG_INFO,"Creating total amount of data to be sent...   ");
	NEGF_FASSERT(this->E1_idx-this->E0_idx+1==this->myNE,"E1_idx=%d, E0_idx=%d, myNE=%d",this->E1_idx,this->E0_idx,this->myNE);
	long buffer_size_needed = 0;
	vector<bool> process_was_computed; 
	process_was_computed.resize(mpi->get_num_procs(), false); 
	while (true)
	{
		vector<bool> receiver; receiver.resize(mpi->get_num_procs(), false);
		vector<bool> sender;   sender.resize(mpi->get_num_procs(), false);
				
		// -------------------------------------
		// set up arrays w/ receiver and sender
		// -------------------------------------
		mpi->determine_senders_receivers(processes_needed, process_was_computed, sender, receiver);
		
		// -------------------------------------------
		// we're done if there are no more receivers
		// -------------------------------------------
		uint num_receivers = 0;
		uint num_senders = 0;
		for (int pp=0; pp < mpi->get_num_procs(); pp++) {
			if (receiver[pp]) num_receivers++;
			if (sender[pp]) num_senders++;
		}
		
		// no screen output! later on...
		if (num_receivers == 0) {
			break;
		}
		
		bool i_am_receiver = receiver[my_rank];
		bool i_am_sender   = sender[my_rank];
		
		if (i_am_sender) {
			NEGF_FASSERT(!i_am_receiver, "p%d is both sender and receiver!",my_rank);
		}
		
		if (i_am_receiver) {
			// is not of our concern right now
		} else 
		{
			if (!i_am_sender) { // this is possible!
			} else {
			
			// PART 1: sender sends his AL's to processes storing energies BELOW his energy interval
			for (uint ee=this->E0_idx; ee<=this->E1_idx; ee++) 
			{			
				// determine receiver(s) BELOW current process - if there is none, continue to next energy
				for (int pp=0; pp<mpi->get_num_procs(); pp++) 
				{
					bool pp_receives = receiver[pp] && pp_E0_plus_hwmin_idx[pp]<=int(ee) && pp_E1_plus_hwmax_idx[pp]>=int(ee);
					if (!pp_receives) continue;
					
					// so the current process needs to send AL[ee-this->E0_idx] plus the size information! (note: 4 Bsends)
					buffer_size_needed += AL_real_comp_size[ee-this->E0_idx] + AL_imag_comp_size[ee-this->E0_idx] + 2*sizeof(unsigned long) + 8*MPI_BSEND_OVERHEAD;
				}			
			}
			
			// PART 2: sender sends his AG's to processes storing energies ABOVE his energy interval
			for (uint ee=this->E0_idx; ee<=this->E1_idx; ee++) 
			{			
				// determine receiver ABOVE current process - if there is none, continue to next energy
				for (int pp=0; pp<mpi->get_num_procs(); pp++) 
				{
					bool pp_receives = receiver[pp] && pp_E0_minus_hwmax_idx[pp]<=int(ee) && pp_E1_minus_hwmin_idx[pp]>=int(ee);
					if (!pp_receives) continue;
					
					// so the current process needs to send AG[ee-this->E0_idx] plus the size information !
					buffer_size_needed += AG_real_comp_size[ee-this->E0_idx] + AG_imag_comp_size[ee-this->E0_idx] + 2*sizeof(unsigned long) + 8*MPI_BSEND_OVERHEAD;
					if (complicated_retarded) {
						buffer_size_needed += AL_real_comp_size[ee-this->E0_idx] + AL_imag_comp_size[ee-this->E0_idx] + 2*sizeof(unsigned long) + 8*MPI_BSEND_OVERHEAD;
						buffer_size_needed += AR_real_comp_size[ee-this->E0_idx] + AR_imag_comp_size[ee-this->E0_idx] + 2*sizeof(unsigned long) + 8*MPI_BSEND_OVERHEAD;
					}
				}	
			}
			} // if(i_am_sender)
		}
	
		// mark receivers as computed
		// needs to be performed in ALL threads (senders, receivers and those which are neither)
		for (int pp=0; pp < mpi->get_num_procs(); pp++) {
			if (receiver[pp]) {
				process_was_computed[pp] = true;
			}
		}
	} // while(true)
	logmsg->emit_all(LOG_INFO_L3,"p%d will need to send %d chars in total.", mpi->get_rank(), buffer_size_needed);
	mpi->synchronize_processes();
	
	// ------------------------------
	// allocate buffer!
	// ------------------------------
	logmsg->emit_noendl(LOG_INFO,"Allocating MPI buffer...   ");
	unsigned long buffersize_long = buffer_size_needed+1000;
	NEGF_ASSERT(buffersize_long < 2147483648, "Buffer size does not fit into an int!!@!");
	int buffersize = int(buffersize_long);
	char * buffer = new char[buffersize];
	int err = MPI_Buffer_attach(buffer, buffersize);
	NEGF_FASSERT(err==0, "MPI_Buffer_attach gave error %d",err);
	mpi->synchronize_processes();
		
	// ------------------------------
	// do MPI sending!
	// ------------------------------
	logmsg->emit(LOG_INFO,"MPI Buffered Send!");
	process_was_computed.clear();
	process_was_computed.resize(mpi->get_num_procs(), false);
	while (true)
	{
		vector<bool> receiver; receiver.resize(mpi->get_num_procs(), false);
		vector<bool> sender;   sender.resize(mpi->get_num_procs(), false);
		
		// -------------------------------------
		// set up arrays w/ receiver and sender
		// -------------------------------------
		mpi->determine_senders_receivers(processes_needed, process_was_computed, sender, receiver);
		
		// -------------------------------------------
		// we're done if there are no more receivers
		// -------------------------------------------
		uint num_receivers = 0;
		uint num_senders = 0;
		for (int pp=0; pp < mpi->get_num_procs(); pp++) {
			if (receiver[pp]) num_receivers++;
			if (sender[pp]) num_senders++;
		}
		
		if (num_receivers == 0) {
			break;
		}
		//mpi->synchronize_processes();
		
		// some screen output
		if (my_rank==root) {
			logmsg->emit_noendl_all(LOG_INFO, "This time we have %d receivers ",num_receivers);
			for (int pp=0; pp < mpi->get_num_procs(); pp++) {
				if (receiver[pp]) logmsg->emit_noendl_all(LOG_INFO_L2, "%d ", pp);
			}
			logmsg->emit_noendl_all(LOG_INFO, " and %d senders ",num_senders);
			for (int pp=0; pp < mpi->get_num_procs(); pp++) {
				if (sender[pp]) logmsg->emit_noendl_all(LOG_INFO_L2, "%d ", pp);
			}
			logmsg->emit_all(LOG_INFO_L2, "");
		}
				
		// ---------------------------------------------------------------------------
		// receiver processes receive all their needed energies from sender processes
		// ---------------------------------------------------------------------------
		bool i_am_receiver = receiver[my_rank];
		bool i_am_sender   = sender[my_rank];
		
		if (i_am_sender) {
			NEGF_FASSERT(!i_am_receiver, "p%d is both sender and receiver!",my_rank);
		}
		
		if (i_am_receiver)
		{
			// receiving will be done later!
		} else 
		{
			if (!i_am_sender) { // this is possible!
			} else {
			
			// PART 1: sender sends his AL's to processes storing energies BELOW his energy interval
			for (uint ee=this->E0_idx; ee<=this->E1_idx; ee++) 
			{			
				// determine receiver(s) BELOW current process - if there is none, continue to next energy
				for (int pp=0; pp<mpi->get_num_procs(); pp++) 
				{
					bool pp_receives = receiver[pp] && pp_E0_plus_hwmin_idx[pp]<=int(ee) && pp_E1_plus_hwmax_idx[pp]>=int(ee);
					if (!pp_receives) continue;
					int receiver_id = pp;
					
					int tag = ee;
					int tag2 = tag+1;
				
					logmsg->emit(LOG_INFO_L3,"p%d sends ee=%d downwards to p%d", my_rank, ee, receiver_id);
					// BUFFERED SEND!
					// mpi->send_antihermitians and mpi->recv_antihermitians send/receive the following 
					// quantities in order: real_char_size, imag_char_size, real_compressed, imag_compressed
					
					mpi->send(AL_norm[ee-this->E0_idx], receiver_id, tag);

					// if all AL is essentially zero for all k-points, skip 
					double ALnrm_max = 0.0;
					for (uint kk=0; kk<Nk; kk++) {
						ALnrm_max = max(ALnrm_max, AL_norm[ee-this->E0_idx][kk]);
					}
					if (ALnrm_max < constants::ALGnorm_neglect) {
						continue;
					}
				
					// send array lengths
					int real_char_size = int(AL_real_comp_size[ee-this->E0_idx]);
					int imag_char_size = int(AL_imag_comp_size[ee-this->E0_idx]);
					err = MPI_Bsend(&real_char_size, 1, MPI_INT, receiver_id, tag, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of real_char_size failed (1): err=%d",err);
					err = MPI_Bsend(&imag_char_size, 1, MPI_INT, receiver_id, tag2, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of imag_char_size failed (1): err=%d",err);
		
					// send arrays
					err = MPI_Bsend(AL_real_compressed[ee-this->E0_idx], real_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of real_compressed failed (1): err=%d",err);
					err = MPI_Bsend(AL_imag_compressed[ee-this->E0_idx], imag_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag2, MPI_COMM_WORLD);
					NEGF_FASSERT(err==0, "MPI send of imag_compressed failed (1): err=%d",err);
				}
			}
			//mpi->synchronize_processes();
			
			// PART 2: sender sends his AG's to processes storing energies ABOVE his energy interval
			for (uint ee=this->E0_idx; ee<=this->E1_idx; ee++) 
			{			
				// determine receiver ABOVE current process - if there is none, continue to next energy
				for (int pp=0; pp<mpi->get_num_procs(); pp++) 
				{
					bool pp_receives = receiver[pp] && pp_E0_minus_hwmax_idx[pp]<=int(ee) && pp_E1_minus_hwmin_idx[pp]>=int(ee);
					if (!pp_receives) continue;
					int receiver_id = pp;
					
					int tag = ee;
					int tag2 = ee+1;
				
					logmsg->emit(LOG_INFO_L3,"p%d sends ee=%d upwards to p%d", my_rank, ee, receiver_id);
					// BUFFERED SEND!
					// mpi->send_antihermitians and mpi->recv_antihermitians send/receive the following 
					// quantities in order: real_char_size, imag_char_size, real_compressed, imag_compressed
					
					mpi->send(AG_norm[ee-this->E0_idx], receiver_id, tag);
					
					// if all AG is essentially zero for all k-points, skip 
					double AGnrm_max = 0.0;
					for (uint kk=0; kk<Nk; kk++) {
						AGnrm_max = max(AGnrm_max, AG_norm[ee-this->E0_idx][kk]);
					}
					if (!(AGnrm_max < constants::ALGnorm_neglect)) {
						// send array lengths
						int real_char_size = int(AG_real_comp_size[ee-this->E0_idx]);
						int imag_char_size = int(AG_imag_comp_size[ee-this->E0_idx]);
						err = MPI_Bsend(&real_char_size, 1, MPI_INT, receiver_id, tag, MPI_COMM_WORLD);
						NEGF_FASSERT(err==0, "MPI send of (2) real_char_size failed: err=%d",err);
						err = MPI_Bsend(&imag_char_size, 1, MPI_INT, receiver_id, tag2, MPI_COMM_WORLD);
						NEGF_FASSERT(err==0, "MPI send of (2) imag_char_size failed: err=%d",err);
			
						// send arrays
						err = MPI_Bsend(AG_real_compressed[ee-this->E0_idx], real_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag, MPI_COMM_WORLD);
						NEGF_FASSERT(err==0, "MPI send of (2) real_compressed failed: err=%d",err);
						err = MPI_Bsend(AG_imag_compressed[ee-this->E0_idx], imag_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag2, MPI_COMM_WORLD);
						NEGF_FASSERT(err==0, "MPI send of (2) imag_compressed failed: err=%d",err);
					}
					
					if (complicated_retarded) 
					{
						// -----------
						// send AL
						// -----------
						mpi->send(AL_norm[ee-this->E0_idx], receiver_id, tag);

						// send array lengths
						int real_char_size = int(AL_real_comp_size[ee-this->E0_idx]);
						int imag_char_size = int(AL_imag_comp_size[ee-this->E0_idx]);
						err = MPI_Bsend(&real_char_size, 1, MPI_INT, receiver_id, tag, MPI_COMM_WORLD);
						NEGF_FASSERT(err==0, "MPI send of (2) real_char_size failed: err=%d",err);
						err = MPI_Bsend(&imag_char_size, 1, MPI_INT, receiver_id, tag2, MPI_COMM_WORLD);
						NEGF_FASSERT(err==0, "MPI send of (2) imag_char_size failed: err=%d",err);
			
						// send arrays
						err = MPI_Bsend(AL_real_compressed[ee-this->E0_idx], real_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag, MPI_COMM_WORLD);
						NEGF_FASSERT(err==0, "MPI send of (2) real_compressed failed: err=%d",err);
						err = MPI_Bsend(AL_imag_compressed[ee-this->E0_idx], imag_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag2, MPI_COMM_WORLD);
						NEGF_FASSERT(err==0, "MPI send of (2) imag_compressed failed: err=%d",err);
						
						// ---------
						// send AR
						// ---------
						mpi->send(AR_norm[ee-this->E0_idx], receiver_id, tag);

						// send array lengths
						real_char_size = int(AR_real_comp_size[ee-this->E0_idx]);
						imag_char_size = int(AR_imag_comp_size[ee-this->E0_idx]);
						err = MPI_Bsend(&real_char_size, 1, MPI_INT, receiver_id, tag, MPI_COMM_WORLD);
						NEGF_FASSERT(err==0, "MPI send of real_char_size failed (1): err=%d",err);
						err = MPI_Bsend(&imag_char_size, 1, MPI_INT, receiver_id, tag2, MPI_COMM_WORLD);
						NEGF_FASSERT(err==0, "MPI send of imag_char_size failed (1): err=%d",err);
			
						// send arrays
						err = MPI_Bsend(AR_real_compressed[ee-this->E0_idx], real_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag, MPI_COMM_WORLD);
						NEGF_FASSERT(err==0, "MPI send of real_compressed failed (1): err=%d",err);
						err = MPI_Bsend(AR_imag_compressed[ee-this->E0_idx], imag_char_size, MPI_UNSIGNED_CHAR, receiver_id, tag2, MPI_COMM_WORLD);
						NEGF_FASSERT(err==0, "MPI send of imag_compressed failed (1): err=%d",err);
					}
					
					logmsg->emit(LOG_INFO_L3,"p%d just sent ee=%d upwards to p%d", my_rank, ee, receiver_id);
				}
			}
			} // if(i_am_sender)
		}
	
		// mark receivers as computed
		// needs to be performed in ALL threads (senders, receivers and those which are neither)
		for (int pp=0; pp < mpi->get_num_procs(); pp++) {
			if (receiver[pp]) {
				process_was_computed[pp] = true;
			}
		}
		//mpi->synchronize_processes();
	} // while(true)
			
	// security check
	for (int pp=0; pp < mpi->get_num_procs(); pp++) {
		NEGF_ASSERT(process_was_computed[pp], "a process was not computed.");
	}
	mpi->synchronize_processes();
	logmsg->emit(LOG_INFO,"");

	// ------------------------------
	// do MPI receiving!
	// ------------------------------
	logmsg->emit(LOG_INFO,"Receive and add to self-energy!");
	process_was_computed.clear();
	process_was_computed.resize(mpi->get_num_procs(), false);
	while (true)
	{
		vector<bool> receiver; receiver.resize(mpi->get_num_procs(), false);
		vector<bool> sender;   sender.resize(mpi->get_num_procs(), false);
		
		// -------------------------------------
		// set up arrays w/ receiver and sender
		// -------------------------------------
		mpi->determine_senders_receivers(processes_needed, process_was_computed, sender, receiver);
		
		// -------------------------------------------
		// we're done if there are no more receivers
		// -------------------------------------------
		uint num_receivers = 0;
		uint num_senders = 0;
		for (int pp=0; pp < mpi->get_num_procs(); pp++) {
			if (receiver[pp]) num_receivers++;
			if (sender[pp]) num_senders++;
		}
		
		if (num_receivers == 0) {
			break;
		}
		//mpi->synchronize_processes();
		
		// some screen output
		if (my_rank==root) {
			logmsg->emit_noendl_all(LOG_INFO, "This time we have %d receivers ",num_receivers);
			for (int pp=0; pp < mpi->get_num_procs(); pp++) {
				if (receiver[pp]) logmsg->emit_noendl_all(LOG_INFO_L2, "%d ", pp);
			}
			logmsg->emit_noendl_all(LOG_INFO, " and %d senders ",num_senders);
			for (int pp=0; pp < mpi->get_num_procs(); pp++) {
				if (sender[pp]) logmsg->emit_noendl_all(LOG_INFO_L2, "%d ", pp);
			}
			logmsg->emit_all(LOG_INFO_L2, "");
		}
				
		// ---------------------------------------------------------------------------
		// receiver processes receive all their needed energies from sender processes
		// ---------------------------------------------------------------------------
		bool i_am_receiver = receiver[my_rank];
		bool i_am_sender   = sender[my_rank];
		
		if (i_am_sender) {
			NEGF_FASSERT(!i_am_receiver, "p%d is both sender and receiver!",my_rank);
		}
		
		if (i_am_receiver)
		{
			SEMat tmp = SPhotMat_create(NxNn);
			vector<SEMat> ALmat; ALmat.resize(Nk, tmp);
			vector<SEMat> AGmat; AGmat.resize(Nk, tmp);
			vector<SEMat> ARmat; if (complicated_retarded) { ARmat.resize(Nk, tmp); }
			
			vector<double> ALnrm; ALnrm.resize(Nk, 0.0);
			vector<double> AGnrm; AGnrm.resize(Nk, 0.0);
			vector<double> ARnrm; ARnrm.resize(Nk, 0.0);
			uint Nk_copy = Nk;

			double mpi_time = 0.0;
			double comp_time = 0.0;
			// PART 1: receiver receives missing AL's ABOVE his energy interval
			for (uint ee=this->E0_plus_hwmin_idx; ee<=this->E1_plus_hwmax_idx; ee++) 
			{
				double  E = energies->get_energy_from_global_idx(ee);
				double dE = energies->get_weight_from_global_idx(ee);
				
				int sender_id = energies->get_process_computing(ee);
				if (sender_id==my_rank) {
					// this is possible, ee could be E1_idx!
					continue;
				}
				NEGF_FASSERT(sender[sender_id] && !receiver[sender_id],"p%d expected sender %d but something went wrong.", my_rank, sender_id);
				int tag = ee;
				
				double t1 = MPI_Wtime();
				logmsg->emit_all(LOG_INFO_L3,"p%d waits for ee=%d from p%d", my_rank, ee, sender_id);
				mpi->recv(ALnrm, Nk_copy, sender_id, tag);	
				
				// if all AL is essentially zero for all k-points, skip 
				double ALnrm_max = 0.0;
				for (uint kk=0; kk<Nk; kk++) {
					ALnrm_max = max(ALnrm_max, ALnrm[kk]);
				}
				if (ALnrm_max < constants::ALGnorm_neglect) continue;
				
				mpi->recv_antihermitians(ALmat, sender_id, tag,
						this->compressed_array_real, this->compressed_array_imag,
						this->decompress_array_real, this->decompress_array_imag); // get stuff from an energy above --> inscattering --> AL, GL
				logmsg->emit_all(LOG_INFO_L3,"p%d just got ee=%d from p%d", my_rank, ee, sender_id);
				double t2 = MPI_Wtime();
				
				// ----------------------------------------------------------------------------
				// ADD TO ALL RELEVANT ENERGIES!
				// ----------------------------------------------------------------------------
				for (uint ee2=0; ee2<this->myNE; ee2++)
				{
					uint   ee_own = energies->get_global_index(ee2);
					double E_own  = energies->get_energy_from_global_idx(ee_own);
					double hw = E - E_own;
					NEGF_ASSERT(hw>0.0, "positive hw expected.");
					
					if (hw<this->hwmin || hw>this->hwmax) {
						continue;
					}
		
					for (uint kk=0; kk<Nk; kk++) {
						// skip if contribution is too small
						if (ALnrm[kk] < constants::ALGnorm_neglect) continue;

						// INscattering --> SL
						SEMat & SL = this->get_lesser(kk,ee_own);
						add(ALmat[kk], dE/(2.0*constants::pi) * hw, SL); // SL += dE/(2.0*constants::pi) * hw * ALmat[kk];
						
						if (complicated_retarded) {
							SEMat & SR = this->get_retarded(kk,ee_own);
							add(ALmat[kk], -0.5 * dE/(2.0*constants::pi) * hw, SR);	// note the -0.5!
						}
					}
				}
				double t3 = MPI_Wtime();
				mpi_time += t2-t1;
				comp_time += t3-t2;
			}
			//mpi->synchronize_processes();
			
			// PART 2: receiver receives missing AG's BELOW his energy interval
			for (uint ee=this->E0_minus_hwmax_idx; ee<=this->E1_minus_hwmin_idx; ee++) 
			{
				double  E = energies->get_energy_from_global_idx(ee);
				double dE = energies->get_weight_from_global_idx(ee);
				
				int sender_id = energies->get_process_computing(ee);
				if (sender_id==my_rank) {
					// this is possible, ee could be E0_idx!
					continue;
				}
				NEGF_FASSERT(sender[sender_id] && !receiver[sender_id],"p%d expected sender %d but something went wrong.", my_rank, sender_id);
				int tag = ee;
								
				double t1 = MPI_Wtime();
				logmsg->emit_all(LOG_INFO_L3,"p%d waits for ee=%d from p%d", my_rank, ee, sender_id);
				mpi->recv(AGnrm, Nk_copy, sender_id, tag);	
				
				// if all AG is essentially zero for all k-points, skip 
				double AGnrm_max = 0.0;
				for (uint kk=0; kk<Nk; kk++) {
					AGnrm_max = max(AGnrm_max, AGnrm[kk]);
				}
				if (!(AGnrm_max < constants::ALGnorm_neglect)) {				
					mpi->recv_antihermitians(AGmat, sender_id, tag,
						this->compressed_array_real, this->compressed_array_imag,
						this->decompress_array_real, this->decompress_array_imag);// get stuff from an energy below --> outscattering --> AG, GG
				}
				
				if (complicated_retarded) {
					mpi->recv(ALnrm, Nk_copy, sender_id, tag);	
					mpi->recv_antihermitians(ALmat, sender_id, tag);// AL is needed for SR
					mpi->recv(ARnrm, Nk_copy, sender_id, tag);	
					mpi->recv(ARmat, sender_id, tag);
				}
				logmsg->emit_all(LOG_INFO_L3,"p%d just got ee=%d from p%d", my_rank, ee, sender_id);
				double t2 = MPI_Wtime();
				
				// ----------------------------------------------------------------------------
				// ADD TO ALL RELEVANT ENERGIES!
				// ----------------------------------------------------------------------------
				for (uint ee2=0; ee2<this->myNE; ee2++)
				{
					uint   ee_own = energies->get_global_index(ee2);
					double E_own  = energies->get_energy_from_global_idx(ee_own);
					double hw = E_own - E;
					NEGF_ASSERT(hw>0.0, "positive hw expected.");
					
					if (hw<this->hwmin || hw>this->hwmax) {
						continue;
					}
		
					for (uint kk=0; kk<Nk; kk++) 
					{
						// OUTscattering --> SG
						SEMat & SG = this->get_greater(kk,ee_own);
						// skip if contribution is too small
						if (!(AGnrm[kk] < constants::ALGnorm_neglect)) {
							add(AGmat[kk], dE/(2.0*constants::pi) * hw, SG); // SG += dE/(2.0*constants::pi) * hw * AGmat[kk];
						}

						if (complicated_retarded) {
							SEMat & SR = this->get_retarded(kk,ee_own);
							if (!(ALnrm[kk] < constants::ALGnorm_neglect)) {
								add(ALmat[kk], 0.5 * dE/(2.0*constants::pi) * hw, SR); // AL, NOT AG! note the +0.5 !
							}
							if (!(ARnrm[kk] < constants::ALGnorm_neglect)) {
								add(ARmat[kk],       dE/(2.0*constants::pi) * hw, SR);
							}
						}
					}
				}
				double t3 = MPI_Wtime();
				mpi_time += t2-t1;
				comp_time += t3-t2;
				
			}
			logmsg->emit_noendl_all(LOG_INFO_L1,"recv p%d: t_MPI=%.3e, t_comput=%.3e           ",mpi->get_rank(),mpi_time,comp_time);

		} else 
		{
			// sending was already done
		}
	
		// mark receivers as computed
		// needs to be performed in ALL threads (senders, receivers and those which are neither)
		for (int pp=0; pp < mpi->get_num_procs(); pp++) {
			if (receiver[pp]) {
				process_was_computed[pp] = true;
			}
		}
		//mpi->synchronize_processes();
	} // while(true)

	// security check
	for (int pp=0; pp < mpi->get_num_procs(); pp++) {
		NEGF_ASSERT(process_was_computed[pp], "a process was not computed.");
	}

	// ----------------------------------------------------------------
	// release buffer
	// ----------------------------------------------------------------
	logmsg->emit(LOG_INFO,"Deallocating MPI buffer");
	err = MPI_Buffer_detach(&buffer, &buffersize);
	delete [] buffer;
	NEGF_FASSERT(err==0, "MPI_Buffer_attach gave error %d",err);
	
	// ----------------------------------------------------------------
	// DO NOT release zipped data! will be released only when 
	// determine_mpi_stuff is computed (which is before each recomputation)
	// the stuff will be used in the spectrum calculation!
	// ----------------------------------------------------------------
	logmsg->emit_noendl_all(LOG_INFO_L2,"p%d   ",mpi->get_rank());
	mpi->synchronize_processes();
);}


void SEPhotonSpontaneous::set_scaling(double new_scaling) throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(new_scaling>=0.0 && new_scaling<=1.0, "."); 
	logmsg->emit(LOG_INFO,"Setting scaling of spontaneous photon scattering to %.4g", new_scaling);
	this->scaling = new_scaling;
);}


const vector<double> & SEPhotonSpontaneous::get_AL_norm(uint ee2) const
{STACK_TRACE(
	NEGF_ASSERT(this->AL_norm.size()>ee2, "ee2 is invalid.");
	return this->AL_norm[ee2];
);}

const vector<double> & SEPhotonSpontaneous::get_AG_norm(uint ee2) const 
{STACK_TRACE(
	NEGF_ASSERT(this->AG_norm.size()>ee2, "ee2 is invalid.");
	return this->AG_norm[ee2];
);}

void SEPhotonSpontaneous::output_debug_info()
{STACK_TRACE(
	logmsg->emit_header("Debug info about GL, GG and SigmaL, SigmaG of el-photon interaction");
	
	// -----------------------------------------
	// compute stuff for own energies
	// -----------------------------------------
	logmsg->emit(LOG_INFO,"Computing...");
	vector<double> my_GL_neg; my_GL_neg.resize(4, 0.0);
	vector<double> my_GL_pos; my_GL_pos.resize(4, 0.0);
	vector<double> my_GG_neg; my_GG_neg.resize(4, 0.0);
	vector<double> my_GG_pos; my_GG_pos.resize(4, 0.0);
	vector<double> my_SL_neg; my_SL_neg.resize(4, 0.0);
	vector<double> my_SL_pos; my_SL_pos.resize(4, 0.0);
	vector<double> my_SG_neg; my_SG_neg.resize(4, 0.0);
	vector<double> my_SG_pos; my_SG_pos.resize(4, 0.0);
	vector<double> my_SLGG_neg; my_SLGG_neg.resize(4, 0.0);
	vector<double> my_SLGG_pos; my_SLGG_pos.resize(4, 0.0);
	vector<double> my_SGGL_neg; my_SGGL_neg.resize(4, 0.0);
	vector<double> my_SGGL_pos; my_SGGL_pos.resize(4, 0.0);
	const double Ethresh = 0.03;
	Matc SLGG(NxNn,NxNn);
	Matc SGGL(NxNn,NxNn);
	for (uint ee2 = 0; ee2 < myNE; ee2++) 
	{
		uint ee = energies->get_global_index(ee2);
		double E = energies->get_energy_from_global_idx(ee);
		for (uint kk=0; kk<Nk; kk++) 
		{
			vector<double> tmp; tmp.resize(4, 0.0);
			
#ifdef USE_BANDED
			Matc GLm(NxNn,NxNn); GLm = gf->get_lesser(kk,ee);
#else
			const Matc & GLm = gf->get_lesser(kk,ee);
#endif
			this->get_CB_VB_norms(GLm, tmp); 
			if (E<Ethresh) { for(uint ii=0; ii<4; ii++) { my_GL_neg[ii] += tmp[ii]; } } else { for(uint ii=0; ii<4; ii++) { my_GL_pos[ii] += tmp[ii]; } }
				
#ifdef USE_BANDED
			Matc GGm(NxNn,NxNn); GGm = gf->get_greater(kk,ee);
#else
			const Matc & GGm = gf->get_greater(kk,ee);
#endif
			this->get_CB_VB_norms(GGm, tmp); 
			if (E<Ethresh) { for(uint ii=0; ii<4; ii++) { my_GG_neg[ii] += tmp[ii]; } } else { for(uint ii=0; ii<4; ii++) { my_GG_pos[ii] += tmp[ii]; } }
			
#ifdef USE_BANDED
			Matc SLm(NxNn,NxNn); SLm = this->get_lesser(kk,ee);
#else
			const Matc & SLm = this->get_lesser(kk,ee);
#endif
			this->get_CB_VB_norms(SLm, tmp); 
			if (E<Ethresh) { for(uint ii=0; ii<4; ii++) { my_SL_neg[ii] += tmp[ii]; } } else { for(uint ii=0; ii<4; ii++) { my_SL_pos[ii] += tmp[ii]; } }
			
#ifdef USE_BANDED
			Matc SGm(NxNn,NxNn); SGm = this->get_greater(kk,ee);
#else
			const Matc & SGm = this->get_greater(kk,ee);
#endif
			this->get_CB_VB_norms(SGm, tmp); 
			if (E<Ethresh) { for(uint ii=0; ii<4; ii++) { my_SG_neg[ii] += tmp[ii]; } } else { for(uint ii=0; ii<4; ii++) { my_SG_pos[ii] += tmp[ii]; } }
			
			mult(SLm, GGm, SLGG); // SLGG = SLm * GGm;
			this->get_CB_VB_norms(SLGG, tmp); 
			if (E<Ethresh) { for(uint ii=0; ii<4; ii++) { my_SLGG_neg[ii] += tmp[ii]; } } else { for(uint ii=0; ii<4; ii++) { my_SLGG_pos[ii] += tmp[ii]; } }
			
			mult(SGm, GLm, SGGL); // SGGL = SGm * GLm;
			this->get_CB_VB_norms(SGGL, tmp); 
			if (E<Ethresh) { for(uint ii=0; ii<4; ii++) { my_SGGL_neg[ii] += tmp[ii]; } } else { for(uint ii=0; ii<4; ii++) { my_SGGL_pos[ii] += tmp[ii]; } }
		}
	}
	mpi->synchronize_processes();
	
	// ------------------------------------------------
	// communicate to master process
	// ------------------------------------------------
	logmsg->emit(LOG_INFO,"Aggregating in master thread...");
	if (mpi->get_rank()==constants::mpi_master_rank) 
	{
		vector<double> total_GL_neg = my_GL_neg;
		vector<double> total_GL_pos = my_GL_pos;
		vector<double> total_GG_neg = my_GG_neg;
		vector<double> total_GG_pos = my_GG_pos;
		vector<double> total_SL_neg = my_SL_neg;
		vector<double> total_SL_pos = my_SL_pos;
		vector<double> total_SG_neg = my_SG_neg;
		vector<double> total_SG_pos = my_SG_pos;
		vector<double> total_SLGG_neg = my_SLGG_neg;
		vector<double> total_SLGG_pos = my_SLGG_pos;
		vector<double> total_SGGL_neg = my_SGGL_neg;
		vector<double> total_SGGL_pos = my_SGGL_pos;
					
		// collect the pieces
		vector<double> tmp; tmp.resize(4, 0.0);
		for (int pp=0; pp<mpi->get_num_procs(); pp++) 
		{
			if (pp==constants::mpi_master_rank) continue;
			
			// receive from other process
			logmsg->emit/*_all*/(LOG_INFO_L3,"Receiving from process %d...",pp);
			int tag = pp;
			uint size = 4;
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_GL_neg[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_GL_pos[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_GG_neg[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_GG_pos[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_SL_neg[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_SL_pos[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_SG_neg[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_SG_pos[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_SLGG_neg[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_SLGG_pos[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_SGGL_neg[ii] += tmp[ii]; }
			mpi->recv(tmp, size, pp, tag); for(uint ii=0; ii<4; ii++) { total_SGGL_pos[ii] += tmp[ii]; }
			logmsg->emit/*_all*/(LOG_INFO_L3,"Receiving done.");
		}
		
		// screen output!
		logmsg->emit(LOG_INFO, "MATRIX NORMS (sum over E,k); E-range splits at %.3g",Ethresh);
		logmsg->emit(LOG_INFO, "quantity   E-range  CB-CB       CB-VB       VB-CB       VB-VB");
		logmsg->emit(LOG_INFO, "------------------------------------------------------------------");
		logmsg->emit(LOG_INFO, "GL         VB       %.3e   %.3e   %.3e   %.3e", total_GL_neg[0], total_GL_neg[1], total_GL_neg[2], total_GL_neg[3]);
		logmsg->emit(LOG_INFO, "GL         CB       %.3e   %.3e   %.3e   %.3e", total_GL_pos[0], total_GL_pos[1], total_GL_pos[2], total_GL_pos[3]);
		logmsg->emit(LOG_INFO, "GG         VB       %.3e   %.3e   %.3e   %.3e", total_GG_neg[0], total_GG_neg[1], total_GG_neg[2], total_GG_neg[3]);
		logmsg->emit(LOG_INFO, "GG         CB       %.3e   %.3e   %.3e   %.3e", total_GG_pos[0], total_GG_pos[1], total_GG_pos[2], total_GG_pos[3]);
		logmsg->emit(LOG_INFO, "SL         VB       %.3e   %.3e   %.3e   %.3e", total_SL_neg[0], total_SL_neg[1], total_SL_neg[2], total_SL_neg[3]);
		logmsg->emit(LOG_INFO, "SL         CB       %.3e   %.3e   %.3e   %.3e", total_SL_pos[0], total_SL_pos[1], total_SL_pos[2], total_SL_pos[3]);
		logmsg->emit(LOG_INFO, "SG         VB       %.3e   %.3e   %.3e   %.3e", total_SG_neg[0], total_SG_neg[1], total_SG_neg[2], total_SG_neg[3]);
		logmsg->emit(LOG_INFO, "SG         CB       %.3e   %.3e   %.3e   %.3e", total_SG_pos[0], total_SG_pos[1], total_SG_pos[2], total_SG_pos[3]);
		logmsg->emit(LOG_INFO, "SL*GG      VB       %.3e   %.3e   %.3e   %.3e", total_SLGG_neg[0], total_SLGG_neg[1], total_SLGG_neg[2], total_SLGG_neg[3]);
		logmsg->emit(LOG_INFO, "SL*GG      CB       %.3e   %.3e   %.3e   %.3e", total_SLGG_pos[0], total_SLGG_pos[1], total_SLGG_pos[2], total_SLGG_pos[3]);
		logmsg->emit(LOG_INFO, "SG*GL      VB       %.3e   %.3e   %.3e   %.3e", total_SGGL_neg[0], total_SGGL_neg[1], total_SGGL_neg[2], total_SGGL_neg[3]);
		logmsg->emit(LOG_INFO, "SG*GL      CB       %.3e   %.3e   %.3e   %.3e", total_SGGL_pos[0], total_SGGL_pos[1], total_SGGL_pos[2], total_SGGL_pos[3]);
	} else {
		// send to master process
		int dest = constants::mpi_master_rank;
		int tag = mpi->get_rank();
		mpi->send(my_GL_neg, dest, tag);
		mpi->send(my_GL_pos, dest, tag);
		mpi->send(my_GG_neg, dest, tag);
		mpi->send(my_GG_pos, dest, tag);
		mpi->send(my_SL_neg, dest, tag);
		mpi->send(my_SL_pos, dest, tag);
		mpi->send(my_SG_neg, dest, tag);
		mpi->send(my_SG_pos, dest, tag);
		mpi->send(my_SLGG_neg, dest, tag);
		mpi->send(my_SLGG_pos, dest, tag);
		mpi->send(my_SGGL_neg, dest, tag);
		mpi->send(my_SGGL_pos, dest, tag);
	}
	mpi->synchronize_processes();
);}


void SEPhotonSpontaneous::get_CB_VB_norms(const Matc & A, vector<double> & result)
{STACK_TRACE(
	NEGF_ASSERT(A.num_rows()==NxNn && A.num_cols()==NxNn && Nn==2, "expected A to be (2*Nx)^2 matrix");
	
	result.assign(4, 0.0);
	for (uint xx = 1; xx <= Nx; xx++)
	{
		for (uint yy = 1; yy <= Nx; yy++)
		{
			cplx CBCB = A((xx-1)*Nn+1, (yy-1)*Nn+1);
			cplx CBVB = A((xx-1)*Nn+1, (yy-1)*Nn+2);
			cplx VBCB = A((xx-1)*Nn+2, (yy-1)*Nn+1);
			cplx VBVB = A((xx-1)*Nn+2, (yy-1)*Nn+2);
			result[0] += std::abs(CBCB*CBCB);
			result[1] += std::abs(CBVB*CBVB);
			result[2] += std::abs(VBCB*VBCB);
			result[3] += std::abs(VBVB*VBVB);
		}
	}
	result[0] = negf_math::sqrt(result[0]);
	result[1] = negf_math::sqrt(result[1]);
	result[2] = negf_math::sqrt(result[2]);
	result[3] = negf_math::sqrt(result[3]);
);}
