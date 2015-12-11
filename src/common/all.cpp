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
#include "all.h"

#include "zlib.h"

//#ifndef NOACML
//	#include "acml_mv.h"
//#endif

// Handle underscore problem of Lapack/Blas routines
#ifdef CSCS
	#define F77_Name(name) name##_
#else 	
	#define F77_Name(name) name
#endif

#ifndef NOACML

// ACML version of zgeev
extern "C" { void zgeev( // ACML does not have an underscore
	char jobvl, // 'N': left eigenvectors of A are not computed;  'V': left eigenvectors are computed.
	char jobvr, // 'N': right eigenvectors of A are not computed; 'V': right eigenvectors are computed.
	int n, 		// order of the matrix A.
	cplx *a, 	// (in/out) array of dimension (LDA,N).  On entry, the matrix A. On exit, A has been overwritten.
	int lda, 	// The leading dimension of the array A.  LDA >= max(1,N).
	cplx * w, 	// (out) will contain the n computed eigenvalues
	cplx *vl, 	// (out) dimension (LDVL,N). If JOBVL = 'N', VL is not referenced.
				// If JOBVL = 'V', the left eigenvectors u(j) are stored one after another in the columns of VL, in the same order
				// as their eigenvalues: u(j) = VL(:,j)
	int ldvl,   // The leading dimension of the array VL.  LDVL >= 1; if JOBVL = 'V', LDVL >= N.
	cplx *vr, 	// same as vl, but with right eigenvalues
	int ldvr, 	// The leading dimension of the array VR.  LDVR >= 1; if JOBVR = 'V', LDVR >= N.
	int *info); // pointer to integer; = 0: success. if INFO = -i, the i-th argument had an illegal value.
}

#else

// "direct" LAPACK version (note the pointers and the work arrays)
// declaration on http://www.netlib.org/clapack/CLAPACK-3.1.1/SRC/zgeev.c
//            and http://www-heller.harvard.edu/people/shaw/programs/zgeev.h:
extern "C" { void F77_Name(zgeev)(char * jobvl, char * jobvr, int * n, cplx * a, int * lda, cplx * w, cplx * vl, int  * ldvl, cplx * vr, int  * ldvr, 
	     cplx * work,  // (workspace/output) array of dimension (max(1,lwork)). On exit, if info=0, work[0] returns the optimal LWORK.   
         int  * lwork, // (input) The dimension of the array work.  LWORK >= max(1,2*N). For good performance, LWORK must generally be larger.   
                       // If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns   
                       // this value as the first entry of the WORK array, and no error message related to LWORK is issued by XERBLA.   
         double * rwork, // (workspace) array of dimension (2*N)
         int  * info);
}

#endif

namespace negf {
	
	// -----------------------------------
	// global variables
	// -----------------------------------
	
	Logger    * logmsg    = new Logger(LOG_INFO_L2);
	Timer     * timer     = new Timer();
	Filenames * fnames    = new Filenames();
	Interrupt * interrupt = new Interrupt();
	myMPI     * mpi       = new myMPI();
	
	Matc work10(1,1); // will be resized to Nn*Nx in constants::prepare_Nx_stuff
	Matc work11(1,1);

	// ------------------------------------
	// constants (descriptions see all.h)
	// ------------------------------------	
		
	uint Nn   = 0;
	//uint Nx   = 0;
	//uint NxNn = 0;
	
	// number of offdiagonals above and below the diagonal (total num. offdiagonals -> *2)
	const uint 		  constants::odSE    			 = /*10*/7;		// this should be the max. of all included scattering mechanisms
	const uint 		  constants::odGL    			 = /*10*/7;		// (dito)
	const uint 		  constants::odOV    			 = 2;		
	const uint 		  constants::odSC    			 = 3;		
	const uint 		  constants::odSAC   			 = 3;		
	const uint 		  constants::odSPOP  			 = 0;       // POP selfenergy w/ offdiagonals does not conserve current as of now
	const uint 		  constants::odSGol  			 = 2;
	const uint 		  constants::odSPhot 			 = /*10*/7;
	const uint 		  constants::odSII 				 = 5;
	
	// simulator constants
	const double	  constants::antiherm_check      = 1e-8;
	const double	  constants::egrid_weight_fL     = 2.0;		// left weight should be higher because it is the injector (NOT FOR PN DIODE)
	const int         constants::mpi_compresslevel   = 9;		// MPI sending takes much more time than decompression!
	const double	  constants::egrid_weight_fR     = 2.0;
	const double	  constants::egrid_weight_res    = 3.0;		// combined weight for all resonances
	//const double	  constants::egrid_weight_fL     = 0.0;	
	//const double	  constants::egrid_weight_fR     = 0.0;
	//const double	  constants::egrid_weight_res    = 0.0;	
	const uint 		  constants::max_inner_iters     = 400;
	const uint 		  constants::max_inner_iters_start=5;
	//const double 	  constants::inner_errcrit       = 1e-3;
	//const double 	  constants::inner_errcrit       = 1e-6;	// SEE CMD.CNF NOW!
	const double 	  constants::inner_ramp_errcrit  = 1e-3;
	const double 	  constants::refdens             = 1e6*1e16;// in m^3
	const double 	  constants::refcurr             = 1e4*1e0;	// in A/m^2
	const uint 		  constants::max_outer_iters     = 20;
	const double 	  constants::outer_errcrit       = /*3e-5*/2e-4; // nanohub...
	const double 	  constants::imag_err            = 1e-8;
	const double      constants::dyson_eta           = 1e-8;	// in eV
	const double      constants::contact_eta         = 1e-3;	// in eV
	const double 	  constants::contact_min_xdecay  = /*0.2*/0.00001;		// in nm (was 1nm) (recently was 0.2 nm - RTD example demanded lower value)
	const bool        constants::old_orthogonal      = false;
	const double      constants::divergence_avoid    = 5e-4;	// in eV; energy must not be in interval [divergence, divergence+Ecrit]
	const double      constants::divergence_distance = 5e-4;	// in eV; energy will be replaced by divergence-dE_div
	const double 	  constants::scattering_decrease = 0.01;
	const double 	  constants::scattering_ramp     = 1.5;
	const double 	  constants::delta_hwmin_Egap    = 0.2;		// in eV
	const double 	  constants::delta_hwmax_Egap    = 0.5;		// in eV
	const double 	  constants::POP_q0              = 2e8;		// in m-1; was 1e7; 1e9=1nm screening length
	const double 	  constants::imp_q0min           = 2e8;		// in m-1; was 1e8; 1e9=1nm screening length
	const double 	  constants::ALGnorm_neglect     = 1e-20;	// simulator units
	const double 	  constants::GXnorm_neglect_ion  = 1e-20;	// simulator units
	const double 	  constants::kT_broad_hack  	 = 2;		// was 5
	const uint 		  constants::ils_max_iters       = 300;
	const int 		  constants::mpi_master_rank     = 0;
	const std::string constants::source_name("Source");
	const std::string constants::doping_fieldname("DopingConcentration");
	const std::string constants::options_name("options");
//	const std::string constants::default_negf_dir("/home/ba01/u107/steiger/src/angel-nanohub");
    const std::string constants::default_negf_dir(NEGF_DIR);
	const double      constants::elem_contain_tol    = 1e-14;	// unitless;
	const double      constants::min_vector_norm     = 1e-10;
	const double      constants::max_vector_norm     = 1e+10;
	const uint		  constants::bm_max_array_size   = 60;		// 50 is OK
	const uint		  constants::bm_max_num_edge_elems = 50;	// 50 is OK
	const uint		  constants::eqn_array_size      = 20000;
	const bool		  constants::eqn_debug           = true;
	const int	 	  constants::eqn_comp_vals_chunk = 50;	 	// was 100
	const int	 	  constants::newton_calc_jac_chunk = 50; 	// was 100
	const int	 	  constants::newton_calc_rhs_chunk = 50; 	// was 100
	const double	  constants::max_pot_change_V    = 0.3;		// 1.0V
	const double      constants::max_kpquasi_change_V= 0.03; 	// kT is not a bad value
	const uint	  	  constants::slc_num_points      = 200; 	// 100-200 energy points should do it
	const double	  constants::gradient_eps        = 1e-3; 	// 1e-3
	const uint		  constants::max_outer_iterations = 20;
	const uint		  constants::min_newton_iterations = 1;
	const uint		  constants::selfconsist_num_inner_iters = 1;
	
	constants Constants;
	constants::constants() {
		
		// ternary materials
		this->ternary_names.push_back("AlGaAs");
		this->ternary_names.push_back("InGaAs");
		this->ternary_names.push_back("InAlAs");
		this->ternary_names.push_back("InGaN");
		this->ternary_names.push_back("AlGaN");
		this->ternary_names.push_back("AlInN");
		this->first_pure_names.push_back("GaAs");	// order must correspond to the one in ternary_names!
		this->first_pure_names.push_back("InAs");
		this->first_pure_names.push_back("InAs");
		this->first_pure_names.push_back("GaN");
		this->first_pure_names.push_back("GaN");
		this->first_pure_names.push_back("InN");
		this->second_pure_names.push_back("AlAs");	// order must correspond to the one in ternary_names!
		this->second_pure_names.push_back("GaAs");
		this->second_pure_names.push_back("AlAs");
		this->second_pure_names.push_back("InN");
		this->second_pure_names.push_back("AlN");
		this->second_pure_names.push_back("AlN");
		
		// oxide names
		this->oxide_names.push_back("SiO2");
		this->oxide_names.push_back("Oxide");
	}
	
	// constants independent of unit systems
	const double constants::pi    = 3.14159265358979323846264;
	const double constants::sqrt2 = 1.414213562373095145474621858738828450441;
	const double constants::sqrt3 = 1.732050807568877193176604123436845839024;
	const double constants::sqrt6 = 2.449489742783177881335632264381274580956;
	const cplx   constants::imag_unit(0.0, 1.0);

	// "master" constants - SI system - never comment this out!
	const double constants::SIhbar  = 1.05459e-34; //  [J s]
	const double constants::SIhbar2 = 1.11216e-68; //  [J^2 s^2]
	const double constants::SIm0    = 9.1095e-31;  //  [kg]
	const double constants::SIec    = 1.60219e-19; //  [C]
	const double constants::SIeps0  = 8.85E-12;    //  [F m^-1 = C V^-1 m^-1 ]
	const double constants::SIkb    = 1.3807e-23;  //  [J K^-1]
	const double constants::SI_speed_of_light = 299792458.0; //  [m/s]

	// handle the choice of unit system via commenting out!

/*	// eV - um - ps - ec;		--> works fine
	const double units::meter   = 1e6;
	const double units::second  = 1e12;
	const double units::coulomb = 1.0/constants::SIec;
	const double units::joule   = 1.0/constants::SIec;
*/
	// eV - nm - ps - ec;
	const double units::meter   = 1e9;
	const double units::second  = 1e12;
	const double units::coulomb = 1.0/constants::SIec;
	const double units::joule   = 1.0/constants::SIec;
/*
	// eV - 10nm - ps - ec;	
	//	--> good for wire/well/bulk structures because 10nm is about the extension of the carriers in the quant dims
	const double units::meter   = 1e8;
	const double units::second  = 1e12;
	const double units::coulomb = 1.0/constants::SIec;
	const double units::joule   = 1.0/constants::SIec;
*//*
	// eV - 100nm - ps - ec;	
	const double units::meter   = 1e7;
	const double units::second  = 1e12;
	const double units::coulomb = 1.0/constants::SIec;
	const double units::joule   = 1.0/constants::SIec;
*//*
	// eV - 10nm - us - ec;	
	const double units::meter   = 1e8;
	const double units::second  = 1e6;
	const double units::coulomb = 1.0/constants::SIec;
	const double units::joule   = 1.0/constants::SIec;
*/
	// derived from above
	const double units::meter2  = units::meter*units::meter;
	const double units::meter3  = units::meter*units::meter*units::meter;
	const double units::second2 = units::second*units::second;


	/** convert any SI values to the chosen unit system.
     *  if an unknown conversion is requested, an exception will be thrown */
	double constants::convert_from_SI(units::UnitType type, double value) {
		switch(type) {
		case units::unitless:
			return value;
		case units::length:
			return value * units::meter;
		case units::area:
			return value * units::meter2;
		case units::volume:
			return value * units::meter3;
		case units::time:
			return value * units::second;
		case units::rate:
			return value / units::second;
		case units::mass:
			return value * units::second2 / units::meter2 * units::joule;
		case units::velocity:
			return value * units::meter / units::second;
		case units::density_1d:
			return value / units::meter;
		case units::density_2d:
			return value / units::meter2;
		case units::density_3d:
			return value / units::meter3;
		case units::charge:
			return value * units::coulomb;
		case units::potential:
			return value * units::joule / units::coulomb;
		case units::temperature:							// special case
			return value;
		case units::energy:
			return value * units::joule;
		case units::force:
			return value * units::joule / units::meter;
		case units::pressure:
			return value * units::joule / units::meter3;
		case units::power:
			return value * units::joule / units::second;
		case units::spectrum:
			return value / units::second;
		case units::hbar:
			return value * units::joule * units::second;
		case units::dielectric:
			return value * units::coulomb * units::coulomb / units::joule / units::meter;
		case units::diffusion:
			return value * units::meter2 / units::second;
		case units::mobility:
			return value * units::meter2 * units::coulomb / units::second / units::joule;
		case units::electric_field:
			return value * units::joule/units::coulomb / units::meter;
		case units::electrical_current:
			return value * units::coulomb / units::second;
		case units::electrical_current_density_3d:
			return value * units::coulomb / units::second / units::meter2;
		case units::electrical_current_density_2d:
			return value * units::coulomb / units::second / units::meter;
		case units::electrical_current_density_1d:
			return value * units::coulomb / units::second;
		case units::particle_current:
			return value / units::second;
		case units::particle_current_density_3d:
			return value / units::second / units::meter2;
		case units::particle_current_density_2d:
			return value / units::second / units::meter;
		case units::particle_current_density_1d:
			return value / units::second;
		case units::recombination_3d:
			return value / units::second / units::meter3;
		case units::recombination_2d:
			return value / units::second / units::meter2;
		case units::recombination_1d:
			return value / units::second / units::meter;
		case units::radcoeff_3d:
			return value / units::second * units::meter3;
		case units::radcoeff_2d:
			return value / units::second * units::meter2;
		case units::radcoeff_1d:
			return value / units::second * units::meter;
		case units::auger_3d:
			return value / units::second * units::meter3*units::meter3;
		case units::auger_2d:
			return value / units::second * units::meter2*units::meter2;
		case units::auger_1d:
			return value / units::second * units::meter*units::meter;
		case units::srhrate:
			return value / units::second;
		case units::unknown:
			logmsg->emit(LOG_WARN,"Warning: Quantity of unknown unit type. Using conversion factor 1.0.");
			return value;
		default:
			logmsg->emit(LOG_ERROR,"Error while trying to convert unit %d",type);
			NEGF_EXCEPTION("Unknown unit type.");
			return 0.0;
		}
	} // constants::convert_from_SI


	vector<double> constants::convert_from_SI(units::UnitType type, vector<double> values) 
	{
		vector<double> result;
		result.resize(values.size(), 0.0);
		for (uint ii = 0; ii < values.size(); ii++) {
			result[ii] = constants::convert_from_SI(type, values[ii]);
		}
		return result;
	}

	double * constants::get_starting_address(vector<double>& v) {
		return &v[0];
	}

	/*void constants::prepare_Nx_stuff(uint Nx_) {
		Nx = Nx_;
		NxNn = Nx*Nn;
#ifdef BANDED
		work10 = Matc(NxNn,NxNn);
		work11 = Matc(NxNn,NxNn);
#endif
		
#ifdef USE_BANDED
		logmsg->emit(LOG_INFO,"Program was compiled with USE_BANDED flag.");
#endif
	}*/
	
	// ------------------------------------------
	// mathematical functions
	// ------------------------------------------
	
	double negf_math::abs(double val) {	return val > 0 ? val:-val; }
	int    negf_math::sign(double a) {	return (a == 0) ? 0 : (a<0 ? -1 : 1); }
	int    negf_math::sign(int a)    { return (a == 0) ? 0 : (a<0 ? -1 : 1); }
	
	double negf_math::exp(double x)   { 
//		#ifdef NOACML
			return std::exp(x);
//		#else
//			return fastexp(x);	// ACML
//		#endif
	}
	double negf_math::log(double x)   {
//		#ifdef NOACML
			return std::log(x); 
//		#else 
//			return fastlog(x);	// ACML
//		#endif
	}
	double negf_math::sqrt(double x)  { return std::sqrt(x); }
	double negf_math::pow(double base, double exponent)   { return std::pow(base, exponent); }
	double negf_math::sin(double x)   { return std::sin(x); }
	double negf_math::cos(double x)   { return std::cos(x); }
	double negf_math::tan(double x)   { return std::tan(x); }
	double negf_math::atan(double x)  { return std::atan(x); }
	double negf_math::tanh(double x)  { return std::tanh(x); }
	
	/** analytical continuation of arccos
	 *  for z=-1...1 acos(z) = +pi...0
	 *  for z<-1, Re(acos(z)) = 0, Im(acos(z)) < 0 
	 *  for z<-1, Re(acos(z)) = 0, Im(acos(z)) > 0 
	 *  note: std::log(cplx z) has the branch cut on the negative real axis (-infty, 0)
	 */
	cplx   negf_math::acos(const cplx & z) {
		return constants::pi/2.0 + constants::imag_unit*std::log(constants::imag_unit*z+std::sqrt(1.0-z*z));
	}
	
	void negf_math::vector_normalize_3d(double vec[3]) {
		double len = negf_math::sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
		vec[0] /= len; 
		vec[1] /= len; 
		vec[2] /= len;
	}

	void negf_math::vector_normalize_2d(double vec[2]) {
		double len = negf_math::sqrt(vec[0]*vec[0]+vec[1]*vec[1]);
		vec[0] /= len;
		vec[1] /= len;
	}
	

	void negf_math::vector_outer_product(const double a[3], const double b[3], double result[3]) {
		result[0] = a[1]*b[2] - a[2]*b[1];
		result[1] = a[2]*b[0] - a[0]*b[2];
		result[2] = a[0]*b[1] - a[1]*b[0];
	}
	
	double negf_math::vector_norm_3d(const double vec[3]){
		return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	}
	double negf_math::vector_norm_2d(const double vec[2]) {
		return sqrt(vec[0]*vec[0] + vec[1]*vec[1]);
	}
	double negf_math::vector_scalar_product_3d(const double a[3], const double b[3]) {
		return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
	}
	double negf_math::vector_scalar_product_2d(const double a[2], const double b[2]) {
		return a[0]*b[0] + a[1]*b[1];
	}
	
	double negf_math::vector_norm(const std::vector<double> & vec) {
		double result = 0.0;
		for (uint ii = 0; ii < vec.size(); ii++)
			result += vec[ii]*vec[ii];
		return negf_math::sqrt(result);
	}
	

/** Compute the (complete) Fermi-Dirac Integral of half-order specified by the user
 *
 *  order must be -1.5, -0.5, 0.5, 1.5, ...
 *  ...or -1, 0 (solution is analytic there)
 * 
 * Macleod, ACM Trans. Math. Softw. 24, 1 (1998) is implemented for order=-0.5, 0.5, 1.5
 * Trellakis et al., Solid-State Electr. 41, 771 (1997) is implemented for order=-1.5
 * Mohankumar and Natarajan, Phys. Stat. Solidi B 188, 635 (1995) is implemented for order >1.5
 */
	double negf_math::fermi_int(double order, double x)
	{

    NEGF_ASSERT((ceil(order+0.5) == floor(order+0.5))
    			|| order==-1.0 || order==0.0, "order must be a half-integer or -1 or 0.");

	if (order<-1.5) {
		NEGF_EXCEPTION("Orders lower than -1.5 arenot implemented.");
		return -1;
	}

	if (order == -1.5)
	{
		// Implementation of Trellakis et al., Solid-State Elec. 41, 771 (1997)
		// Valid for all arguments in [-inf,inf], Relative error is always < 3e-14

		double left[11] = {	-1.414213562373105093, 		-2.808512455255300933,
						-1.917292888003182316, 		-0.51964702858984548964,
						-4.514622493913476685e-2, 	2.703612737394500803e-6,
						3.210663073550497956, 		3.873760373311574755,
						2.152389655957076826, 		0.53418017742142862007,
						4.504399821333562274e-2};

		double center[17] = {	0.38010481260967607091, 	0.22593578285051823702,
							8.661812316913904133e-2, 	2.403125634916255073e-2,
							5.148781019801262283e-3, 	8.608789052829504534e-4,
							1.168557348150129719e-4, 	1.206267113790721295e-5,
							4.072609453165192928e-7, 	0.28217193934711949097,
							0.25532502851212916001, 	5.822115937381929213e-2,
							2.175496311885601000e-2, 	3.477974457807568784e-3,
							6.853605822664495658e-4, 	4.989323625005823784e-5,
							6.558208582427148730e-6};
				
		double rwing[21] = {	0.39153269200636255407, 	0.50064189058861709245,
							0.22896610770798227715, 	-5.24709434796334638e-3,
							1.336144566100145303e-2, 	-3.235697385706716408e-4,
							1.407066727430599809e-4, 	7.735641583414245416e-8,
							3.28791309486191732e-7, 	1.287508599792794304e-8,
							2.344804735751596651e-11, 	1.088865617843124234,
							0.27391832831247438591, 	0.14157238957806381879,
							2.362939847485379411e-2, 	3.859103717605300851e-3,
							2.380195640732865326e-4, 	4.782052236929457305e-5,
							1.724791405963369419e-7, 	1.714439888942051107e-7,
							1.521854237639634267e-9};
				
		double rtail[11] = {	0.69604099973059069850, 	-615.0182719425226928,
							255357.4886948210812, 		-47856043.79944008599,
							3365374936.989075587, 		-30117644029.81147183,
							-893.6701099152226180, 		375615.7301296548115,
							-72319852.35281170832, 		5476463793.538555903,
							-83575260859.85202949};

		double z, sumn, sumd, ratfun;
		if (x<0) {
			z = negf_math::exp(x);
			sumn = left[5];
			for (int ii=4; ii>=0; ii--)
				sumn = sumn * z + left[ii];
			sumd = 0;
			for (int ii=10; ii>=6; ii--)
				sumd = (sumd + left[ii]) * z;
			ratfun = sumn / (1 + sumd);
			return (z * ratfun + 1) * z;
		}
		else if (x<4) {
			sumn = center[8];
			for (int ii=7; ii>=0; ii--)
				sumn = sumn * x + center[ii];
			sumd = 0;
			for (int ii=16; ii>=9; ii--)
				sumd = (sumd + center[ii]) * x;
			ratfun = sumn / (1 + sumd);
			return ratfun;
		}
		else if (x<20) {
			sumn = rwing[10];
			for (int ii=9; ii>=0; ii--)
				sumn = sumn * x + rwing[ii];
			sumd = 0;
			for (int ii=20; ii>=11; ii--)
				sumd = (sumd + rwing[ii]) * x;
			ratfun = sumn / (1 + sumd);
			return ratfun;
		}
		else {
			z = 1/(x*x);
			sumn = rtail[5];
			for (int ii=4; ii>=0; ii--)
				sumn = sumn * z + rtail[ii];
			sumd = 0;
			for (int ii=10; ii>=6; ii--)
				sumd = (sumd + rtail[ii]) * z;
			ratfun = sumn / (1 + sumd);
			return (ratfun * z + 0.564189583547756286) / negf_math::sqrt(x);
		}
	}

	uint nterm1, nterm2, nterm3;
	double XMIN1, XMIN2, XHIGH, XHIGH1, XHIGH2, expx, t, xsq, chv, tmp, twoe;
	if (order == -0.5)
	{
		/*  The function uses Chebyshev expansions which are given to 16 decimal places for x <= 2, 
			but only 10 decimal places for x > 2.
		
			MACHINE-DEPENDENT CONSTANTS:
			NTERMS1 - The number of terms used from the array ARRFD1. The recommended value is such that
									ABS(ARRFD1(NTERMS1)) < EPS/10    subject to 1 <= NTERMS1 <= 14.
			NTERMS2 - The number of terms used from the array ARRFD2. The recommended value is such that
									ABS(ARRFD2(NTERMS2)) < EPS/10    subject to 1 <= NTERMS1 <= 23.
			NTERMS3 - The number of terms used from the array ARRFD3. The recommended value is such that
									ABS(ARRFD3(NTERMS3)) < EPS/10    subject to 1 <= NTERMS3 <= 28.
			XMIN1 - The value of x below which  FDM0P5(x) = exp(x)  to machine precision.
							The recommended value is    LN ( SQRT(2) * EPSNEG )
			XMIN2 - The value of x below which  FDM0P5(x) = 0.0     to machine precision. 
							The recommended value is    LN ( XMIN )
			XHIGH - The value of x above which  FDM0P5(x) = 2 sqrt (x/pi)  to machine precision. 
							The recommended value is    1 / sqrt( 2 * EPSNEG )
		
			For values of EPS, EPSNEG, and XMIN the user should refer to the
			paper by Cody in ACM. Trans. Math. Soft. Vol. 14 (1988) p303-311.
		
			This code is provided with single and double precision values  of the machine-dependent parameters, 
			suitable for machines which satisfy the IEEE floating-point standard.
		
			AUTHOR: DR. ALLAN MACLEOD, DEPT. OF MATHEMATICS AND STATISTICS, UNIVERSITY OF PAISLEY,
										HIGH ST., PAISLEY, SCOTLAND PA1 2BE   (e-mail: macl-ms0@paisley.ac.uk )
			LATEST UPDATE: 20 NOVEMBER, 1996
		*/
	
		const double ARRFD1[15] = {	1.7863596385102264, -0.999372007632333e-1, 0.64144652216054e-2,
				-0.4356415371345e-3, 0.305216700310e-4, -0.2181648110e-5,
				0.158050781e-6, -0.115620570e-7, 0.8525860e-9, -0.632529e-10,
				0.47159e-11, -0.3530e-12, 0.265e-13, -0.20e-14, 0.2e-15};
		
		const double ARRFD2[24] = {	1.6877111526052352, 0.5978360226336983, 0.357226004541669e-1,
				-0.132144786506426e-1, -0.4040134207447e-3, 0.5330011846887e-3,
				-0.148923504863e-4, -0.218863822916e-4, 0.19652084277e-5,
				0.8565830466e-6, -0.1407723133e-6, -0.305175803e-7, 0.83524532e-8,
				0.9025750e-9, -0.4455471e-9, -0.148342e-10,	0.219266e-10, -0.6579e-12,
				-0.10009e-11, 0.936e-13 , 0.420e-13, -0.71e-14, -0.16e-14, 0.4e-15};
		
		const double ARRFD3[59] = {	0.8707195029590563, 0.59833110231733e-2, -0.432670470895746e-1,
				-0.393083681608590e-1, -0.191482688045932e-1, -0.65582880980158e-2,
				-0.22276691516312e-2, -0.8466786936178e-3, -0.2807459489219e-3,
				-0.955575024348e-4, -0.362367662803e-4, -0.109158468869e-4,
				-0.39356701000e-5, -0.13108192725e-5, -0.2468816388e-6,
				-0.1048380311e-6, 0.236181487e-7, 0.227145359e-7, 0.145775174e-7,
				0.153926767e-7, 0.56924772e-8, 0.50623068e-8, 0.23426075e-8, 0.12652275e-8,
				0.8927773e-9, 0.2994501e-9, 0.2822785e-9, 0.910685e-10, 0.696285e-10,
				0.366225e-10, 0.124351e-10, 0.145019e-10, 0.16645e-11, 0.45856e-11,
				0.6092e-12, 0.9331e-12, 0.5238e-12, -0.56e-14, 0.3170e-12, -0.926e-13,
				0.1265e-12, -0.327e-13, 0.276e-13, 0.33e-14, -0.42e-14, 0.101e-13, -0.73e-14,
				0.64e-14, -0.37e-14, 0.23e-14, -0.9e-15, 0.2e-15, 0.2e-15, -0.3e-15, 0.4e-15,
				-0.3e-15, 0.2e-15, -0.1e-15, 0.1e-15};
	
		const double gam1p5 = 0.8862269254527580;
		twoe = 5.4365636569180905;
		nterm1 = 15;
		nterm2 = 24;
		nterm3 = 59;
		XMIN1 = -36.39023;
		XMIN2 = -708.39641;
		XHIGH = 67108864.0;
		
		if (x<-1) {
			if (x>XMIN1) {
				expx = negf_math::exp(x);
				t = twoe * expx -1;
				return expx * cheval(nterm1,ARRFD1,t);
			} else {
				if (x<XMIN2)
					return 0.0;
				else
					return negf_math::exp(x);
			}
		} else if (x<=2) {
				t = (2*x - 1) / 3;
				return cheval(nterm2 , ARRFD2 , t);
		} else {
			tmp = negf_math::sqrt(x) / gam1p5;
			if (x<XHIGH) {
				xsq = x * x;
				t = (50 - xsq) / (42 + xsq);
				chv = cheval(nterm3, ARRFD3, t);
				tmp = tmp * (1 - chv / xsq);
			}
			return tmp;
		}
	}
	
	if (order == +0.5)
	{

		const double ARRFD1[14] = {	1.8862968392734597, -0.543580817644053e-1, 0.23644975439720e-2,
				-0.1216929365880e-3, 0.68695130622e-5, -0.4112076172e-6, 0.256351628e-7,
				-0.16465008e-8, 0.1081948e-9, -0.72392e-11, 0.4915e-12, -0.338e-13, 0.23e-14, -0.2e-15};
		const double ARRFD2[24] = {	2.6982492788170612, 1.2389914141133012, 0.2291439379816278,
				0.90316534687279e-2, -0.25776524691246e-2, -0.583681605388e-4,
				0.693609458725e-4, -0.18061670265e-5, -0.21321530005e-5,
				0.1754983951e-6, 0.665325470e-7, -0.101675977e-7, -0.19637597e-8,
				0.5075769e-9, 0.491469e-10, -0.233737e-10, -0.6645e-12, 0.10115e-11,
				-0.313e-13, -0.412e-13, 0.38e-14, 0.16e-14, -0.3e-15, -0.1e-15};
		const double ARRFD3[54] = {	2.5484384198009122, 0.510439408960652e-1, 0.77493527628294e-2,
				-0.75041656584953e-2, -0.77540826320296e-2, -0.45810844539977e-2,
				-0.23431641587363e-2, -0.11788049513591e-2, -0.5802739359702e-3,
				-0.2825350700537e-3, -0.1388136651799e-3, -0.680695084875e-4,
				-0.335356350608e-4, -0.166533018734e-4, -0.82714908266e-5,
				-0.41425714409e-5, -0.20805255294e-5, -0.10479767478e-5,
				-0.5315273802e-6, -0.2694061178e-6, -0.1374878749e-6, -0.702308887e-7,
				-0.359543942e-7, -0.185106126e-7, -0.95023937e-8, -0.49184811e-8,
				-0.25371950e-8, -0.13151532e-8, -0.6835168e-9, -0.3538244e-9,
				-0.1853182e-9, -0.958983e-10, -0.504083e-10, -0.262238e-10, -0.137255e-10,
				-0.72340e-11, -0.37429e-11, -0.20059e-11, -0.10269e-11, -0.5551e-12,
				-0.2857e-12, -0.1520e-12, -0.811e-13, -0.410e-13, -0.234e-13,
				-0.110e-13, -0.67e-14, -0.30e-14, -0.19e-14, -0.9e-15, -0.5e-15, -0.3e-15, -0.1e-15,
				-0.1e-15};
		const double gam2p5 = 1.329340388179137;
		twoe = 5.4365636569180905;
		nterm1 = 14;
		nterm2 = 24;
		nterm3 = 54;
		XMIN1 = -35.7;
		XMIN2 = -708.394;
		XHIGH1 = 7.45467e7;
		XHIGH2 = 3.8392996e205;
		
		if (x>XHIGH2) {
			NEGF_EXCEPTION("Error, x is too big.");
			return 0.0;
		}
		
		if (x<-1) {
			if (x>XMIN1) {
				expx = negf_math::exp(x);
				t = twoe * expx - 1;
				return expx * cheval(nterm1,ARRFD1,t);
			} else if (x<XMIN2)
				return 0.0;
			else
				return negf_math::exp(x);
		} else if (x<=2) {
			t = (2*x - 1) / 3;
			return cheval(nterm2,ARRFD2,t);
		} else {
			tmp = x * negf_math::sqrt(x) / gam2p5;
			if (x<=XHIGH1) {
				xsq = x * x;
				t = (50 - xsq) / (42 + xsq);
				chv = cheval(nterm3, ARRFD3, t);
				tmp = tmp * (1 + chv / xsq);
			}
			return tmp;
		}
	}
	
	if (order == +1.5)
	{
		const double ARRFD1[13] = {	1.9406549210378650, -0.287867475518043e-1, 0.8509157952313e-3,
				-0.332784525669e-4, 0.15171202058e-5, -0.762200874e-7,
				0.40955489e-8, -0.2311964e-9, 0.135537e-10, -0.8187e-12,
				0.507e-13, -0.32e-14, 0.2e-15};
		const double ARRFD2[23] = {	3.5862251615634306, 1.8518290056265751, 0.4612349102417150,
					0.579303976126881e-1, 0.17043790554875e-2, -0.3970520122496e-3,
				-0.70702491890e-5, 0.76599748792e-5, -0.1857811333e-6, -0.1832237956e-6,
				0.139249495e-7, 0.46702027e-8, -0.6671984e-9, -0.1161292e-9, 0.284438e-10,
				0.24906e-11, -0.11431e-11, -0.279e-13, 0.439e-13, -0.14e-14, -0.16e-14,
				0.1e-15, 0.1e-15};
		const double ARRFD3[56] = {	12.1307581736884627, -0.1547501111287255, -0.739007388850999e-1,
				-0.307235377959258e-1, -0.114548579330328e-1, -0.40567636809539e-2,
				-0.13980158373227e-2, -0.4454901810153e-3, -0.1173946112704e-3,
				-0.148408980093e-4, 0.118895154223e-4, 0.146476958178e-4,
				0.113228741730e-4, 0.75762292948e-5, 0.47120400466e-5, 0.28132628202e-5,
				0.16370517341e-5, 0.9351076272e-6, 0.5278689210e-6, 0.2951079870e-6,
				0.1638600190e-6, 0.905205409e-7, 0.497756975e-7, 0.272955863e-7,
				0.149214585e-7, 0.81420359e-8, 0.44349200e-8, 0.24116032e-8, 0.13105018e-8,
				0.7109736e-9, 0.3856721e-9, 0.2089529e-9, 0.1131735e-9, 0.612785e-10,
				0.331448e-10, 0.179419e-10, 0.96953e-11, 0.52463e-11, 0.28343e-11,
				0.15323e-11, 0.8284e-12, 0.4472e-12, 0.2421e-12, 0.1304e-12, 0.707e-13,
				0.381e-13, 0.206e-13, 0.111e-13, 0.60e-14, 0.33e-14, 0.17e-14, 0.11e-14,
				0.5e-15, 0.3e-15, 0.1e-15, 0.1e-15};
		nterm1 = 13;
		nterm2 = 23;
		nterm3 = 56;
		const double gam3p5 = 3.323350970447843;
		twoe = 5.4365636569180905;
		XMIN1 = -35.004;
		XMIN2 = -708.396418;
		XHIGH1 = 166674733.2;
		XHIGH2 = 3.204467e123;

		if (x>XHIGH2) {
			NEGF_EXCEPTION("Error, x is too large.");
			return 0.0;
		}
	
		if (x<-1) {
			if (x>XMIN1) {
				expx = negf_math::exp(x);
				t = twoe * expx - 1;
				return expx * cheval(nterm1,ARRFD1,t);
			} else if (x<XMIN2)
				return 0.0;
			else
				return negf_math::exp(x);
		} else if (x<=2) {
			t = (2*x - 1) / 3;
			return cheval(nterm2, ARRFD2, t);
		} else {
			tmp = x*x*sqrt(x) / gam3p5;
			if (x<=XHIGH1) {
				xsq = x*x;
				t = (50 - xsq) / (42 + xsq);
				chv = cheval(nterm3, ARRFD3, t);
				tmp = tmp * (1 + chv / xsq);
			}
			return tmp;
		}
	}

	if (order > +1.5) // after Mohankumar and Natarajan, Phys. Stat. Solidi B 188, 635 (1995)
	{
		double u0    = 0.0;
		int N, dh, poles, k_max;
		double h, d, epsilon, eqn_19, temp, rho, theta, delta, lambda, nominator, denominator;
		double pole_contrib;
		std::complex<double> zk(0.0, 0.0);
		const std::complex<double> cunit(0.0, 1.0);
		bool error;
		double pi = constants::pi;

		if (x >= 0)
			u0 = negf_math::sqrt(x+64);
		else
			u0=8;
			
		N  = 16;
		dh = 7;
		h  = u0/N;
		d  = dh*h;
			
		temp = 0;
		for (int n=1; n <=N; n++)
			temp = temp + negf_math::pow(n*h, 2.0*order+1.0) / (1 + negf_math::exp((n*h)*(n*h)-x));
		
		epsilon = 1e-20;
		eqn_19 = h * ( negf_math::pow(epsilon, 2.0*order+1.0) / (1+negf_math::exp(epsilon*epsilon-x)) + 2 * temp);
		
		// determine number of poles of 1/(1+exp(u^2-x))
		poles = 0;
		k_max = 200;
		for (int k=1; k <= k_max; k++) {
			zk = x + cunit*pi*(2.0*k-1.0); // complex!
			zk = std::sqrt(zk);
			if (std::fabs(imag(zk)) < d)
				poles = poles + 1;
		}
			
		pole_contrib = 0.0;
		error = false;
		for (int k=1; k <= poles; k++) {
			rho = negf_math::sqrt(x*x + (2*k-1)*(2*k-1) * pi*pi);

			if (x==0)
				theta = pi / 2;
			else if (x>0)
				theta = negf_math::atan((2*k-1)*pi/x);
			else
				theta = pi + negf_math::atan((2*k-1)*pi/x);

			delta = 2*pi/h * negf_math::sqrt(0.5*(rho + x));
			lambda = negf_math::exp(2*pi/h * negf_math::sqrt(0.5*(rho-x)));
			
			nominator = lambda * negf_math::sin(delta+order*theta) - negf_math::sin(order*theta);
			denominator = 1 + lambda*lambda - 2*lambda*cos(delta);
			if (denominator != 0)
				pole_contrib = pole_contrib + 4*pi * negf_math::pow(rho, order + 0.0) * nominator / denominator;
			else
//				NEGF_EXCEPTION("wanted to divide by zero. Some poles were skipped.");
				logmsg->emit(LOG_INFO,"wanted to divide by zero. Some poles were skipped.");
		}
		if (x > 0)
			pole_contrib = -1 * pole_contrib;
		else
			pole_contrib = +1 * pole_contrib;
		
		eqn_19 = eqn_19 + pole_contrib;
		// logmsg->emit(LOG_INFO, "fermi_int(%d,%e)=%e, d=%e, poles=%d, pole_contrib=%e'
		//                                ,order,x,eqn_19,d,poles,pole_contrib));
		
		return eqn_19 / gamma(order+1);
	}
	
	if (order == 0.0)
	{
		return negf_math::log(1.0 + negf_math::exp(x));
	}
	
	if (order == -1.0)
	{
		return 1.0 / (1.0 + negf_math::exp(-x));
	}
	

	NEGF_EXCEPTION("You should not have got here,");
	return -1.0;
	}  // negf_math::fermi_int


	/** Helper function for fermi_int. This function evaluates a Chebyshev series, using the
     *  Clenshaw method with Reinsch modification, as analysed in the paper by Oliver.
     *  J. Oliver, J.I.M.A., vol. 20, 1977, pp379-391
     */
	double negf_math::cheval(uint n, const double A[], double t)
	{

		double u1 = 0.0;	
		double d1 = 0.0;
		double test = 0.6;
		double tt = 0.0;
		double u0 = 0.0;
		double u2 = 0.0;
		double d2 = 0.0;;

		if (std::fabs(t)<test) {
			u0 = 0;
			tt = t + t;
			for (int ii=n-1; ii >= 0; ii--) {
				u2 = u1;
				u1 = u0;
				u0 = tt * u1 + A[ii] - u2;
			}
			return (u0 - u2) / 2;
		} else if (t>=test) {
			tt = (t-0.5) - 0.5;
			tt = tt + tt;
			for (int ii=n-1; ii>=0; ii--) {
				d2 = d1;
				u2 = u1;
				d1 = tt * u2 + A[ii] + d2;
				u1 = d1 + u2;
			}
			return (d1 + d2) / 2;
		} else { // t<=-test
			tt = (t + 0.5) + 0.5;
			tt = tt + tt;
			for (int ii=n-1; ii>=0; ii--) {
				d2 = d1;
				u2 = u1;
				d1 = tt * u2 + A[ii] - d2;
				u1 = d1 - u2;
			}
			return (d1 - d2) / 2;
		}
	}


	/** Compute the inverse of F_0.5 */
	double negf_math::fermihalf_inverse(double y)
	{
		double v = negf_math::pow(3*sqrt(constants::pi)*y/4, 0.666666666666666666666);
		double initial_guess = 0.0;
		if (std::fabs(y - 1) > 1e-8)
			initial_guess = log(y)/(1-y*y) + v / (1 + 1/((0.24+1.08*v)*(0.24+1.08*v)));
		else
			initial_guess = 0.348747362;
		// This would be the result of Nilsson, phys. stat. sol. A 19, K75 (1973)
		// but the error could be as large as 0.006!
		// therefore we do a Newton iteration to find a better approximation.
	
		double err 		= 1;
		double err_crit = 1e-12;
		double max_iter = 30;
	
		double x = initial_guess;
		for (uint iter=0; iter < max_iter; iter++)
		{
			err = fermi_int(0.5, x) - y;
			if (std::fabs(err)<err_crit) {
				// logmsg->emit(LOG_INFO,5,"Convergence was reached after %d iterations", iter);
				break;
			}
			double jac = fermi_int(-0.5, x);
			x  = x - err/jac;
		}
		if (! (std::fabs(err)<err_crit))
			logmsg->emit(LOG_INFO,"fermihalf_inverse: convergence not reached (error %7.3e).",err);
		return x;
	}
	
	double negf_math::fermim0p5_inverse(const double & x)
	{
		const uint max_iter = 50;						// max. # iterations
	
		/*
		const double eps = 1e-13;						// numerical accuracy
		double y  = aqua_math::exp(-x); 				// initial value
		double fy = negf_math::fermi_int(0.5, y) - x;	// the Newton function, f(y) shall be zero
		double fprime;									// derivative of the Newton function
		uint iter = 0;
		while (iter < max_iter && fabs(fy)>eps) {
			fy     = negf_math::fermi_int(0.5, y) - x;
			fprime = negf_math::fermi_int(-0.5, y);
			NEGF_ASSERT(fabs(fprime)>1e-13,"derivative is too small."); 
			y = y - fy/fprime;
			iter++;
		}
		if (iter==max_iter) NEGF_EXCEPTION("Newton iteration for 1D Fermi energy did not converge.");
		return y;
		*/
		
		double y = negf_math::exp(-x);					// initial value
		double err = 1;
		double err_crit = 1e-12;						// numerical accuracy
		for (uint iter=0; iter < max_iter; iter++)
		{
			err = negf_math::fermi_int(-0.5, y) - x;	// the Newton function, f(y) shall be zero
			if (std::fabs(err)<err_crit) {
				break;
			}
			double jac = negf_math::fermi_int(-1.5, y);	// derivative of the Newton function
			y  = y - err/jac;
		}
		if (! (std::fabs(err)<err_crit))
			logmsg->emit(LOG_INFO,"inverse_m0p5: convergence not reached (error %7.3e).",err);
		return y;
	}

	bool negf_math::solve_linear_equation(double a11, double b1, double& x1)
//	bool solve_linear_equation(double a11, double b1, double& x1)
	{
		if ( a11 == 0 )
			{
			x1 = 0;
			NEGF_EXCEPTION("Cannot solve 1D linear equation");
			return false;
			}
		else
			{
			x1 = b1/a11;
			return true;
			}
		return false;
	}
	
	
	bool negf_math::solve_linear_equation(double a11, double a12, double b1,
//	bool solve_linear_equation(double a11, double a12, double b1,
				double a21, double a22, double b2, double& x1, double& x2)
	{
		double denominator =  ( a11 * a22 - a21 * a12 );
		
		if ( denominator == 0 )	{
			x1 = x2 = 0;
			NEGF_EXCEPTION("Cannot solve 2D linear equation");
			return false;
		}
		else {
			x1 = -( a12 * b2 - b1 * a22 ) / denominator;
			x2 = ( a11 * b2 - a21 * b1 ) / denominator;
			return true;
		}
		return false;
	}
	
	
	bool negf_math::solve_linear_equation( double a11, double a12, double a13, double b1,
								double a21, double a22, double a23, double b2,
								double a31, double a32, double a33, double b3,
								double& x1, double& x2, double& x3)
	{
		double denominator = (  a11*a22*a33 - a11*a32*a23 - a21*a12*a33 + a32*a21*a13 
								- a22*a31*a13 + a31*a12*a23 );
		
		if ( denominator == 0 ) {
			NEGF_EXCEPTION("Cannot solve 3D linear equation");
			x1 = x2 = x3 = 0;
			return false;
		}
		else {
			x1 = -( a12 * b2 * a33 - a12 * a23 * b3 - a13 * a32 * b2 + a13 * a22 * 
				b3 - b1 * a22 * a33 + b1 * a32 * a23 ) / denominator;
			x2 = ( a11 * b2 * a33 - a11 * a23 * b3 + a21 * a13 * b3 + a23 * a31 * b1
				- b2 * a31 * a13 - a21 * b1 * a33 ) / denominator;
			x3 = ( a32 * a21 * b1 - a11 * a32 * b2 + a11 * a22 * b3 - a22 * a31 * b1
				- a21 * a12 * b3 + a31 * a12 * b2 ) / denominator;
			return true;
		}
		return false;
	}
	
	// matrix_norm is implemented in MatFlens.cpp!
	
	void negf_math::zgeev(int n, cplx * matrix, cplx * eigenvalues, cplx * eigenvectors)
	{STACK_TRACE(
		int  info  = 0;
#ifndef NOACML
		::zgeev('N', 'V', n, matrix, n, eigenvalues, NULL, 1, eigenvectors, n, &info); 	// ACML routine
#else
		char     jobvl = 'N';
		char     jobvr = 'V';
		cplx   * vl    = NULL;
		int      ldvl  = 1;
	    int      lwork = 4*n; // (input) The dimension of the array work.  LWORK >= max(1,2*N). For good performance, LWORK must generally be larger.   
                       		// If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns   
                       		// this value as the first entry of the WORK array, and no error message related to LWORK is issued by XERBLA.   
	    cplx   * work = new cplx[lwork];  // (workspace/output) array of dimension (max(1,lwork)). On exit, if info=0, work[0] returns the optimal LWORK. 
	    double * rwork = new double[2*n]; // (workspace) array of dimension (2*N)
		
		::F77_Name(zgeev)(&jobvl, &jobvr, &n, matrix, &n, eigenvalues, vl, &ldvl, eigenvectors, &n, work, &lwork, rwork, &info);
		delete [] work;
		delete [] rwork;
#endif
		NEGF_FASSERT(info==0, "LAPACK zgeev failed: info=%d",info);
	);}
	
	void negf_math::check_compress_error(int err, unsigned long compressed_size, unsigned long num_chars)
	{STACK_TRACE(
		if (err!=Z_OK) { 
			char errstr[1000];
			switch (err) {
				case Z_MEM_ERROR: sprintf(errstr, "out of memory."); break;
				case Z_BUF_ERROR: sprintf(errstr, "output buffer too small (size %ld, num_chars=%ld).", compressed_size, num_chars); break;
				default:          sprintf(errstr, "unknown.");
			}
			NEGF_FEXCEPTION("p%d: error %d occurred during the compression: %s", mpi->get_rank(), err, errstr);
		}
	);}
	
	void negf_math::check_decompress_error(int err, unsigned long decomp_size, unsigned long num_chars)
	{STACK_TRACE(
		if (err!=Z_OK) { 
			char errstr[1000];
			switch (err) {
				case Z_MEM_ERROR: sprintf(errstr, "out of memory."); break;
				case Z_BUF_ERROR: sprintf(errstr, "output buffer too small (size %ld, num_chars=%ld).", decomp_size, num_chars); break;
				default:          sprintf(errstr, "unknown.");
			}
			NEGF_FEXCEPTION("p%d: error %d occurred during the decompression: %s", mpi->get_rank(), err, errstr);
		}
	);}
	
	void negf_math::do_compress_antiherm(const vector<BMatc> & data, double* & real_data, double* & imag_data, unsigned char* & real_char, unsigned char* & imag_char, unsigned long & num_chars, 
						unsigned long & real_comp_size, unsigned long & imag_comp_size, unsigned char* & real_compressed, unsigned char* & imag_compressed)
	{STACK_TRACE(
		unsigned long num_doubles = negf_math::check_and_get_num_doubles(data, true);
		const int max_distance = data[0].num_offdiags;
		// max_distance needs to be the same as in myMPI::recv(vector<BMatc> & data, int & source, int & tag)
		int n = data[0].num_rows();
		int num_matrices = data.size();

		// ----------------------------------------------
		// create double arrays
		// ----------------------------------------------
		unsigned long count = 0;
		real_data = new double[num_doubles];
		imag_data = new double[num_doubles];
		for (int mm=0; mm<num_matrices; mm++) {
			for (int ii=0; ii<n; ii++) {
				for (int jj=max(0,ii-max_distance); jj<=ii; jj++) {	// lower half is sent
					real_data[count] = data[mm](ii+1,jj+1).real();
					imag_data[count] = data[mm](ii+1,jj+1).imag();
					count++;
				}
			}
		}
		NEGF_ASSERT(count==num_doubles, "something went wrong.");
		
		negf_math::do_compress_helper(real_data, imag_data, num_doubles, real_char, imag_char, num_chars, real_comp_size, imag_comp_size, 
				real_compressed, imag_compressed);
	);}
		

	// helper function for compression - real_data, imag_data has already been set up
	// see stuff/compress.cc for a test case
	void negf_math::do_compress_helper(double* & real_data, double* & imag_data, unsigned long & num_doubles, unsigned char* & real_char, unsigned char* & imag_char, unsigned long & num_chars, 
						unsigned long & real_comp_size, unsigned long & imag_comp_size, unsigned char* & real_compressed, unsigned char* & imag_compressed)
	{STACK_TRACE(
		int doublesize = sizeof(double);
		int charsize = sizeof(char);
		NEGF_ASSERT(charsize==1, "class is designed for sizeof(char)=1");
		
		// cast into byte (unsigned int)
		real_char = (unsigned char *)(real_data);
		imag_char = (unsigned char *)(imag_data);
		num_chars = num_doubles*doublesize;
	
		// compress. from the zlib homepage: "...the destination buffer, which must be at least 0.1% larger than sourceLen plus 12 bytes."
		real_comp_size = (unsigned long) (ceil(1.01*double(num_chars))) + 12;
		imag_comp_size = real_comp_size;
		unsigned long * real_compressed_size = &real_comp_size;
		unsigned long * imag_compressed_size = &imag_comp_size;
		real_compressed = new unsigned char[real_comp_size];	
		imag_compressed = new unsigned char[imag_comp_size];
	
		int clevel = constants::mpi_compresslevel;	
		// Z_DEFAULT_COMPRESSION, or between 0 and 9: 1 gives best speed, 9 gives best compression
	
		logmsg->emit_noendl(LOG_INFO_L3,"p%d REAL: uncompressed size:%d, maximum compressed size:%d \n    compressing...", mpi->get_rank(), num_chars, real_comp_size);
		int err = compress2(real_compressed, real_compressed_size, real_char, num_chars, clevel);	// zlib routine
		negf_math::check_compress_error(err, real_comp_size, num_chars);
		logmsg->emit(LOG_INFO_L3, "done. compressed size: %d",(*real_compressed_size));
	
		unsigned long tmp2 = num_chars;
		num_chars = num_doubles*doublesize;
		logmsg->emit_noendl(LOG_INFO_L3,"p%d IMAG: uncompressed size:%d, maximum compressed size:%d \n    compressing...", mpi->get_rank(), num_chars, imag_comp_size);
		err = compress2(imag_compressed, imag_compressed_size, imag_char, num_chars, clevel);	// zlib routine
		negf_math::check_compress_error(err, imag_comp_size, num_chars);
		logmsg->emit(LOG_INFO_L3, "done. compressed size: %d",(*imag_compressed_size));
		NEGF_ASSERT(num_chars==tmp2, "something went wrong (num_chars!=tmp)");
		
		NEGF_ASSERT(sizeof(int)==4, "sizeof(int)==4 expected.");
		unsigned long max_array_size = 2147483648UL; // this gives 268'435'456 doubles or, at an array size of 1000*1000, roughly 268*2>500 possible k-points
		NEGF_ASSERT(*real_compressed_size < max_array_size && *imag_compressed_size < max_array_size, "maximum array size exceeded.");
	);}
	
	void write_matrix(const char * filename, const Matc & matrix, const string & description)
	{STACK_TRACE(
		logmsg->emit_noendl(LOG_INFO_L1,  "Writing %s ...",filename);
		ofstream fout(filename);
		NEGF_FASSERT(fout, "%s could not be opened for output.", filename);
		fout.precision(12);
		fout.setf(ios::right);
		fout << "% NEGF matrix " << matrix.num_rows() << "x" << matrix.num_cols() << "\n";
		if (description!="") fout << "% " << description << "\n";
		fout << "% first " << matrix.num_rows() << " lines are real part, second " << matrix.num_rows() << " lines are imaginary part\n";
		for (uint ii=1; ii <= matrix.num_rows(); ii++) {
			for (uint jj=1; jj <= matrix.num_cols(); jj++) {
				fout << setw(18) << matrix(ii,jj).real();
				if (jj < matrix.num_cols()) fout << "\t";
			}
			fout << "\n";
		}
		for (uint ii=1; ii <= matrix.num_rows(); ii++) {
			for (uint jj=1; jj <= matrix.num_cols(); jj++) {
				fout << setw(18) << matrix(ii,jj).imag();
				if (jj < matrix.num_cols()) fout << "\t";
			}
			fout << "\n";
		}
		fout << "\n";
		fout.close();
		logmsg->emit(LOG_INFO_L1, " Done.");
	);}
	
	
	void write_matrix(const char * filename, const Matd & matrix, const string & description)
	{STACK_TRACE(
		logmsg->emit_noendl(LOG_INFO_L1,  "Writing %s ...",filename);
		ofstream fout(filename);
		NEGF_FASSERT(fout, "%s could not be opened for output.", filename);
		fout.precision(12);
		fout.setf(ios::right);
		fout << "% NEGF matrix " << matrix.num_rows() << "x" << matrix.num_cols() << "\n";
		if (description!="") fout << "% " << description << "\n";
		for (uint ii=1; ii <= matrix.num_rows(); ii++) {
			for (uint jj=1; jj <= matrix.num_cols(); jj++) {
				fout << setw(18) << matrix(ii,jj);
				if (jj < matrix.num_cols()) fout << "\t";
			}
			fout << "\n";
		}
		fout << "\n";
		fout.close();
		logmsg->emit(LOG_INFO_L1, " Done.");
	);}
	
	
	void write_xE_matrix(const char * filename, const Matd & matrix, const vector<double> & xcoord, const vector<double> & energies)
	{STACK_TRACE(
		logmsg->emit_noendl(LOG_INFO_L1,  "Writing %s ...",filename);
		ofstream fout(filename);
		NEGF_FASSERT(fout, "%s could not be opened for output.", filename);
		fout.precision(12);
		fout.setf(ios::right);
		fout << "% NEGF real matrix " << matrix.num_rows() << "x" << matrix.num_cols() << "\n";
		fout << "% first row contains x-coordinates, first column contains energies\n";
		NEGF_ASSERT(matrix.num_cols()==xcoord.size(), "inconsistent x-coordinate vector.");
		NEGF_ASSERT(matrix.num_rows()==energies.size(), "inconsistent energy vector.");
		// write x-coordinates
		fout << 0.0 << "\t";
		for (uint jj=1; jj <= matrix.num_cols(); jj++) {
			fout << setw(18) << xcoord[jj-1];
			if (jj < matrix.num_cols()) fout << "\t";
		}
		fout << "\n";
		// write matrix
		for (uint ii=1; ii <= matrix.num_rows(); ii++) {
			fout << setw(18) << energies[ii-1] << "\t";
			for (uint jj=1; jj <= matrix.num_cols(); jj++) {
				fout << setw(18) << matrix(ii,jj);
				if (jj < matrix.num_cols()) fout << "\t";
			}
			fout << "\n";
		}
		fout << "\n";
		fout.close();
		logmsg->emit(LOG_INFO_L1, " Done.");
	);}
	
	void write_current_matrix(const char * filename, const Matd & matrix, const vector<double> & xcoord, const vector<double> & energies)
	{STACK_TRACE(
		logmsg->emit_noendl(LOG_INFO_L1,  "Writing %s ...",filename);
		ofstream fout(filename);
		NEGF_FASSERT(fout, "%s could not be opened for output.", filename);
		fout.precision(12);
		fout.setf(ios::right);
		fout << "% NEGF real matrix " << matrix.num_rows() << "x" << matrix.num_cols() << "\n";
		fout << "% first row contains x-coordinates, first column contains energies\n";
		NEGF_FASSERT(matrix.num_cols()==xcoord.size()-1, "inconsistent x-coordinate vector: Nx-1=%d, num_cols=%d.",xcoord.size()-1,matrix.num_cols());
		NEGF_FASSERT(matrix.num_rows()==energies.size(), "inconsistent energy vector: NE=%d, num_rows=%d.",energies.size(),matrix.num_rows());
		// write x-coordinates
		fout << 0.0 << "\t";
		for (uint jj=1; jj <= matrix.num_cols(); jj++) {
			fout << setw(18) << 0.5*(xcoord[jj-1]+xcoord[jj]);
			if (jj < matrix.num_cols()) fout << "\t";
		}
		fout << "\n";
		// write matrix
		for (uint ii=1; ii <= matrix.num_rows(); ii++) {
			fout << setw(18) << energies[ii-1] << "\t";
			for (uint jj=1; jj <= matrix.num_cols(); jj++) {
				fout << setw(18) << matrix(ii,jj);
				if (jj < matrix.num_cols()) fout << "\t";
			}
			fout << "\n";
		}
		fout << "\n";
		fout.close();
		logmsg->emit(LOG_INFO_L1, " Done.");
	);}
	
	
	void read_matrix(const char * filename, Matc & matrix)
	{STACK_TRACE(
		logmsg->emit_noendl(LOG_INFO_L1,  "Reading %s ...",filename);
		ifstream fin(filename);
		NEGF_FASSERT(fin, "%s could not be opened.", filename);
		
		char buf[10000];
		// first we expect 2 lines of comments
		NEGF_ASSERT(!fin.eof(), "file is too short!");
		fin.getline(buf,1000);
		NEGF_ASSERT(buf[0]=='%', "expected comment in line 1.");
		NEGF_ASSERT(!fin.eof(), "file is too short!");
		fin.getline(buf,1000);
		NEGF_ASSERT(buf[0]=='%', "expected comment in line 2.");
		
		int N = matrix.num_rows();
		int M = matrix.num_cols();
		double tmp;
		// read real part
		for(int ii = 1; ii <= N; ii++) {
			for(int jj = 1; jj <= M; jj++) {
				NEGF_ASSERT(!fin.eof(), "file is too short!");
				fin >> tmp; 	
				matrix(ii,jj) = tmp;
			}	
		} 
		// read imaginary part
		for(int ii = 1; ii <= N; ii++) {
			for(int jj = 1; jj <= M; jj++) {
				NEGF_ASSERT(!fin.eof(), "file is too short!");
				fin >> tmp; 	
				matrix(ii,jj) = matrix(ii,jj) + constants::imag_unit * tmp;
			}	
		} 
		// expect one more empty line 
		NEGF_ASSERT(!fin.eof(), "file is too short!");
		fin.getline(buf,1000);
		NEGF_ASSERT(!fin.eof(), "file is too short!");
		fin.getline(buf,1000);
		NEGF_ASSERT(fin.eof(), "expected end of file.");
		fin.close();
		logmsg->emit(LOG_INFO_L1, " Done.");
	);}
	
} // namespace negf

