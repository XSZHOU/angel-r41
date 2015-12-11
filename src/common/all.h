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
#ifndef _ALL_H_NEGF
#define _ALL_H_NEGF

/*! \file all.h
    This header should be included by (almost) all other classes and provides accessibility to
    many global constants, classes and mathematical functions.
*/

// -----------------------------------------------------
// compiler includes (for the entire program)
// -----------------------------------------------------
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <complex>
#include <vector>
#include <map>
#include <list>
#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <algorithm>

// -----------------------------------------------------
// type definitions
// -----------------------------------------------------
typedef unsigned int			uint;
typedef unsigned short int		usint;
typedef unsigned long int 		ulint;
typedef std::complex<double> 	cplx;

// ----------------------------------------
// preprocessor defines
// ----------------------------------------
#define FLENS				// include FLENS matrix operations
#define BANDED				// include banded matrix operations
#define USE_BANDED			// use banded matrices for self-energies
//#define DIAGONAL_ONLY 	// Luminescence.cpp 
	
namespace negf {
	// number of degrees of freedom (bands), set by Options class
	extern uint Nn;
	// number of internal vertices, set by main.cpp / python script as soon as grid file is read in
	//extern uint Nx;
	// product
	//extern uint NxNn;
	
	// some forward declarations
	class BMatc;
	class Matc;
	class Matd;
	
	// set here whether self-energies, overlap matrices and propagators are stored in banded form! 
#ifdef USE_BANDED
	typedef BMatc SEMat;   										//!< self-energy matrices
	typedef BMatc GLMat;   										//!< GL, GR, MGLM, MGRM
	typedef BMatc OVMat;   										//!< overlap matrices
	
	#define    SEMat_create(N) BMatc(N,N,constants::odSE)		//!< total self-energy
	#define    GLMat_create(N) BMatc(N,N,constants::odGL)		//!< lesser and greater GF
	#define    OVMat_create(N) BMatc(N,N,constants::odOV)		//!< overlap matrices
	#define    SCMat_create(N) BMatc(N,N,constants::odSC)		//!< contact self-energy
	#define   SACMat_create(N) BMatc(N,N,constants::odSAC)		//!< acoustic phonon self-energy
	#define  SPOPMat_create(N) BMatc(N,N,constants::odSPOP)		//!< polar optical phonon self-energy
	#define  SGolMat_create(N) BMatc(N,N,constants::odSGol)		//!< Golizadeh and Buettiker self-energy
	#define SPhotMat_create(N) BMatc(N,N,constants::odSPhot)	//!< photon self-energy
	#define   SIIMat_create(N) BMatc(N,N,constants::odSII)		//!< ionized impurities self-energy
#else
	typedef Matc SEMat;
	typedef Matc GLMat;
	typedef Matc OVMat;
	
	#define    SEMat_create(N) Matc(N,N)		//!< total self-energy
	#define    GLMat_create(N) Matc(N,N)		//!< lesser and greater GF
	#define    OVMat_create(N) Matc(N,N)		//!< overlap matrices
	#define    SCMat_create(N) Matc(N,N)		//!< contact self-energy
	#define   SACMat_create(N) Matc(N,N)		//!< acoustic phonon self-energy
	#define  SPOPMat_create(N) Matc(N,N)		//!< polar optical phonon self-energy
	#define  SGolMat_create(N) Matc(N,N)		//!< Golizadeh and Buettiker self-energy
	#define SPhotMat_create(N) Matc(N,N)		//!< photon self-energy
	#define   SIIMat_create(N) Matc(N,N)		//!< ionized impurities self-energy
#endif

	
	// some globally accessible matrices used as temporary workhorses s.th. memory does not
	// have to be allocated and deallocated all the time
	extern Matc work10;   // NxNn*NxNn, defined in all.cpp, used in MatBanded.cpp, resized in constants::prepare_Nx_stuff
	extern Matc work11;   // NxNn*NxNn, defined in all.cpp, used in MatBanded.cpp, resized in constants::prepare_Nx_stuff
}

// -----------------------------------------------------------------------
// program includes for globally existing classes which are made available
//    too the vast majority of .h files
// be sure to include these as dependency for each class in your Makefile
// however, include only all.h in the other .h files
// -----------------------------------------------------------------------
#include "Exception.h"
#include "Logger.h" 
#include "Timer.h"
#include "Filenames.h"
#include "Interrupt.h"
// Mat.h and MPI.h are included at the and of this file because it relies on negf_math routines

namespace negf {
		
	// ----------------------------------------------------------------
	// global variables
	// this is only the declaration; creation takes place in all.cpp
	// ----------------------------------------------------------------
	class Logger;	 extern Logger    * logmsg;
	class Timer;	 extern Timer     * timer;
	class Filenames; extern Filenames * fnames;
	class Interrupt; extern Interrupt * interrupt;
	class myMPI;	 extern myMPI 	  * mpi;
	

	// -----------------------------------------------------
	// enumerations
	// -----------------------------------------------------

	/** units for different physical quantities (note: a physical quantity like density may have different units - 1D, 2D, 3D density) */
	class units {
		public:
		enum UnitType {
			unitless,
			length,
			area,
			volume,
			time,
			rate,							//!< 1/time
			mass,
			velocity,
			density_1d,
			density_2d,
			density_3d,
			charge,
			potential,
			temperature,
			energy,
			force,
			pressure,
			power,
			spectrum,						//!< =power/energy =1/time
			hbar,
			dielectric,
			diffusion,
			mobility,
			electric_field,
			electrical_current,				//!< units charge time^-1
			electrical_current_density_3d,  //!< units charge time^-1 length^-2	 (cross-section is an area)
			electrical_current_density_2d,  //!< units charge time^-1 length^-1	 (cross-section is a line)
			electrical_current_density_1d,  //!< units charge time^-1 			 (cross-section is a point)
			particle_current,				//!< units time^-1
			particle_current_density_3d,	//!< units time^-1 length^-2
			particle_current_density_2d,	//!< units time^-1 length^-1
			particle_current_density_1d,	//!< units time^-1
			recombination_3d,				//!< units length^-3 time^-1
			recombination_2d,				//!< units length^-2 time^-1
			recombination_1d,				//!< units length^-1 time^-1
			radcoeff_3d,					//!< B of Bnp, units length^-3 time^-1
			radcoeff_2d,					//!< units length^-2 time^-1
			radcoeff_1d,					//!< units length^-1 time^-1
			auger_3d,						//!< C of (C_n n + C_p p)(np-n_i^2), units length^-6 time^-1
			auger_2d,						//!< units length^-4 time^-1
			auger_1d,						//!< units length^-2 time^-1
			srhrate,						//!< 1/tau in SRH rec formula, units time^-1
			unknown
		};
		
		static const double meter;	 //!< SI meter   expressed in simulator length units  (e.g. eV-um-ps-ec:  1e6)
		static const double second;  //!< SI second  expressed in simulator length units  (e.g. eV-um-ps-ec:  1e12)
		static const double coulomb; //!< SI coulomb expressed in simulator length units  (e.g. eV-um-ps-ec:  ~1e19)
		static const double joule;	 //!< SI joule   expressed in simulator length units  (e.g. eV-um-ps-ec:  ~1e19)
		
		static const double meter2;		// derived quantities
		static const double meter3;
		static const double second2;
		
	}; // end class units

	/** Physical types that an Equation class can have */
	class quantities {
		public:
		enum PhysicalQuantity { 
			length,
			area,
			volume,
			mass,
			time,
			rate,
			electron_density,
			hole_density,
			dopant_density,
			potential,
			force,
			electric_field,
			dielectric_coeff,
			diffusion_coeff,
			mobility,
			temperature,
			recombination,			// units: length^-3 time^-1
			radiative_coefficient,	// B of Bnp
			auger_coefficient,		// C of (C_n n + C_p p)(np-n_i^2), units length^-6 time^-1
			srh_coefficient,		// 1/tau in SRH rec formula, units time^-1
			avalanche_coeff,		// units 1/length
			energy,
			fermilevel,
			electrical_current,
			particle_current,
			electrical_current_density,
			particle_current_density,
			spectrum,
			light_output_power,
			spont_em_rate_per_energy,
			particles,				// an absolute number of particles (unitless), not a density
			mixed
		};
	};

	/** Types of boundary conditions */
	class bndconds {
		public:
		enum BndCond {
			BC_unknown,
			BC_Dirichlet,	//!< fixed value
			BC_Neumann,		//!< derivative is zero
			BC_Ohm,			//!< not used at the moment
			BC_Schottky		//!< not used at the moment
		};
	};

	/** Statistics influencing boundary conditions, etc */
	class statistics {
		public:
		enum StatType {
			fermi,
			boltzmann
		};
	};

	/** Element types */
	class element_type {
		public:
		enum ElementType {
			interval    = 0,	//!< 1D
			triangle    = 1,	//!< 2D
			rectangle   = 2,	//!< 2D - use is discouraged because of problem in vertex ordering
			tetrahedron = 3,	//!< 3D
			pyramid     = 4,	//!< 3D - currently not used
			prism       = 5,	//!< 3D - currently not used
			cuboid      = 6,	//!< 3D - currently not used
			tetrabrick  = 7,	//!< 3D - currently not used
			point       = 8,	//!< 3D - currently not used
			noshape     = 9		//!< erroneous
		};
	};

	/** Face types */
	class face_type {
		public:
		enum FaceType {
			FT_Point     = 0,	//!< currently not used
			FT_Interval  = 1,	//!< currently not used
			FT_Triangle  = 2,	//!< 2D
			FT_Rectangle = 3,	//!< currently not used
			FT_Other     = 4	//!< currently not used
		};
	};

	/** Linear solver types */
	class linear_solver_type {
		public:
		enum LinearSolverType {
			ils     = 0,		//!< Parallel (OpenMP) iterative solver ILS - very stable
			pardiso = 1,		//!< Parallel or sequential direct solver PARDISO - obtainable at www.pardiso-project.org
			superlu = 2,		//!< SuperLU - very slow, never really used
			mumps   = 3,		//!< Parallel (MPI) direct solver MUMPS - free - support is experimental
			umfpack = 4,		//!< Sequential direct solver UMFPACK - free (open source) - very stable, my favorite
			other   = 5
		};
	};

	// ----------------------------------------------------------
	// constants (only declaration; implementation is in all.cpp)
	// ----------------------------------------------------------
	
	class constants {
	public:
		constants();	//!< constructor
		~constants() {}	//!< destructor
			
		// number of stored offdiagonals - only relevant in banded case
		static const uint odSE;					//!< number of offdiagonals in banded total self-energy matrices
		static const uint odGL;					//!< number of offdiagonals in banded Green function matrices, total self-energy matrices
		static const uint odOV;					//!< number of offdiagonals in banded overlap matrices
		static const uint odSC;					//!< number of offdiagonals in banded contact self-energies
		static const uint odSAC;				//!< number of offdiagonals in banded acoustic phonon self-energies
		static const uint odSPOP;				//!< number of offdiagonals in banded optical phonon self-energies
		static const uint odSGol;				//!< number of offdiagonals in banded Golizadeh and Buettiker self-energies
		static const uint odSPhot;				//!< number of offdiagonals in banded photon self-energies
		static const uint odSII;				//!< number of offdiagonals in banded ionized impurities self-energies
		
		// simulator constants
		static const int mpi_master_rank;		//!< will be 0
		static const double antiherm_check;		//!< numerical error criterion for anti-hermiticity checks
		static const int mpi_compresslevel;		//!< level (0=none,1=fastest,9=highest) of zlib compression in MPI communication of self-energies
		static const double egrid_weight_fL;	//!< weight for left fermilevel that enters in monotonic function for determination of energy grid
		static const double egrid_weight_fR;	//!< weight for right fermilevel that enters in monotonic function for determination of energy grid
		static const double egrid_weight_res;	//!< weight for mid-structure resonance that enters in monotonic function for determination of energy grid
		static const uint max_inner_iters;		//!< maximum # inner iterations (for fixed potential)
		static const uint max_inner_iters_start;//!< maximum # inner iterations in first outer iteration and during scattering ramping
		//static const double inner_errcrit;	//!< error criterion used for checking convergence of the inner loop, NOW IN CMD.CNF
		static const double inner_ramp_errcrit;	//!< error criterion for the inner loop during ramping of scattering strength
		static const double refdens;			//!< reference density for error criterion of inner loop, in m^3
		static const double refcurr;			//!< reference current for error criterion of inner loop, in A/m^2
		static const uint max_outer_iters;		//!< maximum # outer iterations (between potential and NEGF loop)
		static const double outer_errcrit;		//!< error criterion for potential (in Volts) used for checking convergence of the outer loop
		static const double imag_err;			//!< maximum imaginary part that real quantities may have
		static const double dyson_eta;			//!< small imaginary part added to energy in Dyson equation (InnerLoop.cpp)
		static const double contact_eta;		//!< default broadening when incoherent contacts are included (SEContacts.cpp)
		static const double contact_min_xdecay;	//!< minimum decay length of contact states (SEContacts.cpp)
		static const bool old_orthogonal;		//!< if true, kpmethod=0 and 3 use old computation w/ unitless unity overlap matrix
		static const double divergence_avoid;	//!< energy points within this range (in eV) of a contact DOS divergence will be replaced
		static const double divergence_distance;//!< new energy point will be at (divergence - this number)
		static const double scattering_decrease;//!< default decrease factor when scattering is ramped in the beginning
		static const double scattering_ramp;	//!< default ramping factor when scattering is ramped in the beginning
		static const double delta_hwmin_Egap;	//!< minimum calculated spontaneous emission energy, energetic distance from bandgap
		static const double delta_hwmax_Egap;	//!< maximum calculated spontaneous emission energy, energetic distance from bandgap
		static const double POP_q0;				//!< inverse screening length in polar optical phonon scattering
		static const double imp_q0min;			//!< minimum inverse screening length in ionized impurity scattering
		static const double ALGnorm_neglect;	//!< norm of AL/AG/AR in SEPhotonSpontaneous below which contribution is nelected
		static const double GXnorm_neglect_ion; //!< norm of MGLM/MGGM/MGRM in SEIonizedImpurities below which contribution is neglected
		static const double kT_broad_hack;		//!< Resonances.cpp: Fermi-Dirac function denominator is x*kT instead of kT in monotonic function
		static const uint ils_max_iters;		//!< maximum # Krylov iterations within ILS linear solver employed in Newton solver
		static const std::string source_name;	//!< name of source contact --> there the potential is fixed to zero
		static const std::string doping_fieldname;	//!< name of the field to read in which contains the doping info
		static const std::string options_name;	//!< name of options identifier in .cmd-file
		static const std::string default_negf_dir;	//!< default location of NEGF directory with conf/, materials/ etc.
		static const double elem_contain_tol;	//!< tolerance if a vertex is contained in an element, see Element.cpp
		static const double min_vector_norm;	//!< minimum vector norm (used e.g. in quantum regions direction vectors, Geometry.cpp)
		static const double max_vector_norm;	//!< maximum vector norm
		static const uint bm_max_array_size;	//!< maximum array size for certain quantities in class BoxMethod
		static const uint bm_max_num_edge_elems;//!< maximum number of elements around edge in BoxMethod::get_ordered_edge_elements()
		static const uint eqn_array_size;		//!< maximum number of nonzeros within a line in the Jacobian of the Newton solution
		static const bool eqn_debug;			//!< if true some error checking is performed (at the expense of a little speed)
		static const int eqn_comp_vals_chunk;   //!< OpenMP chunk size of Equation::compute_values
		static const int newton_calc_jac_chunk; //!< OpenMP chunk size of NewtonSolver::calculate_jacobian
		static const int newton_calc_rhs_chunk; //!< OpenMP chunk size of NewtonSolver::calculate_rhs
		static const double max_pot_change_V;	//!< maximum change of the potential in one Newton step (units: volts)
		static const double max_kpquasi_change_V;//!< maximum change of the KP potential (units: volts)
		static const uint slc_num_points;  		//!< number of points on the energy axis for the TDKP luminescence calculator
		static const double gradient_eps;		//!< factor by which the minimum connected edge length is multiplied for the
												// finite difference spacing of the gradient at a certain vertex
		static const uint max_outer_iterations;	//!< max. number of outer iterations to reach self-consistency
		static const uint selfconsist_num_inner_iters; //!< max. number of iterations for inner (Newton) convergence to mark the criterion of outer loop convergence
		static const uint min_newton_iterations; //!< min. number of iterations before Newton convergence can be reached
		vector<string> ternary_names;			//!< material names which stand for ternary materials (like "AlGaAs")
		vector<string> first_pure_names;		//!< for each ternary material the name of the x=0 material (e.g. "GaAs" is AlGaAs, x=0)
		vector<string> second_pure_names;		//!< for each ternary material the name of the x=1 material (e.g. "AlAs" is AlGaAs, x=1)
		vector<string> oxide_names;				//!< treatment of oxide materials is a little special at times
		
		// mathematical constants independent of unit systems
		static const double pi;
		static const double sqrt2;
		static const double sqrt3;
		static const double sqrt6;
		static const cplx   imag_unit;

		// master constants of nature, in SI systems (constants for other units are derived from this)
		static const double SIhbar;
		static const double SIhbar2;
		static const double SIm0;
		static const double SIec;
		static const double SIeps0;
		static const double SIkb;
		static const double SI_speed_of_light;

		static double convert_from_SI(units::UnitType type, double value);
		static vector<double> convert_from_SI(units::UnitType type, vector<double> values);
		
		// used in linear solver classes s.th. one can work with vector<double> objects in python 
		// and then still create linear solvers from there
		static double * get_starting_address(vector<double>& v);
		
		// NEEDS TO BE CALLED IMMEDIATELY AFTER Nx IS KNOWN!
		//static void prepare_Nx_stuff(const uint Nx_);
	};
	extern constants Constants;	// an object "Constants" is created in all.cpp, needed for vector initializations 
	
	
	// -------------------------------------------------------
	// generally used mathematical functions
	// only decalaration, implemenetation is in all.cpp
	// -------------------------------------------------------
	/** some templates need math functions for cplx/double
	 * unfortunately, the floating point absolute value for double is 
	 * fabs while its abs for complex values. therefore, using 
	 * abs(T) is a template is not possible ... 
	 */		
	class negf_math {
	public:
		static double abs(double val); 						//!< absolute value
		static int    sign(double a);						//!< +1, -1 or 0
		static int    sign(int a);							//!< +1, -1 or 0
		static inline int max(const int & a, const int & b) { return int(std::max(double(a), double(b))); }	//!< the bigger of both
        static inline uint max(const uint & a, const uint & b) { return (uint) negf_math::max((int)a, (int)b); } //!< the bigger of both
        static inline double max(const double & a, const double & b) { return std::max(a, b); } //!< the bigger of both
		static inline int min(const int & a, const int & b) { return int(std::min(double(a), double(b))); }	//!< the smaller of both
        static inline uint min(const uint & a, const uint & b) { return (uint) negf_math::min((int)a, (int)b); } //!< the smaller of both
        static inline double min(const double & a, const double & b) { return std::min(a, b); } //!< the smaller of both
		static double exp(double x);						//!< see all.cpp whether fastexp of ACML is used or not
		static double log(double x);						//!< natural logarithm (base e)
		static double sqrt(double x);						//!< square root
		static double pow(double base, double exponent);	//!< power (inefficient!)
		static double sin(double x);						//!< argument in radians
		static double cos(double x);						//!< argument in radians
		static double tan(double x);						//!< argument in radians
		static double atan(double x);
		static double tanh(double x);						//!< output is in radians
		static cplx   acos(const cplx & z);
		static double fermi_int(double order, double x);
		static double cheval(uint n, const double A[], double t);
		static double fermihalf_inverse(double y);
		static double fermim0p5_inverse(const double & y);
		static double vector_norm_3d(const double vec[3]);	//!< user needs to make sure that pointer really has 3 entries
		static double vector_norm_2d(const double vec[2]);	//!< user needs to make sure that pointer really has 2 entries
		static void   vector_normalize_3d(double vec[3]);	//!< input quantity is overwritten
		static void   vector_normalize_2d(double vec[2]);	//!< input quantity is overwritten
		static void   vector_outer_product(const double a[3], const double b[3], double result[3]);	//!< outer (cross) product result = a x b
		static double vector_scalar_product_3d(const double a[3], const double b[3]);	//!< inner (scalar) product of two 3-vectors
		static double vector_scalar_product_2d(const double a[2], const double b[2]);	//!< inner (scalar) product of two 2-vectors
		static double vector_norm(const std::vector<double> & vec);						//!< 2-norm of a vector of arbitrary length
		static bool   solve_linear_equation(double a11, double b1, double& x1);			//!< solve ax=b for x (scalar equation)
		static bool   solve_linear_equation(double a11, double a12, double b1,
							 				double a21, double a22, double b2, double& x1, double& x2);	//!< solve Ax=b for x (A is 2x2 matrix)
		static bool   solve_linear_equation(double a11, double a12, double a13, double b1,
			 								double a21, double a22, double a23, double b2,
			 								double a31, double a32, double a33, double b3,
			 								double& x1, double& x2, double& x3);	//!< solve Ax=b for x (A is 3x3 matrix)
		static double matrix_norm(const Matc & matrix); 	//!< Frobenius norm sqrt(sum_ij |A_ij|^2)
		static double matrix_norm(const Matd & matrix); 	//!< Frobenius norm sqrt(sum_ij |A_ij|^2)
		static double matrix_norm(const BMatc & matrix); 	//!< Frobenius norm sqrt(sum_ij |A_ij|^2)
		static void   zgeev(int n, cplx * matrix, cplx * eigenvalues, cplx * eigenvectors); 
		
		// used in GreenFunctions, SEOpticalPhonon and SEPhotonSpontaneous!
		// all of those functions allocate heap memory for real_data, imag_data, real_compressed and imag_compressed
		// (except do_compress_helper which works with the already set up real_data, imag_data)
		static void do_compress_antiherm(const vector<BMatc> & data, double* & real_data, double* & imag_data, unsigned char* & real_char, unsigned char* & imag_char, unsigned long & num_chars, 
					unsigned long & real_comp_size, unsigned long & imag_comp_size, unsigned char* & real_compressed, unsigned char* & imag_compressed);
		template<class Mtype>
		static unsigned long check_and_get_num_doubles(const vector<Mtype> & data, bool antiherm);
		template<class Mtype>
		static void do_compress(const vector<Mtype> & data, double* & real_data, double* & imag_data, unsigned char* & real_char, unsigned char* & imag_char, unsigned long & num_chars, 
					unsigned long & real_comp_size, unsigned long & imag_comp_size, unsigned char* & real_compressed, unsigned char* & imag_compressed);
		//static void do_compress(const vector<Matc> & data, double* & real_data, double* & imag_data, unsigned char* & real_char, unsigned char* & imag_char, unsigned long & num_chars, 
		//			unsigned long & real_comp_size, unsigned long & imag_comp_size, unsigned char* & real_compressed, unsigned char* & imag_compressed);
		static void do_compress_helper(double* & real_data, double* & imag_data, unsigned long & num_doubles, unsigned char* & real_char, unsigned char* & imag_char, unsigned long & num_chars, 
						unsigned long & real_comp_size, unsigned long & imag_comp_size, unsigned char* & real_compressed, unsigned char* & imag_compressed);
		static void check_compress_error(int err, unsigned long compressed_size, unsigned long num_chars);
		static void check_decompress_error(int err, unsigned long decomp_size, unsigned long num_chars);
	};
	
	// SAVE TO FILE
	void write_matrix(const char* filename, const Matc & matrix, const string & description = "");	//!< save a full complex matrix to file
	void write_matrix(const char* filename, const Matd & matrix, const string & description = "");	//!< save a full real matrix to file
	void read_matrix(const char* filename, Matc & matrix);											//!< read a full complex matrix from file
	void write_xE_matrix(const char* filename, const Matd & matrix, const vector<double> & xcoord, const vector<double> & energies);	//!< same some quantity f(x,E) to file, including x- and E-vectors
	void write_current_matrix(const char* filename, const Matd & matrix, const vector<double> & xcoord, const vector<double> & energies);	//!< same a quantity f(x,E) to file where there are Nx-1 entries, including x- and E-vectors
}

#include "Mat.h" 
//#include "MPI.h"  commented out by S.Z.
#include "MPIs.h"

using namespace negf;

template<class Mtype>
unsigned long negf_math::check_and_get_num_doubles(const vector<Mtype> & data, bool antiherm)
{STACK_TRACE(
	int num_matrices = data.size();
	NEGF_ASSERT(num_matrices > 0, "there is nothing to send.");
	
	const int max_distance = data[0].num_offdiags;
	// max_distance needs to be the same as in myMPI::recv(vector<Mtype> & data, int & source, int & tag)

	int n = data[0].num_rows();
	for (int ii=0; ii < num_matrices; ii++) {
		NEGF_FASSERT(int(data[ii].num_rows())==n && int(data[ii].num_cols())==n, 
				"expect array of equal size (%d) square matrices; instead, matrix %d is %dx%d.",
				n, ii, data[ii].num_rows(), data[ii].num_cols());
	}

	int n2 = negf_math::max(0, n-max_distance-1);	// size of triagonal slice which is NOT needed
	unsigned long N = 0;
	if (antiherm) {
		N = n*(n+1)/2 - n2*(n2+1)/2;	// only half of the matrix (including diagonal) is needed
	} else {
		N = n*n - n2*(n2+1);			// the whole matrix is needed
	} 
	NEGF_FASSERT(N>0, "N=%d! n=%d, n2=%d",N,n,n2);
	unsigned long num_doubles = N*num_matrices;
	return num_doubles;
);}


template<class Mtype>
void negf_math::do_compress(const vector<Mtype> & data, double* & real_data, double* & imag_data, unsigned char* & real_char, unsigned char* & imag_char, unsigned long & num_chars, 
					unsigned long & real_comp_size, unsigned long & imag_comp_size, unsigned char* & real_compressed, unsigned char* & imag_compressed)
{STACK_TRACE(
	unsigned long num_doubles = negf_math::check_and_get_num_doubles(data, false);
	const int max_distance = data[0].num_offdiags;
	// max_distance needs to be the same as in myMPI::recv(vector<Mtype> & data, int & source, int & tag)
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
			for (int jj=max(0,ii-max_distance); jj<=min(n-1,ii+max_distance); jj++) {
				real_data[count] = data[mm](ii+1,jj+1).real();
				imag_data[count] = data[mm](ii+1,jj+1).imag();
				//cout << "matrix " << mm << " count=" << count << " (i=" << ii+1 << ",j=" << jj+1 << "): data=" << data[mm](ii+1,jj+1) << ", real=" << real_data[count] << ", imag=" << imag_data[count] << endl;
				count++;
			}
		}
	}
	NEGF_ASSERT(count==num_doubles, "something went wrong.");
	
	negf_math::do_compress_helper(real_data, imag_data, num_doubles, real_char, imag_char, num_chars, real_comp_size, imag_comp_size, 
			real_compressed, imag_compressed);
);}


#endif /*_ALL_H_NEGF*/
