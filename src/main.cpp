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
#include "Logger.h"
#include "Timer.h"
#include "Filenames.h"
#include "Interrupt.h"
#include "MPIs.h"

#include "InputParser.h"
#include "MaterialDatabase.h"
#include "PropertyContainer.h"

#include "MaterialDatabase.h"

#include "Vertex.h"
#include "Edge.h"
#include "Element.h"
#include "Region.h"
#include "Contact.h"
#include "Geometry.h"
#include "BoxMethod.h"

#include "Voltages.h"

#include "Energies.h"
#include "Kspace.h"

#include "Options.h"
#include "NEGFObject.h"
#include "GreenFunctions.h"

#include "SelfEnergy.h"
#include "SEContacts.h"
#include "SEBuettiker.h"
#include "SEGolizadeh.h"
#include "SEOpticalPhonon.h"
#include "SEAcousticPhonon.h"
#include "SEPhotonSpontaneous.h"
#include "SelfEnergies.h"

#include "Current.h"
#include "Luminescence.h"
#include "PostProcessing.h"
#include "Hamiltonian.h"
#include "ContactFermilevel.h"
#include "Overlap.h"
#include "InnerLoop.h"
#include "PoissonProblem.h"
#include "OuterLoop.h"

// for mkdir
#include <sys/stat.h>
#include <sys/types.h>

using namespace negf;

int main(int argc, char* argv[])
{
	try {
		
	logmsg->add_listener(std::cout);
	
	for (int ii=0; ii<argc; ii++) {
		logmsg->emit(LOG_WARN,"arg=%d: %s", ii, argv[ii]);
	}
	
	NEGF_ASSERT(argc>1,"you must supply a simulation name as an argument.");
	// 0th argument = run_negf.bin
	// 1st argument = simulation name
	// 2nd argument would be logfile including directory
	// 3rd argument would be output directory
	
	//setenv("NEGFDIR", "/scratch/steiger",1);	// filebase for output!
	char * name = argv[1];
	fnames->init(name);
		
	// log file
	if (mpi->get_rank()==constants::mpi_master_rank) {
		if (argc==2) {
			// create directory for logfiles, if necessary
			logmsg->emit(LOG_INFO,"Creating %s for log files",fnames->get_logfiledirectory().c_str());
			mkdir(fnames->get_logfiledirectory().c_str(),0777);	
		}
		
		string logname = (argc>2) ? string(argv[2]) : fnames->get_logfile(); // 2nd argument would be logfile name including directory
		ofstream  * fout = new ofstream(logname.c_str()/*, ios::app*/);
		NEGF_FASSERT(*fout, "logfile %s could not be created/opened.", logname.c_str());
		logmsg->add_listener(*fout);		
		logmsg->emit(LOG_INFO,"logfile: %s",logname.c_str());
	}
	
	// adjustment of output directory
	//fnames->set_outfile_directory_suffix("_Nk60_NE600");
	if (argc>3) {
		fnames->set_outfile_directory(argv[3]); // 3rd argument would be output directory
	}
	
	// extend path for module searching by current directory. necessary when the given script file is a symbolic link
	//sys.path.append(cvar.fnames.get_path());
	// TDKP configuration
	string tdkpcnfpath(fnames->get_path()); tdkpcnfpath.append("conf");
	setenv("TDKPCONFPATH",tdkpcnfpath.c_str(),1);
	std::cout<<"work started !!!"<<std::endl;	
	// ----------------------------------------------------------------------------------------------
	// master thread only:
	// test existence of necessary directories, create if necessary
	if (mpi->get_rank()==constants::mpi_master_rank) {
		// create base directory for output, if necessary
		logmsg->emit(LOG_INFO,"Creating %s",fnames->get_path().c_str());
		mkdir(fnames->get_path().c_str(),0777);
		// create output directory for specific simulation, if necessary
		logmsg->emit(LOG_INFO,"Creating %s for output",fnames->get_outfiledirectory().c_str());
		mkdir(fnames->get_outfiledirectory().c_str(),0777);
		// no symbolic link creation
	}
	mpi->synchronize_processes();
	
	
	char hostnam[1000]; gethostname(hostnam,1000);
	logmsg->emit_noendl_all(LOG_INFO,"%s (p%d), ", hostnam ,mpi->get_rank());
	logmsg->set_level(LOG_INFO_L1);
	mpi->synchronize_processes();
	logmsg->emit(LOG_INFO,"");
	logmsg->emit_huge_header("*** Welcome to ANGEL (A NEGF solver for LEDs) ***");
	logmsg->emit(LOG_INFO,"Materials directory : %s", fnames->get_materialdirectory().c_str());
	logmsg->emit(LOG_INFO,"Simulation directory: %s", fnames->get_filedirectory().c_str());
		
	// read in options (bandstructure model, which scattiering is included, how many E- and k-points)
	Options * options = new Options();
	std::cout<<"Option reading initialated !!!"<<std::endl;
	// create energies. this also stores which energy is stored in which MPI process
	Energies * energies = new Energies(options);
	
	// read in real-space grid from file (all threads)
	logmsg->emit_header("setting up real-space grid");
	InputParser * parser = new InputParser();

	// NEW --> create Geometry from .cmd-file
	Geometry * xspace = parser->read_grd(options->get_cmdfile()); // includes prepare_molefractions

	// read all relevant material data (all threads) and assign materials to region
	logmsg->emit_header("reading material data");
	MaterialDatabase * material = new MaterialDatabase();
	logmsg->emit(LOG_INFO,"Materials directory is %s", fnames->get_materialdirectory().c_str());
	material->add_search_path(fnames->get_materialdirectory().c_str());
	material->add_search_path(fnames->get_filedirectory().c_str());
	for (uint ii=0; ii<xspace->get_num_regions(); ii++) {
		Region * reg = xspace->get_region(ii);
                string mat_name = reg->get_material_name();
		logmsg->emit(LOG_INFO,"Region %s (material=%s)...",reg->get_name().c_str(), mat_name.c_str());

		if (reg->has_molefraction()) {
			char mf[1000]; sprintf(mf,"%g",reg->get_material_molefraction());
			mat_name.append(mf);
			if (!material->material_exists(mat_name.c_str())) {
				// for some extremely weird reason (likely to be a compiler issue) the following line is necessary to prevent a segfault
				// if included, valgrind shows no errors! if not included, valgrind does show an error. problem is not present in run.py.
				cout << endl; // cout.flush() or cout << "\n"  alone do not work
				
				material->load_ternary_material(reg->get_material_name().c_str(), reg->get_material_molefraction());
			}
			reg->set_material_name(mat_name.c_str());
		} else {
			if (!material->material_exists(mat_name.c_str())) {
				material->load_material(mat_name.c_str());
			}
		}
		reg->set_material(material->get_material(mat_name.c_str()));
	}
	
	// set up k-space automatically (all threads)
	Kspace * kspace = new Kspace(xspace, options, material);
		
	// read in voltages (all threads)
	Voltages * voltages = new Voltages(xspace);
	
	// set up Poisson solver
	// this also reads in the doping
	PoissonProblem * poiss = new PoissonProblem(xspace, material, options);
	//poiss->initial_guess();
	poiss->initial_guess(true, voltages->get_value(0), voltages->get_value(1));
	
	// initialize Hamiltonian (or rather, interface to TDKP Hamilt.)
	Hamiltonian * ham = new Hamiltonian(xspace, kspace, options, material, fnames->get_filename().c_str());
	ham->set_electrostatic_potential(poiss->get_electrostatic_potential());
	ham->set_strain(poiss->get_strain_polarization()); // strain corrections
	//ham->test_flens();
	
	// initialize FEM overlap matrix (also obtained from TDKP)
	Overlap * overlap = new Overlap(ham, xspace);
	
	// initialize Green's functions GR, GL, GA. each thread stores different GF
	bool mangle_overlap_calc = (kspace->get_number_of_points()==1) ? false : true;
	GreenFunctions * gf = new GreenFunctions(options, xspace, kspace, energies, ham, overlap, mangle_overlap_calc);
	
	// initialize self-energies (allocate memory etc.). each thread stores different SE
	SelfEnergies * se = new SelfEnergies(ham, overlap, xspace, kspace, energies, options, gf, material);
	if (se->has_self_energy(SEtype_ion_imp)) {
	    se->get_ionized_impurities_selfenergy()->set_up(poiss->get_doping(), poiss->get_dielectric(), poiss->get_Nc(), poiss->get_Nv());
	 }
	
	// initialize post-processing object (DoS, spectral density, density, spectral current, current)
	PostProcessing * pp = new PostProcessing(ham, overlap, gf, se, options, xspace, kspace, energies);
	if (se->has_self_energy(SEtype_spont_photon)) {
		pp->set_up_luminescence(se->get_spontaneous_photon_selfenergy());
	}
	
	// set up inner loop (self-consistent GF and SE). each thread computes a certain energy window. includes Dyson and Keldysh equations
	InnerLoop * inner = new InnerLoop(ham, overlap, gf, se, pp, material);
	if (options->exists("max_inner_iterations")) {
		uint new_max_inner_iters = (uint) options->get("max_inner_iterations");
		inner->set_max_inner_iterations(new_max_inner_iters);
	}
	
	// set up outer loop (all threads - during Poisson all threads but the master are idle)
	OuterLoop * outer = new OuterLoop(inner, poiss, fnames->get_outfile());
	outer->output_debug(false);
	if (options->exists("outer_loop_monitoring") && options->get("outer_loop_monitoring")==1) {
		outer->output_debug(true);
	}
	if (options->exists("max_outer_iterations")) {
		uint new_max_outer_iters = (uint) options->get("max_outer_iterations");
		outer->set_max_outer_iterations(new_max_outer_iters);
	}
	if (options->exists("outer_convergence_crit")) {
		outer->set_outer_convergence_crit(options->get("outer_convergence_crit"));
	}
	
	// ========================================================================
	// ************************ PERFORM THE SIMULATION ************************
	// ========================================================================		
	logmsg->set_level(LOG_INFO_L1);
	for (uint step=0; step<voltages->get_num_steps(); step++)
	{
		// some screen output
		voltages->emit_header(step);
		
		// update applied voltage for boundary self-energies and Poisson
		logmsg->emit(LOG_INFO,"*** Updating voltages.. ***");
		outer->update_voltages(voltages->get_setting(step));
		if (voltages->get_underrelaxation(step)>=0.0) { // if no specific parameter is set for the step, -1.0 is returned
			outer->update_underrelaxation(voltages->get_underrelaxation(step));
		}
		if (voltages->get_late_underrelax(step)>=0.0) { // if no specific parameter is set for the step, -1.0 is returned
			outer->update_late_underrelax(voltages->get_late_underrelax(step));
		}
		
		// initial guess if this is the first step
		if (step==0) {
			se->initial_guess();  	 // computes contact self-energy
			inner->initial_guess();  // computes retarded and advanced GF
		}
		
		// PERFORM OUTER LOOP (hidden)
		outer->perform();
		
		// output results
		logmsg->emit_header("*** Output... ***");
		if (se->has_self_energy(SEtype_optical_phonon)) {
			pp->compute_scattering_current(SEtype_optical_phonon);
		}
		if (se->has_self_energy(SEtype_spont_photon)) {
			pp->compute_scattering_current(SEtype_spont_photon);
			pp->compute_luminescence();
		}
			
		string output_filename = fnames->get_outfile() + voltages->get_suffix(step);
		if (mpi->get_rank()==constants::mpi_master_rank) {
			outer->save_everything_to_file(output_filename.c_str());
			pp->get_contact_current()->snapshot(voltages->get_ramped_voltage(step));
			pp->get_contact_current()->save_to_file(output_filename.c_str());
			if (se->has_self_energy(SEtype_spont_photon)) {
				pp->get_luminescence2()->snapshot(voltages->get_ramped_voltage(step));
				pp->get_luminescence2()->write_power_to_file(output_filename.c_str());
			}
		}
		pp->output_selfenergies(1.75, output_filename.c_str());
		mpi->synchronize_processes();
	}
		
		
	timer->click("finish");
	double time_in_sec = timer->get_seconds_since_start("finish");
	logmsg->emit(LOG_INFO, "The simulation finished in %7.2f seconds.", time_in_sec);
	logmsg->emit(LOG_INFO, "All of this was written into %s.", fnames->get_logfile().c_str());
	logmsg->emit(LOG_INFO, "Thank you for using NEGF. Please report bugs and comments to steiger@iis.ee.ethz.ch.");
	logmsg->emit(LOG_INFO, "====================================================================================");
	mpi->terminate();
	
	} catch (Exception* e) {
		mpi->synchronize_processes();
		logmsg->emit_all(LOG_INFO,"%s\nterminating\n",e->get_reason().c_str());
		exit(1);
	}
	return 0;
}
 

