import sys
import os
from socket import gethostname 
sys.path.append(os.curdir);
os.putenv("TDKPCONFPATH",os.curdir+"/tdkp_conf");
from negf import *

try:	
	cvar.logmsg.add_device("std::cout");
	if len(sys.argv) != 2 :
		raise NameError, ("you have supplied ",len(sys.argv),"instead of exactly 1 arguments.")
	name = sys.argv[1];
	cvar.fnames.init(name);
	
	#cvar.fnames.set_outfile_directory_suffix("_Nk20_NE800");
	#cvar.fnames.set_outfile_directory_suffix("_golim0.0001");
	#cvar.fnames.set_outfile_directory_suffix("_ballistic_"+str(cvar.mpi.get_num_procs())+"procs");
	#cvar.fnames.set_outfile_directory_suffix("_POP_"+str(cvar.mpi.get_num_procs())+"procs");
	#cvar.fnames.set_outfile_directory_suffix("_phot_"+str(cvar.mpi.get_num_procs())+"procs");
	#cvar.fnames.set_outfile_directory_suffix("_GaN_NOFrey_decay0.2nm");
	
	cvar.fnames.set_outfile_directory_suffix("_POP0.05_pol0.5_AC_II_Goli25e-6");

	# extend path for module searching by current directory.
	# this is necessary when the given script file is a symbolic link
	sys.path.append(cvar.fnames.get_path());
	# TDKP configuration
	os.putenv("TDKPCONFPATH",cvar.fnames.get_path()+"conf");
	
	#if cvar.mpi.get_num_procs()==1:
	#	print "start with mpiexec -n <num_procs> python scripts/run.py structures/...";
	#	raise NameError;
	
	# ----------------------------------------------------------------------------------------------
	# master thread only:
	# test existence of necessary directories, create if necessary
	if cvar.mpi.get_rank()==constants.mpi_master_rank:
		if os.access(cvar.fnames.get_homedirectory(), os.F_OK) == False:
			print "Home directory ",cvar.fnames.get_homedirectory(),"not found. Aborting.";
			raise NameError;
		if os.access(cvar.fnames.get_filedirectory(), os.F_OK) == False:
			print "Directory of input files",cvar.fnames.get_filedirectory(),"not found. Aborting.";
			raise NameError;
		if os.access(cvar.fnames.get_materialdirectory(), os.F_OK) == False:
			print "Directory of material parameters",cvar.fnames.get_materialdirectory(),"not found. Aborting.";
			raise NameError;
		if os.access(cvar.fnames.get_logfiledirectory(), os.F_OK) == False:
			print "Creating directory for logfiles",cvar.fnames.get_logfiledirectory();
			os.mkdir(cvar.fnames.get_logfiledirectory());
		if os.access(cvar.fnames.get_path(), os.F_OK) == False:
			print "Creating base directory for output",cvar.fnames.get_path();
			os.mkdir(cvar.fnames.get_path());
		if os.access(cvar.fnames.get_outfiledirectory(), os.F_OK) == False:
			print "Creating output directory",cvar.fnames.get_outfiledirectory();
			os.mkdir(cvar.fnames.get_outfiledirectory());
		# create a symbolic link to the gridfile in the output directory for convenience
		symbolic_link_file = cvar.fnames.get_outfiledirectory() + cvar.fnames.get_barefile() + ".grd";
		if (os.access(symbolic_link_file, os.F_OK) == False) and (os.name == "posix"):
			print "Creating symbolic link to gridfile in output directory:",symbolic_link_file;
			os.symlink(os.getcwd()+"/"+cvar.fnames.get_filename()+".grd", symbolic_link_file);
	cvar.mpi.synchronize_processes();
	
	if cvar.mpi.get_rank()==constants.mpi_master_rank:
		cvar.logmsg.add_device(cvar.fnames.get_logfile());
	
	cvar.logmsg.emit_all(LOG_INFO,"NEGF is running on %s (p%d)" % (gethostname(),cvar.mpi.get_rank()));
	cvar.logmsg.set_level(LOG_INFO_L1);
	cvar.mpi.synchronize_processes();cvar.mpi.synchronize_processes();
	cvar.logmsg.emit_huge_header("*** Welcome to SEGF (= Sebis NEGF Solver) ***");
	
	# ----------------------------------------------------------------------------------------------
	# read in options (bandstructure model, which scattering is included, how many E- and k-points)
	# sets Nn which is needed for everything else
	# ----------------------------------------------------------------------------------------------
	options = Options();
	
	# --------------------------------------------------------------
	# create energies
	# this also stores which energy is stored in which MPI process
	# --------------------------------------------------------------
	energies = Energies(options);
	
	# --------------------------------------------------------------
	# read in real-space grid from file (all threads)
	# --------------------------------------------------------------
	cvar.logmsg.emit_header('setting up real-space grid');
	parser = InputParser();
	xspace = parser.read_dfise_grd(cvar.fnames.get_filename());
	cvar.logmsg.emit(LOG_INFO,'preparing molefractions...');
	parser.prepare_molefractions(xspace, cvar.fnames.get_filename());	
	
	# --------------------------------------------------------------
	# read all relevant material data (all threads)
	# and assign materials to region
	# --------------------------------------------------------------
	cvar.logmsg.emit_header("reading material data");
	material = MaterialDatabase();
	cvar.logmsg.emit(LOG_INFO,"Materials directory is %s" % cvar.fnames.get_materialdirectory());
	#cvar.logmsg.set_level(LOG_INFO_L3);
	material.add_search_path(cvar.fnames.get_materialdirectory());
	material.add_search_path(cvar.fnames.get_filedirectory());
	#bulk.prepare_molefractions();
	for ii in range(xspace.get_num_regions()):
		reg = xspace.get_region(ii);
		mat_name = reg.get_material_name();
		if (reg.has_molefraction()):
			mat_name = mat_name + '%g' % reg.get_material_molefraction();
			if (not material.material_exists(mat_name)):
				material.load_ternary_material(reg.get_material_name(), reg.get_material_molefraction());
			reg.set_material_name(mat_name);
		else:
			if (not material.material_exists(mat_name)):
				material.load_material(mat_name);
		reg.set_material(material.get_material(mat_name));
	
	# --------------------------------------------------------------
	# set up k-space automatically (all threads)
	# --------------------------------------------------------------
	kspace  = Kspace(xspace, options, material);
		
	# --------------------------------------------------------------
	# read in voltages (all threads)
	# --------------------------------------------------------------
	voltages = Voltages(xspace);
	
	# ---------------------------------------------------
	# set up Poisson solver
	# ---------------------------------------------------
	poiss = PoissonProblem(xspace, material, options);
	
	# --------------------------------------------------------------
	# initialize Hamiltonian (or rather, interface to TDKP Hamilt.)
	# --------------------------------------------------------------
	ham = Hamiltonian(xspace, kspace, options, material, cvar.fnames.get_filename());
	ham.set_electrostatic_potential(poiss.get_poisson_equation().get_values());
	#ham.test_flens();
	
	# --------------------------------------------------------------
	# initialize FEM overlap matrix (also obtained from TDKP)
	# --------------------------------------------------------------
	overlap = Overlap(ham, xspace);
	
	# --------------------------------------------------------------
	# initialize Green's functions GR, GL, GA
	# each thread stores different GF
	# --------------------------------------------------------------
	gf = GreenFunctions(options, xspace, kspace, energies, ham, overlap);
	
	# --------------------------------------------------------------
	# initialize self-energies (allocate memory etc.)
	# each thread stores different SE
	# --------------------------------------------------------------
	se = SelfEnergies(ham, overlap, xspace, kspace, energies, options, gf, material);
	if (se.has_self_energy(SEtype_ion_imp)):
		se.get_ionized_impurities_selfenergy().set_up(poiss.get_doping(), poiss.get_dielectric(), poiss.get_effective_edos(), poiss.get_effective_hdos());
	
	#se.get_optical_phonon_selfenergy().test_operation();

	# -----------------------------------------------------------------------------------------------
	# initialize post-processing object (DoS, spectral density, density, spectral current, current)
	# -----------------------------------------------------------------------------------------------
	pp = PostProcessing(ham, overlap, gf, se, options, xspace, kspace, energies);
	if (se.has_self_energy(SEtype_spont_photon)):
		pp.set_up_luminescence(se.get_spontaneous_photon_selfenergy());
	
	# --------------------------------------------------------------
	# set up inner loop (self-consistent GF and SE) 
	# each thread computes a certain energy window
	# includes Dyson and Keldysh equations
	# --------------------------------------------------------------
	inner = InnerLoop(ham, overlap, gf, se, pp, material);
	
	# ---------------------------------------------------------------------------------------
	# set up outer loop (all threads - during Poisson all threads but the master are idle)
	# ---------------------------------------------------------------------------------------
	outer = OuterLoop(inner, poiss);
	outer.output_debug(True);
	
	# ========================================================================
	# ************************ PERFORM THE SIMULATION ************************
	# ========================================================================		
	cvar.logmsg.set_level(LOG_INFO_L1);
	for step in range(voltages.get_num_steps()):
		cvar.timer.click("voltage"+str(step)+"_start");
		
		# some screen output
		voltages.emit_header(step);
		
		# update applied voltage for boundary self-energies and Poisson
		cvar.logmsg.emit(LOG_INFO,"*** Updating voltages... ***");
		outer.update_voltages(voltages.get_setting(step));
		if voltages.get_underrelaxation(step)>=0.0: # if no specific parameter is set for the step, -1.0 is returned
			outer.update_underrelaxation(voltages.get_underrelaxation(step));
		if voltages.get_late_underrelax(step)>=0.0: # if no specific parameter is set for the step, -1.0 is returned
			outer.update_late_underrelax(voltages.get_late_underrelax(step));
		
		# initial guess if this is the first step
		if step==0:
			se.initial_guess();  	 # computes contact self-energy
			inner.initial_guess();   # computes retarded and advanced GF
		
		# PERFORM OUTER LOOP (hidden)
		outer.perform();
		
		# output results
		cvar.logmsg.emit_header("*** Output... ***");
		pp.compute_scattering_current(SEtype_contact);
		if (se.has_self_energy(SEtype_optical_phonon)):
			pp.compute_scattering_current(SEtype_optical_phonon);
		if (se.has_self_energy(SEtype_spont_photon)):
			pp.compute_scattering_current(SEtype_spont_photon);
			pp.compute_luminescence();
		if (options.exists("Transmission") and options.get("Transmission")==1):
			pp.compute_transmission();
			
		outputdata = poiss.get_output_data();
		output_filename = outputdata.get_filename() + voltages.get_suffix(step);
		if cvar.mpi.get_rank()==constants.mpi_master_rank:
			outer.save_everything_to_file(output_filename);
			pp.get_contact_current().snapshot(voltages.get_ramped_voltage(step));
			pp.get_contact_current().save_to_file(output_filename);
			if (se.has_self_energy(SEtype_spont_photon)):
				pp.get_luminescence2().snapshot(voltages.get_ramped_voltage(step));
				pp.get_luminescence2().write_power_to_file(output_filename);
		
		pp.output_selfenergies(1.75, output_filename);
		cvar.mpi.synchronize_processes();
		
		cvar.timer.click("voltage"+str(step)+"_stop");
		voltage_time = cvar.timer.get_seconds_between("voltage"+str(step)+"_start", "voltage"+str(step)+"_stop");
		cvar.logmsg.emit(LOG_INFO, "The voltage step needed %7.2f seconds." % voltage_time);
		
		
	cvar.timer.click("finish");
	time_in_sec = cvar.timer.get_seconds_since_start("finish");
	cvar.logmsg.emit(LOG_INFO, "The simulation finished in %7.2f seconds." % time_in_sec);
	cvar.logmsg.emit(LOG_INFO, "All of this was written into %s." % cvar.fnames.get_logfile());
	cvar.logmsg.emit(LOG_INFO,"Thank you for using NEGF. Please report bugs and comments to steiger@iis.ee.ethz.ch.");
	cvar.logmsg.emit(LOG_INFO,"====================================================================================");
	cvar.mpi.terminate();
		
except RuntimeError, e:
	print "RuntimeError caught!"
	print e
	cvar.logmsg.emit_all(LOG_INFO, "Thread: %d" % cvar.mpi.get_rank());
	cvar.logmsg.emit_all(LOG_INFO, "Logfile: %s" % cvar.fnames.get_logfile());
	cvar.mpi.terminate();
	sys.exit();

except StandardError, e:
	print e
	cvar.logmsg.emit_all(LOG_INFO, "Thread: %d" % cvar.mpi.get_rank());
	cvar.logmsg.emit_all(LOG_INFO, "Logfile: %s" % cvar.fnames.get_logfile());
	cvar.mpi.terminate();
	sys.exit();
