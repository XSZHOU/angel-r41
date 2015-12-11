import sys
import os
from socket import gethostname 
sys.path.append(os.curdir);
os.putenv("TDKPCONFPATH",os.curdir+"/tdkp_conf");
from negf import *

try:	
	if len(sys.argv) != 2 :
		raise NameError, ("you have supplied ",len(sys.argv),"instead of exactly 1 arguments.")
	name = sys.argv[1];
	cvar.fnames.init(name);
	
	# extend path for module searching by current directory.
	# this is necessary when the given script file is a symbolic link
	sys.path.append(cvar.fnames.get_path());
	# TDKP configuration
	os.putenv("TDKPCONFPATH",cvar.fnames.get_path()+"conf");
	
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
	
	cvar.logmsg.add_device("std::cout");
	if cvar.mpi.get_rank()==constants.mpi_master_rank:
		cvar.logmsg.add_device(cvar.fnames.get_logfile());
	
	cvar.logmsg.emit(LOG_INFO,"NEGF is running on %s (p%d)" % (gethostname(),cvar.mpi.get_rank()));
	cvar.logmsg.set_level(LOG_INFO_L1);
	cvar.mpi.synchronize_processes();
	cvar.logmsg.emit_huge_header("*** Welcome to SEGF (= Sebis NEGF Solver) ***");
		
	# ----------------------------------------------------------------------------------------------
	# read in options (bandstructure model, which scattiering is included, how many E- and k-points)
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
	parser.prepare_molefractions(xspace, cvar.fnames.get_filename());
	
	
	# --------------------------------------------------------------
	# read all relevant material data (all threads)
	# and assign materials to region
	# --------------------------------------------------------------
	cvar.logmsg.emit_header("reading material data");
	material = MaterialDatabase();
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
	
	myNe = energies.get_my_number_of_points();
	Nk   = kspace.get_number_of_points();
	Nx   = xspace.get_num_vertices();
	Nn   = cvar.Nn;
	numNEGF = 9 + 3 + 3; # GF + SEtot + SEcont, each retarded, lesser, greater, and advanced GF, MGM
	cvar.logmsg.emit_header('memory consumption estimate');
	cvar.logmsg.emit(LOG_INFO,'There are %d processes.' % cvar.mpi.get_num_procs());
	cvar.logmsg.emit(LOG_INFO,'Each process will store %d NEGF objects.' % numNEGF);
	cvar.logmsg.emit(LOG_INFO,'Each NEGF object consists of %d*%d full complex matrices of size (%d*%d)^2.' \
					%(myNe, Nk, Nx, Nn));
	GB_per_process = numNEGF * myNe*Nk * (Nx*Nn*Nx*Nn) *16.0 / (1024.*1024.*1024.);
	cvar.logmsg.emit(LOG_INFO,'Memory consumption estimate per process: %.3g GB.' % GB_per_process);
	cvar.logmsg.emit(LOG_INFO,'Total memory used (estimate): %3g GB.' % (GB_per_process * cvar.mpi.get_num_procs()))
	
	# --------------------------------------------------------------
	# read in voltages (all threads)
	# --------------------------------------------------------------
	voltages = Voltages(xspace);
	
	# --------------------------------------------------------------
	# initialize Hamiltonian (or rather, interface to TDKP Hamilt.)
	# --------------------------------------------------------------
	ham = Hamiltonian(xspace, kspace, options, material, cvar.fnames.get_filename());
	#ham.test_flens();
	
	# --------------------------------------------------------------
	# initialize FEM overlap matrix (also obtained from TDKP)
	# --------------------------------------------------------------
	overlap = Overlap(ham.get_tdkp_interface(), xspace);
	
	# --------------------------------------------------------------
	# initialize Green's functions GR, GL, GA
	# each thread stores different GF
	# --------------------------------------------------------------
	gf = GreenFunctions(options, xspace, kspace, energies, overlap);
	
	# --------------------------------------------------------------
	# initialize self-energies (allocate memory etc.)
	# each thread stores different SE
	# --------------------------------------------------------------
	se = SelfEnergies(ham, overlap, xspace, kspace, energies, options, gf, material);
	
	# -----------------------------------------------------------------------------------------------
	# initialize post-processing object (DoS, spectral density, density, spectral current, current)
	# -----------------------------------------------------------------------------------------------
	pp = PostProcessing(ham, overlap, gf, se, options, xspace, kspace, energies);
	if (se.has_self_energy(SEtype_spont_photon)):
		pp.set_up_luminescence(se.get_spontaneous_photon_selfenergy());
	else:
		raise NameError;
	
	Ec = 1.422482;
	Ev = 0.000000;
	#Ec = 1.022482;
	#Ev = -0.400000;
	EFn = Ec + 0.05;
	EFp = Ev + 0.02;
	gf.assign_flat_band_lesser_greater(EFn, EFp, Ec, Ev);
	se.get_spontaneous_photon_selfenergy().calculate();
	# no local DOS because GR, SigmaR were not calculated!
	se.get_spontaneous_photon_selfenergy().output_debug_info();
	pp.compute_spectral_edensity();
	pp.compute_spectral_hdensity();
	pp.compute_edensity();
	pp.compute_hdensity();
	pp.get_luminescence2().calculate();
	#pp.get_luminescence().calculate(); # does not make sense because GR, GA are needed!
	
	if cvar.mpi.get_rank()==constants.mpi_master_rank:
		parser.write_xE_matrix('test_nE',   pp.get_entire_spectral_edensity(), xspace, energies.get_energy_grid());
		parser.write_xE_matrix('test_pE',   pp.get_entire_spectral_hdensity(), xspace, energies.get_energy_grid());
		# no local DOS because GR, SigmaR were not calculated!
		pp.get_luminescence2().write_recombination_to_file('test_radrec');
		pp.get_luminescence2().write_spectrum_to_file('test_spectrum');
		#pp.get_luminescence().write_recombination_to_file('test_radrec1');
		#pp.get_luminescence().write_spectrum_to_file('test_spectrum1');
	
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
