options {

	temperature                   = 300
	
	kp_method                     = 0
	
	# ---------------
	# Numerics
	# ---------------
	
	maximal_k_value               = 1.2 # in 1/nm; 
	# note: start from hbar2*(1e9)^2/(2*m0) = 0.076eV
	# BZ is pi/5A = 6.3nm-1 wide
	num_k_points                  = 150
	
	min_energy                    = 0.8
	max_energy                    = 1.7
#	min_energy                    = 1.3
#	max_energy                    = 2.2
	num_energy_points             = 192
	
	PotentialUnderrelaxation      = 0.4
	PulayMixing                   = 0
	NonzeroBoundaryField          = 0 # (0/1) nonzero electric field for Neumann bnd.cond. for Poisson (like nextnano)
	Resonances                    = 0
	ResonancesAllK                = 0
	PMLResonances                 = 0
	IncludeImaginaryContactStates = 0
	IncoherentContacts            = 1
	IncoherentContactBroadening   = 0.005
	
	KSpaceDiscretization = 1 # 0=equal, 1=sqrt
	
	inner_errcrit                 = 1e-4
	outer_convergence_crit        = 1e-4
	max_outer_iterations          = 20
	max_inner_iterations          = 40
	outer_loop_monitoring         = 1
	
	# ---------------
	# Self Energies
	# ---------------
	ScatteringDecreaseFactor    = 0.125		# S(0) = (this)*S(end)
	ScatteringRampFactor        = 2.0			# S(i) * (this)*S(i-1)
	
	Buettiker                   = 0
	BuettikerParameter          = 0.01 # [eV]
	
	GolizadehMomentumRelaxation = 1
	GolizadehMomentumParameter  = 81e-6 # [eV^2]
	
	GolizadehDephasing          = 0
	GolizadehDephasingParameter = 0.0001 # [eV^2]
	
	OpticalPhonons              = 1 
	LuisierSRpop                = 0
	#PhononCheatFactor          = 0.01
	
	AcousticPhonons             = 1

	IonizedImpurities           = 1
	
	SpontaneousPhotons          = 0
}


# -----------------------
# Voltage ramping
# -----------------------

#experiment_0 {
#	Source_voltage = 0.0
#	Drain_min      = 0.000
#	Drain_max      = -0.099
#	Drain_step     = -0.025
#}
#experiment_1 {
#	Source_voltage = 0.0
#	Drain_min      = -0.100
#	Drain_max      = -0.199
#	Drain_step     = -0.025
#}
#experiment_2 {
#	Source_voltage = 0.0
#	Drain_min      = -0.200
#	Drain_max      = -0.299
#	Drain_step     = -0.02
#}
#experiment_3 {
#	Source_voltage = 0.0
#	Drain_min      = -0.300
#	Drain_max      = -0.399
#	Drain_step     = -0.01
#}
experiment_4 {
	Source_voltage = 0.0
	Drain_min      = -0.400
	Drain_max      = -0.701
	Drain_step     = -0.02
}

# -----------------------
# Structure definition
# -----------------------

regions {
	
	region0_length = 15.0
	region0_dx     = 1.0
	region0_mat    = 0      # GaAs
	region0_molefr = 0.0
	region0_doping = 2e18   # n
	
	region1_length = 3.0
	region1_dx     = 0.5
	region1_mat    = 0      # GaAs
	region1_molefr = 0.3
	region1_doping = 2e18   # n
	
	region2_length = 3.0
	region2_dx     = 0.5
	region2_mat    = 10     # AlGaAs barrier
	region2_molefr = 0.4
	region2_doping = 1e10   # i
	
	region3_length = 3.0
	region3_dx     = 0.5
	region3_mat    = 0      # GaAs well
	region3_molefr = 0.0
	region3_doping = 1e10   # i
	
	region4_length = 3.0
	region4_dx     = 0.5
	region4_mat    = 10     # AlGaAs barrier
	region4_molefr = 0.4
	region4_doping = 1e10   # i
	
	region5_length = 3.0
	region5_dx     = 0.5
	region5_mat    = 0      # GaAs
	region5_molefr = 0.0
	region5_doping = 2e18   # n
	
	region6_length = 15.0
	region6_dx     = 1.0
	region6_mat    = 0      # GaAs
	region6_molefr = 0.0
	region6_doping = 2e18   # n
}
