options {

	temperature                   = 300

	kp_method                     = 0
		
	StrainPolarization 			  = 0
	
	Transmission                  = 0
	QuasiFermilevels              = 0
	
	# ---------------
	# Numerics
	# ---------------
	
	maximal_k_value               = 1.4 # in 1/nm
	num_k_points                  = 40
#	num_k_points                  = 1
	
	num_energy_points             = 300
	
	min_energy                    = 1.3	# for n-doping
	max_energy                    = 2.1
#	min_energy                    = -0.6	# for p-doping
#	max_energy                    = 0.2
# remember to also change the doping file when switching between n- and p-doping!
# and change kpmethod=2!	

	PotentialUnderrelaxation      = 0.1
	PulayMixing                   = 0
	Resonances                    = 0
	PMLResonances                 = 0
	IncludeImaginaryContactStates = 0
	IncoherentContacts            = 1
	IncoherentContactBroadening   = 1e-2
	
	inner_errcrit                 = 1e-8
	outer_convergence_crit        = 1e-4
	max_outer_iterations          = 20
	max_inner_iterations          = 40
	outer_loop_monitoring         = 1
	
	# ---------------
	# Self Energies
	# ---------------
	ScatteringDecreaseFactor    = 0.05	# S(0) = (this)*S(end)
	ScatteringRampFactor        = 2.0	# S(i) * (this)*S(i-1)
	
	Buettiker                   = 0
	BuettikerParameter          = 0.01  # [eV]
	
	GolizadehMomentumRelaxation = 0
	GolizadehMomentumParameter  = 100e-6 # [eV^2]
	
	# Underrelaxation is mandatory for dephasing
	# SelfEnergyUnderrelaxation   = 0.5
	GolizadehDephasing          = 0
	GolizadehDephasingParameter = 0.0001 # [eV^2] sqrt(0.0001)=10meV 
	# ramping doesnt really pay off w/ Golizadeh dephasing
	# self-energy is completely differeent after each Poisson step
	# because it solely smoothes out the DoS
	
	OpticalPhonons              = 1
	LuisierSRpop                = 0
	
	AcousticPhonons             = 0
	
	IonizedImpurities           = 0
}


# -----------------------
# Voltage ramping
# Contact 0 = Source = left
# Contact 1 = Drain = right
# -----------------------

experiment_0 {
	Source_voltage = 0.0
	Drain_min      = 0.0
	Drain_max      = -0.101
	Drain_step     = -0.05
#	Drain_max      = +0.201
#	Drain_step     = +0.05
	PotentialUnderrelaxation = 0.1
	LateUnderrelaxation      = 0.0
}


# -----------------------
# Structure definition
# -----------------------

regions {
	
	region0_length = 12
	region0_dx     = 0.8
	region0_mat    = 0     # 0 = GaAs
	region0_molefr = 0.0
	region0_doping = 1e19
	
	region1_length = 32
	region1_dx     = 0.8
	region1_mat    = 0     # 0 = GaAs
	region1_molefr = 0.0
	region1_doping = 1e18
	
	region2_length = 12
	region2_dx     = 0.8
	region2_mat    = 0     # 0 = GaAs
	region2_molefr = 0.0
	region2_doping = 1e19
}

