options {

	temperature                   = 300
	
	kp_method                     = 3
	
	maximal_k_value               = 1.5 # in 1/nm
	# note: start from hbar2*(1e9)^2/(2*m0) = 0.076eV and multiply w/ kmax^2/m*
	# BZ is pi/5A = 6.3nm-1 wide
	num_k_points                  = 70
	
	min_energy                    = -0.5
	max_energy                    = 4.5
	num_energy_points             = 380
	
	StrainPolarization            = 1
	PolarizationDecreaser         = 0.5

	PotentialUnderrelaxation      = 0.6
	NonzeroBoundaryField          = 0 # (0/1) nonzero electric field for Neumann bnd.cond. for Poisson (like nextnano)
	PulayMixing                   = 0
	Resonances                    = 0
	ResonancesAllK                = 0
	PMLResonances                 = 0
	IncludeImaginaryContactStates = 0
	IncoherentContacts            = 1
	IncoherentContactBroadening   = 1e-2
	FreyModel                     = 0
	
	# ---------------
	# Self Energies
	# ---------------
	ScatteringDecreaseFactor    = 0.01	# S(0) = (this)*S(end)
#	ScatteringDecreaseFactor    = 0.01	# S(0) = (this)*S(end)
	ScatteringRampFactor        = 2.0	# S(i) * (this)*S(i-1)
#	SelfEnergyUnderrelaxation   = 0.3	# S_new = (factor)*S_old + (1-factor)*S_new
	
	Buettiker                   = 0
	BuettikerParameter          = 0.01  # [eV]
	
	GolizadehMomentumRelaxation = 1
	GolizadehMomentumParameter  = 81e-6 # [eV^2]
	GoliQWHack                  = 0
	
	GolizadehDephasing          = 0
	GolizadehDephasingParameter = 9e-6 # [eV^2]
	
	OpticalPhonons              = 0
	LuisierSRpop                = 0
	#PhononCheatFactor           = 1
	
	AcousticPhonons             = 1
	
	IonizedImpurities           = 1
	#IonizedImpurityCheatFactor  = 0.1
	
	SpontaneousPhotons          = 1
	LuisierSRphot               = 0
	LuminescenceRefinementEnergy = 2.85	# hw at which emission is roughly biggest
	LuminescenceRefinementWidth = 0.2
	DeltaEgapHWmin              = 0.05
	DeltaEgapHWmax              = 0.45
	
	inner_errcrit               = 1e-6
	
# be sure to compile the program with as many offdiagonals as there are vertices in the QW!!
}


# ramp in 0.1V-steps until 1.5V, then in 0.025V-steps
# starting at 1.5V works, but not at 1.6V (too much MPI data to be sent)
experiment_0 {
	Source_voltage = 0.0
	Drain_min      = 2.600
#	Drain_min      = 2.800
	Drain_max      = 3.999
	Drain_step     = 0.05
	PotentialUnderrelaxation = 0.6
	LateUnderrelaxation      = 0.2
}
experiment_1 {
	Source_voltage = 0.0
	Drain_min      = 4.000
	Drain_max      = 4.801
	Drain_step     = 0.025
	PotentialUnderrelaxation = 0.4
	LateUnderrelaxation      = 0.2
}


# Nitride LED: GaN - InGaN(x=0.2) - AlGaN(x=0.2)
regions {
	
	region0_length = 14.0
	region0_dx     = 0.4
	region0_mat    = 6      # GaN
	region0_molefr = 0.0
	region0_doping = 1e19   # n
	
	region1_length = 0.2
	region1_dx     = 0.05
	region1_mat    = 6      # GaN side of interface
	region1_molefr = 0.0
	region1_doping = 1e10   # i
	
	region2_length = 0.2
	region2_dx     = 0.05
	region2_mat    = 17     # InGaN side of interface
	region2_molefr = 0.2
	region2_doping = 1e10   # i
	
	region3_length = 3.0
	region3_dx     = 0.25
	region3_mat    = 17     # InGaN
	region3_molefr = 0.2
	region3_doping = 1e10   # i
	
	region4_length = 0.2
	region4_dx     = 0.05
	region4_mat    = 17     # InGaN side of interface
	region4_molefr = 0.2
	region4_doping = 1e10   # i
	
	region5_length = 0.2
	region5_dx     = 0.05
	region5_mat    = 16     # AlGaN side of interface
	region5_molefr = 0.2
	region5_doping = 1e10   # i
	
	region6_length = 10.0
	region6_dx     = 0.25
	region6_mat    = 16     # AlGaN
	region6_molefr = 0.2
	region6_doping = -1e19  # p
}
