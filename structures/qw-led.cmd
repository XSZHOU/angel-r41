options {

	temperature                   = 300
	
	kp_method                     = 3
	
	maximal_k_value               = 1.3 # in 1/nm, was 1.5
	# note: start from hbar2*(1e9)^2/(2*m0) = 0.076eV and multiply w/ kmax^2/m*
	# BZ is pi/5A = 6.3nm-1 wide
	num_k_points                  = 80
	
	min_energy                    = -0.891584
	max_energy                    = 2.0
	num_energy_points             = 300# was 400
	
	#hwLO_AlGaAs = 36.1448meV, :5=3.61448meV
	
	PotentialUnderrelaxation      = 0.6
	NonzeroBoundaryField          = 0 # (0/1) nonzero electric field for Neumann bnd.cond. for Poisson (like nextnano)
	PulayMixing                   = 0
	Resonances                    = 0
	ResonancesAllK                = 1
	PMLResonances                 = 0
	IncludeImaginaryContactStates = 0
	IncoherentContacts            = 1
	IncoherentContactBroadening   = 1e-3
	
	# ---------------
	# Self Energies
	# ---------------
	ScatteringDecreaseFactor    = 0.05	# S(0) = (this)*S(end)
#	ScatteringDecreaseFactor    = 0.01	# S(0) = (this)*S(end)
	ScatteringRampFactor        = 2.0	# S(i) * (this)*S(i-1)
	SelfEnergyUnderrelaxation   = 0.0	# S_new = (factor)*S_old + (1-factor)*S_new
	
	Buettiker                   = 0
	BuettikerParameter          = 0.01  # [eV]
	
	GolizadehMomentumRelaxation = 1
	GolizadehMomentumParameter  = 49e-6 # [eV^2] (7meV)^2
	
	#SelfEnergyUnderrelaxation   = 0.3
	GolizadehDephasing          = 0
	GolizadehDephasingParameter = 0.0001 # [eV^2]
	
	OpticalPhonons              = 1
	LuisierSRpop                = 0
	#PhononCheatFactor          = 0.01
	
	AcousticPhonons             = 1
	
	IonizedImpurities           = 1
	#IonizedImpurityCheatFactor  = 0.1
	
	SpontaneousPhotons          = 1
	LuisierSRphot               = 0
	LuminescenceRefinementEnergy = 1.5	# hw at which emission is roughly biggest
	LuminescenceRefinementWidth = 0.13
	DeltaEgapHWmin              = 0.1
	DeltaEgapHWmax              = 0.4
	
# be sure to compile the program with as many offdiagonals as there are vertices in the QW!!

	# numerical options
	inner_errcrit               = 1e-6    # for 1.1V need 1e-6 or even 1e-7! for >=1.2V 1e-4 is OK. for >1.5V >1e-4 to get convergence
}


# ramp in 0.1V-steps until 1.5V, then in 0.025V-steps
experiment_0 {
	Source_voltage = 0.0
	Drain_min      = 1.300
	Drain_max      = 1.399
	Drain_step     = 0.05
	PotentialUnderrelaxation = 0.6
	LateUnderrelaxation      = 0.3
}
experiment_1 {
	Source_voltage = 0.0
	Drain_min      = 1.450
	Drain_max      = 1.499
	Drain_step     = 0.025
	PotentialUnderrelaxation = 0.5
	LateUnderrelaxation      = 0.3
}
experiment_2 {
	Source_voltage = 0.0
	Drain_min      = 1.550
#	Drain_min      = 1.625
	Drain_max      = 1.801
	Drain_step     = 0.025
	PotentialUnderrelaxation = 0.5
	LateUnderrelaxation      = 0.2
}


# IWCE example
regions {
	
	region0_length = 10.0
	region0_dx     = 0.5
	region0_mat    = 10     # AlGaAs
	region0_molefr = 0.3
	region0_doping = -5e18  # p, was 1e19
	
	region1_length = 6.0    # was 8
	region1_dx     = 0.5
	region1_mat    = 10     # AlGaAs
	region1_molefr = 0.3
	region1_doping = 1e10   # i
	
	region2_length = 0.2
	region2_dx     = 0.1
	region2_mat    = 10     # AlGaAs side of interface
	region2_molefr = 0.3
	region2_doping = 1e10   # i
	
	region3_length = 0.2
	region3_dx     = 0.1
	region3_mat    = 0      # GaAs side of interface
	region3_molefr = 0.0
	region3_doping = 1e10   # i
	
	region4_length = 6.0    # was 8
	region4_dx     = 0.5
	region4_mat    = 0      # GaAs
	region4_molefr = 0.0
	region4_doping = 1e10   # i
	
	region5_length = 0.2
	region5_dx     = 0.1
	region5_mat    = 0      # GaAs side of interface
	region5_molefr = 0.0
	region5_doping = 1e10   # i
	
	region6_length = 0.2
	region6_dx     = 0.1
	region6_mat    = 10     # AlGaAs side of interface
	region6_molefr = 0.3
	region6_doping = 1e10   # i
	
	region7_length = 6.0    # was 8
	region7_dx     = 1.0    # was 0.5
	region7_mat    = 10     # AlGaAs
	region7_molefr = 0.3
	region7_doping = 1e10   # i
	
	region8_length = 10.0
	region8_dx     = 1.0    # was 0.5
	region8_mat    = 10     # AlGaAs
	region8_molefr = 0.3
	region8_doping = 5e18   # n
}
