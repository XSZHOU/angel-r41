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
	Drain_voltage   = 0.0
	Source_min      = 2.600
#	Source_min      = 2.800
	Source_max      = 3.999
	Source_step     = 0.05
	PotentialUnderrelaxation = 0.6
	LateUnderrelaxation      = 0.2
}
experiment_1 {
	Drain_voltage   = 0.0
	Source_min      = 4.000
	Source_max      = 4.801
	Source_step     = 0.025
	PotentialUnderrelaxation = 0.4
	LateUnderrelaxation      = 0.2
}
