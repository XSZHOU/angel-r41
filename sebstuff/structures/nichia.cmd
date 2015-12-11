options {

	temperature                   = 300
	
	kp_method                     = 3
	
	maximal_k_value               = 1.5 # in 1/nm
	# note: start from hbar2*(1e9)^2/(2*m0) = 0.076eV and multiply w/ kmax^2/m*
	# BZ is pi/5A = 6.3nm-1 wide
	num_k_points                  = 60
	
	min_energy                    = -0.5
	max_energy                    = 4.5
	num_energy_points             = 600
	
	StrainPolarization 			  = 1
	PolarizationDecreaser         = 0.5

	PotentialUnderrelaxation      = 0.6
	NonzeroBoundaryField          = 0 # (0/1) nonzero electric field for Neumann bnd.cond. for Poisson (like nextnano)
	PulayMixing                   = 0
	Resonances                    = 0
	ResonancesAllK                = 0
	PMLResonances                 = 0
	IncludeImaginaryContactStates = 0
	IncoherentContacts            = 1
	IncoherentContactBroadening   = 5e-3
	FreyModel                     = 0
	InjectingStatesCutoff         = 0.3  # above/below Ec-cutoff / Ev+cutoff no injecting CB/VB state
	
	# ---------------
	# Self Energies
	# ---------------
	ScatteringDecreaseFactor    = 0.01	# S(0) = (this)*S(end)
	ScatteringRampFactor        = 2.0	# S(i) * (this)*S(i-1)
	SelfEnergyUnderrelaxation   = 0.0	# S_new = (factor)*S_old + (1-factor)*S_new
	
	Buettiker                   = 0
	BuettikerParameter          = 0.01  # [eV]
	
	GolizadehMomentumRelaxation = 1
	GolizadehMomentumParameter  = 1e-6 # [eV^2]
	
	GolizadehDephasing          = 0
	GolizadehDephasingParameter = 0.0001 # [eV^2]
	
	OpticalPhonons              = 1
	LuisierSRpop                = 0
	PhononCheatFactor           = 0.05
	
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
experiment_0 {
	Drain_voltage   = 0.0
	Source_min      = 2.700
	Source_max      = 2.701
	Source_step     = 0.05
	PotentialUnderrelaxation = 0.6
	LateUnderrelaxation      = 0.2
}
