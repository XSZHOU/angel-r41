options {

	temperature                   = 300
	
	kp_method                     = 3
	
	maximal_k_value               = 1.5 # in 1/nm
	# note: start from hbar2*(1e9)^2/(2*m0) = 0.076eV and multiply w/ kmax^2/m*
	# BZ is pi/5A = 6.3nm-1 wide
	num_k_points                  = 80
	
	min_energy                    = -0.891584
	max_energy                    = 2.0
	num_energy_points             = 400
	
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
	Drain_voltage   = 0.0
	Source_min      = -1.200
	Source_max      = -1.399
	Source_step     = -0.05
	PotentialUnderrelaxation = 0.6
	LateUnderrelaxation      = 0.3
}
experiment_1 {
	Drain_voltage   = 0.0
	Source_min      = -1.450
	Source_max      = -1.499
	Source_step     = -0.025
	PotentialUnderrelaxation = 0.5
	LateUnderrelaxation      = 0.3
}
experiment_2 {
	Drain_voltage   = 0.0
	Source_min      = -1.550
#	Source_min      = -1.625
	Source_max      = -1.801
	Source_step     = -0.025
	PotentialUnderrelaxation = 0.5
	LateUnderrelaxation      = 0.2
}
