options {

	temperature                   = 300
	
	kp_method                     = 2
	
	maximal_k_value               = 1.5 # in 1/nm
	# note: start from hbar2*(1e9)^2/(2*m0) = 0.076eV and multiply w/ kmax^2/m*
	# BZ is pi/5A = 6.3nm-1 wide
	num_k_points                  = 40
	
	min_energy                    = -1.0
	max_energy                    = 1.8
	num_energy_points             = 400
	
	PotentialUnderrelaxation      = 0.6
	NonzeroBoundaryField          = 0 # (0/1) nonzero electric field for Neumann bnd.cond. for Poisson (like nextnano)
	PulayMixing                   = 0
	Resonances                    = 0
	PMLResonances                 = 0
	IncludeImaginaryContactStates = 0
	IncoherentContacts            = 1
	IncoherentContactBroadening   = 1e-2
	
	# ---------------
	# Self Energies
	# ---------------
	ScatteringDecreaseFactor    = 0.1	# S(0) = (this)*S(end)
	ScatteringRampFactor        = 2.0	# S(i) * (this)*S(i-1)
	SelfEnergyUnderrelaxation   = 0.0	# S_new = (factor)*S_old + (1-factor)*S_new
	
	Buettiker                   = 0
	BuettikerParameter          = 0.01 # [eV]
	
	GolizadehMomentumRelaxation = 0
	GolizadehMomentumParameter  = 0.002 # [eV^2]
	
	GolizadehDephasing          = 0
	GolizadehDephasingParameter = 0.001 # [eV^2]
	
	OpticalPhonons              = 1
	LuisierSRpop                = 0
	#PhononCheatFactor          = 0.01
	
	AcousticPhonons             = 0
		
	SpontaneousPhotons          = 0
	LuisierSRphot               = 0
#	LuminescenceRefinementEnergy = 1.35
#	LuminescenceRefinementWidth = 0.1
}

#experiment_0 {
#	Drain_voltage   = 0.0
#	Source_min      = 0.0
#	Source_max      = -0.701
#	Source_step     = -0.700
#}
experiment_1 {
	Drain_voltage   = 0.0
	Source_min      = -0.8
	Source_max      = -1.199
	Source_step     = -0.1
	PotentialUnderrelaxation = 0.3
	LateUnderrelaxation = 0.0
}
experiment_2 {
	Drain_voltage   = 0.0
	Source_min      = -1.2
	Source_max      = -1.701
	Source_step     = -0.05
	PotentialUnderrelaxation = 0.3
	LateUnderrelaxation = 0.0
}
