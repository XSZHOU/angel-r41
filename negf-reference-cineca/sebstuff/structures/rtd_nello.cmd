options {

	temperature                   = 300
	
	kp_method                     = 1
	
	maximal_k_value               = 1.5 # in 1/nm; 
	# note: start from hbar2*(1e9)^2/(2*m0) = 0.076eV
	# BZ is pi/5A = 6.3nm-1 wide
	num_k_points                  = 60
	
	min_energy                    = 0.9
	max_energy                    = 1.72
	num_energy_points             = 400
	
	PotentialUnderrelaxation      = 0.6
	PulayMixing                   = 0
	NonzeroBoundaryField          = 0 # (0/1) nonzero electric field for Neumann bnd.cond. for Poisson (like nextnano)
	Resonances                    = 0
	ResonancesAllK                = 1
	PMLResonances                 = 0
	IncludeImaginaryContactStates = 0
	IncoherentContacts            = 0
	IncoherentContactBroadening   = 1e-2
	
	# ---------------
	# Self Energies
	# ---------------
	ScatteringDecreaseFactor    = 0.1		# S(0) = (this)*S(end)
	ScatteringRampFactor        = 2.0			# S(i) * (this)*S(i-1)
	
	Buettiker                   = 0
	BuettikerParameter          = 0.01 # [eV]
	
	GolizadehMomentumRelaxation = 0
	GolizadehMomentumParamete r = 0.002 # [eV^2]
	
	GolizadehDephasing          = 0
	GolizadehDephasingParameter = 0.000001 # [eV^2]
	
	OpticalPhonons              = 0
	LuisierSRpop                = 0
	AcousticPhonons             = 1
	#PhononCheatFactor           = 0.001
	
	SpontaneousPhotons          = 0
}


experiment_0 {
	Drain_voltage   = 0.0
	Source_min      = 0.000
	Source_max      = 0.199
	Source_step     = 0.1
}
experiment_1 {
	Drain_voltage   = 0.0
	Source_min      = 0.200
	Source_max      = 0.299
	Source_step     = 0.02
}
experiment_2 {
	Drain_voltage   = 0.0
	Source_min      = 0.300
	Source_max      = 0.399
	Source_step     = 0.01
}
experiment_3 {
	Drain_voltage   = 0.0
	Source_min      = 0.400
	Source_max      = 0.601
	Source_step     = 0.02
}
