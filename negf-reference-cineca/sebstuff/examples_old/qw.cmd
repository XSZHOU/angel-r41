options {

	temperature                   = 300
	
	kp_method                     = 2
	
	maximal_k_value               = 2.5 # in 1/nm
	# note: start from hbar2*(1e9)^2/(2*m0) = 0.076eV and multiply w/ kmax^2/m*
	# BZ is pi/5A = 6.3nm-1 wide
	num_k_points                  = 100
	
	min_energy                    = -1.0
	max_energy                    = 1.6
	num_energy_points             = 600
	
	StrainPolarization 			  = 1
	
	PotentialUnderrelaxation      = 0.6
	NonzeroBoundaryField          = 0 # (0/1) nonzero electric field for Neumann bnd.cond. for Poisson (like nextnano)
	PulayMixing                   = 0
	Resonances                    = 0
	ResonancesAllK                = 1
	PMLResonances                 = 0
	IncludeImaginaryContactStates = 0
	IncoherentContacts            = 1
	IncoherentContactBroadening   = 1e-2
	
	# ---------------
	# Self Energies
	# ---------------
	ScatteringDecreaseFactor    = 0.05	# S(0) = (this)*S(end)
	ScatteringRampFactor        = 2.0	# S(i) * (this)*S(i-1)
	SelfEnergyUnderrelaxation   = 0.0	# S_new = (factor)*S_old + (1-factor)*S_new
	
	Buettiker                   = 0
	BuettikerParameter          = 0.01  # [eV]
	
	GolizadehMomentumRelaxation = 1
	GolizadehMomentumParameter  = 0.001 # [eV^2]
	
	GolizadehDephasing          = 0
	GolizadehDephasingParameter = 0.001 # [eV^2]
	
	OpticalPhonons              = 0
	LuisierSRpop                = 1
	#PhononCheatFactor          = 0.01
	
	AcousticPhonons             = 1
	
	SpontaneousPhotons          = 1
	LuisierSRphot               = 0
	LuminescenceRefinementEnergy = 0.8	# hw at which emission is roughly biggest
	LuminescenceRefinementWidth = 0.2
	DeltaEgapHWmin              = 0.0	# quantization energy leads to emission strictly above Egap(InAs)
	DeltaEgapHWmax              = 0.8
}


#experiment_0 {
#	Drain_voltage   = 0.0
#	Source_min      = -0.5
#	Source_max      = -0.699
#	Source_step     = -0.1
#}
experiment_1 {
	Drain_voltage   = 0.0
#	Source_min      = -0.700
#	Source_min      = -0.850
	Source_min      = -0.950
	Source_max      = -1.501
	Source_step     = -0.05
	PotentialUnderrelaxation = 0.6
	LateUnderrelaxation      = 0.2
}
