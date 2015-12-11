options {

	temperature                   = 300

	kp_method                     = 0
	maximal_k_value               = 2.5 # in 1/nm
	num_k_points                  = 40
	
	num_energy_points             = 300
	
	min_energy                    = 0.9	# for n-doping
	max_energy                    = 1.7
#	min_energy                    = -0.5	# for p-doping
#	max_energy                    = 0.2
# remember to also change the doping file when switching between n- and p-doping!
# and change kpmethod=2!	
		
	StrainPolarization 			  = 0

	PotentialUnderrelaxation      = 0.1
	PulayMixing                   = 0
	Resonances                    = 0
	PMLResonances                 = 0
	IncludeImaginaryContactStates = 0
	IncoherentContacts            = 1
	IncoherentContactBroadening   = 1e-2
	FreyModel                     = 2
	
	Transmission                 = 1
	
	# ---------------
	# Self Energies
	# ---------------
	ScatteringDecreaseFactor    = 0.01	# S(0) = (this)*S(end)
	ScatteringRampFactor        = 2.0	# S(i) * (this)*S(i-1)
	
	Buettiker                   = 0
	BuettikerParameter          = 0.01  # [eV]
	
	GolizadehMomentumRelaxation = 0
	GolizadehMomentumParameter  = 0.002 # [eV^2]
	
	# Underrelaxation is mandatory for dephasing
	#SelfEnergyUnderrelaxation   = 0.3
	GolizadehDephasing          = 0
	GolizadehDephasingParameter = 0.0004 # [eV^2] sqrt(0.0004)=20meV 
	# ramping doesnt really pay off w/ Golizadeh dephasing
	# self-energy is completely differeent after each Poisson step
	# because it solely smoothes out the DoS
	
	OpticalPhonons              = 0
	LuisierSRpop                = 0
	
	AcousticPhonons             = 1
}


experiment_0 {
	Drain_voltage   = 0.0
	Source_min      = 0.0
	Source_max      = 0.251
	Source_step     = 0.05
#	Source_max      = -0.251
#	Source_step     = -0.05
	PotentialUnderrelaxation = 0.1
	LateUnderrelaxation      = 0.0
}
