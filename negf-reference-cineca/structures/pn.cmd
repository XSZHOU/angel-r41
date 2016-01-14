# --------------------------------------------------
# Simple pin-diode (left: p, middle: i, right: n)
#
# Working configurations:
# ballistic, Nk=40, NE=375, Emin=-0.5, Emax=2.1, V=0.8V; mh=0.2, dx=1.0
# Acoustic phonons, otherwise same as ball.
# POP, otherwise same as ball.
# --------------------------------------------------

options {

	temperature                   = 300
	
	kp_method                     = 3
	
	Transmission                  = 0
	QuasiFermilevels              = 0
	
	# ---------------
	# Numerics
	# ---------------
	
	maximal_k_value               = 1.5 # in 1/nm
	# note: start from hbar2*(1e9)^2/(2*m0) = 0.076eV and multiply w/ kmax^2/m*
	# BZ is pi/5A = 6.3nm-1 wide
	num_k_points                  = 40
	
	# equilibrium settings - both Fermilevels near E=0
	#min_energy                    = -1.0
	#max_energy                    = 1.0
	
	# V_Source=-0.8V settings - EF_left~0, EF_right~+0.8V
	#min_energy                    = -0.5
	#max_energy                    = 1.5
	
	# V_Source=-1.2V settings
	#min_energy                    = -0.5
	#max_energy                    = 1.7
	
	# V_Source=-1.4V settings
	min_energy                    = -0.5
	max_energy                    = 2.0
	
	num_energy_points             = 501
	
	PotentialUnderrelaxation      = 0.4
	NonzeroBoundaryField          = 0 # (0/1) nonzero electric field for Neumann bnd.cond. for Poisson (like nextnano)
	PulayMixing                   = 0
	Resonances                    = 0
	PMLResonances                 = 0
	IncludeImaginaryContactStates = 0
	IncoherentContacts            = 1
	IncoherentContactBroadening   = 5e-3
	FreyModel                     = 0
	
	inner_errcrit                 = 1e-5
	outer_convergence_crit        = 1e-4
	max_outer_iterations          = 20
	max_inner_iterations          = 20
	outer_loop_monitoring         = 1

	# ---------------
	# Self Energies
	# ---------------
	
	ScatteringDecreaseFactor    = 0.01	# S(0) = (this)*S(end)
	ScatteringRampFactor        = 4.0	# S(i) * (this)*S(i-1)
	SelfEnergyUnderrelaxation   = 0.0	# S_new = (factor)*S_old + (1-factor)*S_new
	
	Buettiker                   = 0
	BuettikerParameter          = 0.01 # [eV]
	
	GolizadehMomentumRelaxation = 0
	GolizadehMomentumParameter  = 100e-6 # [eV^2]
	
	GolizadehDephasing          = 0
	GolizadehDephasingParameter = 0.0001 # [eV^2]
	
	OpticalPhonons              = 0
	LuisierSRpop                = 0
	#PhononCheatFactor          = 0.01
	
	AcousticPhonons             = 0
	
	IonizedImpurities           = 1
	#IonizedImpurityCheatFactor  = 0.1
		
	SpontaneousPhotons          = 0
	LuisierSRphot               = 0
#	LuminescenceRefinementEnergy = 1.35
#	LuminescenceRefinementWidth = 0.1	
#	DeltaEgapHWmin              = 0.05
#	DeltaEgapHWmax              = 0.45
	
}


# -----------------------
# Voltage ramping
# Contact 0 = Source = left
# Contact 1 = Drain = right
# -----------------------

#experiment_0 {
#	Source_voltage = 0.0
#	Drain_min      = 0.0
#	Drain_max      = 0.701
#	Drain_step     = 0.700
#}
#experiment_1 {
#	Source_voltage = 0.0
#	Drain_min      = 0.8       # this is roughly when current can be seen numerically
#	Drain_max      = 1.001
#	Drain_step     = 0.1
#	PotentialUnderrelaxation = 0.3
#	LateUnderrelaxation = 0.0
#}
experiment_2 {
	Source_voltage = 0.0
#	Drain_min      = 1.2
	Drain_min      = 1.3
	Drain_max      = 1.701
	Drain_step     = 0.05
	PotentialUnderrelaxation = 0.3
	LateUnderrelaxation = 0.0
}


# -----------------------
# Structure definition
# -----------------------

regions {
	
	region0_length = 15
	region0_dx     = 1.0
	region0_mat    = 0     # 0 = GaAs
	region0_molefr = 0.0
	region0_doping = -1e19 # p-doping
	
	region1_length = 32
	region1_dx     = 1.0
	region1_mat    = 0     # 0 = GaAs
	region1_molefr = 0.0
	region1_doping = 1e10  # intrinsic region
	
	region2_length = 15
	region2_dx     = 1.0
	region2_mat    = 0     # 0 = GaAs
	region2_molefr = 0.0
	region2_doping = 1e19  # n-doping
}

