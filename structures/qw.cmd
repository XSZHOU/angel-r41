options {

	temperature                   = 300
	
	kp_method                     = 3
	
	Transmission                  = 0
	QuasiFermilevels              = 1
	
	# ---------------
	# Numerics
	# ---------------
	
	maximal_k_value               = 1.7 # in 1/nm
	# note: start from hbar2*(1e9)^2/(2*m0) = 0.076eV and multiply w/ kmax^2/m*
	# BZ is pi/5A = 6.3nm-1 wide
	num_k_points                  = 80
	
	min_energy                    = -0.5
	max_energy                    = 1.65
	num_energy_points             = 301
	
	StrainPolarization 			  = 1
	
	PotentialUnderrelaxation      = 0.6
	NonzeroBoundaryField          = 0 # (0/1) nonzero electric field for Neumann bnd.cond. for Poisson (like nextnano)
	PulayMixing                   = 0
	Resonances                    = 0
	ResonancesAllK                = 1
	PMLResonances                 = 0
	IncludeImaginaryContactStates = 0
	InjectingStatesCutoff         = 0.050 # set contact self-energy to zero below min(Ec)-cutoff (CB), above max(Ev)+cutoff (VB)
	IncoherentContacts            = 1
	IncoherentContactBroadening   = 1e-2
	
	inner_errcrit                 = 1e-5
	outer_convergence_crit        = 1e-4
	max_outer_iterations          = 20
	max_inner_iterations          = 50
	outer_loop_monitoring         = 1
	
	# ---------------
	# Self Energies
	# ---------------
	ScatteringDecreaseFactor    = 0.05	# S(0) = (this)*S(end)
	ScatteringRampFactor        = 2.0	# S(i) * (this)*S(i-1)
	SelfEnergyUnderrelaxation   = 0.0	# S_new = (factor)*S_old + (1-factor)*S_new
	
	Buettiker                   = 0
	BuettikerParameter          = 0.01  # [eV]
	
	GolizadehMomentumRelaxation = 1
	GolizadehMomentumParameter  = 225e-6 # [eV^2]
	GoliQWOnly  = 1
	GoliQWLeft  = 17
	GoliQWRight = 30 
	
	GolizadehDephasing          = 0
	GolizadehDephasingParameter = 0.001 # [eV^2]
	
	OpticalPhonons              = 1
	LuisierSRpop                = 0
	#PhononCheatFactor          = 0.01
	
	AcousticPhonons             = 1
	
	IonizedImpurities           = 1
	#IonizedImpurityCheatFactor  = 0.1
	
	SpontaneousPhotons          = 1
	LuisierSRphot               = 0
	LuminescenceRefinementEnergy = 0.8	# hw at which emission is roughly biggest
	LuminescenceRefinementWidth = 0.2
	DeltaEgapHWmin              = 0.0	# quantization energy leads to emission strictly above Egap(InAs)
	DeltaEgapHWmax              = 0.8
}


# -----------------------
# Voltage ramping
# -----------------------

#experiment_0 {
#	Source_voltage = 0.0
#	Drain_min      = 0.5
#	Drain_max      = 0.699
#	Drain_step     = 0.1
#}
experiment_1 {
	Source_voltage = 0.0
#	Drain_min      = 0.850
	Drain_min      = 1.000
	Drain_max      = 1.501
	Drain_step     = 0.05
	PotentialUnderrelaxation = 0.6
	LateUnderrelaxation      = 0.2
}


# -----------------------
# Structure definition
# -----------------------

# GaAs-InGaAs
regions {
	
	region0_length = 10.0
	region0_dx     = 0.4
	region0_mat    = 0      # GaAs
	region0_molefr = 0.0
	region0_doping = -2e18  # p
	
	region1_length = 4.0
	region1_dx     = 0.4
	region1_mat    = 0      # GaAs
	region1_molefr = 0.0
	region1_doping = 1e10   # i
	
	region2_length = 0.2
	region2_dx     = 0.1
	region2_mat    = 0      # GaAs
	region2_molefr = 0.0
	region2_doping = 1e10   # i
	
	region3_length = 0.2
	region3_dx     = 0.1
	region3_mat    = 11     # InGaAs
	region3_molefr = 0.8    # Ga molefraction!
	region3_doping = 1e10   # i
	
	region4_length = 5.0
	region4_dx     = 0.5
	region4_mat    = 11     # InGaAs
	region4_molefr = 0.8    # Ga molefraction!
	region4_doping = 1e10   # i
	
	region5_length = 0.2
	region5_dx     = 0.1
	region5_mat    = 11     # InGaAs
	region5_molefr = 0.8    # Ga molefraction!
	region5_doping = 1e10   # i
	
	region6_length = 0.2
	region6_dx     = 0.1
	region6_mat    = 0      # GaAs
	region6_molefr = 0.0
	region6_doping = 1e10   # i
	
	region7_length = 4.0
	region7_dx     = 0.5
	region7_mat    = 0      # GaAs
	region7_molefr = 0.0
	region7_doping = 1e10   # i
	
	region8_length = 10.0
	region8_dx     = 1.0
	region8_mat    = 0      # GaAs
	region8_molefr = 0.0
	region8_doping = 2e18   # n
}
