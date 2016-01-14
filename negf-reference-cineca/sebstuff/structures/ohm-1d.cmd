options {

	temperature = 300
	
	kp_method = 1
	maximal_k_value = 3.0
	num_k_points = 100
	
	min_energy = 1.2
	max_energy = 1.9
	num_energy_points = 400
	
	PotentialUnderrelaxation = 0.2
	
	# ---------------
	# Self Energies
	# ---------------
	ScatteringDecreaseFactor = 0.01		# S(0) = (this)*S(end)
	ScatteringRampFactor = 2.0			# S(i) * (this)*S(i-1)
	
	Buettiker = 0
	BuettikerParameter = 0.01 # [eV]
	
	GolizadehMomentumRelaxation = 1
	GolizadehMomentumParameter = 0.002 # [eV^2]
	
	GolizadehDephasing = 0
	GolizadehDephasingParameter = 0.0001 # [eV^2]
	
	OpticalPhonons = 0
}


experiment_0 {
	Drain_voltage   = 0.0
	Source_min      = 0.0
	Source_max      = 0.301
	Source_step     = 0.1
}
