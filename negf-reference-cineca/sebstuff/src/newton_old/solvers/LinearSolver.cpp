/*
Copyright (c) 2009 Sebastian Steiger, Integrated Systems Laboratory, ETH Zurich.
Comments, suggestions, criticism or bug reports are welcome: steiger@purdue.edu. 

This file is part of ANGEL, a simulator for LEDs based on the NEGF formalism.
The software is distributed under the Lesser GNU General Public License (LGPL).
ANGEL is free software: you can redistribute it and/or modify it under the terms 
of the Lesser GNU General Public License v3 or later. ANGEL is distributed
without any warranty; without even the implied warranty of merchantability or 
fitness for a particular purpose. See also <http://www.gnu.org/licenses/>.
*/
#include "LinearSolver.h"
using namespace negf;

#ifdef _OPENMP
	#include <omp.h>
#endif
	
LinearSolver::LinearSolver(CSRMatrix<double>* matrix_, double* rhs_, double* solution_):
	matrix(matrix_),
	rhs(rhs_),
	solution(solution_)
{STACK_TRACE(
	#ifdef _OPENMP
		this->nthreads = omp_get_max_threads();
	#else
		this->nthreads = 1;
	#endif
	this->fidx = matrix->get_prow()[0];
);}

