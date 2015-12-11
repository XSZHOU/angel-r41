#ifndef LINEARSOLVERILS_H_
#define LINEARSOLVERILS_H_

#include "all.h"

#include "CSRMatrix.h"
#include "LinearSolver.h"
#include "src/ils.h"	// $(INC_PATH)/...  must contain ILS path!

using namespace std;

namespace negf {

	/** Implements the parallel iterative solver ILS-2.0 of S. Roellin */
	class LinearSolverILS: public LinearSolver {

	public:
		
		LinearSolverILS(CSRMatrix<double>* matrix_, double * rhs_, double * solution_, const char * ilsr_filename);
		~LinearSolverILS();

		// set the maximum number of Krylov iteration (if not specified or set to zero, the .ilsrc value is used)
		void set_max_num_iterations(uint max_iterations_) { this->max_iterations = max_iterations_; }
		
		void solve();
	
	protected:
	
		void throw_error(int error_code);
	
		// some options, set in constructor
		const bool scale_lines_individually;
		const uint recompute_nonsym_perm;   /*  0: do not recompute nonsymmetric permutation and scalings, even if values have changed
 		  									    1: recompute nonsymmetric permutation (choose this one for normal usage!) */
		const uint recompute_preconditioner; /* 0: reuse preconditioner
                                    		    1: reuse pattern of preconditioner
                                    		    2: recompute preconditioner */
		int ils_handle;
		bool solved_at_least_once;
		
		int max_iterations;
		
		int * ils_prow;
		int * ils_icol;
	};

} // end namespace

#endif /*LINEARSOLVERILS_H_*/
