#include "LinearSolverILS.h"
using namespace negf;

LinearSolverILS::LinearSolverILS(CSRMatrix<double>* matrix_, double * rhs_, double * solution_, const char * ilsrc_filename):
	LinearSolver(matrix_, rhs_, solution_),
	scale_lines_individually(false),
	recompute_nonsym_perm(1),
	recompute_preconditioner(2),	// 1 didnt work with qwr-kapon3
	ils_handle(-1),
	solved_at_least_once(false),
	max_iterations(0)
{STACK_TRACE(
	// set environment variable ILSRC to point to the file .ilsrc containing the settings for the linear solver
	#ifdef AMD
	setenv("ILSRC",ilsrc_filename, 1);  // not overwritten if already existing
	//if (nthreads > 1) setenv("TIME_MPS","999", 1);	// hack from S. Roellin to bypass a malfunctioning section of parallel ILS
	#endif
	#ifdef SUN
	putenv("ILSRC",ilsrc_filename, 1);
	//if (nthreads > 1) putenv("TIME_MPS","999", 1);
	#endif
	#ifdef IA64
	setenv("ILSRC",ilsrc_filename, 1); 
	//if (nthreads > 1) setenv("TIME_MPS","999", 1);
	#endif
	#ifdef PPC
	setenv("ILSRC",ilsrc_filename, 1); 
	//if (nthreads > 1) setenv("TIME_MPS","999", 1);
	#endif

	// modify prow, icol s.th. C indices are used (SHOULD not be necessary)
	const uint n = matrix->get_size();
	const uint nnz = matrix->get_num_nonzeros();
	if (this->fidx!=0) {
		int * prow = matrix->get_prow();
		int * icol = matrix->get_icol();
		this->ils_prow = new int[n+1];
		for (uint ii = 0; ii < n+1; ii++) {
			ils_prow[ii] = prow[ii] - fidx;
		}
		this->ils_icol = new int[nnz];
		for (uint ii = 0; ii < nnz; ii++) {
			ils_icol[ii] = icol[ii] - fidx;
		}
	} else {
		this->ils_prow = matrix->get_prow();
		this->ils_icol = matrix->get_icol();
	}
	

	// ------------------------------------------------------
	// initialize ILS-matrix using parameter set 1
	// ------------------------------------------------------
	if (nthreads==1) {
		logmsg->emit(LOG_INFO_L3, "      Calling ils_init in unparallelized mode...");
		this->ils_handle = ils_init( n, ils_prow, ils_icol, matrix->get_nonzeros(), ILS_EXTERN, 1 );
	} else {
		logmsg->emit(LOG_INFO_L3, "      Calling ils_init in parallelized mode...");
		this->ils_handle = ils_init( n, ils_prow, ils_icol, matrix->get_nonzeros(), ILS_EXTERN | ILS_PARALLEL, 1 );
	}
	logmsg->emit(LOG_INFO_L3, "      ils_init finished.");
	
	// the handle must be >=0, otherwise its an error
	if (this->ils_handle < 0) this->throw_error(ils_handle);
);}


LinearSolverILS::~LinearSolverILS()
{STACK_TRACE(
	if (this->fidx!=0) {
		delete [] this->ils_prow;
		delete [] this->ils_icol;
	}
	logmsg->emit(LOG_INFO_L3, "      Calling ils_free...");
	ils_free(this->ils_handle);
	logmsg->emit(LOG_INFO_L3, "      ils_free finished.");
);}


void LinearSolverILS::solve()
{STACK_TRACE(
	int iterations = 20;
	
	// debug
	//logmsg->emit(LOG_INFO_L3, "      Calling ils_test_matrix...");
	//int test_code = ils_test_matrix( this->ils_handle );
	//logmsg->emit(LOG_INFO_L3,  "      ils_test_matrix returned code %d.", test_code);
	int error_code = 0;
	
	if (scale_lines_individually)
	{
		// ils_prow, ils_icol have C indices
		double * vals = matrix->get_nonzeros();
		int      n    = matrix->get_size();
		
		// ------------------------------------------------------------------------
		// instead of Ax=b solve DAx=Db, D a diagonal matrix s.th. (DA)_ii=1.0
		// therefore D_ii = 1.0/A_ii, (Db)_i = b_i / A_ii
		// ------------------------------------------------------------------------
		// 1. create a vector with scalings for each line
		logmsg->emit(LOG_INFO_L3, "      Determining scaling...");
		double * scalings = new double[n];
		for (int ii=0; ii<n; ii++) {
			int jdiag = -1;	// will store position of diagonal element of row i in icol and vals vectors
			double maxelem = 0.0;
			for (int jj=ils_prow[ii]; jj < ils_prow[ii+1]; jj++) {
				if (fabs(vals[jj]) > maxelem) {
					jdiag = jj;
					maxelem = fabs(vals[jdiag]);
				}
			}
			NEGF_ASSERT(jdiag!=-1, "Seems like an entire matrix row was zero.");
			scalings[ii] = maxelem;
			if (scalings[ii]>1e10) scalings[ii] = 1e10;
			//logmsg->emit(LOG_INFO,"ii=%d: setting scaling to 1.0/%e", ii, scalings[ii]);
		}
		// 2. scale jacobian and RHS
		logmsg->emit(LOG_INFO_L3, "      Scaling...");
		for (int ii=0; ii<n; ii++) {
			rhs[ii] /= scalings[ii];
			for (int jj=ils_prow[ii]; jj<ils_prow[ii+1]; jj++) {
				vals[jj] /= scalings[ii]; // vals[jj] stores A[ii][icol[jj]]
			}
		}
		// save matrix to file
		//matrix->save_to_file("scaled_matrix.csr");
		
		// 4. solve
		logmsg->emit(LOG_INFO_L3, "      Solving...");
		if (!solved_at_least_once) {
			// 4a. reinitialize matrix
			ils_free(this->ils_handle);
			if (nthreads==1) {
				this->ils_handle = ils_init( n, ils_prow, ils_icol, vals, ILS_EXTERN, 1 );
			} else {
				this->ils_handle = ils_init( n, ils_prow, ils_icol, vals, ILS_EXTERN | ILS_PARALLEL, 1 );
			}
			if (this->ils_handle < 0) this->throw_error(ils_handle);
			
			// 4b. solve
			error_code = ils_solve( this->ils_handle, this->rhs, this->solution, &iterations );
		} else {
			error_code = ils_solve_with_new_values( this->ils_handle, vals, this->rhs, this->solution,
           						&iterations, this->recompute_nonsym_perm, this->recompute_preconditioner );
		}
		logmsg->emit_noendl(LOG_INFO_L1,"ils needed %d iterations...",iterations);
		
		logmsg->emit(LOG_INFO_L3, "      Scaling back...");
		// 5. scale back jacobian and RHS
		for (int ii=0; ii<n; ii++) {
			rhs[ii] *= scalings[ii];
			for (int jj=ils_prow[ii]; jj<ils_prow[ii+1]; jj++) {
				vals[jj] *= scalings[ii];
			}
		}
		
		delete [] scalings;
	} else {
		if (!solved_at_least_once) {
			logmsg->emit(LOG_INFO_L3, "      Calling ils_solve...");
			if (this->max_iterations != 0) {
				ils_set_maxit( this->ils_handle, this->max_iterations);
			}
			error_code = ils_solve( this->ils_handle, this->rhs, this->solution, &iterations );
		} else {
			logmsg->emit(LOG_INFO_L3, "      Calling ils_solve_with_new_values...");
			if (this->max_iterations != 0) {
				ils_set_maxit( this->ils_handle, this->max_iterations);
			}
			error_code = ils_solve_with_new_values( this->ils_handle, matrix->get_nonzeros(), this->rhs, this->solution,
           						&iterations, this->recompute_nonsym_perm, this->recompute_preconditioner );
			if (iterations==this->max_iterations) {
				logmsg->emit(LOG_WARN,"Warning: ils_solve_with_new_values needed the maximum # iterations."); 
				if (error_code != 0) this->throw_error(error_code);
				
				logmsg->emit_noendl(LOG_WARN,"ils_free...");
				ils_free(this->ils_handle);
				logmsg->emit_noendl(LOG_WARN," ils_init...");
				if (nthreads==1) {
					this->ils_handle = ils_init( this->matrix->get_size(), this->ils_prow, this->ils_icol, this->matrix->get_nonzeros(), ILS_EXTERN, 1 );
				} else {
					this->ils_handle = ils_init( this->matrix->get_size(), this->ils_prow, this->ils_icol, this->matrix->get_nonzeros(), ILS_EXTERN | ILS_PARALLEL, 1 );
				}
				logmsg->emit_noendl(LOG_WARN," ils_solve...");
				error_code = ils_solve( this->ils_handle, this->rhs, this->solution, &iterations );
			}
		}
		logmsg->emit_noendl(LOG_INFO_L1," needed %d iterations...",iterations);
	}
	solved_at_least_once = true;
	if (error_code != 0) this->throw_error(error_code);
);}


void LinearSolverILS::throw_error(int error_code)
{STACK_TRACE(
	string error_msg;
	switch(error_code) {
	case  -1: error_msg = "OUT_OF_MEMORY"; 					break;
	case  -2: error_msg = "FILE_OPEN"; 						break;
	case  -3: error_msg = "FILE_WRITE"; 					break;
	case  -4: error_msg = "WRONG_INPUT_DATA"; 				break;
	case  -5: error_msg = "WRONG_HANDLE"; 					break;
	case  -6: error_msg = "NO_FREE_HANDLE"; 				break;
	case  -7: error_msg = "FILE_READ"; 						break;
	case -10: error_msg = "LANCZOS_BREAKDOWN"; 				break;
	case -11: error_msg = "PIVOT_BREAKDOWN"; 				break;
	case -12: error_msg = "MINIMIZATION_FAILED"; 			break;
	case -20: error_msg = "PARDISO_ZERO_PIVOT"; 			break;
	case -21: error_msg = "PARDISO_OUT_OF_MEMORY"; 			break;
	case -30: error_msg = "WRONG_MATRIX_SIZE"; 				break;
	case -31: error_msg = "MATRIX_INCONSISTENT_IA"; 		break;
	case -32: error_msg = "MATRIX_INCONSISTENT_JA"; 		break;
	case -33: error_msg = "DIAG_DOESNT_EXIST"; 				break;
	case -34: error_msg = "INDICES_NOT_SORTED"; 			break;
	case -35: error_msg = "MATRIX_STRUCTURALLY_SINGULAR"; 	break;
	case -41: error_msg = "ZERO_PIVOT"; 					break;
	case -51: error_msg = "METHOD_NOT_AVAILABLE"; 			break;
	case -52: error_msg = "INCONSISTENT_PARAMETERS"; 		break;
	case -61: error_msg = "ARPACK"; 						break;
	case -71: error_msg = "UNKNOWN"; 						break;
	default:  error_msg = "not specified"; 					break;
	}
	NEGF_FEXCEPTION("ILS returned error code %d: %s", error_code, error_msg.c_str());
);}