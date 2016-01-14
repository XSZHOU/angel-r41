#ifndef MATH_DEFINITION_H_
#define MATH_DEFINITION_H_

// ------------------------------------------
// standard includes
// ------------------------------------------
#include "all.h"
#include "cxxtest/TestSuite.h"

// ------------------------------------------
// class includes
// ------------------------------------------
// flens.h is included in all.h
typedef flens::GeMatrix<flens::FullStorage<cplx, flens::ColMajor> > 	CGEMatrix; 
typedef flens::DenseVector<flens::Array<cplx> > 						CDEVector; 
typedef flens::GeMatrix<flens::FullStorage<double, flens::ColMajor> > 	DGEMatrix; 
typedef flens::DenseVector<flens::Array<double> > 						DDEVector; 

// LAPACK
typedef std::complex<double> doublecomplex;
extern "C" { 
void zgeev(
	char jobvl, // 'N': left eigenvectors of A are not computed;  'V': left eigenvectors are computed.
	char jobvr, // 'N': right eigenvectors of A are not computed; 'V': right eigenvectors are computed.
	int n, 		// order of the matrix A.
	doublecomplex *a, 	// (in/out) array of dimension (LDA,N).  On entry, the matrix A. On exit, A has been overwritten.
	int lda, 	// The leading dimension of the array A.  LDA >= max(1,N).
	doublecomplex * w, 	// (out) will contain the n computed eigenvalues
	doublecomplex *vl, 	// (out) dimension (LDVL,N). If JOBVL = 'N', VL is not referenced.
				// If JOBVL = 'V', the left eigenvectors u(j) are stored one after another in the columns of VL, in the same order
				// as their eigenvalues: u(j) = VL(:,j)
	int ldvl,   // The leading dimension of the array VL.  LDVL >= 1; if JOBVL = 'V', LDVL >= N.
	doublecomplex *vr, 	// same as vl, but with right eigenvalues
	int ldvr, 	// The leading dimension of the array VR.  LDVR >= 1; if JOBVR = 'V', LDVR >= N.
	int *info); // pointer to integer; = 0: success. if INFO = -i, the i-th argument had an illegal value.
				// if INFO = i, the QR algorithm failed to compute all the  eigenvalues, and no eigenvectors have been computed;
				// elements and i+1:N of W contain eigenvalues which have converged.
}

using namespace negf;
using namespace flens;

namespace negf {

class MathTest : public CxxTest::TestSuite 
{
public:
	
	int N;
	int M;
	double eps;
	std::fstream fout;	
	void setUp() {
		N = 20;
		M = 30;
		eps = 1e-13;	
	}
	void tearDown() {
	}
		
	/** Test the following operations:
	 *  zgeev - eigenvalues and eigenvectors of a complex unsymmetric matrix
	 *  log - complex logarithm (tested against Mathematica results)
	 *  acos - complex arc cosine (tested against Mathematica results)
	 */
	void test_operation()
	{
		fout.open("./negf/Math_output.log", std::ios::out); 
		// current directory is Makefile-directory
		logmsg->add_listener(fout);
		
		const double eps = 1e-13;
		
		// NOTE: all numerical results were computed using built-in Mathematica functions
		
		// ---------------------
		// test std::log
		// ---------------------
		
		cplx a =  1.0 + constants::imag_unit * 1.0;
		cplx b =  1.0 - constants::imag_unit * 1.0;
		cplx c = -1.0 + constants::imag_unit * 1.0;
		cplx d = -1.0 - constants::imag_unit * 1.0;
		cplx e = -2.0 + constants::imag_unit * 1e-14;
		cplx f = -2.0 - constants::imag_unit * 1e-14;
		
		cplx log_a = std::log(a);
		TS_ASSERT(fabs(log_a.real() - 0.34657359027997) < eps);
		TS_ASSERT(fabs(log_a.imag() - 0.78539816339745) < eps);
		
		cplx log_b = std::log(b);
		TS_ASSERT(fabs(log_b.real() - 0.34657359027997) < eps);
		TS_ASSERT(fabs(log_b.imag() + 0.78539816339745) < eps);
		
		cplx log_c = std::log(c);
		TS_ASSERT(fabs(log_c.real() - 0.3465735902800) < eps);
		TS_ASSERT(fabs(log_c.imag() - 2.3561944901923) < eps);
		      
		cplx log_d = std::log(d);
		TS_ASSERT(fabs(log_d.real() - 0.3465735902800) < eps);
		TS_ASSERT(fabs(log_d.imag() + 2.3561944901923) < eps);
		
		cplx log_e = std::log(e);
		TS_ASSERT(fabs(log_e.real() - 0.6931471805599) < eps);
		TS_ASSERT(fabs(log_e.imag() - 3.1415926535898) < eps);
		
		cplx log_f = std::log(f);
		TS_ASSERT(fabs(log_f.real() - 0.6931471805599) < eps);
		TS_ASSERT(fabs(log_f.imag() + 3.1415926535898) < eps);
		
		// ---------------------
		// test negf_math::acos
		// ---------------------
		
		cplx acos_a = negf_math::acos(a);
		TS_ASSERT(fabs(acos_a.real() - 0.90455689430238) < eps);
		TS_ASSERT(fabs(acos_a.imag() + 1.06127506190504) < eps);
		
		cplx acos_b = negf_math::acos(b);
		TS_ASSERT(fabs(acos_b.real() - 0.90455689430238) < eps);
		TS_ASSERT(fabs(acos_b.imag() - 1.06127506190504) < eps);
		
		cplx acos_c = negf_math::acos(c);
		TS_ASSERT(fabs(acos_c.real() - 2.2370357592874) < eps);
		TS_ASSERT(fabs(acos_c.imag() + 1.0612750619050) < eps);
		
		cplx acos_d = negf_math::acos(d);
		TS_ASSERT(fabs(acos_d.real() - 2.2370357592874) < eps);
		TS_ASSERT(fabs(acos_d.imag() - 1.0612750619050) < eps);
		
		cplx acos_e = negf_math::acos(e);
		TS_ASSERT(fabs(acos_e.real() - 3.1415926535898) < eps);
		TS_ASSERT(fabs(acos_e.imag() + 1.3169578969248) < eps);
		
		cplx acos_f = negf_math::acos(f);
		TS_ASSERT(fabs(acos_f.real() - 3.1415926535898) < eps);
		TS_ASSERT(fabs(acos_f.imag() - 1.3169578969248) < eps);
		
		// ---------------------------------------
		// test zgeev with a diagonal matrix
		// ---------------------------------------
		CGEMatrix A(6,6);
		A(1,1) = a;
		A(2,2) = b;
		A(3,3) = c;
		A(4,4) = d;
		A(5,5) = e;
		A(6,6) = f;
		int n = A.numRows();
		
		// convert FLENS matrix to array
		cplx AA[n*n];				// SOOOOOOO INEFFICIENT
		for (int ii=1; ii <= n; ii++) {
			for (int jj=1; jj <= n; jj++) {
				AA[(jj-1)*n+(ii-1)] = A(ii,jj);
			}
		}
		
		// solve w/ LAPACK
		int info = -88;
		cplx lambda[n];	// eigenvalues
		cplx vr[n*n];	// eigenvectors
		zgeev('N', 'V', n, AA, n, lambda, NULL, 1, vr, n, &info); 
		TS_ASSERT(info==0);
		TS_ASSERT(abs(lambda[0] - a) < eps);
		TS_ASSERT(abs(lambda[1] - b) < eps);
		TS_ASSERT(abs(lambda[2] - c) < eps);
		TS_ASSERT(abs(lambda[3] - d) < eps);
		TS_ASSERT(abs(lambda[4] - e) < eps);
		TS_ASSERT(abs(lambda[5] - f) < eps);
		
		vector< vector<cplx> > psi;
		psi.resize(n);
		for (uint ii=0; ii<n; ii++) {
			psi[ii].resize(n);			// psi[ii] stores eigenvector ii
			for (uint jj=0; jj<n; jj++) {
				psi[ii][jj] = vr[ii*n+jj];	// vr(jj,ii)
				if (ii!=jj) {
					TS_ASSERT(abs(psi[ii][jj]) < eps);
				} else {
					TS_ASSERT(fabs(abs(psi[ii][jj]) - 1.0) < eps);
				}
			}
		}
		
		fout.close();
	}
	
};

} // end namespace negf

#endif /*MATH_DEFINITION_H_*/
