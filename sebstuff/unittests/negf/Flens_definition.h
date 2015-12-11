#ifndef FLENS_DEFINITION_H_
#define FLENS_DEFINITION_H_

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


using namespace negf;
using namespace flens;

namespace negf {

class FlensTest : public CxxTest::TestSuite 
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
	
	/* MAIN OBSTACLES:
	 * 1. Multiplication with itself, i.e. something like A = B*A, will yield incorrect results
	 * 2. G = -A * F does not work correctly for complex matrices; use { G=A*F; G=-G; } instead
	 * 3. B = A + B will fail at runtime. B = B + A works.
	 */
	
	
	/** Tested operations:
	 *  +, -, +=, -=, *
	 *  transpose
	 *  inversion (trf,tri)
	 *  conjugateTranspose (=transpose)
	 */
	void test_double_operation()
	{
		fout.open("./negf/Flens_output.log", std::ios::out); 
		// current directory is Makefile-directory
		logmsg->add_listener(fout);
		
		logmsg->emit(LOG_INFO,"*********** TESTING DOUBLE OPERATION ***********");
		
		DGEMatrix A(N,M), B(N,M);
		TS_ASSERT( A.numRows()==N );
		TS_ASSERT( A.numCols()==M );
		TS_ASSERT( A.firstRow()==1 );
		TS_ASSERT( A.lastRow()==N );
		TS_ASSERT( A.firstCol()==1 );
		TS_ASSERT( A.lastCol()==M );
		
		double tmp;
		TS_ASSERT_THROWS_NOTHING( tmp = A(N,M) );
		
		// create random NxM matrices
		srand48(888);
		double test1[N][M], test2[N][M];
		for (int ii = 0; ii < N; ii++) {
			for (int jj = 0; jj < M; jj++) {
				test1[ii][jj] = drand48();	// between 0 and 1
				test2[ii][jj] = drand48();	
			}
		}
		
		// fill A and B with random numbers
		for (int ii = A.firstRow(); ii <= A.lastRow(); ii++) {
			for (int jj = A.firstCol(); jj <= A.lastCol(); jj++) {
				TS_ASSERT_THROWS_NOTHING( A(ii,jj) = test1[ii-1][jj-1] );
				TS_ASSERT_THROWS_NOTHING( B(ii,jj) = test2[ii-1][jj-1] );
			}
		}
		
		// test C=A+B
		DGEMatrix C(N,M);
		C = A + B;
		for (int ii = C.firstRow(); ii <= C.lastRow(); ii++) {
			for (int jj = C.firstCol(); jj <= C.lastCol(); jj++) {
				TS_ASSERT( fabs(C(ii,jj) - (A(ii,jj)  + B(ii,jj))) < eps );
				TS_ASSERT( fabs(C(ii,jj) - (test1[ii-1][jj-1]  + test2[ii-1][jj-1])) < eps );
			}
		}
		
		// test C=-A-B
		C = -A-B;
		for (int ii = A.firstRow(); ii <= A.lastRow(); ii++) {
			for (int jj = A.firstCol(); jj <= A.lastCol(); jj++) {
				TS_ASSERT( fabs(C(ii,jj) + (test1[ii-1][jj-1]  + test2[ii-1][jj-1])) < eps );
			}
		}
		
		// test C=A; C+=B
		C = A;
		C += B;
		for (int ii = A.firstRow(); ii <= A.lastRow(); ii++) {
			for (int jj = A.firstCol(); jj <= A.lastCol(); jj++) {
				TS_ASSERT( fabs(C(ii,jj) - (test1[ii-1][jj-1]  + test2[ii-1][jj-1])) < eps );
			}
		}
		
		// test C=-A; C-=B
		C = -A;
		C -= B;
		for (int ii = A.firstRow(); ii <= A.lastRow(); ii++) {
			for (int jj = A.firstCol(); jj <= A.lastCol(); jj++) {
				TS_ASSERT( fabs(C(ii,jj) + (test1[ii-1][jj-1]  + test2[ii-1][jj-1])) < eps );
			}
		}
		cout << "C=C+A" << endl;
		// test C=C+A(=-B);
		C = C + A;
		for (int ii = A.firstRow(); ii <= A.lastRow(); ii++) {
			for (int jj = A.firstCol(); jj <= A.lastCol(); jj++) {
				TS_ASSERT( fabs(C(ii,jj) + (test2[ii-1][jj-1])) < eps );
			}
		}
		// C=A+C woudl fail!!!
		
		// test resizing of smaller matrix
		DGEMatrix D(1,1);
		D = A;
		TS_ASSERT( D.numRows()==N );
		TS_ASSERT( D.numCols()==M );
		
		// test resizing of bigger matrix
		DGEMatrix E(N+10,M+10);
		E = A;
		TS_ASSERT( E.numRows()==N );
		TS_ASSERT( E.numCols()==M );
		
		// test matrix transposititon (transpose and conjugateTranspose)
		DGEMatrix F(N,M);
		//std::cout << "A(1,2)=" << A(1,2)  << ", A(2,1)=" << A(2,1) << endl;
		F = transpose(A);
		// F = A;	transpose(F); does not alter F!
		//std::cout << "F(1,2)=" << F(1,2) << ", F(2,1)=" << F(2,1) << endl;
		TS_ASSERT( F.numRows()==M );
		TS_ASSERT( F.numCols()==N );
		for (int ii = F.firstRow(); ii <= F.lastRow(); ii++) {
			for (int jj = F.firstCol(); jj <= F.lastCol(); jj++) {
				TS_ASSERT( fabs(F(ii,jj) - A(jj,ii)) < eps );
			}
		}
		F = conjugateTranspose(A);
		TS_ASSERT( F.numRows()==M );
		TS_ASSERT( F.numCols()==N );
		for (int ii = F.firstRow(); ii <= F.lastRow(); ii++) {
			for (int jj = F.firstCol(); jj <= F.lastCol(); jj++) {
				TS_ASSERT( fabs(F(ii,jj) - A(jj,ii)) < eps );
			}
		}	
		
		// test matrix multiplication
		DGEMatrix G;
		timer->click("d1");
		G = A * F;
		timer->click("d2");
		TS_ASSERT( G.numRows()==N );
		TS_ASSERT( G.numCols()==N );
		for (int ii = G.firstRow(); ii <= G.lastRow(); ii++) {
			for (int jj = G.firstCol(); jj <= G.lastCol(); jj++) {
				double result = 0.0;
				for (int kk = A.firstCol(); kk <= A.lastCol(); kk++) {
					result += A(ii,kk) * F(kk,jj);
				}
				TS_ASSERT( fabs(G(ii,jj) - result) < eps );	// this also implies that A and F are left unaltered.
			}
		}
		timer->click("d3");
		for (int ii = G.firstRow(); ii <= G.lastRow(); ii++) {
			for (int jj = G.firstCol(); jj <= G.lastCol(); jj++) {
				G(ii,jj) = 0.0;
				for (int kk = A.firstCol(); kk <= A.lastCol(); kk++) {
					G(ii,jj) += A(ii,kk) * F(kk,jj);
				}
			}
		}	
		timer->click("d4");	
		double t1 = timer->get_click("d2")-timer->get_click("d1");
		double t2 = timer->get_click("d4")-timer->get_click("d3");
		logmsg->emit(LOG_INFO,"Multiplication of %dx%d matrix with %dx%d matrix (double):",N,M,M,N);
		logmsg->emit(LOG_INFO,"   FLENS time used: %6.3e, Naive approach: %6.3e",t1,t2);
		logmsg->emit(LOG_INFO,"   Improvement: factor %6.3g",t2/t1);
		
		// test G = -A * F;
		G = -A * F;
		for (int ii = G.firstRow(); ii <= G.lastRow(); ii++) {
			for (int jj = G.firstCol(); jj <= G.lastCol(); jj++) {
				double result = 0.0;
				for (int kk = A.firstCol(); kk <= A.lastCol(); kk++) {
					result += -A(ii,kk) * F(kk,jj);
				}
				TS_ASSERT_DELTA( G(ii,jj), result, eps );	// this also implies that A and F are left unaltered.
			}
		}
		
		// note: double matrix multiplication is not supported
		
		DGEMatrix G2;
		G2 = G * G;
		for (int ii = G.firstRow(); ii <= G.lastRow(); ii++) {
			for (int jj = G.firstCol(); jj <= G.lastCol(); jj++) {
				double result = 0.0;
				for (int kk = G	.firstCol(); kk <= G.lastCol(); kk++) {
					result += G(ii,kk) * G(kk,jj);
				}
				TS_ASSERT_DELTA(result, G2(ii,jj), eps);
			}
		}	
		// test multiplication with itself
		// TEST WILL FAIL!
		/*
		G = G * G;
		TS_ASSERT( G.numRows()==N );
		TS_ASSERT( G.numCols()==N );
		for (int ii = G.firstRow(); ii <= G.lastRow(); ii++) {
			for (int jj = G.firstCol(); jj <= G.lastCol(); jj++) {
				TS_ASSERT( G(ii,jj) == G2(ii,jj) );
			}
		}
		*/	
		/*
		DGEMatrix G3;
		G3 = G;
		G = G3 * G;
		for (int ii = G.firstRow(); ii <= G.lastRow(); ii++) {
			for (int jj = G.firstCol(); jj <= G.lastCol(); jj++) {
				TS_ASSERT_DELTA( G(ii,jj), G2(ii,jj), 1e-12 );
			}
		}
		*/
		
		/*
		G = G3;
		G = G * G3;
		for (int ii = G.firstRow(); ii <= G.lastRow(); ii++) {
			for (int jj = G.firstCol(); jj <= G.lastCol(); jj++) {
				TS_ASSERT_DELTA( G(ii,jj), G2(ii,jj), 1e-12 );
			}
		}
		*/
		
		// --------------------------------------------------------------------------------
		// test matrix inversion of 100x100-matrix
		// the matrix is read in from a file, and the result is compared to the result 
		// obtained previously by matlab
		// --------------------------------------------------------------------------------
		const uint Ninv = 100;
		
		// read in matrix
		//cout << "\nreading matrix..." << endl;
		DGEMatrix H(Ninv,Ninv);
		std::fstream fin2("./negf/mat100.mat", std::ios::in);
		TS_ASSERT(fin2 != 0);
		for (int ii = H.firstRow(); ii <= H.lastRow(); ii++) {
			for (int jj = H.firstCol(); jj <= H.lastCol(); jj++) {
				fin2 >> H(ii,jj);
			}
		}	
		double dump = 0.0; fin2 >> dump;  // dont know why I need one more..
		TS_ASSERT(fin2.eof());
		fin2.close();
		
		//cout << "performing inversion..." << endl;
		// perform inversion:
		flens::DenseVector<flens::Array<int> > p(Ninv);		// initialization of right length is mandatory!
		timer->click("d5");	
		int error = trf(H,p);	// LU factorization --> stored in H
		TS_ASSERT(error==0);
		error = tri(H,p);	// backward substitution --> result also stored in H
		TS_ASSERT(error==0);
		timer->click("d6");
		
		// read in matlab result
		//cout << "reading matlab result..." << endl;
		DGEMatrix I(Ninv,Ninv);
		std::fstream fin3("./negf/mat100inv.mat", std::ios::in);
		TS_ASSERT(fin3 != 0);
		for (int ii = I.firstRow(); ii <= I.lastRow(); ii++) {
			for (int jj = I.firstCol(); jj <= I.lastCol(); jj++) {
				fin3 >> I(ii,jj);
			}
		}	
		fin3 >> dump;
		TS_ASSERT(fin3.eof());
		fin3.close();
		
		// compare!
		//cout << "comparing..." << endl;
		for (int ii = H.firstRow(); ii <= H.lastRow(); ii++) {
			for (int jj = H.firstCol(); jj <= H.lastCol(); jj++) {
				if (fabs(H(ii,jj) - I(ii,jj)) > eps)
					cout << "ii=" << ii << " jj=" << jj << ": H=" << H(ii,jj) <<", I=" << I(ii,jj) << endl;
				TS_ASSERT( fabs(H(ii,jj) - I(ii,jj)) < eps );
			}
		}
		
		logmsg->emit(LOG_INFO,"Time needed for inversion of a %dx%d matrix: %e",
					Ninv,Ninv,timer->get_click("d6")-timer->get_click("d5"));
		
		// test that inversion of a singular matrix throws an error
		DGEMatrix S(3,3);
		S =  1,  0,  0, 
			-2,  0,  0, 
			 4,  6,  1;
		flens::DenseVector<flens::Array<int> > p2(3);
		int error2 = trf(S,p2);	// LU factorization
		TS_ASSERT(error2!=0);
		int error3 = tri(S,p2);
		TS_ASSERT(error3!=0);
		logmsg->emit(LOG_INFO,"flens::trf threw error %d when used with a singular matrix, subsequently flens::tri threw error %d.",error2,error3);
	}
	

	/** tested operations:
	 *  +, *, +=, -=, transpose, conjugateTranspose
	 *  inv !!!!!
	 * */
	void test_cplx_operation()
	{
		logmsg->emit(LOG_INFO,"*********** TESTING COMPLEX OPERATION ***********");
		
		const cplx i(0.,1.);
		
		CGEMatrix A(N,M), B(N,M);
		TS_ASSERT( A.numRows()==N );
		TS_ASSERT( A.numCols()==M );
		TS_ASSERT( A.firstRow()==1 );
		TS_ASSERT( A.lastRow()==N );
		TS_ASSERT( A.firstCol()==1 );
		TS_ASSERT( A.lastCol()==M );
		
		cplx tmp;
		TS_ASSERT_THROWS_NOTHING( tmp = A(N,M) );
		// unfortunaltely FLENS does not throw an Exception but just assert(...)'s. This is not catchable.
		//TS_ASSERT_THROWS_ANYTHING( tmp = A(N+1,M) ); 
		//TS_ASSERT_THROWS_ANYTHING( tmp = A(N,M+1) ); 
		//TS_ASSERT_THROWS_ANYTHING( tmp = A(0,1) );
		//TS_ASSERT_THROWS_ANYTHING( tmp = A(1,0) );
		
		// create random complex NxM matrices
		cplx test1[N][M], test2[N][M];
		srand48(999);
		for (int ii = 0; ii < N; ii++) {
			for (int jj = 0; jj < M; jj++) {
				test1[ii][jj] = cplx(drand48()) + i*cplx(drand48());
				test2[ii][jj] = cplx(drand48()) + i*cplx(drand48());
			}
		}
		
		// fill A and B with random numbers
		for (int ii = A.firstRow(); ii <= A.lastRow(); ii++) {
			for (int jj = A.firstCol(); jj <= A.lastCol(); jj++) {
				TS_ASSERT_THROWS_NOTHING( A(ii,jj) = test1[ii-1][jj-1] );
				TS_ASSERT_THROWS_NOTHING( B(ii,jj) = test2[ii-1][jj-1] );
			}
		}		
	
		// test C = A+B
		CGEMatrix C;
		C = A+B;
		TS_ASSERT( C.numRows()==A.numRows() && C.numCols()==A.numCols() );
		for (int ii = A.firstRow(); ii <= A.lastRow(); ii++) {
			for (int jj = A.firstCol(); jj <= A.lastCol(); jj++) {
				TS_ASSERT_THROWS_NOTHING( abs(C(ii,jj) - (test1[ii-1][jj-1]  + test2[ii-1][jj-1])) < eps );
			}
		}
		
		// test C=A; C+=B
		C = A;
		C += B;
		for (int ii = C.firstRow(); ii <= C.lastRow(); ii++) {
			for (int jj = C.firstCol(); jj <= C.lastCol(); jj++) {
				if (abs(C(ii,jj) - (test1[ii-1][jj-1]  + test2[ii-1][jj-1])) > eps) {
					cout << "C(" << ii << "," << jj << ")=" << C(ii,jj) << ", (A+B)(" << ii << "," << jj << ")=" << test1[ii-1][jj-1]  + test2[ii-1][jj-1] << endl;
				}
				TS_ASSERT( abs(C(ii,jj) - (test1[ii-1][jj-1]  + test2[ii-1][jj-1])) < eps );
			}
		}
		
		// test C=-A; C-=B
		C = -A;
		C -= B;
		for (int ii = C.firstRow(); ii <= C.lastRow(); ii++) {
			for (int jj = C.firstCol(); jj <= C.lastCol(); jj++) {
				if (abs(C(ii,jj) + (test1[ii-1][jj-1]  + test2[ii-1][jj-1])) > eps) {
					cout << "C(" << ii << "," << jj << ")=" << C(ii,jj) << ", -(A+B)(" << ii << "," << jj << ")=" << -(test1[ii-1][jj-1]  + test2[ii-1][jj-1]) << endl;
				}
				TS_ASSERT( abs(C(ii,jj) + (test1[ii-1][jj-1]  + test2[ii-1][jj-1])) < eps );
			}
		}
		
		
		// test matrix transposititon (transpose and conjugateTranspose)
		CGEMatrix F;
		F = transpose(A);
		TS_ASSERT( F.numRows()==M );
		TS_ASSERT( F.numCols()==N );
		for (int ii = F.firstRow(); ii <= F.lastRow(); ii++) {
			for (int jj = F.firstCol(); jj <= F.lastCol(); jj++) {
				TS_ASSERT( abs(F(ii,jj) - A(jj,ii)) < eps );
			}
		}
		
		F = conjugateTranspose(A);
		TS_ASSERT( F.numRows()==M );
		TS_ASSERT( F.numCols()==N );
		for (int ii = F.firstRow(); ii <= F.lastRow(); ii++) {
			for (int jj = F.firstCol(); jj <= F.lastCol(); jj++) {
				TS_ASSERT_THROWS_NOTHING( abs(F(ii,jj) - std::conj(A(jj,ii))) < eps );
			}
		}
		
		
		
		// test matrix multiplication
		CGEMatrix G(N,N);
		timer->click("c1");
		G = A * F;		// DOES NOT WORK
		timer->click("c2");
		TS_ASSERT( G.numRows()==N );
		TS_ASSERT( G.numCols()==N );
		for (int ii = G.firstRow(); ii <= G.lastRow(); ii++) {
			for (int jj = G.firstCol(); jj <= G.lastCol(); jj++) {
				complex<double> result = 0.0;
				for (int kk = A.firstCol(); kk <= A.lastCol(); kk++) {
					result = result + A(ii,kk) * F(kk,jj);
				}
				TS_ASSERT( abs(G(ii,jj) - result) < eps );
			}
		}	
		timer->click("c3");
		for (int ii = G.firstRow(); ii <= G.lastRow(); ii++) {
			for (int jj = G.firstCol(); jj <= G.lastCol(); jj++) {
				G(ii,jj) = 0.0;
				for (int kk = A.firstCol(); kk <= A.lastCol(); kk++) {
					G(ii,jj) += A(ii,kk) * F(kk,jj);
				}
			}
		}	
		timer->click("c4");	
		double t1 = timer->get_click("c2")-timer->get_click("c1");
		double t2 = timer->get_click("c4")-timer->get_click("c3");
		logmsg->emit(LOG_INFO,"Multiplication of %dx%d matrix with %dx%d matrix (cplx):",N,M,M,N);
		logmsg->emit(LOG_INFO,"   FLENS time used: %6.3e, Naive approach: %6.3e",t1,t2);
		logmsg->emit(LOG_INFO,"   Improvement: factor %6.3g",t2/t1);
		
		// test G = -A * F; DOES NOT WORK!!!
		/*G = -A * F;
		for (int ii = G.firstRow(); ii <= G.lastRow(); ii++) {
			for (int jj = G.firstCol(); jj <= G.lastCol(); jj++) {
				cplx result = 0.0;
				for (int kk = A.firstCol(); kk <= A.lastCol(); kk++) {
					result += -A(ii,kk) * F(kk,jj);
				}
				if (abs(G(ii,jj)-result) >= eps) {
					logmsg->emit(LOG_ERROR, "G(ii,jj)=(%.3e,%.3e), result=(%.3e,%.3e)", 
						G(ii,jj).real(), G(ii,jj).imag(), result.real(), result.imag());
				}
				TS_ASSERT( abs(G(ii,jj)-result) < eps );	// this also implies that A and F are left unaltered.
			}
		}*/
		// test G = A * F; G = -G;
		G = A * F;
		G = -G;
		for (int ii = G.firstRow(); ii <= G.lastRow(); ii++) {
			for (int jj = G.firstCol(); jj <= G.lastCol(); jj++) {
				cplx result = 0.0;
				for (int kk = A.firstCol(); kk <= A.lastCol(); kk++) {
					result += -A(ii,kk) * F(kk,jj);
				}
				if (abs(G(ii,jj)-result) >= eps) {
					logmsg->emit(LOG_ERROR, "G(ii,jj)=(%.3e,%.3e), result=(%.3e,%.3e)", 
						G(ii,jj).real(), G(ii,jj).imag(), result.real(), result.imag());
				}
				TS_ASSERT( abs(G(ii,jj)-result) < eps );	// this also implies that A and F are left unaltered.
			}
		}
		
		// test multiplication with itself
		// TEST WILL FAIL
		/*
		CGEMatrix G2;
		G2 = G * G;
		G = G * G;
		TS_ASSERT( G.numRows()==N );
		TS_ASSERT( G.numCols()==N );
		for (int ii = G.firstRow(); ii <= G.lastRow(); ii++) {
			for (int jj = G.firstCol(); jj <= G.lastCol(); jj++) {
				TS_ASSERT( G(ii,jj) == G2(ii,jj) );
			}
		}
		*/	
		// test tmp = G; G2 = tmp * G;
		CGEMatrix tmpG;
		tmpG = G;
		CGEMatrix G2;
		G2 = G * tmpG;
		for (int ii = G.firstRow(); ii <= G.lastRow(); ii++) {
			for (int jj = G.firstCol(); jj <= G.lastCol(); jj++) {
				cplx result = 0.0;
				for (int kk = G.firstCol(); kk <= G.lastCol(); kk++) {
					result += G(ii,kk) * G(kk,jj);
				}
				TS_ASSERT( abs(G2(ii,jj) - result) < eps );
			}
		}
		
		// --------------------------------------------------------------------------------
		// test the speed of multiplication of two 1x1-matrices compared to complex numbers
		// --------------------------------------------------------------------------------
		uint num_times = 100000;
		CGEMatrix S1(1,1);	S1(1,1) = 3.0 + constants::imag_unit * 5.0;
		CGEMatrix S2(1,1);	S2(1,1) = 1.0 - constants::imag_unit * 2.0;
		CGEMatrix S3(1,1);
		timer->click("s1");
		for (uint ii=0; ii<num_times; ii++) {
			S3 = S1*S2;
		}
		timer->click("s2");
		cplx s1 = 3.0 + constants::imag_unit * 5.0;
		cplx s2 = 1.0 - constants::imag_unit * 2.0;
		cplx s3;
		timer->click("s3");
		for (uint ii=0; ii<num_times; ii++) {
			s3 = s1*s2;
		}
		timer->click("s4");
		logmsg->emit(LOG_INFO,"Time needed for multiplication of %d 1x1-matrices: %.3g[s] compared to complex numbers: %.3g[s]",
					num_times, timer->get_seconds_between("s2", "s1"), timer->get_seconds_between("s4", "s3"));
		
		// --------------------------------------------------------------------------------
		// test matrix inversion of a typical retarded Greens function
		// the matrix is read in from a file, and the result is compared to the result 
		// obtained previously by matlab
		// --------------------------------------------------------------------------------
		const uint nGR = 63;
		
		// read in matrix
		CGEMatrix H(nGR,nGR);
		TS_ASSERT_THROWS_NOTHING( flens_io::read_matrix("./negf/GR_inv.mat", H); );
				
		//cout << "performing inversion..." << endl;
		// perform inversion:
		flens::DenseVector<flens::Array<int> > p(nGR);		// initialization of right length is mandatory!
		timer->click("d7");	
		int error = trf(H,p);	// LU factorization --> stored in H
		TS_ASSERT(error==0);
		error = tri(H,p);	// backward substitution --> result also stored in H
		TS_ASSERT(error==0);
		timer->click("d8");
		
		// read in matlab result
		CGEMatrix I(nGR,nGR);
		flens_io::read_matrix("./negf/GR.mat", I);
		
		// compare!
		//cout << "comparing..." << endl;
		for (int ii = H.firstRow(); ii <= H.lastRow(); ii++) {
			for (int jj = H.firstCol(); jj <= H.lastCol(); jj++) {
				if (abs(H(ii,jj) - I(ii,jj)) > eps)
					cout << "ii=" << ii << " jj=" << jj << ": H=" << H(ii,jj) <<", I=" << I(ii,jj) << endl;
				TS_ASSERT( abs(H(ii,jj) - I(ii,jj)) < eps );
			}
		}
		
		logmsg->emit(LOG_INFO,"Time needed for inversion of a complex %dx%d matrix: %e",
					nGR,nGR,timer->get_click("d8")-timer->get_click("d7"));
		
		// test that inversion of a singular matrix throws an error
		CGEMatrix S(3,3);
		S =  1,  0,  0, 
			-2,  0,  0, 
			 4,  6,  1;
		flens::DenseVector<flens::Array<int> > p2(3);
		int error2 = trf(S,p2);	// LU factorization
		TS_ASSERT(error2!=0);
		int error3 = tri(S,p2);
		TS_ASSERT(error3!=0);
		logmsg->emit(LOG_INFO,"flens::trf threw error %d when used with a singular matrix, subsequently flens::tri threw error %d.",error2,error3);
	}

	/** tested operations:
	 *  - DGEMatrix::View yields correct matrix
	 *  - copy constructor of DGEMatrix::View object is correct (both to DGEMatrix and DGEMatrix::View)
	 *  - matrix multiplication of Views and ConstViews does not change the Views or ConstViews
	 */
	void test_views()
	{	
		logmsg->emit(LOG_INFO,"*********** TESTING VIEWS ***********");
		
		TS_ASSERT(N>=5);
		DGEMatrix A(N,N);
		for (int ii=1; ii<N; ii++) {
			for (int jj=1; jj<=N; jj++) {
				A(ii,jj) = ii*jj;
			}
		}
		
		DGEMatrix::View Aview = A(_(2,4),_(2,4));
		TS_ASSERT(Aview.numCols()==3 && Aview.numRows()==3);
		for (uint ii=1; ii<3; ii++) {
			for (uint jj=1; jj<3; jj++) {
				TS_ASSERT(Aview(ii,jj)==A(ii+1,jj+1));
			}
		}
		
		DGEMatrix::View B = Aview; // GEMatrix alleine geht nicht
		for (uint ii=1; ii<=3; ii++) {
			for (uint jj=1; jj<=3; jj++) {
				TS_ASSERT(B(ii,jj)==Aview(ii,jj));
			}
		}
		
		//DGEMatrix::NoView C = Aview; // geht nicht
		//DGEMatrix C = A(_(1,3),_(1,3)); // geht nicht
		//DGEMatrix::NoView C = A(_(1,3),_(1,3)); // geht nicht
		DGEMatrix C(Aview.numRows(), Aview.numCols());
		C = Aview;
		for (uint ii=1; ii<=3; ii++) {
			for (uint jj=1; jj<=3; jj++) {
				C(ii,jj) = 0.555;
			}
		}
		// check that A, Aview and B did not change
		for (uint ii=1; ii<=3; ii++) {
			for (uint jj=1; jj<=3; jj++) {
				TS_ASSERT(A(ii+1,jj+1)==(ii+1)*(jj+1));
				TS_ASSERT(Aview(ii,jj)==(ii+1)*(jj+1));
				TS_ASSERT(B(ii,jj)==(ii+1)*(jj+1));
			}
		}
		
		DGEMatrix::ConstView AAA  = A(_(3,5),_(3,5));
		Aview = AAA;
		//cout << A;
		
		// -------------------------------------------------------------------
		// test that matrix multiplication of views does not change the views
		// -------------------------------------------------------------------
		CGEMatrix X(N,N);
		CGEMatrix Y(N,N);
		for (int ii=1; ii<=N; ii++) {
			for (int jj=1; jj<=N; jj++) {
				X(ii,jj) = ii*jj/(ii+jj);
				Y(ii,jj) = ii*jj/(ii+jj*ii);
			}
		}
		const CGEMatrix::ConstView Xview = X(_(3,5),_(3,5));
		const CGEMatrix::ConstView Yview = Y(_(3,5),_(3,5));
		CGEMatrix Z(3,3);
		Z = Xview * Yview;
		for (uint ii=1; ii<=3; ii++) {
			for (uint jj=1; jj<=3; jj++) {
				uint i2 = ii+2;
				uint j2 = jj + 2;
				TS_ASSERT(Xview(ii,jj)==i2*j2/(i2+j2));
				TS_ASSERT(Yview(ii,jj)==i2*j2/(i2+j2*i2));
			}
		}
		CGEMatrix::View Xview2 = X(_(3,5),_(3,5));
		CGEMatrix::View Yview2 = Y(_(3,5),_(3,5));
		Z = Xview2 * Yview2;
		for (uint ii=1; ii<=3; ii++) {
			for (uint jj=1; jj<=3; jj++) {
				uint i2 = ii+2;
				uint j2 = jj + 2;
				TS_ASSERT(Xview2(ii,jj)==i2*j2/(i2+j2));
				TS_ASSERT(Yview2(ii,jj)==i2*j2/(i2+j2*i2));
			}
		}
		
		// -------------------------------------------------------------------
		// test the speed of creating a view
		// -------------------------------------------------------------------
		const uint num_views = 10000;
		timer->click("v1");
		for (uint ii=0; ii < num_views; ii++) {
			const CGEMatrix::ConstView Xview = X(_(3,5),_(3,5));
		}
		timer->click("v2");
		for (uint ii=0; ii < num_views; ii++) {
			CGEMatrix::View Xview = X(_(3,5),_(3,5));
		}
		timer->click("v3");
		for (uint ii=0; ii < num_views; ii++) {
			CGEMatrix::View Xview = X(_(3,3),_(3,3));
		}
		timer->click("v4");
		for (uint ii=0; ii < num_views; ii++) {
			CGEMatrix::View Xview = X(_(1,5),_(1,5));
		}
		timer->click("v5");
		for (uint ii=0; ii < num_views; ii++) {
			CGEMatrix::View Xview = this->get_view(X, 3, 5, 1);
		}
		timer->click("v6");
		for (uint ii=0; ii < num_views; ii++) {
			const CGEMatrix::ConstView Xview = this->get_const_view(X, 3, 5, 1);
		}
		timer->click("v7");
		logmsg->emit(LOG_INFO,"time needed for creating %d ConstViews of size 3: %.3g[s]",num_views,timer->get_seconds_between("v2", "v1"));
		logmsg->emit(LOG_INFO,"time needed for creating %d Views of size 3: %.3g[s]",num_views,timer->get_seconds_between("v3", "v2"));
		logmsg->emit(LOG_INFO,"time needed for creating %d Views of size 1: %.3g[s]",num_views,timer->get_seconds_between("v4", "v3"));
		logmsg->emit(LOG_INFO,"time needed for creating %d Views of size 5: %.3g[s]",num_views,timer->get_seconds_between("v5", "v4"));
		logmsg->emit(LOG_INFO,"time needed for creating %d indirect Views of size 3: %.3g[s]",num_views,timer->get_seconds_between("v6", "v5"));
		logmsg->emit(LOG_INFO,"time needed for creating %d indirect ConstViews of size 3: %.3g[s]",num_views,timer->get_seconds_between("v7", "v6"));
		
		// --------------------------------------------------
		// test the speed of a buettiker probe code extract
		// --------------------------------------------------
		uint Nx = 50;
		CGEMatrix BuettH(Nx,Nx);
		CGEMatrix BuettGRM(Nx,Nx);
		CGEMatrix HGRM(1,Nx);	 // re-init to zero
		CGEMatrix tmp(1,1);
		timer->click("v8");
		for (uint ii=0; ii < num_views; ii++) {
			for (uint jj=1; jj<=Nx; jj++) 
			{
				GEMatrix::View HGRM_ij = HGRM(_(1,1),_(jj,jj)); //this->get_view(HGRM, 1, jj, 1);
				
				for (uint ll=1; ll<=Nx; ll++) 
				{					
					uint gi = 3;
					uint gl = ll;
					if (fabs(gi-gl) > 1 + 0.0001) continue;
					
					const GEMatrix::ConstView H_il = BuettH(_(gi,gi),_(gl,gl));
					
					const GEMatrix::ConstView GRM_lj = BuettGRM(_(ll,ll),_(jj,jj)); //this->get_const_view(BuettGRM,ll,jj,1);
					HGRM_ij += H_il*GRM_lj;
				}
			}
		}
		timer->click("v9");
		for (uint ii=0; ii < num_views; ii++) {
			for (uint jj=1; jj<=Nx; jj++) 
			{
				cplx & HGRM_ij = HGRM(1,jj);
				
				for (uint ll=1; ll<=Nx; ll++) 
				{					
					uint gi = 3;
					uint gl = ll;
					if (fabs(gi-gl) > 1 + 0.0001) continue;
					
					const cplx & H_il = BuettH(gi,gl);
					
					const cplx & GRM_lj = BuettGRM(ll,jj);
					
					HGRM_ij += H_il*GRM_lj;
				}
			}
		}
		timer->click("v10");
		double time1 = timer->get_seconds_between("v9", "v8");
		double time2 = timer->get_seconds_between("v10", "v9");
		logmsg->emit(LOG_INFO,"%d*Buettiker time (Nx=%d) needed w/ views: %.3g[s], w/ cplx: %.3g[s], factor=%.2g",
				num_views,Nx,time1,time2,time1/time2);
		
		// ---------------------------------------------------------
		// test Views and conjugateTranspose
		// ---------------------------------------------------------
		uint Nn = 10;
		cplx imag_unit(0.,1.0);
		CGEMatrix L1(Nn,Nn);
		for (uint ii=1; ii<=Nn; ii++) {
			for (uint jj=1; jj<=Nn; jj++) {
				L1(ii,jj) = double(ii) + imag_unit * double(jj);
			}
		}
		const CGEMatrix::ConstView L1_24 = L1(_(2,4),_(2,4));
		CGEMatrix L2(3,3);
		L2 = conjugateTranspose(L1_24);
		cout << "L2=" << L2 << endl;
		for (uint ii=1; ii<=3; ii++) {
			for (uint jj=1; jj<=3; jj++) {
				TS_ASSERT(   fabs(L2(ii,jj).real() - L1(1+jj,1+ii).real()) < 1e-14 
						  && fabs(L2(ii,jj).imag() + L1(1+jj,1+ii).imag()) < 1e-14);
			}
		}
		CGEMatrix::View L1_68 = L1(_(6,8),_(6,8));
		L1_68 = conjugateTranspose(L1_24);;
		for (uint ii=1; ii<=3; ii++) {
			for (uint jj=1; jj<=3; jj++) {
				TS_ASSERT(   fabs(L1_68(ii,jj).real() - L1(1+jj,1+ii).real()) < 1e-14 
						  && fabs(L1_68(ii,jj).imag() + L1(1+jj,1+ii).imag()) < 1e-14);
			}
		}
		L1_68 = -conjugateTranspose(L1_24);
		for (uint ii=1; ii<=3; ii++) {
			for (uint jj=1; jj<=3; jj++) {
				TS_ASSERT(   fabs(L1_68(ii,jj).real() + L1(1+jj,1+ii).real()) < 1e-14 
						  && fabs(L1_68(ii,jj).imag() - L1(1+jj,1+ii).imag()) < 1e-14);
			}
		}
		
		CGEMatrix & L1_ref = L1;
		const CGEMatrix::ConstView L1_22 = L1_ref(_(2,2),_(2,2));
		CGEMatrix L3(1,1);
		L3 = -conjugateTranspose(L1_22);
		cout << "L3=" << L3 << endl;
		TS_ASSERT(fabs(L3(1,1).real() + L1(2,2).real()) < 1e-14 
			   && fabs(L3(1,1).imag() - L1(2,2).imag()) < 1e-14);
		
		fout.close();	
	}
	
	      CGEMatrix::View      get_view      (      CGEMatrix & X, uint start, uint end, uint Nn);
	const CGEMatrix::ConstView get_const_view(const CGEMatrix & X, uint start, uint end, uint Nn);
	
	
};

CGEMatrix::View FlensTest::get_view(CGEMatrix & X, uint xx, uint yy, uint Nn) {
	NEGF_ASSERT(xx>=1 && yy>=1, "FLENS indices start with 1!");
	NEGF_ASSERT(X.numRows()>=(int)xx && X.numCols()>=(int)yy, "index out of range.");
		return X(_((xx-1)*Nn+1,(xx-1)*Nn+Nn),_((yy-1)*Nn+1,(yy-1)*Nn+Nn));
}

const CGEMatrix::ConstView FlensTest::get_const_view(const CGEMatrix & X, uint xx, uint yy, uint Nn) {
	NEGF_ASSERT(xx>=1 && yy>=1, "FLENS indices start with 1!");
	NEGF_ASSERT(X.numRows()>=(int)xx && X.numCols()>=(int)yy, "index out of range.");
		return X(_((xx-1)*Nn+1,(xx-1)*Nn+Nn),_((yy-1)*Nn+1,(yy-1)*Nn+Nn));
}

} // end namespace negf

#endif /*FLENS_DEFINITION_H_*/
