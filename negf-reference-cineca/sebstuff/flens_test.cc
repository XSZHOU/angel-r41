/* compile using:
g++4 -c -Wall -fPIC -g -m64 -I/home/steiger/src/release/amd64/flens/ -I/home/steiger/src/release/amd64/gsl/cblas/ -o flens_test.o flens_test.cc
g++4 -Wall -fPIC -g -m64 -L/home/steiger/src/release/amd64/lib/ -lflens -lacml -lacml_mv -lgfortran -o flens_test flens_test.o
*/

#include <complex>
typedef std::complex<double> cplx;

#include <flens/flens.h>
typedef flens::GeMatrix<flens::FullStorage<cplx, flens::ColMajor> > 	GEMatrix; 
typedef flens::DenseVector<flens::Array<cplx> > 						DEVector; 

int main(int argc, char* argv[])
{
	GEMatrix A(4,4), B(4,4), C(4,4);
	A = 1, 2, 3, 4,
		5, 6, 7, 8,
		9, 8, 7, 6,
		5, 4, 3, 20;
	B = 8, 9, 6, 5,
		3, 1, 7, 3,
		6, 0, 3, 1,
		3, 6, 7, 88;
	
	C = A - B;
	C = A - conjugateTranspose(B);
	C = A - conjugateTranspose(A);
	C = A * B;
	
	DEVector::ConstView col = A(2,flens::_);
	double vec_norm = flens::nrm2(col);
	
	flens::DenseVector<flens::Array<int> > p(4);
	flens::trf(A,p);	// LU factorization
	flens::tri(A,p);	// backward substitution column by column
	
}
