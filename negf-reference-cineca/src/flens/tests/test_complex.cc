#include <tests/evaltest.h>
#include <flens/flens.h>
#include <complex>

typedef std::complex<double> cplx;
typedef flens::GeMatrix<flens::FullStorage<cplx, flens::ColMajor> > 	GEMatrix; 
typedef flens::DenseVector<flens::Array<cplx> > 						DEVector;

void 
test1()
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
	cerr << "C = A - B = " << C << endl;										
	C = A + B;
	cerr << "C = A + B = " << C << endl;	
	C = A * B;
	cerr << "C = A * B = " << C << endl;	
}

void
test2()
{
	GEMatrix B(4,4), C(4,4);
	B = 8, 9, 6, 5,
		3, 1, 7, 3,
		6, 0, 3, 1,
		3, 6, 7, 88;

	B = 2.0 * B;
	cerr << "B = 2.0 * B = " << B << endl;		
	B = cplx(1.0, 1.0) * B;
	cerr << "B = (1.0 + i) * B = " << B << endl;
	
}

void
test3()
{	
	GEMatrix B(4,4), C(4,4);
	B = 8, 9, 6, 5,
		3, 1, 7, 3,
		6, 0, 3, 1,
		3, 6, 7, 88;
	B = cplx(1.0, 1.0) * B;
	C = transpose(B);
	cerr << "C = B^t = " << C << endl;	
	C = conjugateTranspose(B);
	cerr << "C = conj(B^t) = " << C << endl;
	DEVector::ConstView col = C(2,flens::_);
	cplx vec_norm = flens::nrm2(col);
	cerr << "norm2(C(2,_)) = " << vec_norm << endl;
}

void
test4()
{	
	GEMatrix A(4,4), B(4,4), C(4,4);
	A = 8, 9, 6, 5,
		3, 1, 7, 3,
		6, 0, 3, 1,
		3, 6, 7, 88;
	B = A;
	flens::DenseVector<flens::Array<int> > p(4);
	flens::trf(A,p);	// LU factorization
	flens::tri(A,p);	// backward substitution column by column
	
	C = A * B;
	cerr << "C = inv(A) * A = " << C;		
} 


int main(int argc, char* argv[])
{	
	run("test1", test1, true);
	run("test2", test2, true);
	run("test3", test3, true);
	run("test4", test4, true);
}
