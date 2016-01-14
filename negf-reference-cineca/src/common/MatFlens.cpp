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
#include "Mat.h"
using namespace negf;
 
#ifdef FLENS

GEMatrix negf::work1(1,1); // will be resized to Nn in Options.cpp (as soon as Nn is known)
GEMatrix negf::work2(1,1);
GEMatrix negf::work3(1,1);

Matc::Matc() {
	this->flens_matrix = GEMatrix(1,1);
	this->num_offdiags = 0;
}
Matc::Matc(int m, int n) {
	this->flens_matrix = GEMatrix(m,n);
	this->num_offdiags = m-1;
}
Matc::Matc(const Matc & copy):
	flens_matrix(copy.flens_matrix) {
	this->num_offdiags = copy.num_offdiags;
}
Matc::Matc(const GEMatrix & flens_matrix_):
	flens_matrix(flens_matrix_) {
	this->num_offdiags = flens_matrix.numRows()-1;
}
Matc::Matc(const GEMatrix::ConstView flens_matrix_) {
	this->flens_matrix = flens_matrix_;
	this->num_offdiags = flens_matrix.numRows()-1;
}
Matc::~Matc() {}

void Matc::get_submatrix(uint x1, uint x2, uint y1, uint y2, Matc & result) const {
	const GEMatrix::ConstView v = this->flens_matrix(flens::_(x1,x2),flens::_(y1,y2));
	result = v;
}

Matd::Matd() {
	this->flens_matrix = DGEMatrix(1,1);
}
Matd::Matd(int m, int n) {
	this->flens_matrix = DGEMatrix(m,n);
}
Matd::Matd(const DGEMatrix & flens_matrix_):
	flens_matrix(flens_matrix_) {
}
Matd::Matd(const DGEMatrix::ConstView flens_matrix_) {
	this->flens_matrix = flens_matrix_;
}
#ifdef CSCS
Matd::Matd(const Matd & copy) {
	this->flens_matrix = copy.flens_matrix;
}
#endif
Matd::~Matd() {}


Vecd::Vecd(int n) {
	this->flens_vector = DDEVector(n);
}
Vecd::Vecd(const DDEVector & flens_vector_):
	flens_vector(flens_vector_) {
}
Vecd::~Vecd() {}


void negf::solve_linear_problem(/*const*/ Matd & A, Vecd & rhs) {
	flens::DenseVector<flens::Array<int> > pivots(rhs.length()+1);
	flens::sv(A.flens_matrix, pivots, rhs.flens_vector);
}

void negf::mult(const Matd & A, const Vecd & v, Vecd & rhs) {
	rhs.flens_vector = A.flens_matrix * v.flens_vector;
}

void negf::invert(Matc & A) 
{STACK_TRACE(
	flens::DenseVector<flens::Array<int> > p(A.num_cols());	// initialization of right length is mandatory!
	int error = flens::trf(A.flens_matrix,p);				// LU factorization
	NEGF_FASSERT(error==0, "error %d encountered when LU-factorising a complex matrix of size %d.",error,A.num_cols());
	error = flens::tri(A.flens_matrix,p);					// backward substitution column by column
	NEGF_FASSERT(error==0, "error %d encountered when backward-substituting a complex matrix of size %d.",error,A.num_cols());
);}
	
void negf::invert(Matd & A) 
{STACK_TRACE(
	flens::DenseVector<flens::Array<int> > p(A.num_cols());	// initialization of right length is mandatory!
	int error = flens::trf(A.flens_matrix,p);				// LU factorization
	NEGF_FASSERT(error==0, "error %d encountered when LU-factorising a real matrix of size %d.",error,A.num_cols());
	error = flens::tri(A.flens_matrix,p);					// backward substitution column by column
	NEGF_FASSERT(error==0, "error %d encountered when backward-substituting a double matrix of size %d.",error,A.num_cols());
);}

double negf_math::matrix_norm(const Matc & matrix)
{
	double result = 0.0;
	for (uint ii = matrix.first_row(); ii <= matrix.last_row(); ii++)
	{
		DEVector::ConstView col = matrix.flens_matrix(ii,flens::_);
		cplx vec_norm = flens::nrm2(col); // double does not work
		result += std::abs(vec_norm*vec_norm);
	}
	result = negf_math::sqrt(result);
	return result;
}

double negf_math::matrix_norm(const Matd & matrix)
{
	double result = 0.0;
	for (uint ii = matrix.first_row(); ii <= matrix.last_row(); ii++)
	{
		DDEVector::ConstView col = matrix.flens_matrix(ii,flens::_);
		double vec_norm = flens::nrm2(col);
		result += std::fabs(vec_norm*vec_norm);
	}
	result = negf_math::sqrt(result);
	return result;
}

// (xx-1)*Nn+nn
// (nn-1)*Nx+xx

void Matc::fill_block(uint xx, uint yy, const Matc & block_matrix, const uint & Nx) 
{STACK_TRACE(
	NEGF_ASSERT(this->num_rows()==Nx*Nn, "number of rows is not Nx*Nn");
	NEGF_ASSERT(block_matrix.num_rows()==Nn, "number of rows of block_matrix is not Nn");
	for (uint mm=1; mm<=Nn; mm++) {
		for (uint nn=1; nn<=Nn; nn++) {
			NEGF_ASSERT(flens_matrix.numRows()>=int(get_mat_idx(xx,mm,Nx)) && flens_matrix.numCols()>=int(get_mat_idx(yy,nn,Nx)), "fill_block1 (1)");
			NEGF_ASSERT(block_matrix.num_rows()>=mm           && block_matrix.num_cols()>=nn,           "fill_block1 (2)");
			this->flens_matrix(get_mat_idx(xx,mm,Nx),get_mat_idx(yy,nn,Nx)) = block_matrix(mm,nn);
		}
	}
);}

void Matc::fill_block(uint xx, uint yy, const Matc & whole_matrix, uint x2, uint y2, const uint & Nx) 
{STACK_TRACE(
	NEGF_ASSERT(this->num_rows()==Nx*Nn, "number of rows is not Nx*Nn");
	NEGF_ASSERT(whole_matrix.num_rows()==Nx*Nn, "number of rows of whole_matrix is not Nx*Nn");
	for (uint mm=1; mm<=Nn; mm++) {
		for (uint nn=1; nn<=Nn; nn++) {
			NEGF_ASSERT(flens_matrix.numRows()>=get_mat_idx(xx,mm,Nx) && flens_matrix.numCols()>=get_mat_idx(yy,nn,Nx), "fill_block2 (1)");
			NEGF_ASSERT(int(whole_matrix.num_rows())>=get_mat_idx(x2,mm,Nx) && int(whole_matrix.num_cols())>=get_mat_idx(y2,nn,Nx), "fill_block2 (2)");
			this->flens_matrix(get_mat_idx(xx,mm,Nx),get_mat_idx(yy,nn,Nx)) = whole_matrix(get_mat_idx(x2,mm,Nx),get_mat_idx(y2,nn,Nx));
		}
	}
);}

/*
void Matc::fill_block(uint xx, uint yy, const Matc & whole_matrix, uint x2, uint y2, uint Nvert) 
{STACK_TRACE(
	NEGF_ASSERT(this->num_rows()==NxNn, "number of rows is not Nx*Nn");
	for (uint mm=1; mm<=Nn; mm++) {
		for (uint nn=1; nn<=Nn; nn++) {
			NEGF_ASSERT(flens_matrix.numRows()>=(mm-1)*Nvert+xx && flens_matrix.numCols()>=(nn-1)*Nvert+yy, "fill_block2 (1)");
			NEGF_ASSERT(whole_matrix.num_rows()>=(mm-1)*Nx+x2 && whole_matrix.num_cols()>=(nn-1)*Nx+y2, "fill_block2 (2)");
			this->flens_matrix((mm-1)*Nvert+xx,(nn-1)*Nvert+yy) = whole_matrix(get_mat_idx(x2,mm),get_mat_idx(y2,nn));
		}
	}
);}*/

void Matd::fill_block(uint xx, uint yy, const Matd & block_matrix, const uint & Nx) 
{STACK_TRACE(
	NEGF_ASSERT(this->num_rows()==Nx*Nn, "number of rows is not Nx*Nn");
	NEGF_ASSERT(block_matrix.num_rows()==Nn, "number of rows of block_matrix is not Nn");
	for (uint mm=1; mm<=Nn; mm++) {
		for (uint nn=1; nn<=Nn; nn++) {
			NEGF_ASSERT(flens_matrix.numRows()>=int(get_mat_idx(xx,mm,Nx)) && flens_matrix.numCols()>=int(get_mat_idx(yy,nn,Nx)), "fill_block3 (1)");
			NEGF_ASSERT(block_matrix.num_rows()>=mm          && block_matrix.num_cols()>=nn,          "fill_block3 (2)");
			this->flens_matrix(get_mat_idx(xx,mm,Nx),get_mat_idx(yy,nn,Nx)) = block_matrix(mm,nn);
		}
	}
);}

void Matc::subtract_block(uint xx, uint yy, const Matc & whole_matrix, uint x2, uint y2, const uint & Nx) 
{STACK_TRACE(
	NEGF_ASSERT(this->num_rows()==Nx*Nn, "number of rows is not Nx*Nn");
	NEGF_ASSERT(whole_matrix.num_rows()==Nx*Nn, "number of rows of whole_matrix is not Nx*Nn");
	for (uint mm=1; mm<=Nn; mm++) {
		for (uint nn=1; nn<=Nn; nn++) {
			NEGF_ASSERT(flens_matrix.numRows()>=get_mat_idx(xx,mm,Nx) && flens_matrix.numCols()>=get_mat_idx(yy,nn,Nx), "subtract_block (1)");
			NEGF_ASSERT(int(whole_matrix.num_rows())>=get_mat_idx(x2,mm,Nx) && int(whole_matrix.num_cols())>=get_mat_idx(y2,nn,Nx), "subtract_block (2)");
			this->flens_matrix(get_mat_idx(xx,mm,Nx),get_mat_idx(yy,nn,Nx)) -= whole_matrix(get_mat_idx(x2,mm,Nx),get_mat_idx(y2,nn,Nx));
		}
	}
);}

void Matc::get_block(uint xx, uint yy, Matc & result, const uint & Nx) const 
{STACK_TRACE(
	NEGF_ASSERT(this->num_rows()==Nx*Nn, "number of rows is not Nx*Nn");
	result = Matc(Nn,Nn);
	for (uint mm=1; mm<=Nn; mm++) {
		for (uint nn=1; nn<=Nn; nn++) {
			NEGF_ASSERT(flens_matrix.numRows()>=int(get_mat_idx(xx,mm,Nx)) && flens_matrix.numCols()>=int(get_mat_idx(yy,nn,Nx)), "get_block");
			result(mm,nn) = this->flens_matrix(get_mat_idx(xx,mm,Nx),get_mat_idx(yy,nn,Nx));
		}
	}
);}
		
void Matc::multiply_block(uint xx, uint yy, const cplx & alpha, const uint & Nx) 
{STACK_TRACE(
	NEGF_ASSERT(this->num_rows()==Nx*Nn, "number of rows is not Nx*Nn");
	for (uint mm=1; mm<=Nn; mm++) {
		for (uint nn=1; nn<=Nn; nn++) {
			NEGF_ASSERT(flens_matrix.numRows()>=int(get_mat_idx(xx,mm,Nx)) && flens_matrix.numCols()>=int(get_mat_idx(yy,nn,Nx)), "multiply_block");
			this->flens_matrix(get_mat_idx(xx,mm,Nx),get_mat_idx(yy,nn,Nx)) *= alpha;
		}
	}
);}

void negf::multiply_blocks(const Matc & A, uint xa, uint ya, const Matc & B, uint xb, uint yb, Matc & C, const uint & Nx) 
{STACK_TRACE(	
	NEGF_ASSERT(A.num_rows()==Nx*Nn, "number of rows of A is not Nx*Nn");
	NEGF_ASSERT(B.num_rows()==Nx*Nn, "number of rows of B is not Nx*Nn");
	NEGF_ASSERT(C.num_rows()==Nn, "number of rows of C is not Nn");
	for (uint mm=1; mm<=Nn; mm++) {
		for (uint nn=1; nn<=Nn; nn++) {
			NEGF_ASSERT(A.flens_matrix.numRows()>=int(get_mat_idx(xa,mm,Nx)) && A.flens_matrix.numCols()>=int(get_mat_idx(ya,nn,Nx)), "multiply_blocks1 (1)");
			NEGF_ASSERT(B.flens_matrix.numRows()>=int(get_mat_idx(xb,mm,Nx)) && B.flens_matrix.numCols()>=int(get_mat_idx(yb,nn,Nx)), "multiply_blocks1 (2)");
			work1(mm,nn) = A.flens_matrix(get_mat_idx(xa,mm,Nx),get_mat_idx(ya,nn,Nx));
			work2(mm,nn) = B.flens_matrix(get_mat_idx(xb,mm,Nx),get_mat_idx(yb,nn,Nx));
		}
	}
	C.flens_matrix = work1*work2;
);}

void negf::multiply_blocks(const Matc & A, uint xa, uint ya, const Matc & B, uint xb, uint yb, Matc & C, uint xc, uint yc, const uint & Nx) 
{STACK_TRACE(	
	NEGF_ASSERT(A.num_rows()==Nx*Nn, "number of rows of A is not Nx*Nn");
	NEGF_ASSERT(B.num_rows()==Nx*Nn, "number of rows of B is not Nx*Nn");
	NEGF_ASSERT(C.num_rows()==Nx*Nn, "number of rows of C is not Nx*Nn");
	for (uint mm=1; mm<=Nn; mm++) {
		for (uint nn=1; nn<=Nn; nn++) {
			work1(mm,nn) = A.flens_matrix(get_mat_idx(xa,mm,Nx),get_mat_idx(ya,nn,Nx));
			work2(mm,nn) = B.flens_matrix(get_mat_idx(xb,mm,Nx),get_mat_idx(yb,nn,Nx));
		}
	}
	
	work3 = work1*work2;
	
	for (uint mm=1; mm<=Nn; mm++) {
		for (uint nn=1; nn<=Nn; nn++) {
			C.flens_matrix(get_mat_idx(xc,mm,Nx),get_mat_idx(yc,mm,Nx)) = work3(mm,nn);
		}
	}
);}

#endif // FLENS
