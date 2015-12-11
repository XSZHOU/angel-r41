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

#ifdef BANDED

BMatc::BMatc(int m_, int n_, int num_offdiags_) {
// note: REORDER flag puts band=band'-terms near diagonal - we would not expect NEGF to work with banded
// matrices when this flag is not set because x'=x+-1 would otherwise be far-offdiagonal and not stored
	
	NEGF_ASSERT(m_==n_, "m==n required for banded matrices");
	NEGF_ASSERT(num_offdiags_>=0, "nonnegative num_nonzeros expected!");
	this->n    = m_;
	this->info = 0;
	this->zero = 0;
	this->nrhs = n;
	
	this->num_offdiags = num_offdiags_;
	this->KL   = this->num_offdiags;
	this->KU   = this->num_offdiags;
	this->lda  = KL + KU + 1;
	
	this->num_nonzeros = n*lda;
	this->data = new cplx[num_nonzeros];
}

	
BMatc::~BMatc() {
	delete [] data;
}

BMatc::BMatc(const BMatc & copy) {
	this->n            = copy.n;
	this->num_offdiags = copy.num_offdiags;
	this->KL           = copy.KL;
	this->KU           = copy.KU;
	this->nrhs         = copy.nrhs;
	this->lda          = copy.lda;
	this->info         = 0;
	this->zero         = 0;
	this->num_nonzeros = copy.num_nonzeros;
	this->data         = new cplx[num_nonzeros];
	for (uint i=0; i<num_nonzeros; i++) {
		this->data[i] = copy.data[i];
	}
}


// ======================================================================================================
// BANDED STORAGE SCHEME:
//
// A NxN banded matrix with KL sub- and KU super-diagonals is stored in a LDAxN-matrix where LDA=KL+KU+1.
// 
// This matrix is stored in an array according to A(ib,j) = data[N*ib+j], 0<=j<=N, 0<=ib<=LDA
// 
// A(ib,j) corresponds to column j and row j-KU-1+ib of the original matrix
//
// Conversely, row i and column j of the original matrix can be found in data[N*ib+j], where ib=i-j+KU+1.
//
// Finally, entry k of the data-array corresponds to A(i,j) where j = k mod N and i = k div N + j - KU - 1
//
// ======================================================================================================


// get_submatrix is irrespective of how bands/positions are arranged within the matrix
void BMatc::get_submatrix(uint x1, uint x2, uint y1, uint y2, Matc & result) const {
	if (result.num_rows()!=x2-x1+1 || result.num_cols()!=y2-y1+1) {
		result = Matc(x2-x1+1, y2-y1+1);
	}
	for (uint ii=x1; ii<=x2; ii++) {
		for (uint jj=y1; jj<=y2; jj++) {
			if (this->stores(ii,jj)) {
				result(ii-x1+1,jj-y1+1) = (*this)(ii,jj);
			}
		}
	}
}

double negf_math::matrix_norm(const BMatc & matrix)
{
	double result = 0.0;
	for (uint ii = 0; ii < matrix.num_nonzeros; ii++) {
		result += std::abs(matrix.data[ii]*matrix.data[ii]);
	}
	result = negf_math::sqrt(result);
	return result;
}

// ZGBSV computes the solution to a complex system of linear equations A*X=B, where A is a band matrix of order N with 
// KL subdiagonals and KU superdiagonals, and X and B are N-by-NRHS matrices.
// zgbsv_(&n, &KL, &KU, &nrhs, A, &ldab, ipivot, B, &ldb, &info);


// BLOCKS: bands (m,n) position (x,y) (both one-based) is stored in (xx-1)*Nx+xx, (nn-1)*Nx+yy
// these indices are accessed via get_mat_idx!

void BMatc::fill_block(uint xx, uint yy, const Matc & block_matrix, const uint & Nx) {
	for (uint mm=1; mm<=Nn; mm++) {
		for (uint nn=1; nn<=Nn; nn++) {
			uint i = get_mat_idx(xx,mm,Nx) - 1;	// i, j are 0-based; xx,yy are 1-based
			uint j = get_mat_idx(yy,nn,Nx) - 1;
			if (this->stores(i+1,j+1)) {
				data[this->n*get_ib(i,j)+j] = block_matrix(mm,nn);
			}
		}
	}
}

void BMatc::fill_block(uint xx, uint yy, const Matc & whole_matrix, uint x2, uint y2, const uint & Nx) {
	for (uint mm=1; mm<=Nn; mm++) {
		for (uint nn=1; nn<=Nn; nn++) {
			uint i2 = get_mat_idx(x2,mm,Nx); // 1-based
			uint j2 = get_mat_idx(y2,nn,Nx);
			
			uint i = get_mat_idx(xx,mm,Nx) - 1;	// i, j are 0-based; xx,yy are 1-based
			uint j = get_mat_idx(yy,nn,Nx) - 1;
			if (this->stores(i+1,j+1)) {
				data[this->n*get_ib(i,j)+j] = whole_matrix(i2,j2);
			}
		}
	}
}

void BMatc::fill_block(uint xx, uint yy, const BMatc & whole_matrix, uint x2, uint y2, const uint & Nx) {
	for (uint mm=1; mm<=Nn; mm++) {
		for (uint nn=1; nn<=Nn; nn++) {
			uint i2 = get_mat_idx(x2,mm,Nx); // 1-based
			uint j2 = get_mat_idx(y2,nn,Nx);
			
			uint i = get_mat_idx(xx,mm,Nx) - 1;	// i, j are 0-based; xx,yy are 1-based
			uint j = get_mat_idx(yy,nn,Nx) - 1;
			if (this->stores(i+1,j+1)) {
				if (whole_matrix.stores(i2,j2)) {
					data[this->n*get_ib(i,j)+j] = whole_matrix(i2,j2);
				} else {
					data[this->n*get_ib(i,j)+j] = 0.0;
				}
			}
		}
	}
}

void BMatc::subtract_block(uint xx, uint yy, const Matc & whole_matrix, uint x2, uint y2, const uint & Nx) {
	for (uint mm=1; mm<=Nn; mm++) {
		for (uint nn=1; nn<=Nn; nn++) {
			uint i2 = get_mat_idx(x2,mm,Nx); // 1-based
			uint j2 = get_mat_idx(y2,nn,Nx);
			
			uint i = get_mat_idx(xx,mm,Nx) - 1;	// i, j are 0-based; xx,yy are 1-based
			uint j = get_mat_idx(yy,nn,Nx) - 1;
			if (this->stores(i+1,j+1)) {
				data[this->n*get_ib(i,j)+j] -= whole_matrix(i2,j2);
			}
		}
	}
}

void BMatc::subtract_block(uint xx, uint yy, const BMatc & whole_matrix, uint x2, uint y2, const uint & Nx) {
	for (uint mm=1; mm<=Nn; mm++) {
		for (uint nn=1; nn<=Nn; nn++) {
			uint i2 = get_mat_idx(x2,mm,Nx); // 1-based
			uint j2 = get_mat_idx(y2,nn,Nx);
			
			uint i = get_mat_idx(xx,mm,Nx) - 1;	// i, j are 0-based; xx,yy are 1-based
			uint j = get_mat_idx(yy,nn,Nx) - 1;
			if (this->stores(i+1,j+1)) {
				if (whole_matrix.stores(i2,j2)) {
					data[this->n*get_ib(i,j)+j] -= whole_matrix(i2,j2);
				} else {
					// subtract nothing
				}
			}
		}
	}
}
		
void BMatc::multiply_block(uint xx, uint yy, const cplx & alpha, const uint & Nx) {
	for (uint mm=1; mm<=Nn; mm++) {
		for (uint nn=1; nn<=Nn; nn++) {
			uint i = get_mat_idx(xx,mm,Nx) - 1;	// i, j are 0-based; xx,yy are 1-based
			uint j = get_mat_idx(yy,nn,Nx) - 1;
			if (this->stores(i+1,j+1)) {
				data[n*get_ib(i,j)+j] *= alpha;
			}
		}
	}
}

void BMatc::get_block(uint xx, uint yy, Matc & result, const uint & Nx) const {
	if (result.num_rows()!=Nn || result.num_cols()!=Nn) {
		result = Matc(Nn,Nn);
	}
	for (uint mm=1; mm<=Nn; mm++) {
		for (uint nn=1; nn<=Nn; nn++) {
			uint i = get_mat_idx(xx,mm,Nx) - 1;	// i, j are 0-based; xx,yy are 1-based
			uint j = get_mat_idx(yy,nn,Nx) - 1;
			if (this->stores(i+1,j+1)) {
				result(mm,nn) = data[n*get_ib(i,j)+j];
			}
		}
	}
}

#endif
