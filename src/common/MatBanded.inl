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
#ifdef BANDED
	
inline const uint BMatc::num_rows()  const { return this->n; }
inline const uint BMatc::num_cols()  const { return this->n; }
inline const uint BMatc::first_row() const { return 1; }
inline const uint BMatc::last_row()  const { return this->n; }
inline const uint BMatc::first_col() const { return 1; }
inline const uint BMatc::last_col()  const { return this->n; }

// ib will be 0-based. note that (i-1)-(j-1)=i-j, so it does not matter if i,j are both 0- or 1-based!
inline int BMatc::get_ib(const uint & i, const uint & j) const { return i-j+KU/*+1*/; }

// assumption: i, j are 1-based!
inline bool BMatc::stores(uint ii, uint jj) const {
	return (ii >= max(1.0, double(jj)-this->KU) && ii <= min(double(this->n), double(jj)+this->KL));
}

// assumption: i, j are 1-based!
inline const cplx   & BMatc::operator()(const int & i, const int & j) const { 
	if (this->stores(i,j)) {
		return data[n*get_ib(i,j)+(j-1)];
	} else {
		return this->zero;
	}
}

inline const cplx   & BMatc::nocheck(const int & i, const int & j) const { 
	return data[n*get_ib(i,j)+(j-1)];
}

// assumption: i, j are 1-based!
inline       cplx   & BMatc::operator()(const int & i, const int & j)       {
	if (this->stores(i,j)) {
		return data[n*get_ib(i,j)+(j-1)];
	} else {
		NEGF_FEXCEPTION("cannot reference non-constant banded matrix entry outside of band! (i=%d, j=%d, num_offdiags=%d)",i,j,num_offdiags);
	}
}

inline BMatc & BMatc::operator*=(const cplx & alpha)   { 
	for (uint ii=0; ii<num_nonzeros; ii++) {
		this->data[ii] *= alpha;
	}
	return *this;
}

inline BMatc & BMatc::operator+=(const cplx & alpha)   { 
	for (uint ii=0; ii<num_nonzeros; ii++) {
		this->data[ii] += alpha;
	}
	return *this; 
}

inline BMatc & BMatc::operator-=(const cplx & alpha)   { 
	for (uint ii=0; ii<num_nonzeros; ii++) {
		this->data[ii] -= alpha;
	}
	return *this; 
}

inline BMatc & BMatc::operator+=(const BMatc & A)   	 { 
	NEGF_ASSERT(A.n==this->n, "different matrix sizes (BMatc::operator+=(const BMatc & A)).");
	if (A.num_offdiags == this->num_offdiags) {
		for (uint ii=0; ii<num_nonzeros; ii++) {
			this->data[ii] += A.data[ii];
		}
	} else if (A.num_offdiags > this->num_offdiags) {
		for (int jj=0; jj<n; jj++) { 											// 0-based
			for (int ii=max(0,jj-this->KU); ii<min(n,jj+this->KL+1); ii++) {	// 0-based
				NEGF_FASSERT(n*get_ib(ii,jj)+jj < int(num_nonzeros), "ii=%d jj=%d ib=%d n=%d: n*ib+jj=%d, num_nonzeros=%d", ii, jj, get_ib(ii,jj), n, n*get_ib(ii,jj)+jj, num_nonzeros);
				this->data[n*get_ib(ii,jj)+jj] += A(ii+1, jj+1);				// A(.,.) expects 1-based indices
	      	}
		}
	} else { // A.num_offdiags < this->num_offdiags
		for (int jj=0; jj<n; jj++) { 											// 0-based
			for (int ii=max(0,jj-A.KU); ii<min(n,jj+A.KL+1); ii++) {			// 0-based
				NEGF_FASSERT(n*get_ib(ii,jj)+jj < int(num_nonzeros), "ii=%d jj=%d ib=%d n=%d: n*ib+jj=%d, num_nonzeros=%d", ii, jj, get_ib(ii,jj), n, n*get_ib(ii,jj)+jj, num_nonzeros);
				this->data[n*get_ib(ii,jj)+jj] += A(ii+1, jj+1);				// A(.,.) expects 1-based indices
	      	}
		}
	}
	return *this;
}

inline BMatc & BMatc::operator-=(const BMatc & A)   	 { 
	NEGF_ASSERT(A.n==this->n, "different matrix sizes (BMatc::operator-=(const BMatc & A)).");
	if (A.num_offdiags == this->num_offdiags) {
		for (uint ii=0; ii<num_nonzeros; ii++) {
			this->data[ii] -= A.data[ii];
		}
	} else if (A.num_offdiags > this->num_offdiags) {
		for (int jj=0; jj<n; jj++) { 											// 0-based
			for (int ii=max(0,jj-this->KU); ii<min(n,jj+this->KL+1); ii++) {	// 0-based
				NEGF_FASSERT(n*get_ib(ii,jj)+jj < int(num_nonzeros), "ii=%d jj=%d ib=%d n=%d: n*ib+jj=%d, num_nonzeros=%d", ii, jj, get_ib(ii,jj), n, n*get_ib(ii,jj)+jj, num_nonzeros);
				this->data[n*get_ib(ii,jj)+jj] -= A(ii+1, jj+1);				// A(.,.) expects 1-based indices
	      	}
		}
	} else { // A.num_offdiags < this->num_offdiags
		for (int jj=0; jj<n; jj++) { 											// 0-based
			for (int ii=max(0,jj-A.KU); ii<min(n,jj+A.KL+1); ii++) {			// 0-based
				NEGF_FASSERT(n*get_ib(ii,jj)+jj < int(num_nonzeros), "ii=%d jj=%d ib=%d n=%d: n*ib+jj=%d, num_nonzeros=%d", ii, jj, get_ib(ii,jj), n, n*get_ib(ii,jj)+jj, num_nonzeros);
				this->data[n*get_ib(ii,jj)+jj] -= A(ii+1, jj+1);				// A(.,.) expects 1-based indices
	      	}
		}
	}
	return *this; 
}

inline BMatc & BMatc::operator= (const BMatc & A)   	{ 
	NEGF_ASSERT(A.n==this->n, "different matrix sizes (BMatc::operator= (const BMatc & A)).");
	if (A.num_offdiags == this->num_offdiags) {
		for (uint ii=0; ii<num_nonzeros; ii++) {
			this->data[ii] = A.data[ii];
		}
	} else {
		for (int jj=0; jj<n; jj++) { 											// 0-based
			for (int ii=max(0,jj-this->KU); ii<min(n,jj+this->KL+1); ii++) {	// 0-based
				NEGF_FASSERT(n*get_ib(ii,jj)+jj < int(num_nonzeros), "ii=%d jj=%d ib=%d n=%d: n*ib+jj=%d, num_nonzeros=%d", ii, jj, get_ib(ii,jj), n, n*get_ib(ii,jj)+jj, num_nonzeros);
				this->data[n*get_ib(ii,jj)+jj] = A(ii+1, jj+1);					// A(.,.) expects 1-based indices
				// note that A(i,j) will just be zero when it is outside the band
	      	}
		}
	} 
	return *this; 
} 

// discard far-offdiagonal terms!
inline BMatc & BMatc::operator= (const Matc & A)   	{ 
	NEGF_FASSERT(A.num_rows()==this->num_rows() && A.num_cols()==this->num_cols(), "different matrix sizes (BMatc::operator=(const Matc & A)): A.num_rows()=%d, A.num_cols()=%d, B.num_rows()=%d, B.num_cols()=%d.",A.num_rows(),A.num_cols(),this->num_rows(),this->num_cols());
	for (int jj=0; jj<n; jj++) { 								// 0-based
		for (int ii=max(0,jj-KU); ii<min(n,jj+KL+1); ii++) {	// 0-based
			NEGF_FASSERT(n*get_ib(ii,jj)+jj < int(num_nonzeros), "ii=%d jj=%d ib=%d n=%d: n*ib+jj=%d, num_nonzeros=%d", ii, jj, get_ib(ii,jj), n, n*get_ib(ii,jj)+jj, num_nonzeros);
			this->data[n*get_ib(ii,jj)+jj] = A(ii+1, jj+1);		// A(.,.) expects 1-based indices
      	}
	}
	return *this; 
} 

// discard far-offdiagonal terms!
inline BMatc & BMatc::operator+=(const Matc & A)   	{ 
	NEGF_ASSERT(A.num_rows()==this->num_rows() && A.num_cols()==this->num_cols(), "different matrix sizes (BMatc::operator+=(const Matc & A)).");
	for (int jj=0; jj<n; jj++) { 								// 0-based
		for (int ii=max(0,jj-KU); ii<min(n,jj+KL+1); ii++) {	// 0-based
			NEGF_FASSERT(n*get_ib(ii,jj)+jj < int(num_nonzeros), "ii=%d jj=%d ib=%d n=%d: n*ib+jj=%d, num_nonzeros=%d",	ii, jj, get_ib(ii,jj), n, n*get_ib(ii,jj)+jj, num_nonzeros);
			this->data[n*get_ib(ii,jj)+jj] += A(ii+1, jj+1);				// A(.,.) expects 1-based indices
      	}
	}
	return *this; 
} 

// discard far-offdiagonal terms!
inline BMatc & BMatc::operator-=(const Matc & A)   	{ 
	NEGF_ASSERT(A.num_rows()==this->num_rows() && A.num_cols()==this->num_cols(), "different matrix sizes (BMatc::operator-=(const Matc & A)).");
	for (int jj=0; jj<n; jj++) { 								// 0-based
		for (int ii=max(0,jj-KU); ii<min(n,jj+KL+1); ii++) {	// 0-based
			NEGF_FASSERT(n*get_ib(ii,jj)+jj < int(num_nonzeros), "ii=%d jj=%d ib=%d n=%d: n*ib+jj=%d, num_nonzeros=%d", ii, jj, get_ib(ii,jj), n, n*get_ib(ii,jj)+jj, num_nonzeros);
			this->data[n*get_ib(ii,jj)+jj] -= A(ii+1, jj+1);				// A(.,.) expects 1-based indices
      	}
	}
	return *this; 
} 

// --------------------------------------------
// MATC DEFINITIONS INVOLVING BANDED MATRICES
// --------------------------------------------

inline Matc & Matc::operator= (const BMatc & A)   	{ 
	NEGF_ASSERT(A.num_rows()==this->num_rows() && A.num_cols()==this->num_cols(), "different matrix sizes (Matc::operator=(const BMatc & A)).");
	
	for (int ii=1; ii<=A.n; ii++) {
		int m1 = negf_math::max(1,ii-A.KL);
		int m2 = negf_math::min(A.n, ii+A.KU);		// ii, m1, m2 are 1-based
		for (int jj=1; jj<=m1-1; jj++) {
			(*this)(ii,jj) = 0.0;
		}
		for (int jj=m1; jj<=m2; jj++) {
			(*this)(ii,jj) = A(ii,jj);
		}
		for (int jj=m2+1; jj<=A.n; jj++) {
			(*this)(ii,jj) = 0.0;
		}
	}	
	return *this; 
} 

inline Matc & Matc::operator+=(const BMatc & A)   	{ 
	NEGF_ASSERT(A.num_rows()==this->num_rows() && A.num_cols()==this->num_cols(), "different matrix sizes (Matc::operator+=(const BMatc & A)).");
	for (int ii=1; ii<=A.n; ii++) {
		for (int jj=negf_math::max(1,ii-A.KL); jj<=negf_math::min(A.n,ii+A.KU); jj++) { // ii,jj are 1-based
			(*this)(ii,jj) += A(ii,jj);
		}
	}
	return *this; 
} 

inline Matc & Matc::operator-=(const BMatc & A)   	{ 
	NEGF_ASSERT(A.num_rows()==this->num_rows() && A.num_cols()==this->num_cols(), "different matrix sizes (Matc::operator-=(const BMatc & A)).");
	for (int ii=1; ii<=A.n; ii++) {
		for (int jj=negf_math::max(1,ii-A.KL); jj<=negf_math::min(A.n,ii+A.KU); jj++) { // ii,jj are 1-based
			(*this)(ii,jj) -= A(ii,jj);
		}
	}
	return *this; 
}

// ----------------------------------------
// COMPLEX BANDED OPERATIONS
// ----------------------------------------

// C = A+B
inline void add (const BMatc & A, const BMatc & B, BMatc & C) {
	NEGF_FASSERT(A.n==B.n, "different matrix sizes (add(BMatc,BMatc,BMatc)) (A.n=%d,B.n=%d).",A.n,B.n);
	NEGF_FASSERT(A.n==C.n, "different matrix sizes (add(BMatc,BMatc,BMatc)) (A.n=%d,C.n=%d).",A.n,C.n);
	if (C.num_offdiags==A.num_offdiags && B.num_offdiags==A.num_offdiags) {
		for (uint ii=0; ii<C.num_nonzeros; ii++) {
			C.data[ii] = A.data[ii] + B.data[ii];
		}
	} else {
		for (int jj=0; jj<C.n; jj++) { 								    // 0-based
			for (int ii=max(0,jj-C.KU); ii<min(C.n,jj+C.KL+1); ii++) {	// 0-based
				C(ii+1,jj+1) = A(ii+1,jj+1) + B(ii+1,jj+1);				// A(i,j), B(i,j), C(i,j) expect 1-based indices
	      	}
		}
	}
}

// C = A - B
inline void sub (const BMatc & A, const BMatc & B, BMatc & C) {
	NEGF_FASSERT(A.n==B.n, "different matrix sizes (sub(BMatc,BMatc,BMatc)) (A.n=%d,B.n=%d).",A.n,B.n);
	NEGF_FASSERT(A.n==C.n, "different matrix sizes (sub(BMatc,BMatc,BMatc)) (A.n=%d,C.n=%d).",A.n,C.n);
	if (C.num_offdiags==A.num_offdiags && B.num_offdiags==A.num_offdiags) {
		for (uint ii=0; ii<C.num_nonzeros; ii++) {
			C.data[ii] = A.data[ii] - B.data[ii];
		}
	} else {
		for (int jj=0; jj<C.n; jj++) { 								    // 0-based
			for (int ii=max(0,jj-C.KU); ii<min(C.n,jj+C.KL+1); ii++) {	// 0-based
				C(ii+1,jj+1) = A(ii+1,jj+1) - B(ii+1,jj+1);				// A(i,j), B(i,j), C(i,j) expect 1-based indices
	      	}
		}
	}
}

// B += a*A
inline void add (const BMatc & A, const  cplx & a, BMatc & B) {
	NEGF_ASSERT(A.n==B.n, "different matrix sizes (A,B).");
	if (A.num_offdiags==B.num_offdiags) {
		for (uint ii=0; ii<B.num_nonzeros; ii++) {
			B.data[ii] += a * A.data[ii];
		}
	} else {
		for (int ii=1; ii<=B.n; ii++) {
			for (int jj=negf_math::max(1,ii-B.KL); jj<=negf_math::min(B.n,ii+B.KU); jj++) { // ii,jj are 1-based
				B(ii,jj) += a*A(ii,jj);
			}
		}
	}
}

// B += a*A
inline void add (const BMatc & A, const  cplx & a, Matc & B) {
	NEGF_ASSERT(A.num_rows()==B.num_rows() && A.num_cols()==B.num_cols(), "different matrix sizes (A,B).");
	for (int ii=1; ii<=A.n; ii++) {
		for (int jj=negf_math::max(1,ii-A.KL); jj<=negf_math::min(A.n,ii+A.KU); jj++) { // ii,jj are 1-based
			B(ii,jj) += a*A(ii,jj);
		}
	}
}

// B += a*A, far-offdiagonal entries are not added
inline void add (const  Matc & A, const  cplx & a, BMatc & B) {
	NEGF_ASSERT(A.num_rows()==B.num_rows() && A.num_cols()==B.num_cols(), "different matrix sizes.");
	for (int ii=1; ii<=B.n; ii++) {
		for (int jj=negf_math::max(1,ii-B.KL); jj<=negf_math::min(B.n,ii+B.KU); jj++) { // ii,jj are 1-based
			B(ii,jj) += a*A(ii,jj);
		}
	}
}

// C = A - B, far-offdiagonal entries are discarded
inline void sub (const  Matc & A, const  Matc & B, BMatc & C) {
	NEGF_ASSERT(A.num_rows()==B.num_rows() && A.num_cols()==B.num_cols(), "different matrix sizes (A,B).");
	NEGF_ASSERT(A.num_rows()==C.num_rows() && A.num_cols()==C.num_cols(), "different matrix sizes (A,C).");
	for (int ii=1; ii<=C.n; ii++) {
		for (int jj=negf_math::max(1,ii-C.KL); jj<=negf_math::min(C.n,ii+C.KU); jj++) { // ii,jj are 1-based
			C(ii,jj) = A(ii,jj) - B(ii,jj);
		}
	}
}

// C = A * B, far-offdiagonal entries are discarded
inline void mult(const  Matc & A, const  Matc & B, BMatc & C) {
	for (int ii=1; ii<=C.n; ii++) {
		for (int jj=negf_math::max(1,ii-C.KL); jj<=negf_math::min(C.n,ii+C.KU); jj++) {
			// compute C(ii,jj)
			C(ii,jj) = A(ii,1)*B(1,jj);
			for (int kk=2; kk<=C.n; kk++) {
				C(ii,jj) += A(ii,kk)*B(kk,jj);
			}
		}
	}		
}

// C = A * B, far-offdiagonal entries are discarded
inline void mult(const BMatc & A, const  Matc & B, BMatc & C) {
	for (int ii=1; ii<=C.n; ii++) {
		for (int jj=negf_math::max(1,ii-C.KL); jj<=negf_math::min(C.n,ii+C.KU); jj++) {
			// compute C(ii,jj)
			C(ii,jj) = 0.0;
			for (int kk=1; kk<=C.n; kk++) {
				if (A.stores(ii,kk)) {
					C(ii,jj) += A(ii,kk)*B(kk,jj);
				}
			}
		}
	}
}

// C = A * B
inline void mult(const  Matc & A, const BMatc & B,  Matc & C) {
	if (work10.num_rows()!=B.num_rows()) {
		work10 = Matc(B.num_rows(), B.num_rows());
	}
	work10 = B;				// full matrix of size NxNn
	mult(A, work10, C);
}

// C = A * B
inline void mult(const BMatc & A, const  Matc & B,  Matc & C) {
	if (work10.num_rows()!=A.num_rows()) {
		work10 = Matc(A.num_rows(), A.num_rows());
	}
	work10 = A;				// full matrix of size NxNn
	mult(work10, B, C);
}

// C = A * B
inline void mult(const BMatc & A, const BMatc & B,  Matc & C) {
	C = Matc(C.num_rows(),C.num_cols());
	for (int ii=1; ii<=A.n; ii++) {
		for (int kk=negf_math::max(1,ii-A.KL); kk<=negf_math::min(A.n, ii+A.KU); kk++) {
			for (int jj=negf_math::max(1,kk-B.KU); jj<=negf_math::min(B.n,kk+B.KL); jj++) {
				C(ii,jj) += A(ii,kk)*B(kk,jj);
			}
		}
	}
}

// C = A * B
inline void mult(const BMatc & A, const BMatc & B, BMatc & C) {
	NEGF_EXCEPTION("This there.");
}

// B = a * A
inline void mult(const BMatc & A, const  cplx & a, BMatc & B) {
	NEGF_ASSERT(A.n==B.n && A.lda==B.lda, "different matrix sizes (A,B).");
	for (int ii=1; ii<=B.n; ii++) {
		for (int jj=negf_math::max(1,ii-B.KL); jj<=negf_math::min(B.n,ii+B.KU); jj++) {
			B(ii,jj) = a*A(ii,jj);
		}
	}		
}

// B = a * A
inline void mult(const BMatc & A, const  cplx & a,  Matc & B) {
	NEGF_ASSERT(A.num_rows()==B.num_rows() && A.num_cols()==B.num_cols(), "different matrix sizes.");
	//B = Matc(A.num_rows(), A.num_cols());
	for (int ii=1; ii<=A.n; ii++) {
		for (int jj=1; jj<negf_math::max(1,ii-A.KL); jj++) {
			B(ii,jj) = 0.0;
		}
		for (int jj=negf_math::max(1,ii-A.KL); jj<=negf_math::min(A.n,ii+A.KU); jj++) {
			cplx res = a*A(ii,jj);
			B(ii,jj) = res;
		}
		for (int jj=negf_math::min(A.n, ii+A.KU)+1; jj<=A.n; jj++) {
			B(ii,jj) = 0.0;
		}
	}
}

// B = a * A, far-offdiagonal entries are discarded
inline void mult(const  Matc & A, const  cplx & a, BMatc & B) {
	NEGF_ASSERT(A.num_rows()==B.num_rows() && A.num_cols()==B.num_cols(), "different matrix sizes.");
	for (int ii=1; ii<=B.n; ii++) {
		for (int jj=negf_math::max(1,ii-B.KL); jj<=negf_math::min(B.n,ii+B.KU); jj++) {
			B(ii,jj) = a*A(ii,jj);
		}
	}	
}

// B = A^T
inline void trans(const BMatc & A, BMatc & B) {
	NEGF_ASSERT(A.n==B.n && A.lda==B.lda, "different matrix sizes (A,B).");
	for (int ii=1; ii<=B.n; ii++) {
		for (int jj=negf_math::max(1,ii-B.KL); jj<=negf_math::min(B.n,ii+B.KU); jj++) {
			B(ii,jj) = A(jj,ii);
		}
	}
}

// B = A^+
inline void conjtrans(const BMatc & A, BMatc & B) {
	NEGF_ASSERT(A.n==B.n && A.lda==B.lda, "different matrix sizes (A,B).");
	cplx iu(0.0, 1.0);
	for (int ii=1; ii<=B.n; ii++) {
		for (int jj=negf_math::max(1,ii-B.KL); jj<=negf_math::min(B.n,ii+B.KU); jj++) {
			B(ii,jj) = A(jj,ii).real() - iu * A(jj,ii).imag();
		}
	}
}

// calculate diagonal of A*B
inline void get_diag(const BMatc & A, const Matc & B, vector<cplx> & diag) {
	diag.assign(A.n, 0.0);
	for (int ii=1; ii <= A.n; ii++) {
		for (int jj=negf_math::max(1,ii-A.KL); jj<=negf_math::min(A.n,ii+A.KU); jj++) {
			diag[ii-1] += A(ii,jj) * B(jj,ii);
		}
	}
}

// calculate diagonal of A*B
inline void get_diag(const Matc & A, const BMatc & B, vector<cplx> & diag) {
	diag.assign(B.n, 0.0);
	for (int ii=1; ii <= B.n; ii++) {
		for (int jj=negf_math::max(1,ii-B.KU); jj<=negf_math::min(B.n,ii+B.KL); jj++) {
			diag[ii-1] += A(ii,jj) * B(jj,ii);
		}
	}
}

// calculate diagonal of A*B
inline void get_diag(const BMatc & A, const BMatc & B, vector<cplx> & diag) {
	NEGF_ASSERT(A.n==B.n, "expected same matrix sizes.");
	diag.assign(B.n, 0.0);
	for (int ii=1; ii <= B.n; ii++) {
		for (int jj=negf_math::max(1,ii-B.KU); jj<=negf_math::min(B.n, ii+B.KL); jj++) {
			diag[ii-1] += A(ii,jj) * B(jj,ii);
		}
	}
}

#endif
