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
#ifdef FLENS
	
// will be Nn*Nn global matrices used as workhorses, saves reallocation of memory each time
// defined in MatFlens.cpp
extern GEMatrix work1;
extern GEMatrix work2;
extern GEMatrix work3;
	
inline const uint Matc::num_rows()  const { return this->flens_matrix.numRows(); }
inline const uint Matc::num_cols()  const { return this->flens_matrix.numCols(); }
inline const uint Matc::first_row() const { return this->flens_matrix.firstRow(); }
inline const uint Matc::last_row()  const { return this->flens_matrix.lastRow(); }
inline const uint Matc::first_col() const { return this->flens_matrix.firstCol(); }
inline const uint Matc::last_col()  const { return this->flens_matrix.lastCol(); }
inline const uint Matd::num_rows()  const { return this->flens_matrix.numRows(); }
inline const uint Matd::num_cols()  const { return this->flens_matrix.numCols(); }
inline const uint Matd::first_row() const { return this->flens_matrix.firstRow(); }
inline const uint Matd::last_row()  const { return this->flens_matrix.lastRow(); }
inline const uint Matd::first_col() const { return this->flens_matrix.firstCol(); }
inline const uint Matd::last_col()  const { return this->flens_matrix.lastCol(); }

inline const cplx   & Matc::operator()(const int & i, const int & j) const {return this->flens_matrix(i,j); }
inline       cplx   & Matc::operator()(const int & i, const int & j)       {return this->flens_matrix(i,j); }
inline const double & Matd::operator()(const int & i, const int & j) const {return this->flens_matrix(i,j); }
inline       double & Matd::operator()(const int & i, const int & j)       {return this->flens_matrix(i,j); }

inline Matc & Matc::operator*=(const cplx & alpha)   { this->flens_matrix *= alpha; return *this; }
inline Matc & Matc::operator+=(const cplx & alpha)   { this->flens_matrix += alpha; return *this; }
inline Matc & Matc::operator-=(const cplx & alpha)   { this->flens_matrix -= alpha; return *this; }
inline Matc & Matc::operator+=(const Matc & A)   	 { this->flens_matrix += A.flens_matrix; return *this; }
inline Matc & Matc::operator-=(const Matc & A)   	 { this->flens_matrix -= A.flens_matrix; return *this; }
inline Matc & Matc::operator= (const Matc & A)   	 { this->flens_matrix =  A.flens_matrix; return *this; }
inline Matd & Matd::operator*=(const double & alpha) { this->flens_matrix *= alpha; return *this; }
inline Matd & Matd::operator+=(const double & alpha) { this->flens_matrix += alpha; return *this; }
inline Matd & Matd::operator-=(const double & alpha) { this->flens_matrix -= alpha; return *this; }
inline Matd & Matd::operator+=(const Matd & A)   	 { this->flens_matrix += A.flens_matrix; return *this; }
inline Matd & Matd::operator-=(const Matd & A)   	 { this->flens_matrix -= A.flens_matrix; return *this; }
inline Matd & Matd::operator=(const Matd & A)   	 { this->flens_matrix =  A.flens_matrix; return *this; }

// COMPLEX OPERATIONS
inline void add (const Matc & A, const Matc & B, Matc & C)	{ C.flens_matrix = A.flens_matrix + B.flens_matrix; } // C = A + B
inline void add (const Matc & A, const cplx & a, Matc & B)  { B.flens_matrix += a*A.flens_matrix; }				  // B += a*A
inline void sub (const Matc & A, const Matc & B, Matc & C)	{ C.flens_matrix = A.flens_matrix - B.flens_matrix; } // C = A - B
inline void mult(const Matc & A, const Matc & B, Matc & C)	{ C.flens_matrix = A.flens_matrix * B.flens_matrix; } // C = A * B
inline void mult(const Matc & A, const cplx & a, Matc & B)	{ B.flens_matrix = a*A.flens_matrix; } 				  // B = a * A
inline void mult(Matc & A, const cplx & alpha) 				{ A.flens_matrix *= alpha; }						  // A *= alpha
inline void trans(const Matc & A, Matc & B)					{ B.flens_matrix = transpose(A.flens_matrix); }		  // B = A^T
inline void conjtrans(const Matc & A, Matc & B)				{ B.flens_matrix = conjugateTranspose(A.flens_matrix); } // B = A^+
	
// calculate diagonal of A*B
inline void get_diag(const Matc & A, const Matc & B, vector<cplx> & diag) {
	diag.assign(A.num_rows(), 0.0);
	for (uint ii=1; ii <= A.num_rows(); ii++) {
		for (uint jj=1; jj<=A.num_rows(); jj++) {
			diag[ii-1] += A(ii,jj)*B(jj,ii);
		}
	}
}

// DOUBLE OPERATIONS 	
inline void add (const Matd & A, const Matd & B, Matd & C)	{ C.flens_matrix = A.flens_matrix + B.flens_matrix; }	// C = A + B
inline void sub (const Matd & A, const Matd & B, Matd & C)	{ C.flens_matrix = A.flens_matrix - B.flens_matrix; }	// C = A - B
inline void mult(const Matd & A, const Matd & B, Matd & C)	{ C.flens_matrix = A.flens_matrix * B.flens_matrix; }	// C = A * B
inline void mult(const Matd & A, const double&a, Matd & B)	{ B.flens_matrix = a*A.flens_matrix; } 				  	// B = a * A
inline void mult(Matd & A, const double & alpha) 			{ A.flens_matrix *= alpha; }							// A *= alpha
inline void trans(const Matd & A, Matd & B) 				{ B.flens_matrix = transpose(A.flens_matrix); }			// B = A^T

inline const double & Vecd::operator()(const int & i) const { return this->flens_vector(i); }
inline       double & Vecd::operator()(const int & i)       { return this->flens_vector(i); }
inline const uint     Vecd::length()				  const { return this->flens_vector.length(); }

#endif	// ifdef FLENS
