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
#ifndef MAT_H_NEGF
#define MAT_H_NEGF

// see all.h for preprocessor defines
#include "all.h"

#ifdef FLENS
#include <flens/flens.h>
typedef flens::GeMatrix<flens::FullStorage<std::complex<double>, flens::ColMajor> > GEMatrix; 
typedef flens::GeMatrix<flens::FullStorage<double, flens::ColMajor> > 				DGEMatrix; 
typedef flens::DenseVector<flens::Array<std::complex<double> > > 					DEVector; 
typedef flens::DenseVector<flens::Array<double> > 									DDEVector; 
#endif

using namespace std;

// a+=b is the same as a.operator+=(b)

namespace negf {
	
	class BMatc;
	
	/** Wrapper for COMPLEX matrix storage used and matrix operations */
	class Matc {
	public:
		Matc(); 				 //!< default constructor --> 1x1-matrix
		Matc(int n, int m);		 //!< n=rows, m=columns
		~Matc();
		Matc(const Matc & copy); //!< copy constructor
		
		/** Get the Nn*Nn-block-matrix part of the matrix corresponding to positions xx,yy and overwrite it with something */
		void fill_block(uint xx, uint yy, const Matc & block_matrix, const uint & Nx);
		
		/** Get the Nn*Nn-block-matrix part of the matrix corresponding to positions xx,yy and overwrite it with
		  * the Nn*Nn-matrix corresponding to position x2,y2 in the matrix "whole_matrix" */
		void fill_block(uint xx, uint yy, const Matc & whole_matrix, uint x2, uint y2, const uint & Nx);
		
		/** Get the Nn*Nn-block-matrix part of the matrix corresponding to positions xx,yy OUT OF NVERT IN TOTAL and overwrite it with
		  * the Nn*Nn-matrix corresponding to position x2,y2 in the matrix "whole_matrix" */
		//void fill_block(uint xx, uint yy, const Matc & whole_matrix, uint x2, uint y2, uint Nvert);
		
		/** Get the Nn*Nn-block-matrix part of the matrix corresponding to positions xx,yy and subtract
		  * the Nn*Nn-matrix corresponding to position x2,y2 in the matrix "whole_matrix" */
		void subtract_block(uint xx, uint yy, const Matc & whole_matrix, uint x2, uint y2, const uint & Nx);
		
		/** Multiply the Nn*Nn-block corresponding to xx,yy by some scalar */
		void multiply_block(uint xx, uint yy, const cplx & alpha, const uint & Nx);
		
		/** Return a Nn*Nn-block (copy, not reference!) corresponding to the place xx,yy */
		void get_block(uint xx, uint yy, Matc & result, const uint & Nx) const;
		
		/** Return a submatrix (copy, not reference!) */
		void get_submatrix(uint x1, uint x2, uint y1, uint y2, Matc & result) const;
		
		inline const uint num_rows()  const;
		inline const uint num_cols()  const;
		inline const uint first_row() const;
		inline const uint last_row()  const;
		inline const uint first_col() const;
		inline const uint last_col()  const;
		
		inline const cplx & operator()(const int & i, const int & j) const;
		inline       cplx & operator()(const int & i, const int & j);
		inline 		 Matc & operator*=(const cplx & alpha);
		inline 		 Matc & operator+=(const cplx & alpha);
		inline 		 Matc & operator-=(const cplx & alpha);
		inline 		 Matc & operator+=(const Matc & A);
		inline 		 Matc & operator-=(const Matc & A);
		inline 		 Matc & operator= (const Matc & A);
		inline 		 Matc & operator+=(const BMatc & A);
		inline 		 Matc & operator-=(const BMatc & A);
		inline 		 Matc & operator= (const BMatc & A);
		
		int num_offdiags;						//!< will always be n-1. necessary for do_compress, MPI::recv(vector<..> data)
		template<int offdiags>
		static Matc & create(int n_, int m_);   //!< need it for NEGF constructor. offdiags is irrelevant
#ifdef FLENS
		Matc(const GEMatrix & flens_matrix_);
		Matc(const GEMatrix::ConstView flens_matrix_);

		GEMatrix flens_matrix;
#endif

	protected:
		//Matc(const Matc & copy) {} // protected implementation of copy constructor prevents erroneous calling of default c.c.
		Matc(const BMatc & copy) {}
	};	
	
	/** Wrapper for DOUBLE matrix storage used and matrix operations */
	class Matd {
	public:
		Matd(); 			//!< 1x1-matrix
		Matd(int m, int n);	//!< m=rows, n=columns
		~Matd();
#ifdef CSCS
		Matd(const Matd & copy); //!< for some reason need public copy constructor when compiling on palu.cscs.ch
#endif
				
		/** get the Nn*Nn-block-matrix part of the matrix corresponding to positions xx,yy and overwrite it with something */
		void fill_block(uint xx, uint yy, const Matd & block_matrix, const uint & Nx);
		
		inline const double & operator()(const int & i, const int & j) const;
		inline       double & operator()(const int & i, const int & j);
		inline 		 Matd & operator*=(const double & alpha);
		inline 		 Matd & operator+=(const double & alpha);
		inline 		 Matd & operator-=(const double & alpha);
		inline 		 Matd & operator+=(const Matd & A);
		inline 		 Matd & operator-=(const Matd & A);
		inline 		 Matd & operator= (const Matd & A);
		
		inline const uint num_rows()  const;
		inline const uint num_cols()  const;
		inline const uint first_row() const;
		inline const uint last_row()  const;
		inline const uint first_col() const;
		inline const uint last_col()  const;
		
#ifdef FLENS
		Matd(const DGEMatrix & flens_matrix_);
		Matd(const DGEMatrix::ConstView flens_matrix_);
		
		DGEMatrix flens_matrix;
#endif
	protected:
#ifndef CSCS
		Matd(const Matd & copy) {} //!< protected implementation of copy constructor prevents erroneous calling of default c.c.
#endif
		Matd(const Matc & copy) {}
		Matd(const BMatc & copy) {}
	};	
	
	/** Wrapper for COMPLEX BANDED matrix storage used and matrix operations */
	class BMatc {
	public:
		BMatc(int m_, int n_, int num_offdiags_);
		
		virtual ~BMatc();
		BMatc(const BMatc & copy);
				
		void     fill_block(uint xx, uint yy, const Matc  & block_matrix, const uint & Nx);
		void     fill_block(uint xx, uint yy, const Matc  & whole_matrix, uint x2, uint y2, const uint & Nx);
		void     fill_block(uint xx, uint yy, const BMatc & whole_matrix, uint x2, uint y2, const uint & Nx);
		void subtract_block(uint xx, uint yy, const Matc  & whole_matrix, uint x2, uint y2, const uint & Nx);
		void subtract_block(uint xx, uint yy, const BMatc & whole_matrix, uint x2, uint y2, const uint & Nx);
		
		void multiply_block(uint xx, uint yy, const cplx & alpha, const uint & Nx);
		void      get_block(uint xx, uint yy, Matc & result, const uint & Nx) const;  //!< copy, not reference!
		void  get_submatrix(uint x1, uint x2, uint y1, uint y2, Matc & result) const; //!< copy, not reference!
		
		inline const uint num_rows()  const;
		inline const uint num_cols()  const;
		inline const uint first_row() const;
		inline const uint last_row()  const;
		inline const uint first_col() const;
		inline const uint last_col()  const;
		
		inline const cplx & operator()(const int & i, const int & j) const;
		inline       cplx & operator()(const int & i, const int & j);
		inline const cplx & nocheck(const int & i, const int & j) const;			//!< fast version of operator()
		inline 		 BMatc & operator*=(const cplx & alpha);
		inline 		 BMatc & operator+=(const cplx & alpha);
		inline 		 BMatc & operator-=(const cplx & alpha);
		inline 		 BMatc & operator+=(const BMatc & A);
		inline 		 BMatc & operator+=(const Matc  & A);	//!< skips all far-offdiagonal entries!
		inline 		 BMatc & operator-=(const BMatc & A);
		inline 		 BMatc & operator-=(const Matc  & A);	//!< skips all far-offdiagonal entries!
		inline 		 BMatc & operator= (const BMatc & A);
		inline 		 BMatc & operator= (const Matc  & A);	//!< skips all far-offdiagonal entries!
		
#ifdef BANDED
		uint 		 num_offdiags;
		cplx 		 zero;
		
		// LAPACK stuff
		int 		 n;
		int 		 KL;
		int 		 KU;
		int 		 nrhs;
		int 		 lda;
		int 		 info;
		
		cplx * 		 data;
		uint 		 num_nonzeros;
		
		inline bool  stores(uint i, uint j) const;
		inline int   get_ib(const uint & i, const uint & j) const;
#endif
		
	protected:
		BMatc(const Matc & copy) { NEGF_EXCEPTION("Do not call!"); }
		BMatc(const Matd & copy) { NEGF_EXCEPTION("Do not call!");}
		
	};	

	// get the index of the matrix entry corresponding to position xx, band nn
	// assumption: input AND result are 1-based
	// it's important to return an int because get_mat_idx(x,m)-get_mat_idx(y,n) will be used!
	inline int get_mat_idx(const uint & xx, const uint & nn, const uint & Nx) {
		return (nn-1)*Nx+xx;
	}
	
	// get the index of the matrix entry corresponding to position xx, band nn, given that the matrix has Nvert "vertices"
	// so the matrix is assumed to have a size (Nn*Nvert)^2
	// Nvert could be e.g. the total number of vertices, including the internal contact vertices, for accessing the Hamiltonian
	// assumption: input AND result are 1-based
	//inline int get_special_mat_idx(const uint & Nvert, const uint & xx, const uint & nn) {
	//	return (nn-1)*Nvert+xx;
	//}
	
	// COMPLEX OPERATIONS
	inline void add (const Matc & A, const Matc & B, Matc & C); //!< C = A + B
	inline void add (const Matc & A, const cplx & a, Matc & B); //!< B += a*A
	inline void sub (const Matc & A, const Matc & B, Matc & C); //!< C = A - B
	inline void mult(const Matc & A, const Matc & B, Matc & C); //!< C = A * B
	inline void mult(const Matc & A, const cplx & a, Matc & B); //!< B = a * A
	inline void mult(Matc & A, const cplx & alpha);				//!< A *= alpha
	inline void trans(const Matc & A, Matc & B);		  		//!< B = A^T
	inline void conjtrans(const Matc & A, Matc & B); 			//!< B = A^+
	void 		invert(Matc & A);								//!< A --> inv(A)
	
	/** calculate the diagonal of A*B */
	inline void get_diag(const Matc & A, const Matc & B, vector<cplx> & diag);
	
	/** get Nn*Nn-submatrices at position xa,ya (A) and xb,yb (B), multiply (A*B) and store in the Nn*Nn-matrix C */
	void multiply_blocks(const Matc & A, uint xa, uint ya, const Matc & B, uint xb, uint yb, Matc & C, const uint & Nx);
		
	/** get Nn*Nn-submatrices at position xa,ya (A) and xb,yb (B), multiply (A*B) and store in the Nn*Nn-submatrix position xc,yc in C */
	void multiply_blocks(const Matc & A, uint xa, uint ya, const Matc & B, uint xb, uint yb, Matc & C, uint xc, uint yc, const uint & Nx);
	
	// DOUBLE OPERATIONS 	
	inline void add (const Matd & A, const Matd & B, Matd & C);	//!< C = A + B
	inline void sub (const Matd & A, const Matd & B, Matd & C);	//!< C = A - B
	inline void mult(const Matd & A, const Matd & B, Matd & C);	//!< C = A * B
	inline void mult(const Matd & A, const double&a, Matd & B); //!< B = a * A
	inline void mult(Matd & A, const double & alpha);			//!< A *= alpha
	inline void trans(const Matd & A, Matd & B);				//!< B = A^T
	void 		invert(Matd & A);								//!< A --> inv(A)

	// COMPLEX OPERATIONS w/ BANDED MATRICES
	inline void add (const BMatc & A, const BMatc & B, BMatc & C); //!< C = A + B
	inline void add (const BMatc & A, const  cplx & a, BMatc & B); //!< B += a*A
	inline void add (const BMatc & A, const  cplx & a,  Matc & B); //!< B += a*A
	inline void add (const  Matc & A, const  cplx & a, BMatc & B); //!< B += a*A, far-offdiagonal entries are not added
	inline void sub (const BMatc & A, const BMatc & B, BMatc & C); //!< C = A - B
	inline void sub (const  Matc & A, const  Matc & B, BMatc & C); //!< C = A - B, far-offdiagonal entries are discarded
	inline void mult(const  Matc & A, const  Matc & B, BMatc & C); //!< C = A * B, far-offdiagonal entries are discarded
	inline void mult(const BMatc & A, const  Matc & B, BMatc & C); //!< C = A * B, far-offdiagonal entries are discarded
	inline void mult(const  Matc & A, const BMatc & B,  Matc & C); //!< C = A * B
	inline void mult(const BMatc & A, const  Matc & B,  Matc & C); //!< C = A * B
	inline void mult(const BMatc & A, const BMatc & B,  Matc & C); //!< C = A * B
	inline void mult(const BMatc & A, const BMatc & B, BMatc & C); //!< C = A * B
	inline void mult(const BMatc & A, const  cplx & a, BMatc & B); //!< B = a * A
	inline void mult(const BMatc & A, const  cplx & a,  Matc & B); //!< B = a * A
	inline void mult(const  Matc & A, const  cplx & a, BMatc & B); //!< B = a * A, far-offdiagonal entries are discarded
	inline void trans(const BMatc & A, BMatc & B);		  		   //!< B = A^T
	inline void conjtrans(const BMatc & A, BMatc & B); 			   //!< B = A^+
	inline void get_diag(const  Matc & A, const BMatc & B, vector<cplx> & diag);
	inline void get_diag(const BMatc & A, const  Matc & B, vector<cplx> & diag);
	inline void get_diag(const BMatc & A, const BMatc & B, vector<cplx> & diag);
	
	/* Wrapper for DOUBLE VECTORS */
	class Vecd {
	public:
		Vecd(int n);
		~Vecd();
		
		inline const double & operator()(const int & i) const;
		inline       double & operator()(const int & i);
		
		inline const uint length() const; 
		
#ifdef FLENS
		Vecd(const DDEVector & flens_vector_);
		DDEVector flens_vector;
#endif
	protected:
	};
	
	void solve_linear_problem(/*const*/ Matd & A, Vecd & rhs); 	//!< rhs is overwritten w/ solution
	void mult(const Matd & A, const Vecd & v, Vecd & rhs); 		//!< rhs = A*v;
	
	// =====================================================================
	// IMPLEMENTATIONS OF INLINE FUNCTIONS
	// note that these are independent on how the Nn bands are stored
	// =====================================================================
#ifdef FLENS
	#include "MatFlens.inl"
#endif
#ifdef BANDED
	#include "MatBanded.inl"
#endif
		
} // end namespace 
		
#endif /*MAT_H_NEGF*/
