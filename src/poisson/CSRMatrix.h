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
#ifndef CSRMATRIX_H_NEGF
#define CSRMATRIX_H_NEGF

#include "all.h"

#define MAX_INT_SIZE 2147483637

namespace negf {
		
	enum MatrixType { symmetric_matrix, nonsymmetric_matrix };
	
	/** Home-made sparse matrix class. Was originally templated by Ratko, but only the double version is used in AQUA. <BR>
	 *  The matrix is stored in the compressed sparse row (CSR) format. In this format, three vectors describe the matrix: <BR>
	 *  icol[ii] gives the position in prow/nonzeros of the first nonzero of line ii. Hence the number of nonzeros in line ii is icol[ii+1]-icol[ii]. <BR>
	 *  prow[ii] gives the column number of a nonzero. <BR>
	 *  nonzeros[ii] gives the value of that matrix entry. <BR>
	 *  Note that the length of the vectors is n+1, num_nonzeros and num_nonzeros. <BR>
	 *  There is an additional complication with the index base. Most linear solvers use fortran where indices are 1-based. <BR>
	 *  The index base is controlled by fidx. */
	template<class T>
	class CSRMatrix {
	
		friend class CSRMatrix<cplx>;
		friend class CSRMatrix<double>;

	public:
		CSRMatrix(uint size) throw (Exception *);							//!< creates symmetric CSR matrix
		CSRMatrix(uint size, MatrixType matrix_type) throw (Exception *);	//!< creates symmetric or nonsymmetric CSR matrix
		CSRMatrix(const CSRMatrix<T> &copy) throw (Exception *);			//!< 1:1 copy
		CSRMatrix(const CSRMatrix<T> &copy, MatrixType matrix_type) throw (Exception *); //!< can create nonsymmetric from symmetric matrix
		CSRMatrix(const char* filename) throw (Exception *);				//!< reads from file
		~CSRMatrix() throw (Exception *);

		// ---------------------------------------------------------
		// structure
		// ---------------------------------------------------------
		void clear_all();							//!< deletes icol, prow, nonzeros, num_nz
		void clear_but_keep_structure();			//!< sets every entry of nonzeros to 0
		void announce(const int &ii,const int &jj); //!< announces a nonzero entry; zero-based ii,jj!
		void set_structure();						//!< finalize the matrix structure. afterwards: set, add, get
		void print_sparse_structure() const;		//!< prints out icol, prow, nonzeros, num_nz
		bool check_singular_structure() const; 		//!< check if a row/column contains no nonzero entry
		bool check_matrix() const;					//!< check order of prows, order and bounds of icols
		bool symmetric() const;						//!< true if is_symmetric=true
		bool symmetric_by_value() const;			//!< true also if values are symmetric
		bool structure_is_set() const { return this->structure_set; }

		// ---------------------------------------------------------
		// access
		// ---------------------------------------------------------
		inline void  set(int ii, int jj, T val);	//!< 0-based ii,jj!
		inline void  add(int ii, int jj, T val);	//!< 0-based ii,jj!
		inline T     get(int ii, int jj) const;		//!< 0-based ii,jj!
		uint 		 get_size()        const { return size; }	//!< get the matrix size (this class has only square matrices)

		uint		get_num_nonzeros() const { return this->num_nz; }		//!< get the number of nonzero matrix entries
		T*   		get_nonzeros()     const { return this->nonzeros; }		//!< get the nonzeros array (CSR formatting)
		int* 		get_icol()         const { return this->icol; }			//!< get the icol array (CSR formatting). 1-based indices by default!
		int* 		get_prow()         const { return this->prow; }			//!< get the prow array (CSR formatting). 1-based indices by default!

		T 			get_nonzeros_entry(uint idx) const throw (Exception *); //!< Get value at a certain position in the nonzeros array
		int 		get_icol_entry(uint idx) const throw (Exception *);		//!< Get value at a certain position in the icol array
		int 		get_prow_entry(uint idx) const throw (Exception *);		//!< Get value at a certain position in the prow array

		int 		get_index_base() const { return fidx; }	// 0 or 1
	
		// ---------------------------------------------------------
		// I/O
		// ---------------------------------------------------------
		void 			print_stats() const;						// rather useless stats about avg diag/offdiag values
		T*   			create_full_matrix() const;
		CSRMatrix<cplx>* create_complex_matrix() const;				// converts this matrix into a complex one
		void 			save_to_file(const char* filename) const;	//!< in text format
		void 			read_from_file(const char* filename);		//!< file should have been created by save_to_file() or have same format

		// ---------------------------------------------------------
		// operations
		// ---------------------------------------------------------
		template<class B>
		void mult_vec(const B* in, B* out) const;					// calculates out=A*in, A = this matrix

		template<class B>
		void mult_vec_multiple(const B* in, B* out, int kpn) const;
	
		template<class B>
		void ILU_precondition(B* x, const B*b) const;
		void ILU_factorize();
		bool factorized() const { return this->is_factorized; }

	protected:

		class announce_list {
			public:
				announce_list() { col = -1; next = 0; prev = 0; }
				int insert(const int &col_);
				int col;
				announce_list* next;
				announce_list* prev;
		};
		void clear_announced();  //!< Destroy the announced structure
		void init_announced();	//!< Creates the announce structure on the diagonal
		void copy_raw(const CSRMatrix<T> &copy);
		void init();			//!< Init all member variables except size to zero
		T    internal_conj(T z) const; //!< <ss> dunno what this is

		uint size;				//!< the size of the square matrix

		int* prow;				//!< array of length size + 1, gives position of row ii in nonzeros
		int* icol;  			//!< array of length nonzeros, gives column index of element ii at nonzeros(ii)
		T*   nonzeros;			//!< nonzeros ...
		int  num_nz;    		//!< size of nonzeros
		int  fidx; 				//!< fortran index -> if = 1, prow/icol use fortran indices

		bool structure_set;		//!< initially false, will be set to true by set_structure()
		bool is_symmetric;		//!< always false in ANGEL
		bool is_factorized;		//!< irrelevant in ANGEL
		
		announce_list** announced;

	};


/** create SYMMETRIC compressed sparse row matrix
 * 
 * to allow switching between CSCMatrix used by arpack++ and symmetric CSR matrices 
 * used by JDQZ and ARPACK + Pardiso, we have to provide the class with a constructor
 * which only takes the size as input value. therefore we default here to a symmetric
 * matrix to allow simple switching via redefinition of type GMatrix.
 * @param size the size of the symmetric matrix
 */
template<class T>
CSRMatrix<T>::CSRMatrix(unsigned int size_)  throw (Exception *)
{STACK_TRACE(
	if(size_ == 0) NEGF_EXCEPTION("invalid matrix size(0) set");
	this->init(); // sets all to zero except size
	this->is_symmetric = true;
	this->size         = size_;
	this->init_announced();
	this->clear_all();	
);}


/** create compressed sparse row matrix
 * 
 * @param size_     size of matrix
 * @param symmetric boolean to indicate whether matrix is symmetric or not
 */
template<class T>
CSRMatrix<T>::CSRMatrix(unsigned int size_, MatrixType matrix_type) throw (Exception *): 
	size(size_) 
{STACK_TRACE(	
	if(size_ == 0) NEGF_EXCEPTION("invalid matrix size(0) set");
	this->init(); // sets all to zero except size
	this->is_symmetric = (matrix_type == symmetric_matrix ? true : false);
	this->init_announced();
	this->clear_all();
);}


/** create matrix and read data from file
 * 
 * @param matrix file where matrix was saved with @see CSRMatrix<T>::save_to_file
 */
template<class T>
CSRMatrix<T>::CSRMatrix(const char* filename) throw (Exception *)
{STACK_TRACE(
	this->init();
	this->clear_all();		
	this->read_from_file(filename);
);}


/** copy constructor to create a deep copy of a CSRMatrix. calls copy_raw.
 * 
 * only matrices with initialized structures may be copied. if
 * copy has no structure set, we throw an exception
 */
template<class T>
CSRMatrix<T>::CSRMatrix(const CSRMatrix<T> &copy) throw (Exception *)
{STACK_TRACE(
	this->init();
	this->copy_raw(copy);
);}


/** copy constructor to create a deep copy of a CSRMatrix, but allowing to create a nonsymmetric matrix from a symmetric one
 * 
 * only matrices with initialized structures may be copied. if
 * copy has no structure set, we throw an exception
 * also if symmetric == true && copy.symmetric == false, an exception is thrown
 * 
 * @param copy matrix to copy
 * @param symmetric boolean to indicate if new matrix should be reagarded as symmetric or not
 */
template<class T>
CSRMatrix<T>::CSRMatrix(const CSRMatrix<T> &copy, MatrixType matrix_type) throw (Exception *)
{STACK_TRACE(
	this->init();
	if(copy.icol && copy.prow && copy.nonzeros) {
		if(matrix_type == symmetric_matrix) {
			// symmetric to symmetric -> copy raw
			if(copy.is_symmetric) {
				this->copy_raw(copy);	
			} else {
				NEGF_EXCEPTION("cannot create symmetric matrix from unsymmetric matrix");	
			}
		} else {
			// nonsymmetric to nonsymmetric -> copy raw
			if(!copy.is_symmetric) {
				this->copy_raw(copy);
			} else {				
				this->prow 	 	    = 0;
				this->icol  	 	= 0;
				this->nonzeros 	    = 0;		
				this->size          = copy.size;
				this->structure_set = false;	
				this->is_factorized = false;		
				this->is_symmetric  = false; // structure is not symmetric
				this->init_announced();
		
				this->clear_all();				
				int until, from, colidx;
				T tmp;	
				// announce		
				for(int ii = 0; ii < (signed)copy.size; ii++) {
					until    = copy.prow[ii + 1] - copy.fidx;
					from     = copy.prow[ii] - copy.fidx;
					this->announce(ii,ii);
					for(int jj = from + 1; jj < until; jj++) {
						colidx = copy.icol[jj] - copy.fidx;
						this->announce(ii,colidx);
						this->announce(colidx,ii);
					}					
				}		
				// create matrix
				this->set_structure();
				// fill in
				for(int ii = 0; ii < (signed)copy.size; ii++) {
					until    = copy.prow[ii + 1] - copy.fidx;
					from     = copy.prow[ii] - copy.fidx;
					tmp      = copy.nonzeros[from];
					this->set(ii,ii, tmp);	
					NEGF_ASSERT((copy.icol[from] - copy.fidx) == ii, "(copy.icol[from] - copy.fidx) == ii");			
					for(int jj = from + 1; jj < until; jj++) {
						tmp = copy.nonzeros[jj];
						colidx = copy.icol[jj] - copy.fidx;				
						this->set(ii,colidx, tmp);
						// set cplx conjugate
						this->set(colidx,ii, internal_conj(tmp));
					}					
				}		
			}
		}
	} else {
		NEGF_EXCEPTION("cannot create CSRMatrix from uninitialized one");	
	}									
);}


/** destructor
 * 
 * calls clear_all and clear_announced to kill allocated data in heap
 */
template<class T>
CSRMatrix<T>::~CSRMatrix()  throw (Exception *)
{STACK_TRACE(
	this->clear_all();
	this->clear_announced();
);}


/** copy function
 * 
 * performs a deep copy of sparse data
 * @param copy CSRMatrix to copy
 */
template<class T>
void CSRMatrix<T>::copy_raw(const CSRMatrix<T> &copy) 
{STACK_TRACE(
	// abort if object already initialized
	if(this->announced || this->prow || this->icol || this->nonzeros) {
		NEGF_EXCEPTION("copy_raw may not be used on initialized objects");	
	}
	
	this->size 	   		= copy.size;
	this->fidx          = copy.fidx;
	this->num_nz   		= copy.num_nz;
	this->prow 	   		= 0;
 	this->icol 	   		= 0;
	this->nonzeros 		= 0;	
	this->announced 	= 0;	
	this->structure_set = copy.structure_set;
	this->is_symmetric  = copy.is_symmetric;
	this->is_factorized = copy.is_factorized;
			
	// check if we copy the matrix itself
	if(copy.icol && copy.prow && copy.nonzeros) {
		logmsg->emit(LOG_INFO_L2, "copy values");
		// create copy of full csc
		icol     = new int[num_nz];
		prow     = new int[size + 1];				
		nonzeros = new T[num_nz];
		if(!icol || !prow || ! nonzeros) {
			NEGF_EXCEPTION("cannot allocate memory for matrix");	
		}
		for(int ii = 0; ii < num_nz; ii++) {
			this->icol[ii]     = copy.icol[ii];	
			this->nonzeros[ii] = copy.nonzeros[ii];
		}
		for(unsigned int ii = 0; ii <= size; ii++) {
			this->prow[ii] = copy.prow[ii];	
		}
	} else {
		// structure not set, throw exception
		NEGF_EXCEPTION("cannot copy matrix with empty structure");
	}
);}


template<class T>
CSRMatrix<cplx>* CSRMatrix<T>::create_complex_matrix() const 
{STACK_TRACE(
	CSRMatrix<cplx>* ret = new CSRMatrix<cplx>(this->get_size());
	ret->clear_announced();
	ret->size 	   	   = this->size;
	ret->fidx          = this->fidx;
	ret->num_nz   	   = this->num_nz;
	ret->prow 	   	   = new int[this->size + 1];
 	ret->icol 	   	   = new int[this->num_nz];
	ret->nonzeros 	   = new cplx[this->num_nz];	
	ret->announced 	   = 0;	
	ret->structure_set = this->structure_set;
	ret->is_symmetric  = this->is_symmetric;
	ret->is_factorized = this->is_factorized;
	if(!ret->icol || !ret->prow || !ret->nonzeros) {
		NEGF_EXCEPTION("cannot allocate memory for matrix");	
	}
	for(int ii = 0; ii < this->num_nz; ii++) {
		ret->icol[ii]     = this->icol[ii];	
		ret->nonzeros[ii] = this->nonzeros[ii];
	}
	for(unsigned int ii = 0; ii <= this->size; ii++) {
		ret->prow[ii] = this->prow[ii];	
	}
	return ret;	
);}


/** init all member variables except size to zero 
 */
template<class T>
void CSRMatrix<T>::init() 
{STACK_TRACE(
	this->fidx          = 1; // prow/icol have fortran indices (for paradiso)
//	this->fidx          = 0;
	this->prow 	 	    = 0;
	this->icol  	 	= 0;
	this->nonzeros 	    = 0;		
	this->structure_set = false;	
	this->is_factorized = false;
	this->is_symmetric  = false;
	this->announced     = 0;
	this->num_nz        = 0;	
);}


/** set value at (ii,jj)
 * 
 * if matrix was initialized as symmetric, only the upper triangle matrix may be accessed. if lower triangle is accessed, an exception will be thrown
 * also, setting of an nonannouced position results in an exception
 * 
 * @param ii row index (0-based)
 * @param jj col index (0-based)
 */
template<class T>
inline void CSRMatrix<T>::set(int ii, int jj, T val) 
{STACK_TRACE(
	NEGF_ASSERT(ii >= 0 && ii < (signed)size && jj >= 0 && jj < (signed)size, "matrix index out of bounds in set()");

	if(this->is_symmetric && ii > jj) {
		NEGF_EXCEPTION("setting of lower triangle matrix for a symmetric matrix is forbidden");	
	}
	
	NEGF_ASSERT(prow[ii]   >= 0 + fidx && prow[jj]   <= num_nz + fidx, "prow[ii]   >= 0 + fidx && prow[jj]   <= num_nz + fidx");
	NEGF_ASSERT(prow[ii+1] >= 0 + fidx && prow[ii+1] <= num_nz + fidx, "prow[ii+1] >= 0 + fidx && prow[ii+1] <= num_nz + fidx");		
	
	for(int rr = prow[ii]; rr < prow[ii + 1]; rr++) {
		if(icol[rr - fidx] == (jj + fidx)) {
			nonzeros[rr - fidx] = val;
			return;	
		}	
	}
		
	NEGF_FEXCEPTION("tried to set value to not announced pos: %d/%d", ii, jj);
);}


/** add to value at (ii,jj)
 * 
 * if matrix was initialized as symmetric, only the upper triangle matrix may be accessed. if lower triangle is accessed, an exception will be thrown
 * also, setting of an nonannouced position results in an exception
 * 0-based indices ii,jj!
 * 
 * @param ii row index
 * @param jj col index 
 */
template<class T>
inline void CSRMatrix<T>::add(int ii, int jj, T val) 
{STACK_TRACE(
	NEGF_ASSERT(ii >= 0 && ii < (signed)size && jj >= 0 && jj < (signed)size, "matrix index out of bounds in add()");
	if(this->is_symmetric && ii > jj) {
		NEGF_EXCEPTION("setting of lower triangle matrix for a symmetric matrix is forbidden");	
	}	
	NEGF_ASSERT(prow[ii]   >= 0 + fidx && prow[jj]   <= num_nz + fidx, "prow[ii]   >= 0 + fidx && prow[jj]   <= num_nz + fidx");
	NEGF_ASSERT(prow[ii+1] >= 0 + fidx && prow[ii+1] <= num_nz + fidx, "prow[ii+1] >= 0 + fidx && prow[ii+1] <= num_nz + fidx");		
	
	for(int rr = prow[ii]; rr < prow[ii + 1]; rr++) {
		if(icol[rr - fidx] == jj + fidx) {
			nonzeros[rr - fidx] += val;		
			return;
		}	
	}
	NEGF_FEXCEPTION("tried to set value to not announced pos: %d/%d", ii, jj);
);}


/** returns value at (ii,jj)
 *
 * @param ii row index [0, size - 1], 0-based!
 * @param jj col index [0, size - 1], 0-based!
 */
template<class T>
inline T CSRMatrix<T>::get(int ii, int jj) const 
{STACK_TRACE(
	int tmp;	
	bool return_conjugate = false;
	NEGF_ASSERT(ii >= 0 && ii < (signed)size && jj >= 0 && jj < (signed)size, "matrix index out of bounds in get()");
	
	if(this->is_symmetric && ii > jj) {
		tmp = jj;
		jj  = ii;
		ii  = tmp;	
		return_conjugate = true;
	}	
	NEGF_ASSERT(prow[ii]   >= 0 + fidx && prow[jj]   <= num_nz + fidx, "prow[ii]   >= 0 + fidx && prow[jj]   <= num_nz + fidx");
	NEGF_ASSERT(prow[ii+1] >= 0 + fidx && prow[ii+1] <= num_nz + fidx, "prow[ii+1] >= 0 + fidx && prow[ii+1] <= num_nz + fidx");			
	
	for(int rr = prow[ii]; rr < prow[ii + 1]; rr++) {
		if(icol[rr - fidx] == jj + fidx) {
			if(return_conjugate) {
				return internal_conj(nonzeros[rr - fidx]);
			} else {
				return nonzeros[rr - fidx];
			}
		}	
	}	
	return (T)0;
);} 


/** returns value at a certain position in the nonzeros array */
template<class T>
T CSRMatrix<T>::get_nonzeros_entry(uint idx) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->nonzeros != NULL, "nonzeros array was not yet set up.");
	NEGF_ASSERT(idx < this->get_num_nonzeros(), "invalid index.");
	return this->nonzeros[idx];
);}

/** returns value at a certain position in the icol array */
template<class T>
int CSRMatrix<T>::get_icol_entry(uint idx) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->icol != NULL, "icol array was not yet set up.");
	NEGF_ASSERT(idx < this->get_num_nonzeros(), "invalid index.");
	return this->icol[idx];
);}

/** returns value at a certain position in the prow array */
template<class T>
int CSRMatrix<T>::get_prow_entry(uint idx) const throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(this->prow != NULL, "icol array was not yet set up.");
	NEGF_ASSERT(idx < this->get_size() + 1, "invalid index.");
	return this->prow[idx];
);}


/** delete sparse matrix data and structure
 * 
 * attention: used for destructors, does not reinit the matrix
 */
template<class T>
void CSRMatrix<T>::clear_all() 
{STACK_TRACE(
	if(icol)     { delete[] icol; icol = 0; }
	if(prow)     { delete[] prow; prow = 0; }
	if(nonzeros) { delete[] nonzeros; nonzeros = 0;}	
	num_nz 		  = 0;
	structure_set = false;
);}


/** creates the announce structure on the diagonal 
 */
template<class T>
void CSRMatrix<T>::init_announced() 
{STACK_TRACE(
	this->announced     = new announce_list*[this->size];
	for(unsigned int ii = 0; ii < this->size; ii++) {
		this->announced[ii] = new announce_list();
		this->announced[ii]->col = ii; // set the diagonal	
	}		
);}


/** destroy the announced structure
 */
template<class T>
void CSRMatrix<T>::clear_announced() 
{STACK_TRACE(
	announce_list* ptr;
	announce_list* next;
	announce_list* prev;
	
	// if announced array is initialized
	if(this->announced) {
		// loop over announced array
		for(unsigned int ii = 0; ii < size; ii++) {
			if(this->announced[ii]) {			
				// if element is set, descend tree and delete remaining list entries			
				// delete prev
				ptr  = 0; 
				prev = this->announced[ii]->prev;
				while(prev) {
					ptr  = prev;
					prev = ptr->prev;
					delete ptr;	
				}
				// delete next			
				ptr  = 0;
				next = this->announced[ii];
				while(next) {
					ptr = next;
					next = ptr->next;
					delete ptr;	
				}
			}
		}	
		delete[] this->announced;
		this->announced = 0;
	}	
);}


/** insert col_ into list tree 
 * 
 * @param col some value (here: column index)
 * @return 1 if col was inserted, 0 if col is already in tree (therefore allowing to count number of nonzeros)
 */ 
template<class T>
int CSRMatrix<T>::announce_list::insert(const int &col_) 
{STACK_TRACE(
	if(col_ == this->col) {
		return 0;
	} else if(col_ > this->col) {
		// pass to next? 
		if(this->next && this->next->col <= col_) {
			return this->next->insert(col_);	
		// no next avialble, add new at end
		} else if(!next) {
			this->next       = new announce_list();
			this->next->col  = col_;
			this->next->next = 0;
			this->next->prev = this;
			return 1;
		// value is between this and next -> insert new one
		} else {
			announce_list* ptr = new announce_list();
			ptr->col         = col_;
			ptr->next        = this->next;
			ptr->prev  		 = this; // me is prev of next
			this->next->prev = ptr;  // old next has now new next as prev instead of me
			this->next = ptr;
			return 1;
		}
	} else if(col_ < this->col) {		
		// previous available, pass to prev?
		if(this->prev && this->prev->col >= col_) {
			return this->prev->insert(col_);
		// no prev available, add new at beginning
		} else if(!this->prev) {
			announce_list* ptr = new announce_list();
			ptr->col 		   = col_;
			ptr->next          = this;
			ptr->prev          = 0;
			this->prev         = ptr;
			return 1;
		// add in between prev and this
		} else {			
			announce_list* ptr = new announce_list();
			ptr->col         = col_;
			ptr->next        = this;
			ptr->prev        = this->prev;
			this->prev->next = ptr;
			this->prev       = ptr;	
			return 1;
		} 	
	} else {
		NEGF_EXCEPTION("should not reach this point");
	}
);}


/** set values to zero
 */
template<class T>
void CSRMatrix<T>::clear_but_keep_structure() 
{STACK_TRACE(
	if(this->nonzeros) {
		for(int ii = 0; ii < num_nz; ii++) {
			this->nonzeros[ii] = (T)0.0;	
		}
	}
);}


/** announce value at position (ii,jj)
 * 
 * announce that there is going to be a nonzero value at position (ii,jj)
 * it is mandatory that all nonzero values are announced first, then 
 * the matrix structure is build with set_structure
 * after that, the functions set, add, and get can be used
 *
 * @param ii row index [0, size - 1]
 * @param jj col index [0, size - 1]
 */
template<class T>
void CSRMatrix<T>::announce(const int &ii,const int &jj) 
{STACK_TRACE(		
	if(this->is_symmetric && ii > jj) {
		NEGF_EXCEPTION("announcing lower triangle matrix elements for a symmetric matrix is forbidden");	
	}	
	if(!this->announced) {
		NEGF_EXCEPTION("announced not initialized");	
	}	
	if(!this->announced[ii]) {
		this->announced[ii] = new announce_list();
		this->announced[ii]->col = jj;	
		this->num_nz++;
	} else {
		int tt = this->announced[ii]->insert(jj);
		this->num_nz += tt;
		NEGF_ASSERT(this->num_nz < MAX_INT_SIZE, "this>num_nz < MAX_INT_SIZE");		
//		this->num_nz +=	this->announced[ii]->insert(jj);
	}
);}

/** sets announced structure to matrix and deletes announcement tree
 */
template<class T>
void CSRMatrix<T>::set_structure()
{STACK_TRACE(		
	// we need to start at a fresh point
	if(icol || prow || nonzeros) {
		NEGF_EXCEPTION("set structure called twice.");	
	}			
	// we need an announcement tree
	if(!this->announced) {
		NEGF_EXCEPTION("announcment tree not build");	
	}	
	
	// structure must have diagonal elements (requested by pardiso), 
	// therefore we have to add them to the counter
	this->num_nz += this->size;
	
	// we need at least size  nonzero elements
	if((signed)this->size > this->num_nz) {
		NEGF_EXCEPTION("there are less nonzero elements than matrix size ... matrix is therefore surely singular :-(");	
	}

	logmsg->emit(LOG_INFO_L2, "building matrix structure with %d nonzero entries, matrix size: %d, filling of %5.5f %%"
					, num_nz, size, 100.0 * ((double)num_nz / ((double)(size * size))));

	// allocate 
	nonzeros = new   T[num_nz];
	icol     = new int[num_nz];
	prow     = new int[size + 1];  
	
	if(!nonzeros || !icol || !prow) {
		NEGF_EXCEPTION("could not allocate memory for matrix");	
	}
					
	announce_list* ptr;
	announce_list* ptr_del;
	
	int ii = 0;
	// loop over rows and create structure
	for(unsigned int rr = 0; rr < size; rr++) {
		// stop if row was not announced
		if(!this->announced[rr]) {
			std::ostringstream sout;
			sout << "error in row " << rr << ". no element announced";
			NEGF_EXCEPTION(sout.str().c_str());
		} else {
			// set current pos to prow
			prow[rr] = ii + fidx; // include fortran indices				
			ptr = this->announced[rr];
			// get first entry
			while(ptr->prev) {
				ptr = ptr->prev;	
			}
			// loop over entries, get colums and delete linked list
			while(ptr) {
				this->nonzeros[ii] = (T)0.0;
			 	this->icol[ii]     = ptr->col + fidx;
			 	ii++;
			 	ptr_del = ptr;
				ptr     = ptr->next;
				delete ptr_del;				
			}			
			this->announced[rr] = 0;
		}	
	}
	// set last row pointer
	prow[size] 		= ii + fidx;
	NEGF_ASSERT(prow[size] == num_nz + fidx, "prow[size] == num_nz + fidx");	
		
	// delete announced
	delete[] this->announced;
	this->announced 	= 0;	
	this->structure_set = true;
	logmsg->emit(LOG_INFO_L2, "finished creating matrix structure");
);}


/** creates full matrix 
 * 
 * @return pointer to newly created array of length size * size where element (ii,jj) is at ret[ii * size + jj]
 */
template<class T>
T* CSRMatrix<T>::create_full_matrix() const 
{STACK_TRACE(
	T* full = new T[this->size * this->size];
	for(unsigned int ii = 0; ii < this->size; ii++) {
		for(unsigned int jj = 0; jj < this->size; jj++) {
			full[ii * this->size + jj] = this->get(ii,jj);	
		}	
	}	
	return full;
);}


/** check whether the matrix structure has errors */
template<class T>
bool CSRMatrix<T>::check_matrix() const 
{STACK_TRACE(
  int i, j, k;

  // Checking if prow is in ascending order.
  i = 0;
  while ((i!=(signed)this->size)&&(prow[i]<=prow[i+1])) i++;
  if (i!=(signed)this->size) {
  	 logmsg->emit(LOG_INFO_L2, "prow not in ascendig order at i = %d: %d next %d, size: %d, num_nz: %d", i, prow[i], prow[i+1],this->size, this->num_nz) ;
  	 return false;
  }

  // Checking if icol components are in order and within bounds.
  for (i=0; i!=(signed)this->size; i++) {
    j = prow[i]     - fidx;
    k = prow[i+1]-1 - fidx;
    if (j<=k) {
      if ((icol[j] < (0 + fidx)) || (icol[k] >=((signed)this->size + fidx))) {
      	logmsg->emit(LOG_INFO_L2, "icol[k] %d >= this->n %d,  k: %d", icol[k], this->size, k);
      	return false;
      }
      while ((j!=k)&&(icol[j]<icol[j+1])) j++;
      if (j!=k) {
      	logmsg->emit(LOG_INFO_L2, "i != k %d %d", i, k);
      	return false;      
      }
    }
  }
  return true;
);}


/** check if a row or column of the matrix contains no nonzero entry */
template<class T>
bool CSRMatrix<T>::check_singular_structure() const 
{STACK_TRACE(
	if(prow && icol) {
		logmsg->init_progress_bar(LOG_INFO_L2,"checking whether matrix structure is singular", size);
		int  next = 0;
		bool have_col_ii;
		for(int ii = 0; ii < (signed)size; ii++) {
			have_col_ii = false;
			if(prow[ii] == prow[ii + 1]) {
				logmsg->end_progress_bar();				
				logmsg->emit(LOG_ERROR, "matrix is singular for row: %d", ii);
				return true;	
			}	
			for(int jj = 0; jj < num_nz; jj++) {
				if(icol[jj] == ii + fidx) {
					have_col_ii = true;	
					break;
				}	
			}
			if(!have_col_ii) {
				logmsg->end_progress_bar();				
				logmsg->emit(LOG_ERROR, "matrix is singular for col: %d", ii);
				return true;	
			}
			if(ii == next) {
				next = logmsg->set_progress_bar(ii, size);	
			}
		}
		logmsg->end_progress_bar();		
	} else {
		logmsg->emit(LOG_WARN, "check_singular_structure called while structure not set yet");	
	}
	return false;
);}


/** returns whether the matrix has symmetric structure
 */
template<class T>
bool CSRMatrix<T>::symmetric() const 
{STACK_TRACE(
	return this->is_symmetric;
);}


/** returns whether the matrix has symmetric strcture OR the values of a nonsymmetric structure are symmetric
 */
template<class T>
bool CSRMatrix<T>::symmetric_by_value() const 
{STACK_TRACE(	
	logmsg->emit(LOG_INFO, "checking for symmetry");	
	if(this->is_symmetric) return true;
	if(this->size < 10000) {
		for(int ii = 0; ii < (signed)this->size; ii++) {			
			for(int jj = ii; jj < (signed)this->size; jj++) {
				if(std::fabs(this->get(ii,jj) - internal_conj(this->get(jj,ii))) > 1.0e-9) {
					return false;	
				}
			}
		}	
		return true;	
	} else {
		logmsg->emit(LOG_WARN, "not checking whether matrix is symmetric");
		return true;
	}
);}


template<class T>
void CSRMatrix<T>::print_sparse_structure() const 
{STACK_TRACE(
	std::cout << "icol: ";
	for(int ii = 0; ii < num_nz; ii++) {
		std::cout << icol[ii] << "  ";
	}
	std::cout << "\nprow: ";
	for(uint ii = 0; ii <= size; ii++) {
		std::cout << prow[ii] << "  ";
	}
	std::cout << " / size would be: " << num_nz << "\nnz: ";
	for(int ii = 0; ii < num_nz; ii++) {
		std::cout << nonzeros[ii] << "  ";
	}		
	std::cout << "\n";
);}	

template<class T>
std::ostream& operator<<(std::ostream& stream,const CSRMatrix<T> &mat) 
{STACK_TRACE(
	int width = 14;
	stream.precision(width - 5);	
	stream.setf(std::ios::left);

	for(int ii = 0; ii < (signed)mat.get_size(); ii++) {
		for(int jj = 0; jj < (signed)mat.get_size(); jj++) {
			stream << std::setw(width) << mat.get(ii,jj) << " ";		
		}	
		stream << "\n";
	}
	return stream;	
);}

template<class T>
void CSRMatrix<T>::print_stats() const 
{STACK_TRACE(
	T offdiag = 0;
	T diag    = 0;
	int num_offdiag = 0;
	for(int ii = 0; ii < (signed)this->get_size(); ii++) {
		diag += abs(this->get(ii,ii));
		for(int jj = ii + 1; jj < (signed)this->get_size(); jj++) {
			offdiag += abs(this->get(ii,jj));		
			num_offdiag++;
		}				
	}
	std::cout << "diagonal: num = " << this->get_size() << " val: " << diag / (double)this->get_size() << "  "
	          << "offdiag:  num = " << 2 * num_offdiag << " val: " << offdiag / (double)num_offdiag << "\n";
);}
	

template<> template<> void CSRMatrix<cplx>::mult_vec(const cplx* in, cplx* out) const;


/** multiply a vector with matrix Ain = out
 * 
 * @param in  input vector
 * @param out output vector
 */
template<class T> template<class B>
void CSRMatrix<T>::mult_vec(const B* in, B* out) const 
{STACK_TRACE(
	if(this->prow && this->icol	&& this->nonzeros) {		
						
		int until, from;	
		// --------------------------------------------------
		// distinguish between symmetric and nonsymmetric
		// --------------------------------------------------		
		if(this->is_symmetric) {					
			// clean up out
			for(int ii = 0; ii < (signed)this->size; ii++) {
				out[ii] = (T)0.0;	
			}				
			// perform symmetric matrix vector product
			for(int ii = 0; ii < (signed)this->size; ii++) {
				until    = this->prow[ii + 1] - fidx;
				from     = this->prow[ii] - fidx;
				out[ii] += this->nonzeros[from] * in[ii]; // diagonal	
		//		NEGF_ASSERT(from >= 0 && from < this->num_nz && until >= 0 && until <= this->num_nz, "from >= 0 && from < this->size && until >= 0 && until <= this->size");
				for(int jj = from + 1; jj < until; jj++) {
					out[ii] += this->nonzeros[jj] * in[(this->icol[jj] - fidx)];
					out[(this->icol[jj] - fidx)] += this->nonzeros[jj] * in[ii];
				}		
			}			
		} else {
			// perform nonsymmetric matrix vector product
			for(int ii = 0; ii < (signed)this->size; ii++) {
				until    = this->prow[ii + 1] - this->fidx;
				from     = this->prow[ii] - this->fidx;
				out[ii]  = 0.0;			
				for(int jj = from; jj < until; jj++) {
					out[ii] += this->nonzeros[jj] * in[(this->icol[jj] - this->fidx)];
				}		
			}			
		}		
	} else {
		NEGF_EXCEPTION("cannot multiply with noninitialized matrix");	
	}		
);}


/** multiply 'multiple' vectors with matrix
 * 
 * in FEM formulation of multiple partial differential equations on a grid
 * as the kp stuff here, the rhs matrix depends on the mesh and leads to 
 * a 'block' sparse matrix, where each block is a diagonal identity matrix.
 * this leads to an independence between the single solutions, therefore allowing
 * only to save this part and to perform the multiplication at the same time
 * 
 * now, suppose that the matrix B is a (n * kpn) x (n * kpn) matrix, where
 * kpn denotes the number of partial differential equations solved at the same time
 * 
 * @param in input vector of length n * kpn, values stored as [x1_0, x1_1, ..., x1_kpn, x2_0, x2_1, ...]
 * @param out output vector of length n * kpn
 * @param kpn number of partial differential equataions > 1
 */
template<class T> template<class B>
void CSRMatrix<T>::mult_vec_multiple(const B* in, B* out, int kpn) const 
{STACK_TRACE(
	if(kpn <= 1) {
		NEGF_EXCEPTION("kpn must be greater than one");
	}
	int msize    = (signed)this->size;
	int vec_size = msize * kpn;
	
	if(this->prow && this->icol	&& this->nonzeros) {
		int until, from, off;					
		// --------------------------------------------------
		// distinguish between symmetric and nonsymmetric
		// --------------------------------------------------		
		if(this->is_symmetric) {					
			// clean up out
			for(int ii = 0; ii < vec_size; ii++) {
				out[ii] = (T)0.0;	
			}				
			// perform symmetric matrix vector product
			for(int ii = 0; ii < msize; ii++) {
				until    = this->prow[ii + 1] - fidx;
				from     = this->prow[ii] - fidx;
				off      = ii * kpn;
				for(int nn = 0; nn < kpn; nn++) {
					out[off + nn] += this->nonzeros[from] * in[off + nn]; // diagonal	
				}
				NEGF_ASSERT(from >= 0 && from < this->num_nz && until >= 0 && until <= this->num_nz, "from >= 0 && from < this->size && until >= 0 && until <= this->size");
				for(int jj = from + 1; jj < until; jj++) {
					for(int nn = 0; nn < kpn; nn++) {
						out[off + nn] += this->nonzeros[jj] * in[(this->icol[jj] - fidx) * kpn + nn];
						out[(this->icol[jj] - fidx) * kpn + nn] += this->nonzeros[jj] * in[off + nn];
					}
				}		
			}			
		} else {
			// perform nonsymmetric matrix vector product
			for(int ii = 0; ii < (signed)this->size; ii++) {
				until    = this->prow[ii + 1] - this->fidx;
				from     = this->prow[ii] - this->fidx;
				for(int nn = 0; nn < kpn; nn++) {
					out[ii * kpn + nn]  = 0.0;	
				}		
				for(int jj = from; jj < until; jj++) {
					for(int nn = 0; nn < kpn; nn++) {
						out[ii * kpn + nn] += this->nonzeros[jj] * in[(this->icol[jj] - this->fidx) * kpn + nn];
					}
				}		
			}			
		}		
	} else {
		NEGF_EXCEPTION("cannot multiply with noninitialized matrix");	
	}		
);}


template<class T>	
void CSRMatrix<T>::save_to_file(const char* filename) const 
{STACK_TRACE(
	if(prow && icol	&& nonzeros) {
		std::fstream fout(filename, std::ios::out);
		if(fout) {
			fout << this->size << "\n" << this->num_nz << "\n";
			for(uint ii = 0; ii < size + 1; ii++) {
				fout << prow[ii] << "\n";	
			}
			for(int ii = 0; ii < num_nz; ii++) {
				fout << icol[ii] << "\n";
			}
			fout.precision(20);
			for(int ii = 0; ii < num_nz; ii++) {
				fout << nonzeros[ii] << "\n";	
			}		
			fout.close();
		} else {
			logmsg->emit(LOG_ERROR, "error while trying to read %s:",filename);
			NEGF_FEXCEPTION("cannot open file %s",filename);	
		}			
	} else {
		NEGF_EXCEPTION("cannot save noninitialized matrix");	
	}
);}


/** read matrix from file
 *
 * the file is expected to consist of one entry per row, entries in the following order:
 * 1. matrix size
 * 2. number of nonzeros
 * 3. boolean if matrix is symmetric
 * 4. boolean if matrix is factorized
 * 5. prow, N+1 lines where N is the matrix dimension
 * 6. icol, NZ  lines where NZ is the number of nonzeros
 * 7. nonzeros, NZ lines
 * 8. the control number 101
 * @param filename of matrix which was stored by save_to_file
 */
template<class T>
void CSRMatrix<T>::read_from_file(const char* filename) 
{STACK_TRACE(
	if(this->prow || this->icol || this->nonzeros) {
		this->clear_all();
	}
	std::fstream fin(filename, std::ios::in);
	if(fin) {
		fin >> this->size;
		fin >> this->num_nz;
//		fin >> this->is_symmetric;
//		fin >> this->is_factorized;
		if(this->size <= 0) {
			NEGF_EXCEPTION("invalid matrix size (<=0)");
		}	
		if(this->num_nz <= 0) {
			NEGF_EXCEPTION("invalid number of nonzeros (<=0)");	
		}
		this->prow     = new int[this->size + 1];
		this->icol     = new int[this->num_nz];
		this->nonzeros = new   T[this->num_nz];
		NEGF_ASSERT(this->prow && this->icol && this->nonzeros, "memory allocation failed");
		for(uint ii = 0; ii < size + 1; ii++) {
			NEGF_ASSERT(!fin.eof(), "end of file reached too early.");
			fin >> this->prow[ii];	
		}
		for(int ii = 0; ii < num_nz; ii++) {
			NEGF_ASSERT(!fin.eof(), "end of file reached too early.");
			fin >> this->icol[ii];	
		} 
		for(int ii = 0; ii < num_nz; ii++) {
			NEGF_ASSERT(!fin.eof(), "end of file reached too early.");
			fin >> this->nonzeros[ii]; 
		}
		double tmp = -88888.88888;
		fin >> tmp;		// last line should be empty, expect to read in nothing
		//if (fin.eof() || tmp!=-88888.88888) {
		//	NEGF_FEXCEPTION("reached end of file %s or found strange value of %g",tmp);
		//}
		fin >> tmp;		// expect to read in nothing one more time at the end before reaching EOF (dunno why)
		//if (fin.eof() && tmp!=-88888.88888) {
		//	NEGF_FEXCEPTION("reached end of file %s but found strange value of %g before",tmp);
		//}
		if (!fin.eof()) {
			tmp = -88888.88888;
			fin >> tmp;
			char buf[1000];
			if (tmp==-88888.88888) {
				sprintf(buf,"<nothing>");
			} else {
				sprintf(buf,"%g",tmp);
			}
			if (fin.eof()) {
				NEGF_FEXCEPTION("reached end of file %s one line too late (found %s before)",filename,buf);
			}
			NEGF_FEXCEPTION("should have reached end of file %s now (after %d lines), found %s instead.",
					filename, size+1 + num_nz + num_nz, buf);
		}
		fin.close();
		this->structure_set = true;
		return;
	} else {
		NEGF_FEXCEPTION("cannot open file %s",filename);
	}	
);}
		
		
template<class T>
void CSRMatrix<T>::ILU_factorize() 
{STACK_TRACE(	
	if(this->is_symmetric) {
		NEGF_EXCEPTION("ILU factorization of symmetric matrix is not possible");	
	}
	
	if(this->prow && this->icol	&& this->nonzeros) {

		if(this->factorized()) {
			NEGF_EXCEPTION("matrix is already ILU factorized");	
		}

		int until, from, mid, kk;
		int next_bar_update = 0;
		logmsg->init_progress_bar(LOG_INFO_L1,"performing ilu factorization", size);
		// for each row
		for(int ii = 0; ii < (signed)this->size; ii++) {			
			until    = this->prow[ii + 1] - this->fidx;
			from     = this->prow[ii] - this->fidx;			
			// find mid element
			mid = -1;
			for(int jj = from; jj < until; jj++) {
				if((this->icol[jj] - this->fidx) == ii) {
					mid = jj;
					break;	
				} 		
			}						
			if(mid == -1) {
				NEGF_EXCEPTION("unable to find mid element");	
			}
			// for left part of row
			for(int skk = from; skk < mid; skk++) {
				// get 'real' matrix index
				kk = this->icol[skk] - this->fidx;
				// divide element by diagonal
				this->nonzeros[skk] /= this->get(kk,kk);				
				// for rest of row
				for(int sjj = skk + 1; sjj < until; sjj++) {
					this->nonzeros[sjj] -= this->nonzeros[skk] * this->get(kk, this->icol[sjj] - this->fidx);		
				} 			
			}					
			if(ii == next_bar_update) {
				next_bar_update = logmsg->set_progress_bar(ii, size);	
			}				
		}
		logmsg->end_progress_bar();
		this->is_factorized = true;
	} else {
		NEGF_EXCEPTION("cannot perform ILU on noninitialized matrix");	
	}				
);}

/** use ILU factorization to precondition x
 * 
 *  means, we solve approximately \f$\mathbf{A}x = b\f$ by using the incomplete ILU decomposition 
 * 
 * @param x array of length size where the result will be stored
 * @param b array of length size, giving the rhs of the equation system
 */
template<class T> template<class B>
void CSRMatrix<T>::ILU_precondition(B* x, const B*b) const 
{STACK_TRACE(
	if(!this->factorized()) {
		NEGF_EXCEPTION("ILU_precondition called on not factorized matrix");	
	}	
	
	B* y = new T[this->size];
	for(int ii = 0; ii < (signed)this->size; ii++) {
		x[ii] = 0.0;	
		y[ii] = 0.0;
	}	
	int until, from, mid, kk;

	// fwd substitution
	for(int ii = 0; ii < (signed)this->size; ii++) {			
		until    = this->prow[ii + 1] - this->fidx;
		from     = this->prow[ii] - this->fidx;				
		// find mid element
		mid = -1;
		for(int jj = from; jj < until; jj++) {
			if((this->icol[jj] - this->fidx) == ii) {
				mid = jj;
				break;	
			} 		
		}						
		if(mid == -1) {
			NEGF_EXCEPTION("unable to find mid element");	
		}
		y[ii] = b[ii];
		// lower triangular with ones on diagonal
		for(int jj = from; jj < mid; jj++) {
			kk = this->icol[jj] - this->fidx;
			NEGF_ASSERT(kk >= 0 && kk < ii, "kk >= 0 && kk < ii");
			y[ii] -= this->nonzeros[jj] * y[kk];
		}	
	}
	// backward substitution
	for(int ii = this->size - 1; ii >= 0; ii--) {
		until    = this->prow[ii + 1] - this->fidx;
		from     = this->prow[ii] - this->fidx;				
		// find mid element
		mid = -1;
		for(int jj = from; jj < until; jj++) {
			if((this->icol[jj] - this->fidx) == ii) {
				mid = jj;
				break;	
			} 		
		}						
		if(mid == -1) {
			NEGF_EXCEPTION("unable to find mid element");	
		}
		x[ii] = y[ii];
		// upper triangular part			
		// lower triangular with ones on diagonal
		for(int jj = mid + 1; jj < until; jj++) {
			kk = this->icol[jj] - this->fidx;
			NEGF_ASSERT(kk > ii && kk < (signed)this->size, "kk > ii && kk < (signed)this->size");
			x[ii] -= this->nonzeros[jj] * x[kk];
		}		
		NEGF_ASSERT(this->nonzeros[mid] != 0.0, "this->nonzeros[mid] != 0.0");
		x[ii] /= this->nonzeros[mid];
	}	
);}
		
		
} // end namespace

#endif /*CSRMATRIX_H_NEGF*/
