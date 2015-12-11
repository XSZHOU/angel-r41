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
#ifndef MPIS_H_NEGF
#define MPIS_H_NEGF

#include "all.h"
//#include "mpi.h"  Commented out by S.Z.
#include <mpi.h>
#include "zlib.h"

using namespace std;

namespace negf {
	
	/** Interface for all functions involving MPI calls. Ideally only methods from this class should be used, and no direct MPI methods. <BR>
	 *  Make sure that there is only one instance of this class created at the very beginning
	 *  and destroyed at the very end of the simulation! */
	class myMPI
	{
	public:
		myMPI() throw (Exception *);
		~myMPI();
		
		int get_rank();		 // >=0
		int get_num_procs(); // >=1
		
		// these functions do not change the arguments, but MPI consists of C code which does not have the const-concept
		void send(double & data, 		int & destination, int & tag);
		void send(int & data, 			int & destination, int & tag);
		void send(bool & data, 			int & destination, int & tag);
		void send(const vector<double> & data,int & destination, int & tag);
		void send(const vector<int> & data, 	int & destination, int & tag);
		void send(const vector<bool> & data, 	int & destination, int & tag);
		void send(string & data, 		int & destination, int & tag);
		void send(Matc & data, int & destination, int & tag);
		template<class Mtype>
		void send(const vector<Mtype> & data, int & destination, int & tag) throw (Exception *);
		void send_antihermitian(const BMatc & data, int & destination, int & tag) throw (Exception *);
		void send_antihermitians(const vector<BMatc> & data, int & destination, int & tag) throw (Exception *);
		void send(Matd & data, int & destination, int & tag);
		void broadcast(double & data, int & root);
		void broadcast(int & data, int & root);
		void broadcast(bool & data, int & root);
		void broadcast(vector<double> & data, int & root);
		void broadcast(vector<int> & data, int & root);
		void broadcast(vector<bool> & data, int & root);
		
		// these functions only change the first argument (output data)
		void recv(double & data, 		int & source, int & tag) throw (Exception *);
		void recv(int & data, 			int & source, int & tag) throw (Exception *);
		void recv(bool & data, 			int & source, int & tag) throw (Exception *);
		void recv(vector<double> & data, uint & size, int & source, int & tag) throw (Exception *);
		void recv(vector<int> & data, 	uint & size,  int & source, int & tag) throw (Exception *);
		void recv(vector<bool> & data, 	uint & size,  int & source, int & tag) throw (Exception *);
		void recv(string & data, 		uint & size,  int & source, int & tag) throw (Exception *);
		void recv(Matc & data, int & source, int & tag) throw (Exception *);
		template<class Mtype>
		void recv(vector<Mtype> & data, int & source, int & tag) throw (Exception *);
		void recv_antihermitian(BMatc & data, int & source, int & tag) throw (Exception *);
		void recv_antihermitians(vector<BMatc> & data, int & source, int & tag, 
				unsigned char * real_compressed = 0, unsigned char * imag_compressed = 0,
				unsigned char * real_uncompressed = 0, unsigned char * imag_uncompressed = 0) throw (Exception *);
		void recv(Matd & data, int & source, int & tag) throw (Exception *);
		
		void synchronize_processes() throw (Exception *);	//!< synchronize all MPI processes
		
		/** Determine which processes receive and which send from a list of needed processes for every processes and a list which
		 *  process was already previously computed */
		void determine_senders_receivers(const vector< vector<int> > & processes_needed, const vector<bool> & process_was_computed, 
										vector<bool> & senders, vector<bool> & receivers);
		
		void terminate();
		
		const string & get_hostname() const { return this->hostname; }
	
	protected:
	
		string hostname;
	};
	
	
	/** Receive an array of complex matrices from a specific process
	 *  the transmission is done in a compressed (zlib) format */
	template<class Mtype>
	void myMPI::recv(vector<Mtype> & data, int & source, int & tag) throw (Exception *)
	{try { 			// SWIG doesn't like STAKC_TRACE( macro
		MPI_Status status; 
		status.MPI_ERROR = 0;
		
		int num_matrices = data.size();
		NEGF_ASSERT(num_matrices > 0, "matrix vector must have size >0");
	
		//const int max_distance = data[0].num_rows()-1;
		const int max_distance = data[0].num_offdiags; // for Matc type, this will be num_rows()-1
		
		int n = data[0].num_rows();
		for (int ii=0; ii < num_matrices; ii++) {
			NEGF_ASSERT(int(data[ii].num_rows())==n && int(data[ii].num_cols())==n, "expect array of equal size square matrices.");
		}
		
		int n2 = max(0, n-max_distance-1);
		
		int N = n*n - n2*(n2+1);
		unsigned long num_doubles = num_matrices * N;
		
		// --------------------------------------------
		// receive length of char arrays
		// --------------------------------------------
		int real_char_size = -1;
		int imag_char_size = -1;
		int error = MPI_Recv(&real_char_size, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
		NEGF_FASSERT(error==0, "An error (%d) occurred while receiving an integer via MPI.",error);
		NEGF_ASSERT(status.MPI_ERROR == 0, "An error occurred while receiving an integer via MPI.");
		int tag2 = tag + 1;
		error = MPI_Recv(&imag_char_size, 1, MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
		NEGF_FASSERT(error==0, "An error (%d) occurred while receiving an integer via MPI.",error);
		NEGF_ASSERT(status.MPI_ERROR == 0, "An error occurred while receiving an integer via MPI.");
	
		// --------------------------------------------
		// receive char arrays
		// --------------------------------------------	
	#ifdef HEAPARRAYS
		unsigned char * real_data = new unsigned char[real_char_size];
		unsigned char * imag_data = new unsigned char[imag_char_size];
	#else
		unsigned char real_data[real_char_size];
		unsigned char imag_data[imag_char_size];
	#endif	
		
		error = MPI_Recv(real_data, real_char_size, MPI_UNSIGNED_CHAR, source, tag, MPI_COMM_WORLD, &status);
		NEGF_FASSERT(error==0, "An error (%d) occurred while receiving the real part of a matrix via MPI.",error);
		NEGF_ASSERT(status.MPI_ERROR == 0, "An error occurred while receiving the real part of a complex matrix via MPI.");
		
		error = MPI_Recv(imag_data, imag_char_size, MPI_UNSIGNED_CHAR, source, tag2, MPI_COMM_WORLD, &status);
		NEGF_FASSERT(error==0, "An error (%d) occurred while receiving the imag part of a matrix via MPI.",error);
		NEGF_ASSERT(status.MPI_ERROR == 0, "An error occurred while receiving the imag part of a complex matrix via MPI.");
		
		// ------------------------------------
		// decompress
		// from the zlib homepage: "...the destination buffer, which must be large enough to hold the entire uncompressed data."
		// ------------------------------------
		int doublesize = sizeof(double);
		int charsize = sizeof(char);
		NEGF_ASSERT(charsize==1, "class is designed for sizeof(char)=1");
		unsigned long num_chars = num_doubles * doublesize;
		unsigned long real_decomp_size = num_chars + 100; // 100 should not even be necessary
		unsigned long imag_decomp_size = real_decomp_size; 
	#ifdef HEAPARRAYS
		unsigned char * real_uncompressed = new unsigned char[real_decomp_size];
		unsigned char * imag_uncompressed = new unsigned char[imag_decomp_size];
	#else
		unsigned char real_uncompressed[real_decomp_size];
		unsigned char imag_uncompressed[imag_decomp_size];
	#endif	
		unsigned long * real_decompressed_size = &real_decomp_size;
		unsigned long * imag_decompressed_size = &imag_decomp_size;
		
		error = uncompress(real_uncompressed, real_decompressed_size, real_data, real_char_size);
		negf_math::check_decompress_error(error, real_decomp_size, num_chars);
		
		error = uncompress(imag_uncompressed, imag_decompressed_size, imag_data, imag_char_size);
		negf_math::check_decompress_error(error, imag_decomp_size, num_chars);
		
		// ------------------------------------------------
		// cast into doubles 
		// ------------------------------------------------
		double * reals = (double *)(real_uncompressed);
		double * imags = (double *)(imag_uncompressed);
		
		// --------------------------------------------------
		// fill matrices
		// --------------------------------------------------
		int count = 0;
		for (int mm=0; mm < num_matrices; mm++) {
			for (int ii=0; ii<n; ii++) {
				for (int jj=max(0,ii-max_distance); jj<=min(n-1,ii+max_distance); jj++) {
					data[mm](ii+1,jj+1) = reals[count] + constants::imag_unit * imags[count];
					//cout << "data["<<mm<<"]("<<ii+1<<","<<jj+1<<")="<<data[mm](ii+1,jj+1)<<endl;
					count++;
				}
#ifdef USE_BANDED
				// no more entries stored
				continue;
#endif
				for (int jj=0; jj<max(0,ii-max_distance); jj++) {
					data[mm](ii+1,jj+1) = 0.0;
				}
				for (int jj=ii+max_distance+1; jj<n; jj++) {
					data[mm](ii+1,jj+1) = 0.0;
				}
			}
		}
		NEGF_ASSERT((unsigned long)(count)==num_doubles, "count==num_doubles");
		
		// --------------------------------------------------
		// clean up heap arrays
		// --------------------------------------------------
	#ifdef HEAPARRAYS
		delete [] real_data;
		delete [] imag_data;
		delete [] real_uncompressed;
		delete [] imag_uncompressed;
	#endif	
	} catch (Exception * e) { e->append_info(__LINE__,__FILE__,__DATE__,__TIME__, __func__); throw e; };
	}
	
	/** Send an an array of complex matrices to some other process 
	 *  This is done by creating real vectors and compressing them */
	template<class Mtype>
	void myMPI::send(const vector<Mtype> & data, int & destination, int & tag) throw (Exception *)
	{try { // SWIG doesn't like STACK_TRACE()-macro
		//double t1 = MPI_Wtime();
		
		unsigned long num_doubles = negf_math::check_and_get_num_doubles(data,false);
		const int max_distance = data[0].num_offdiags;	// will be n-1 for Matc type
		int n = data[0].num_rows();
		int num_matrices = data.size();
		
		// ----------------------------------------------
		// create double arrays
		// ----------------------------------------------
		unsigned long count = 0;
	#ifdef HEAPARRAYS
		double * real_data = new double[num_doubles];
		double * imag_data = new double[num_doubles];
	#else
		double real_data[num_doubles];
		double imag_data[num_doubles];
	#endif
		for (int mm=0; mm<num_matrices; mm++) {
			for (int ii=0; ii<n; ii++) {
				for (int jj=max(0,ii-max_distance); jj<=min(n-1,ii+max_distance); jj++) {
					real_data[count] = data[mm](ii+1,jj+1).real();
					imag_data[count] = data[mm](ii+1,jj+1).imag();
					count++;
				}
			}
		}
		NEGF_ASSERT(count==num_doubles, "something went wrong.");
			
		// ----------------------------------------------
		// compress those arrays
		// see stuff/compress.cc for a test case
		// ----------------------------------------------
		int doublesize = sizeof(double);
		int charsize = sizeof(char);
		NEGF_ASSERT(charsize==1, "class is designed for sizeof(char)=1");
		
		// cast into byte (unsigned int)
		unsigned char * real_char = (unsigned char *)(real_data);
		unsigned char * imag_char = (unsigned char *)(imag_data);
		unsigned long  num_chars = num_doubles*doublesize;
		
		// compress. from the zlib homepage: "...the destination buffer, which must be at least 0.1% larger than sourceLen plus 12 bytes."
		unsigned long real_comp_size = (unsigned long) (ceil(1.01*double(num_chars))) + 12;
		unsigned long imag_comp_size = real_comp_size;
		unsigned long * real_compressed_size = &real_comp_size;
		unsigned long * imag_compressed_size = &imag_comp_size;
	#ifdef HEAPARRAYS
		unsigned char * real_compressed = new unsigned char[real_comp_size];	
		unsigned char * imag_compressed = new unsigned char[imag_comp_size];
	#else
		unsigned char real_compressed[real_comp_size];
		unsigned char imag_compressed[imag_comp_size];
	#endif
		
		int clevel = constants::mpi_compresslevel;	
		// Z_DEFAULT_COMPRESSION, or between 0 and 9: 1 gives best speed, 9 gives best compression
		
		int err = compress2(real_compressed, real_compressed_size, real_char, num_chars, clevel);	// zlib routine
		negf_math::check_compress_error(err, real_comp_size, num_chars);
		
		err = compress2(imag_compressed, imag_compressed_size, imag_char, num_chars, clevel);	// zlib routine
		negf_math::check_compress_error(err, imag_comp_size, num_chars);
		
		// ----------------------------------------------
		// send array lengths
		// ----------------------------------------------
		// integer is 4 bytes --> max. array size is 2^(4*8-1)
		NEGF_ASSERT(sizeof(int)==4, "sizeof(int)==4 expected.");
		unsigned long max_array_size = 2147483648UL; // this gives 268'435'456 doubles or, at an array size of 1000*1000, roughly 268*2>500 possible k-points
		NEGF_ASSERT(*real_compressed_size < max_array_size && *imag_compressed_size < max_array_size, "maximum array size exceeded.");
		int real_char_size = int(*real_compressed_size);
		int imag_char_size = int(*imag_compressed_size);
		
		//double t2 = MPI_Wtime();
		
		MPI_Send(&real_char_size, 1, MPI_INT, destination, tag, MPI_COMM_WORLD);
		int tag2 = tag+1;
		MPI_Send(&imag_char_size, 1, MPI_INT, destination, tag2, MPI_COMM_WORLD);
		
		//double t3 = MPI_Wtime();
			
		// ----------------------------------------------
		// send arrays
		// ----------------------------------------------
		MPI_Send(real_compressed, real_char_size, MPI_UNSIGNED_CHAR, destination, tag, MPI_COMM_WORLD);
		MPI_Send(imag_compressed, imag_char_size, MPI_UNSIGNED_CHAR, destination, tag2, MPI_COMM_WORLD);
			
		//double t4 = MPI_Wtime();
		
		// --------------------------------------------------
		// clean up heap arrays
		// --------------------------------------------------
	#ifdef HEAPARRAYS
		delete [] real_data;
		delete [] imag_data;
		delete [] real_compressed;
		delete [] imag_compressed;
	#endif	
		
		//double t5 = MPI_Wtime();
		
		/*double own_time = t2-t1 + t5-t4;
		double mixed_time = t3-t2;
		double mpi_time = t4-t3;
		if (this->get_rank()==11 || this->get_rank()==10) {
			char buf[1000]; sprintf(buf,"p%d(%s)->p%d: own+mpi=%.5g, wait=%.5g, (%d reals, %d imags)", 
					this->get_rank(), this->get_hostname().c_str(), destination, own_time+mpi_time, mixed_time, real_char_size, imag_char_size);
			//logmsg->emit(LOG_INFO_L3,buf);
			cout << buf << endl;
		}*/
	} catch (Exception * e) { e->append_info(__LINE__,__FILE__,__DATE__,__TIME__, __func__); throw e; };
	}

} // end namespace 
	
#endif /*MPI_H_NEGF*/
