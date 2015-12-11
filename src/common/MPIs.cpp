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
#include "MPIs.h"
using namespace negf;

#include <unistd.h> // for gethostname

//#include "zlib.h"

/** Start up MPI */
myMPI::myMPI() throw (Exception *)
{STACK_TRACE(
	int argc = 0;
	char** argv = /*NULL*/ new char*;
	MPI_Init(&argc, &argv);
	
	char name[32] ;
	gethostname(name,32);
	NEGF_ASSERT(name!=NULL, "gethostname failed.");
	this->hostname = name;
	
	if (argv) delete argv;
);}

/** Destroy the myMPI class. MPI itself is shut down in terminate() */
myMPI::~myMPI()
{STACK_TRACE(
	//MPI_Finalize();
);}

/** Shut down MPI. */
void myMPI::terminate()
{STACK_TRACE(
	MPI_Finalize();
);}


// ---------------------------------
// SEND ROUTINES
// ---------------------------------

/** Send a double to some other process
 * @param data the double to be sent
 * @param destination the receiver process rank
 * @param tag a tag to give the receiver the chance of doublechecking if the right msg was received
 */
void myMPI::send(double & data, int & destination, int & tag)
{STACK_TRACE(
	MPI_Send(&data, 1, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);
);}

/** Send an integer to some other process
 * @param data the integer to be sent
 * @param destination the receiver process rank
 * @param tag a tag to give the receiver the chance of doublechecking if the right msg was received
 */
void myMPI::send(int & data, int & destination, int & tag)
{STACK_TRACE(
	MPI_Send(&data, 1, MPI_INT, destination, tag, MPI_COMM_WORLD);
);}

/** Send a boolean to some other process
 * @param data the bool to be sent
 * @param destination the receiver process rank
 * @param tag a tag to give the receiver the chance of doublechecking if the right msg was received
 */
void myMPI::send(bool & data, int & destination, int & tag)
{STACK_TRACE(
	int data_int = (data) ? 1 : 0;
	MPI_Send(&data_int, 1, MPI_INT, destination, tag, MPI_COMM_WORLD);
);}

/** Send an array of doubles to some other process
 * @param data the double array to be sent
 * @param destination the receiver process rank
 * @param tag a tag to give the receiver the chance of doublechecking if the right msg was received
 */
void myMPI::send(const vector<double> & data, int & destination, int & tag)
{STACK_TRACE(
	double data_copy[data.size()];
	for (uint ii = 0; ii < data.size(); ii++) {
		data_copy[ii] = data[ii];
	}
	MPI_Send(&data_copy, data.size(), MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);
);}

/** Send an array of integers to some other process
 * @param data the integer array to be sent
 * @param destination the receiver process rank
 * @param tag a tag to give the receiver the chance of doublechecking if the right msg was received
 */
void myMPI::send(const vector<int> & data, int & destination, int & tag)
{STACK_TRACE(
	int data_copy[data.size()];
	for (uint ii = 0; ii < data.size(); ii++) {
		data_copy[ii] = data[ii];
	}
	MPI_Send(&data_copy, data.size(), MPI_INT, destination, tag, MPI_COMM_WORLD);
);}

/** Send an array of booleans to some other process
 * @param data the boolean array to be sent
 * @param destination the receiver process rank
 * @param tag a tag to give the receiver the chance of doublechecking if the right msg was received
 */
void myMPI::send(const vector<bool> & data, int & destination, int & tag)
{STACK_TRACE(
	int data_copy[data.size()];
	for (uint ii = 0; ii < data.size(); ii++) {
		data_copy[ii] = (data[ii]) ? 1 : 0;
	}
	MPI_Send(&data_copy, data.size(), MPI_INT, destination, tag, MPI_COMM_WORLD);
);}

/** Send a string to some other process
 * @param data the string to be sent
 * @param destination the receiver process rank
 * @param tag a tag to give the receiver the chance of doublechecking if the right msg was received
 */
void myMPI::send(string & data, int & destination, int & tag)
{STACK_TRACE(
	char data_char[data.length()];
	sprintf(data_char,"%s",data.c_str());
	//const char * data_char = data.c_str();
	MPI_Send(&data_char, data.length(), MPI_CHAR, destination, tag, MPI_COMM_WORLD);
);}

/** Send a full complex matrix to some other process
 * @param data the matrix to be sent
 * @param destination the receiver process rank
 * @param tag a tag to give the receiver the chance of doublechecking if the right msg was received
 */
void myMPI::send(Matc & data, int & destination, int & tag)
{STACK_TRACE(
	uint n = data.num_rows();
	uint m = data.num_cols();
	double real_data[n*m];
	double imag_data[n*m];
	for (uint ii = 0; ii < n; ii++) {
		for (uint jj = 0; jj < m; jj++) {
			real_data[ii*m+jj] = data(ii+1,jj+1).real();
			imag_data[ii*m+jj] = data(ii+1,jj+1).imag();
		}
	}
	MPI_Send(&real_data, n*m, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);
	int tag2 = tag+1;
	MPI_Send(&imag_data, n*m, MPI_DOUBLE, destination, tag2, MPI_COMM_WORLD);
);}


/** Send an antihermitian complex matrix (A=-A+) to some other process 
 *  Only the lower part (including the diagonal) is sent over the network,
 *  the rest will be reconstructed in recv_antihermitian
 */
void myMPI::send_antihermitian(const BMatc & data, int & destination, int & tag) throw (Exception *)
{STACK_TRACE(
	NEGF_ASSERT(data.num_rows()==data.num_cols(), "only possible w/ square matrix");
	int N = data.num_rows()*(data.num_rows()+1) / 2;
	int count = 0;
	double real_data[N];
	double imag_data[N];
	for (uint ii = 0; ii < data.num_rows(); ii++) {
		for (uint jj = 0; jj <= ii; jj++) {
			real_data[count] = data(ii+1,jj+1).real();
			imag_data[count] = data(ii+1,jj+1).imag();
			count++;
		}
	}
	NEGF_ASSERT(count==N, "something went wrong.");
	MPI_Send(&real_data, N, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);
	int tag2 = tag+1;
	MPI_Send(&imag_data, N, MPI_DOUBLE, destination, tag2, MPI_COMM_WORLD);
);}


/** Send a full double matrix to some other process
 * @param data the matrix to be sent
 * @param destination the receiver process rank
 * @param tag a tag to give the receiver the chance of doublechecking if the right msg was received
 */
void myMPI::send(Matd & data, int & destination, int & tag)
{STACK_TRACE(
	uint n = data.num_rows();
	uint m = data.num_cols();
	double real_data[n*m];
	for (uint ii = 0; ii < n; ii++) {
		for (uint jj = 0; jj < m; jj++) {
			real_data[ii*m+jj] = data(ii+1,jj+1);
		}
	}
	MPI_Send(&real_data, n*m, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);
);}


/** Send a double to every process 
 * @param data the double to be sent
 * @param root broadcast() is called by all processes. the process with number root sends data, the rest receives
 */
void myMPI::broadcast(double & data, int & root)
{STACK_TRACE(
	MPI_Bcast(&data, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
);}


/** Send an integer to every process 
 * @param data the integer to be sent
 * @param root broadcast() is called by all processes. the process with number root sends data, the rest receives
 */
void myMPI::broadcast(int & data, int & root)
{STACK_TRACE(
	MPI_Bcast(&data, 1, MPI_INT, root, MPI_COMM_WORLD);
);}


/** Send a boolean to every process 
 * @param data the integer to be sent
 * @param root broadcast() is called by all processes. the process with number root sends data, the rest receives
 */
void myMPI::broadcast(bool & data, int & root)
{STACK_TRACE(
	int data_int = data;
	MPI_Bcast(&data_int, 1, MPI_INT, root, MPI_COMM_WORLD);
	if (this->get_rank()!=root) { // receivers
		data = (data_int==1) ? true : false;
	}
);}


/** Send a double array to every process 
 * @param data the double array to be sent. The array must have the correct length in all processes
 * @param root broadcast() is called by all processes. the process with number root sends data, the rest receives
 */
void myMPI::broadcast(vector<double> & data, int & root)
{STACK_TRACE(
	double data_copy[data.size()];
	for (uint ii = 0; ii < data.size(); ii++) {
		data_copy[ii] = data[ii];
	}
	MPI_Bcast(&data_copy, data.size(), MPI_DOUBLE, root, MPI_COMM_WORLD);
	for (uint ii = 0; ii < data.size(); ii++) {
		data[ii] = data_copy[ii];
	}
);}

/** Send an integer array to every process 
 * @param data the integer array to be sent. The array must have the correct length in all processes
 * @param root broadcast() is called by all processes. the process with number root sends data, the rest receives
 */
void myMPI::broadcast(vector<int> & data, int & root)
{STACK_TRACE(
	int data_copy[data.size()];
	for (uint ii = 0; ii < data.size(); ii++) {
		data_copy[ii] = data[ii];
	}
	MPI_Bcast(&data_copy, data.size(), MPI_INT, root, MPI_COMM_WORLD);
	for (uint ii = 0; ii < data.size(); ii++) {
		data[ii] = data_copy[ii];
	}
);}


/** Send a boolean array to every process 
 * @param data the boolean array to be sent. The array must have the correct length in all processes
 * @param root broadcast() is called by all processes. the process with number root sends data, the rest receives
 */
void myMPI::broadcast(vector<bool> & data, int & root)
{STACK_TRACE(
	int data_copy[data.size()];
	for (uint ii = 0; ii < data.size(); ii++) {
		data_copy[ii] = (data[ii]==true) ? 1 : 0;
	}
	MPI_Bcast(&data_copy, data.size(), MPI_INT, root, MPI_COMM_WORLD);
	for (uint ii = 0; ii < data.size(); ii++) {
		data[ii] = (data_copy[ii]==1) ? true : false;
	}
);}


// ---------------------------------
// RECEIVE ROUTINES
// ---------------------------------

/** Receive a double from a specific process
 * @param data the double in which the data will be stored
 * @param source the process rank from which the data is expected
 * @param tag a tag to double check if the received data was the right one
 */
void myMPI::recv(double & data, int & source, int & tag) throw (Exception *)
{STACK_TRACE(
	MPI_Status status; 
	status.MPI_ERROR = 0;
	
	int error = MPI_Recv(&data, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
	
	NEGF_FASSERT(error==0, "An error (%d) occurred while receiving a double via MPI.",error);
	NEGF_ASSERT(status.MPI_ERROR == 0, "An error occurred while receiving a double via MPI.");
);}


/** Receive an integer from a specific process
 * @param data the integer in which the data will be stored
 * @param source the process rank from which the data is expected
 * @param tag a tag to double check if the received data was the right one
 */
void myMPI::recv(int & data, int & source, int & tag) throw (Exception *)
{STACK_TRACE(
	MPI_Status status; 
	status.MPI_ERROR = 0;
	
	int error = MPI_Recv(&data, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
	
	NEGF_FASSERT(error==0, "An error (%d) occurred while receiving an integer via MPI.",error);
	NEGF_ASSERT(status.MPI_ERROR == 0, "An error occurred while receiving an integer via MPI.");
);}


/** Receive a boolean from a specific process
 * @param data the boolean in which the data will be stored
 * @param source the process rank from which the data is expected
 * @param tag a tag to double check if the received data was the right one
 */
void myMPI::recv(bool & data, int & source, int & tag) throw (Exception *)
{STACK_TRACE(
	MPI_Status status; 
	status.MPI_ERROR = 0;
	
	int data_int = -1;
	int error = MPI_Recv(&data_int, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
	
	NEGF_FASSERT(error==0, "An error (%d) occurred while receiving an integer via MPI.",error);
	NEGF_ASSERT(status.MPI_ERROR == 0, "An error occurred while receiving an integer via MPI.");
	NEGF_ASSERT(data_int==0 || data_int==1, "expected 0 or 1.");
	data = (data_int==1) ? true : false;
);}


/** Receive a double array from a specific process
 * @param data the vector in which the array will be stored
 * @param size the expected size of the vector
 * @param source the process rank from which the data is expected
 * @param tag a tag to double check if the received data was the right one
 */
void myMPI::recv(vector<double> & data, uint & size, int & source, int & tag) throw (Exception *)
{STACK_TRACE(
	MPI_Status status; 
	status.MPI_ERROR = 0;
	double data_copy[size + 1];
	data_copy[size] = 1e100;
	
	int error = MPI_Recv(&data_copy, size, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
	
	NEGF_FASSERT(error==0, "An error (%d) occurred while receiving a vector<double> via MPI.",error);
	NEGF_ASSERT(status.MPI_ERROR == 0, "An error occurred while receiving a vector<double> via MPI.");
	NEGF_ASSERT(data_copy[size]==1e100, "inconsistent array.");
	data.resize(size, 0.0);
	for (uint ii = 0; ii < size; ii++) {
		data[ii] = data_copy[ii];
	}
);} 


/** Receive an integer array from a specific process
 * @param data the vector in which the array will be stored
 * @param size the expected size of the vector
 * @param source the process rank from which the data is expected
 * @param tag a tag to double check if the received data was the right one
 */
void myMPI::recv(vector<int> & data, uint & size, int & source, int & tag) throw (Exception *)
{STACK_TRACE( 
	MPI_Status status; 
	status.MPI_ERROR = 0;
	int data_copy[size + 1];
	data_copy[size] = 100000;
	
	int error = MPI_Recv(&data_copy, size, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
	
	NEGF_FASSERT(error==0, "An error (%d) occurred while receiving a vector<int> via MPI.",error);
	NEGF_ASSERT(status.MPI_ERROR == 0, "An error occurred while receiving a vector<int> via MPI.");
	NEGF_ASSERT(data_copy[size]==100000, "inconsistent array.");
	data.resize(size, 0);
	for (uint ii = 0; ii < size; ii++) {
		data[ii] = data_copy[ii];
	}
);}


/** Receive an integer array from a specific process
 * @param data the vector in which the array will be stored
 * @param size the expected size of the vector
 * @param source the process rank from which the data is expected
 * @param tag a tag to double check if the received data was the right one
 */
void myMPI::recv(vector<bool> & data, uint & size, int & source, int & tag) throw (Exception *)
{STACK_TRACE( 
	MPI_Status status; 
	status.MPI_ERROR = 0;
	int data_copy[size + 1];
	data_copy[size] = 100000;
	
	int error = MPI_Recv(&data_copy, size, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
	
	NEGF_FASSERT(error==0, "An error (%d) occurred while receiving a vector<int> via MPI.",error);
	NEGF_ASSERT(status.MPI_ERROR == 0, "An error occurred while receiving a vector<int> via MPI.");
	NEGF_ASSERT(data_copy[size]==100000, "inconsistent array.");
	data.resize(size, false);
	for (uint ii = 0; ii < size; ii++) {
		NEGF_ASSERT(data_copy[ii]==0 || data_copy[ii]==1, "expected 0 or 1.");
		data[ii] = (data_copy[ii]==1) ? true : false;
	}
);}


/** Receive a string from a specific process
 * @param data the string in which the data will be stored
 * @param size the expected length of the string
 * @param source the process rank from which the data is expected
 * @param tag a tag to double check if the received data was the right one
 */
void myMPI::recv(string & data, uint & size, int & source, int & tag) throw (Exception *)
{STACK_TRACE(
	MPI_Status status; 
	status.MPI_ERROR = 0;
	char data_copy[size + 1];
	data_copy[size] = 'X';
	
	int error = MPI_Recv(&data_copy, size, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
	
	NEGF_FASSERT(error==0, "An error (%d) occurred while receiving a string via MPI.",error);
	NEGF_ASSERT(status.MPI_ERROR == 0, "An error occurred while receiving a string via MPI.");
	NEGF_ASSERT(data_copy[size]=='X', "inconsistent array.");
	data_copy[size]='\0';
	data = data_copy;
	NEGF_ASSERT(data.length()==size, "error.");
);}


/** Receive a full complex matrix from a specific process
 * note: the given matrix must already have the expected size
 * @param data the matrix in which the data will be stored
 * @param source the process rank from which the data is expected
 * @param tag a tag to double check if the received data was the right one
 */
void myMPI::recv(Matc & data, int & source, int & tag) throw (Exception *)
{STACK_TRACE(
	MPI_Status status; 
	status.MPI_ERROR = 0;
	uint n = data.num_rows();
	uint m = data.num_cols();
	
	double real_data[n*m];
	double imag_data[n*m];
	int error = MPI_Recv(&real_data, n*m, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
	
	NEGF_FASSERT(error==0, "An error (%d) occurred while receiving the real part of a matrix via MPI.",error);
	NEGF_ASSERT(status.MPI_ERROR == 0, "An error occurred while receiving the real part of a complex matrix via MPI.");
	
	uint tag2 = tag + 1;
	error = MPI_Recv(&imag_data, n*m, MPI_DOUBLE, source, tag2, MPI_COMM_WORLD, &status);
	NEGF_FASSERT(error==0, "An error (%d) occurred while receiving the imag part of a matrix via MPI.",error);
	NEGF_ASSERT(status.MPI_ERROR == 0, "An error occurred while receiving the imag part of a complex matrix via MPI.");
	
	for (uint ii=0; ii<n; ii++) {
		for (uint jj=0; jj<m; jj++) {
			data(ii+1,jj+1) = real_data[ii*m+jj] + constants::imag_unit * imag_data[ii*m+jj];
		}
	}
);}


/** Receive an antihermitian complex matrix from a specific process (A=-A+)
 *  Only the lower part is expected from the sender, the rest is reconstructed from
 *  A_ij = -A_ji^* = -(Re[A_ji] - i*Im[A_ji])
 */
void myMPI::recv_antihermitian(BMatc & data, int & source, int & tag) throw (Exception *)
{STACK_TRACE(
	MPI_Status status; 
	status.MPI_ERROR = 0;
	NEGF_ASSERT(data.num_rows()==data.num_cols(), "receiver matrix must be square.");
	int n = data.num_rows();
	int N = n*(n+1) / 2;
	
	double real_data[N];
	double imag_data[N];
	int error = MPI_Recv(&real_data, N, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
	NEGF_FASSERT(error==0, "An error (%d) occurred while receiving the real part of a matrix via MPI.",error);
	NEGF_ASSERT(status.MPI_ERROR == 0, "An error occurred while receiving the real part of a complex matrix via MPI.");
	
	uint tag2 = tag + 1;
	error = MPI_Recv(&imag_data, N, MPI_DOUBLE, source, tag2, MPI_COMM_WORLD, &status);
	NEGF_FASSERT(error==0, "An error (%d) occurred while receiving the imag part of a matrix via MPI.",error);
	NEGF_ASSERT(status.MPI_ERROR == 0, "An error occurred while receiving the imag part of a complex matrix via MPI.");
	
	// fill lower half (including diagonal)
	int count = 0;
	for (int ii=0; ii<n; ii++) {
		for (int jj=0; jj<=ii; jj++) {
			data(ii+1,jj+1) = real_data[count] + constants::imag_unit * imag_data[count];
			count++;
		}
	}
	NEGF_ASSERT(count==N, "count==N");
	
	// reconstruct upper half
	for (int ii=1; ii<=n; ii++) {
		for (int jj=ii+1; jj<=n; jj++) {
			data(ii,jj) = -data(jj,ii).real() + constants::imag_unit * data(jj,ii).imag();
		}
	}
);}



/** Receive a full double matrix from a specific process
 * note: the given matrix must already have the expected size
 * @param data the matrix in which the data will be stored
 * @param source the process rank from which the data is expected
 * @param tag a tag to double check if the received data was the right one
 */
void myMPI::recv(Matd & data, int & source, int & tag) throw (Exception *)
{STACK_TRACE(
	MPI_Status status; 
	status.MPI_ERROR = 0;
	uint n = data.num_rows();
	uint m = data.num_cols();
	
	double real_data[n*m];
	int error = MPI_Recv(&real_data, n*m, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
	NEGF_FASSERT(error==0, "An error (%d) occurred while receiving the real part of a matrix via MPI.",error);
	NEGF_ASSERT(status.MPI_ERROR == 0, "An error occurred while receiving the a real matrix via MPI.");
	
	for (uint ii=0; ii<n; ii++) {
		for (uint jj=0; jj<m; jj++) {
			data(ii+1,jj+1) = real_data[ii*m+jj];
		}
	}
);}


// --------------------------------------
// OTHER ROUTINES
// --------------------------------------

/** Get the rank
 * @return the communication rank of the calling process in the group MPI_COMM_WORLD */
int myMPI::get_rank()
{STACK_TRACE(
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	return my_rank;
);}

/** Get the number of processes
 * @return the number of processes in the group MPI_COMM_WORLD */
int myMPI::get_num_procs()
{STACK_TRACE(
	int num_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	return num_procs;
);}

/** Insert this routine as a blocker that waits until all other processes have reached it as well */
void myMPI::synchronize_processes() throw (Exception *)
{STACK_TRACE(
    MPI_Barrier(MPI_COMM_WORLD);

/*	int OK = 88;
	int tag = 88;
	int OK_plus_one = OK+1;
	int master = constants::mpi_master_rank;
	if(mpi->get_rank()==master) 
	{
		int check = 0;
		for (int cpu = 0; cpu < mpi->get_num_procs(); cpu++) {
			if (cpu==master) continue;
			mpi->recv(check, cpu, tag);
			NEGF_FASSERT(check==OK, "received %d instead of %d from process %d.", check, OK, cpu);
		}
		// re-broadcast that everybody is finished.
		mpi->broadcast(OK_plus_one, master);
	} else {
		mpi->send(OK, master, tag);
		int check = 0;
		mpi->broadcast(check, master);
		NEGF_FASSERT(check==OK+1, "received %d instead of %d from process 0.", check, OK+1)
	}
	logmsg->emit(LOG_INFO_L3, "Synchronized MPI processes.");*/
);}


/** Send an an array of antihermitian complex matrices (A=-A+) to some other process 
 *  This is done by creating real vectors and compressing them
 */
void myMPI::send_antihermitians(const vector<BMatc> & data, int & destination, int & tag) throw (Exception *)
{STACK_TRACE(
	double t1 = MPI_Wtime();

	unsigned long num_doubles = negf_math::check_and_get_num_doubles(data, true);
	const int max_distance = data[0].num_offdiags;
	int n = data[0].num_rows();
	int num_matrices = data.size();
	
	// ----------------------------------------------
	// create double arrays
	// ----------------------------------------------
	unsigned long count = 0;
	double * real_data = new double[num_doubles];
	double * imag_data = new double[num_doubles];
	for (int mm=0; mm<num_matrices; mm++) {
		for (int ii=0; ii<n; ii++) {
			for (int jj=max(0,ii-max_distance); jj<=ii; jj++) {
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
	unsigned char * real_compressed = new unsigned char[real_comp_size];	
	unsigned char * imag_compressed = new unsigned char[imag_comp_size];
	
	int clevel = 9;	// Z_DEFAULT_COMPRESSION, or between 0 and 9: 1 gives best speed, 9 gives best compression
	
	logmsg->emit_noendl(LOG_INFO_L3,"p%d REAL: uncompressed size:%d, maximum compressed size:%d \n    compressing...", this->get_rank(), num_chars, real_comp_size);
	int err = compress2(real_compressed, real_compressed_size, real_char, num_chars, clevel);	// zlib routine
	negf_math::check_compress_error(err, real_comp_size, num_chars);
	logmsg->emit(LOG_INFO_L3, "done. compressed size: %d",(*real_compressed_size));
	
	logmsg->emit_noendl(LOG_INFO_L3,"p%d IMAG: uncompressed size:%d, maximum compressed size:%d \n    compressing...", this->get_rank(), num_chars, imag_comp_size);
	err = compress2(imag_compressed, imag_compressed_size, imag_char, num_chars, clevel);	// zlib routine
	negf_math::check_compress_error(err, imag_comp_size, num_chars);
	logmsg->emit(LOG_INFO_L3, "done. compressed size: %d",(*imag_compressed_size));
	
	// ----------------------------------------------
	// send array lengths
	// ----------------------------------------------
	// integer is 4 bytes --> max. array size is 2^(4*8-1)
	NEGF_ASSERT(sizeof(int)==4, "sizeof(int)==4 expected.");
	unsigned long max_array_size = 2147483648UL; // this gives 268'435'456 doubles or, at an array size of 1000*1000, roughly 268*2>500 possible k-points
	NEGF_ASSERT(*real_compressed_size < max_array_size && *imag_compressed_size < max_array_size, "maximum array size exceeded.");
	int real_char_size = int(*real_compressed_size);
	int imag_char_size = int(*imag_compressed_size);
	
	double t2 = MPI_Wtime();
	
	MPI_Send(&real_char_size, 1, MPI_INT, destination, tag, MPI_COMM_WORLD);
	int tag2 = tag+1;
	MPI_Send(&imag_char_size, 1, MPI_INT, destination, tag2, MPI_COMM_WORLD);
	
	double t3 = MPI_Wtime();
		
	// ----------------------------------------------
	// send arrays
	// ----------------------------------------------
	MPI_Send(real_compressed, real_char_size, MPI_UNSIGNED_CHAR, destination, tag, MPI_COMM_WORLD);
	MPI_Send(imag_compressed, imag_char_size, MPI_UNSIGNED_CHAR, destination, tag2, MPI_COMM_WORLD);
		
	double t4 = MPI_Wtime();
	
	// --------------------------------------------------
	// clean up heap arrays
	// --------------------------------------------------
	delete [] real_data;
	delete [] imag_data;
	delete [] real_compressed;
	delete [] imag_compressed;
	
	double t5 = MPI_Wtime();
	
	double own_time = t2-t1 + t5-t4;
	double mixed_time = t3-t2;
	double mpi_time = t4-t3;
	if (this->get_rank()==11 || this->get_rank()==10) {
		logmsg->emit(LOG_INFO_L3,"p%d(%s)->p%d: own+mpi=%.5g, wait=%.5g, (%d reals, %d imags)", 
				this->get_rank(), this->get_hostname().c_str(), destination, own_time+mpi_time, mixed_time, real_char_size, imag_char_size);
	}
);}


/** Receive an array of antihermitian (A=-A+) complex matrices from a specific process
 *  Only the lower part is expected from the sender, the rest is reconstructed from
 *  A_ij = -A_ji^* = -(Re[A_ji] - i*Im[A_ji])
 *  the transmission is done in a compressed (zlib) format
 */
void myMPI::recv_antihermitians(vector<BMatc> & data, int & source, int & tag, 
				unsigned char * real_compressed,   unsigned char * imag_compressed,
				unsigned char * real_uncompressed, unsigned char * imag_uncompressed) throw (Exception *)
{STACK_TRACE(
	double t1 = MPI_Wtime();
	MPI_Status status; 
	status.MPI_ERROR = 0;
	int num_matrices = data.size();
	NEGF_ASSERT(num_matrices > 0, "matrix vector must have size >0");
	
	unsigned long num_doubles = negf_math::check_and_get_num_doubles(data, true);
	const int max_distance = data[0].num_offdiags;
	int n = data[0].num_rows();
	
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
	bool allocate_compressed_arrays = false;
	if (real_compressed==0 || imag_compressed==0) {
		allocate_compressed_arrays = true;
		real_compressed = new unsigned char[real_char_size];
		imag_compressed = new unsigned char[imag_char_size];
	}

	double t2 = MPI_Wtime();
	error = MPI_Recv(real_compressed, real_char_size, MPI_UNSIGNED_CHAR, source, tag, MPI_COMM_WORLD, &status);
	NEGF_FASSERT(error==0, "An error (%d) occurred while receiving the real part of a matrix via MPI.",error);
	NEGF_ASSERT(status.MPI_ERROR == 0, "An error occurred while receiving the real part of a complex matrix via MPI.");
	
	double t2b = MPI_Wtime();
	error = MPI_Recv(imag_compressed, imag_char_size, MPI_UNSIGNED_CHAR, source, tag2, MPI_COMM_WORLD, &status);
	NEGF_FASSERT(error==0, "An error (%d) occurred while receiving the imag part of a matrix via MPI.",error);
	NEGF_ASSERT(status.MPI_ERROR == 0, "An error occurred while receiving the imag part of a complex matrix via MPI.");
	
	double t3 = MPI_Wtime();
	// ------------------------------------
	// decompress
	// from the zlib homepage: "...the destination buffer, which must be large enough to hold the entire uncompressed data."
	// ------------------------------------
	int doublesize = sizeof(double);
	int charsize = sizeof(char);
	NEGF_ASSERT(charsize==1, "class is designed for sizeof(char)=1");
	unsigned long num_chars = num_doubles * doublesize;
	unsigned long real_decomp_size = num_chars/* + 100*/; // 100 should not even be necessary
	unsigned long imag_decomp_size = real_decomp_size; 

	bool allocate_uncompressed_arrays = false;
	if (real_uncompressed==0 || imag_uncompressed==0) {
		allocate_uncompressed_arrays = true;
		real_uncompressed = new unsigned char[real_decomp_size];
		imag_uncompressed = new unsigned char[imag_decomp_size];
	}
	unsigned long * real_decompressed_size = &real_decomp_size;
	unsigned long * imag_decompressed_size = &imag_decomp_size;
	
	logmsg->emit_noendl(LOG_INFO_L3,"decompressing REAL part...   ");
	error = uncompress(real_uncompressed, real_decompressed_size, real_compressed, real_char_size);
	negf_math::check_decompress_error(error, real_decomp_size, num_chars);
	logmsg->emit(LOG_INFO_L3, "done. decompressed size: %d",(*real_decompressed_size));
	
	logmsg->emit_noendl(LOG_INFO_L3,"decompressing IMAG part...   ");
	error = uncompress(imag_uncompressed, imag_decompressed_size, imag_compressed, imag_char_size);
	negf_math::check_decompress_error(error, imag_decomp_size, num_chars);
	logmsg->emit(LOG_INFO_L3, "done. decompressed size: %d",(*imag_decompressed_size));
	
	double t4 = MPI_Wtime();
	// ------------------------------------------------
	// cast into doubles 
	// ------------------------------------------------
	double * reals = (double *)(real_uncompressed);
	double * imags = (double *)(imag_uncompressed);
	
	// --------------------------------------------------
	// fill lower half (including diagonal)
	// --------------------------------------------------
	int count = 0;
	for (int mm=0; mm < num_matrices; mm++) {
		for (int ii=0; ii<n; ii++) {
			for (int jj=max(0,ii-max_distance); jj<=ii; jj++) {
//#ifdef USE_BANDED
//				// m!=n will probably be outside band
//				if (fabs(ii-jj) > data[mm].num_offdiags+1e-8) continue;
//#endif
				data[mm](ii+1,jj+1) = reals[count] + constants::imag_unit * imags[count];
				count++;
			}
/*			for (int jj=0; jj<max(0,ii-max_distance); jj++) {
#ifdef USE_BANDED
				// m!=n will probably be outside band
				if (fabs(ii-jj) > data[mm].num_offdiags+1e-8) continue;
#endif
				data[mm](ii+1,jj+1) = 0.0;
			}*/
		}
	}
	NEGF_ASSERT((unsigned long)(count)==num_doubles, "count==num_doubles");
		
	// --------------------------------------------------
	// reconstruct upper half
	// --------------------------------------------------
	for (int mm=0; mm < num_matrices; mm++) {
		for (int ii=1; ii<=n; ii++) {
/*			for (int jj=ii+1; jj<=n; jj++) {
#ifdef USE_BANDED
				// m!=n will probably be outside band
				if (fabs(ii-jj) > data[mm].num_offdiags+1e-8) continue;
#endif*/
			for (int jj=ii+1; jj<=min(n,ii+max_distance); jj++) {
				data[mm](ii,jj) = -data[mm](jj,ii).real() + constants::imag_unit * data[mm](jj,ii).imag();
			}
		}
	}
	
	// --------------------------------------------------
	// clean up heap arrays
	// --------------------------------------------------
	if (allocate_compressed_arrays) {
		delete [] real_compressed;
		delete [] imag_compressed;
	}
	if (allocate_uncompressed_arrays) {
		delete [] real_uncompressed;
		delete [] imag_uncompressed;
	}
	double t5 = MPI_Wtime();

	double       pre_time = t2-t1;
	double recv_time_real = t2b-t2;
	double recv_time_imag = t3-t2b;
	double    decomp_time = t4-t3;
	double     alloc_time = t5-t4;
	if (mpi->get_rank()==100) {
		logmsg->emit_all(LOG_INFO_L3,"\np100 recv_antihermitians timing: MPI_Recv (%ld=%ld+%ld chars) %.2e (%.2e+%.2e), rest %.2e", 
					real_char_size+imag_char_size, real_char_size, imag_char_size, 
					recv_time_real+recv_time_imag, recv_time_real, recv_time_imag, pre_time+decomp_time+alloc_time);
	}
);}


void myMPI::determine_senders_receivers(const vector< vector<int> > & processes_needed, const vector<bool> & process_was_computed, 
										vector<bool> & senders, vector<bool> & receivers)
{STACK_TRACE(
	NEGF_ASSERT(int(processes_needed.size())==this->get_num_procs() && int(process_was_computed.size())==this->get_num_procs()
			&& int(senders.size())==this->get_num_procs() && int(receivers.size())==this->get_num_procs(), "array sizes are wrong.");
	for (int pp=0; pp < this->get_num_procs(); pp++) 
	{
		// check if pp is not already engaged in communication and was not yet computed
		// if so, it's pp's turn to be computed and we add all necessary processes as sender
		if (!receivers[pp] && !senders[pp] && !process_was_computed[pp]) {
			// it is also possible that one of the needed processes is marked as a receiver already
			// (when such a needed process has a smaller rank than the currently investigated process)
			// we also need to exclude this case because a process cannot be receiver and sender at the same time
			bool one_of_the_needed_processes_is_receiver = false;
			for (uint ii=0; ii < processes_needed[pp].size(); ii++) {
				if (receivers[processes_needed[pp][ii]]) {
					one_of_the_needed_processes_is_receiver = true;
				}
			}
			if (one_of_the_needed_processes_is_receiver) {
				// do not mark as a receiver
				continue;
			}
			
			receivers[pp] = true;
			for (uint ii=0; ii < processes_needed[pp].size(); ii++) {
				senders[processes_needed[pp][ii]] = true;
			}
		}
	}
);}
		
